
# coding: utf-8

# Import Statements

# In[46]:


from threading import Thread
from subprocess import Popen
from time import sleep, time
import signal
import os
import re
import numpy as np

from bokeh.plotting import figure
from bokeh.io import show, output_notebook, push_notebook
from bokeh.models import ColumnDataSource, CDSView
from bokeh.transform import factor_cmap
from bokeh.palettes import brewer
from collections import OrderedDict

output_notebook()


# # Setup

# In[47]:


#set paths/parameters
fastq_dir = "/home/anthony/siRNA_profiling/FASTQ/"
output_dir = "/home/anthony/in-frame/pipeline_testing/"
available_cores = 48 #TODO:deal with situation if threads are fewer than number of samples.
star_ncRNA_dir = "/home/anthony/reference_genome/boris_profiling_annotations/hg38_Anno_BZ/ncRNA_STAR"
star_genome_dir = "/home/anthony/reference_genome/STAR_index/"
star_reporter_dir = "/home/anthony/siRNA_profiling/reporter_fasta/reporter_STAR_index/" # if no reporter, set to False

fastq_dir += "/" if not fastq_dir.endswith("/") else ""
output_dir += "/" if not output_dir.endswith("/") else ""


# In[48]:


#rename fasta files
barcode_name_mapping = {
    "GCCAAT":"ASCC3",
    "CAGATC":"ASCC2-TRIP4",
    "ACTTGA":"ASCC3-ASCC2-TRIP4-ZNF598",
    "GATCAG":"Scrambled",
    "GGCTAC":"Mock",
}

files = os.listdir(fastq_dir)
for file in files:
    match = re.search("(.+)(\.fastq|\.fastq\.gz)$", file)
    if match:
        prefix = match.group(1)
        suffix = match.group(2)
        for barcode in barcode_name_mapping:
            if barcode in prefix:
                os.rename(fastq_dir+"/"+file, fastq_dir+"/"+barcode_name_mapping[barcode]+suffix)
                break


# In[49]:


#set up folder structure for output
if os.access(output_dir, os.F_OK):
    if not os.access(output_dir+"output", os.F_OK):
        os.chdir(output_dir)
        os.mkdir("output")
        os.mkdir(output_dir+"/output/logs")
        os.mkdir(output_dir+"/output/logs/pipeline_completion")
        os.mkdir(output_dir+"/output/deduplicated")
        os.mkdir(output_dir+"/output/trimmed")
        os.mkdir(output_dir+"/output/ncRNA_aligned")
        if star_reporter_dir:
            os.mkdir(output_dir+"/output/reporter_aligned")
        os.mkdir(output_dir+"/output/genome_aligned")
    output_dir += "output/"
else:
    print("WARNING:", output_dir, "does not exist.")
if not os.access(output_dir, os.F_OK):
    print("WARNING:", fastq_dir, "does not exist.")

#check file endings for correct data processing.
filenames = os.listdir(fastq_dir)
fastq_files = []
for filename in filenames:
    if filename.endswith(".fastq") or filename.endswith(".fastq.gz"):
        fastq_files.append(filename)

#calculate the number of cores by dividing available cores by number of samples.
numcores = max(available_cores//len(fastq_files), 1)

print("Ready to process files:")
[print("\t\t"+file) for file in fastq_files]
print("\nUsing", str(numcores)+"/"+str(available_cores), "cores for each sample.")
print("Data will be output to", output_dir)


# # Main

# In[56]:


#defines process handler for running one sample through the pipeline
class run_one_sample():
    '''handles the running of one sample through the entire pipeline. 
    This class should be called by an individual thread. It will output stdout/stderr 
    readouts from the individual steps in the .../output/logs/ folder.
    
    It also uses a simple text log to pick up where the sample left off if the pipeline 
    fails. This log is stored in .../output/logs/pipeline_completion/
    Steps are defined in self.steps
    
    The self.check method will terminate the run and kill all running process groups if
    one of the processes in this thread returns an exit code other than 0, or if any
    other thread changes the global kill_pipeline variable to True. Default behavior is 
    to set kill_pipeline to True whenever any thread fails, thus halting processing of 
    all samples.
    
    On completion, each thread will add one to the global samples_done counter, to 
    let the main thread know when all are finished.
    
    TODO:use queues/messages instead of changing global variables
    '''
    def __init__(self, filename):
        global kill_pipeline
        global samples_done
        global star_ncRNA_dir
        global star_genome_dir
        global star_reporter_dir
        global fastq_dir
        global output_dir
        global numcores
        global failed_sample
        global thread_tracker
        self.star_ncRNA_dir = star_ncRNA_dir
        self.star_genome_dir = star_genome_dir
        self.star_reporter_dir = star_reporter_dir
        self.fastq_dir = fastq_dir
        self.filename = filename
        self.output_dir = output_dir
        self.numcores = numcores
        self.killflag = False
        match = re.search("(.+)(\.fastq|\.fastq\.gz)$", self.filename)
        self.prefix = match.group(1)
        self.suffix = match.group(2)
        if self.prefix not in thread_tracker:
            thread_tracker[self.prefix] = "initiate"
        self.steps = self.define_steps()
        
        with open(output_dir+"/logs/"+self.prefix+"_stdout.txt", "a") as self.stdout,             open(output_dir+"/logs/pipeline_completion/"+self.prefix+"_completion_log.txt", "a") as self.completion:
            try:
                for step_string in self.steps.keys():
                    if self.killflag == False and kill_pipeline == False:
                        function_name, function_args = self.steps.pop(step_string)
                        self.process = function_name(*function_args)
                        self.check(self.process)
                        self.completion.write(step_string+"\n")
                        self.completion.flush()
                        thread_tracker[self.prefix] = step_string
                if len(self.steps) == 0:
                    samples_done += 1
                        
            except:
                self.killflag = True
                kill_pipeline = True
                samples_done += 1 
                raise
                
    
    def define_steps(self):
        steps = OrderedDict([
            ("deduplicate", (self.deduplicate, (fastq_dir, output_dir, self.filename, numcores))),
            ("trim", (self.trim_reads, (output_dir, self.filename, numcores))),
            ("align_ncRNA", (self.align_ncRNA, (output_dir, self.filename, numcores))),
            ("align_reporter", (self.align_reporter, (output_dir, self.filename, numcores))),
            ("remove_reverse_reads", (self.samtools_remove_rv_reads, ("reporter_aligned", numcores))),
            ("sort_reporter_aligned", (self.samtools_sort, ("reporter_aligned", numcores))),
            ("index_reporter_aligned", (self.samtools_index, ("reporter_aligned", numcores))),
            ("align_genome", (self.align_genome, (output_dir, self.filename, numcores))),
            ("sort_genome_aligned", (self.samtools_sort, ("genome_aligned", numcores))),
            ("index_genome_aligned", (self.samtools_index, ("genome_aligned", numcores))),
        ])
        if not self.star_reporter_dir:    
            [steps.pop(step) for step in ["align_reporter", 
                                          "remove_reverse_reads", 
                                          "sort_reporter_aligned", 
                                          "index_reporter_aligned"]]
            
        if os.access(output_dir+"logs/pipeline_completion/"+self.prefix+"_completion_log.txt", os.F_OK):
            with open(output_dir+"logs/pipeline_completion/"+self.prefix+"_completion_log.txt", "r") as self.completion:
                for line in self.completion.readlines():
                    steps.pop(line.strip())
        return steps
    
    
    def deduplicate(self, fastq_dir, output_dir, filename, numcores):
        '''deduplicate ribosome profiling reads using dedupe.sh from 
        the BBTools suite.
        '''
        self.stdout.write("Deduplicate\n"+"".join(["-"]*20)+"\n\n")
        self.stdout.flush()
        dedupe_dir = output_dir+"deduplicated/"+self.prefix+"/"
        if not os.access(dedupe_dir, os.F_OK):
            os.mkdir(dedupe_dir)
        return Popen([
            "dedupe.sh",
            "in="+fastq_dir+filename,
            "out="+dedupe_dir+self.prefix+".deduped"+self.suffix,
            "absorbmatch=t", #absorb exact matches of contigs
            "absorbcontainment=f", #do not absorb full containments of contigs
            "absorbrc=f", #do not absorb reverse-compliments
            "threads="+str(numcores),
            "overwrite=t",
        ], stdout=self.stdout, stderr=self.stdout, preexec_fn=os.setsid)
            
            
    
    def trim_reads(self, output_dir, filename, numcores):
        '''Trim adapters and low quality regions from reads using bbduk.sh
        from the BBTools suite.
        '''
        self.stdout.write("\n\nTrim Reads\n"+"".join(["-"]*20)+"\n\n")
        self.stdout.flush()
        dedupe_dir = output_dir+"deduplicated/"+self.prefix+"/"
        trimmed_dir = output_dir+"trimmed/"+self.prefix+"/"
        if not os.access(trimmed_dir, os.F_OK):
            os.mkdir(trimmed_dir)
            os.mkdir(trimmed_dir+"failedQC")
        return Popen([
            "bbduk.sh",
            "in="+dedupe_dir+self.prefix+".deduped"+self.suffix,
            "out="+trimmed_dir+self.prefix+".trimmed.fastq",
            "outm="+trimmed_dir+"failedQC/"+self.prefix+".failedQC"+self.suffix,
            "rpkm="+trimmed_dir+"rpkm.txt",
            "refstats="+trimmed_dir+"trimming_stats.txt",
            "literal=NNNNNNCACTCGGGCACCAAGGAC",
            "k=24", # this parameter sets the minimum kmer being trimmed. 
                                  #Longer = more specific, shorter = more sensitive
            "mink=8", #includes truncations of the kmers down to 8
            "mm=f", #do not ignore middle base mismatch of kmer
            "rcomp=f", #do not allow reverse complement kmer matches
            "copyundefined=t",
            "ktrim=r",
            "forcetrimleft=4", #removes random barcode on left of reads.
            "minavgquality=10",
            "minlength=10",
            "threads="+str(numcores),
            "overwrite=t",
        ],
        stdout=self.stdout, stderr=self.stdout, preexec_fn=os.setsid)
            

    def align_ncRNA(self, output_dir, filename, numcores):
        '''Align reads to ncRNA using STAR. ncRNA fasta sequences from Boris.
        Output unaligned reads.
        '''
        self.stdout.write("\n\nAlign to ncRNA\n"+"".join(["-"]*20)+"\n\n")
        self.stdout.flush()
        trimmed_dir = output_dir+"trimmed/"+self.prefix+"/"
        ncRNA_aligned_dir = output_dir+"ncRNA_aligned/"+self.prefix+"/"
        if not os.access(ncRNA_aligned_dir, os.F_OK):
            os.mkdir(ncRNA_aligned_dir)
        command = [
            "STAR",
            "--runThreadN", str(numcores),
            "--genomeDir", self.star_ncRNA_dir,
            "--readFilesIn", trimmed_dir+self.prefix+".trimmed.fastq",
            "--outFileNamePrefix", ncRNA_aligned_dir+self.prefix+"_",
            "--outSAMtype", "BAM", "Unsorted",
            "--outReadsUnmapped", "Fastx",
            "--alignSJDBoverhangMin", "1",
            "--alignSJoverhangMin", "8",
            "--outFilterMultimapNmax", "20",
            "--outFilterType", "BySJout",
        ]
        return Popen(command, stderr=self.stdout, stdout=self.stdout, 
                     preexec_fn=os.setsid)
            
    
    def align_reporter(self, output_dir, filename, numcores):
        '''Align reads to reporter sequence using STAR.
        Output unaligned reads.
        '''
        if self.star_reporter_dir:
            self.stdout.write("\n\nAlign to reporter\n"+"".join(["-"]*20)+"\n\n")
            self.stdout.flush()
            ncRNA_aligned_dir = output_dir+"ncRNA_aligned/"+self.prefix+"/"
            reporter_aligned_dir = output_dir+"reporter_aligned/"+self.prefix+"/"
            if not os.access(reporter_aligned_dir, os.F_OK):
                os.mkdir(reporter_aligned_dir)
            command = [
                "STAR",
                "--runThreadN", str(numcores),
                "--genomeDir", self.star_reporter_dir,
                "--readFilesIn", ncRNA_aligned_dir+self.prefix+"_Unmapped.out.mate1",
                "--outFileNamePrefix", reporter_aligned_dir+self.prefix+"_",
                "--outSAMtype", "BAM", "Unsorted",
                "--outReadsUnmapped", "Fastx",
                "--alignSJDBoverhangMin", "1",
                "--alignSJoverhangMin", "8",
                "--outFilterMultimapNmax", "1",
                "--outFilterType", "BySJout",
            ]
            return Popen(command, stderr=self.stdout, stdout=self.stdout, 
                         preexec_fn=os.setsid)
                
    
    
    def align_genome(self, output_dir, filename, numcores):
        '''Align remaining reads to genome.
        '''
        self.stdout.write("\n\nAlign to genome\n"+"".join(["-"]*20)+"\n\n")
        self.stdout.flush()
        if self.star_reporter_dir:
            previous_aligned_dir = output_dir+"reporter_aligned/"+self.prefix+"/"
        else:
            previous_aligned_dir = output_dir+"ncRNA_aligned/"+self.prefix+"/"
        tx_aligned_dir = output_dir+"genome_aligned/"+self.prefix+"/"
        if not os.access(tx_aligned_dir, os.F_OK):
            os.mkdir(tx_aligned_dir)
        command = [
            "STAR",
            "--runThreadN", str(self.numcores),
            "--genomeDir", self.star_genome_dir,
            "--readFilesIn", previous_aligned_dir+self.prefix+"_Unmapped.out.mate1",
            "--outFileNamePrefix", tx_aligned_dir+self.prefix+"_",
            "--outSAMtype", "BAM", "Unsorted",
            "--outReadsUnmapped", "Fastx",
            "--alignSJDBoverhangMin", "1",
            "--alignSJoverhangMin", "8",
            "--outFilterMultimapNmax", "1", #how many multimap sites allowed for read
            "--outSAMmultNmax", "1", #how many map sites to write to output for each read
            "--outMultimapperOrder", "Random", #assign read to random alignment if multimapper
            "--outFilterType", "BySJout",
        ]
        return Popen(command, stderr=self.stdout, stdout=self.stdout,
                     preexec_fn=os.setsid)
            
    
    def samtools_sort(self, input_dir, numcores):
        '''Sort BAM file from STAR ouput.
        '''
        self.stdout.write("\n\nSort "+input_dir+" BAM file\n"+"".join(["-"]*20)+"\n\n")
        self.stdout.flush()
        if input_dir == "reporter_aligned":
            aligned_suffix = "_Aligned.out.plusstrand"
        else:
            aligned_suffix = "_Aligned.out"
        return Popen([
            "samtools",
            "sort",
            "-@", str(numcores),
            self.output_dir+input_dir+"/"+self.prefix+"/"+self.prefix+aligned_suffix+".bam",
            "-o", self.output_dir+input_dir+"/"+self.prefix+"/"+self.prefix+aligned_suffix+".sorted.bam"
        ], stderr=self.stdout, stdout=self.stdout, preexec_fn=os.setsid)
            
    
    def samtools_index(self, input_dir, numcores):
        '''Index BAM file from STAR output
        '''
        self.stdout.write("\n\nIndex "+input_dir+" BAM file\n"+"".join(["-"]*20)+"\n\n")
        self.stdout.flush()
        if input_dir == "reporter_aligned":
            aligned_suffix = "_Aligned.out.plusstrand"
        else:
            aligned_suffix = "_Aligned.out"
        try:
            os.remove(self.output_dir+input_dir+"/"+self.prefix+"/"+self.prefix+aligned_suffix+".bam",)
        except FileNotFoundError:
            pass
        return Popen([
            "samtools",
            "index",
            "-@", str(numcores),
            self.output_dir+input_dir+"/"+self.prefix+"/"+self.prefix+aligned_suffix+".sorted.bam"
        ], stderr=self.stdout, stdout=self.stdout, preexec_fn=os.setsid)
            
        
    def samtools_remove_rv_reads(self, input_dir, numcores):
        '''Remove reverse reads aligned to the reporter minus strand with samtools
        '''
        if self.star_reporter_dir:
            self.stdout.write("\n\nRemove reporter reverse strand reads\n"+"".join(["-"]*20)+"\n\n")
            self.stdout.flush()
            return Popen([
                "samtools",
                "view",
                "-@", str(numcores),
                "-F", "0x10", #only include reads with this flag, which means plus strand alignment
                "-o",
                self.output_dir+input_dir+"/"+self.prefix+"/"+self.prefix+"_Aligned.out.plusstrand.bam",
                self.output_dir+input_dir+"/"+self.prefix+"/"+self.prefix+"_Aligned.out.bam"
            ], stderr=self.stdout, stdout=self.stdout, preexec_fn=os.setsid)
                
            
        
    def check(self, proc):
        '''Poll Popen processes returned by each method to determine if a nonzero
        error code was returned. If so, kill process group and set global kill_pipeline
        to True, signalling to other processing threads to shut down as well.
        Polling happens every 1 second.
        '''
        global kill_pipeline
        global thread_tracker
        exit_code = proc.poll()
        while exit_code == None:
            sleep(1)
            if kill_pipeline == True:
                self.killflag = True
            if self.killflag == True:
                os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            exit_code = proc.poll()
        else:
            if exit_code != 0:
                self.killflag = True
                kill_pipeline = True


# In[58]:


start_time = time()

kill_pipeline = False
samples_done = 0
sample_runs = {}
thread_tracker = OrderedDict()

#start threads for running each sample through the pipeline
for filename in fastq_files:
    sample_runs[filename] = Thread(target=run_one_sample, args=(filename,))
    sample_runs[filename].start()
    
print("Pipeline Running...")

if star_reporter_dir:
    steps_mapper = OrderedDict([ #offset by one because if it's recorded in the log, it's already done
        ("initiate", "Deduplicating"), 
        ("deduplicate", "Trimming Adapter"),
        ("trim", "Aligning to ncRNA"),
        ("align_ncRNA", "Aligning to Reporter"),
        ("align_reporter", "Filtering reporter reverse reads"),
        ("remove_reverse_reads", "Sorting reporter-aligned reads"),
        ("sort_reporter_aligned", "Indexing reporter-aligned reads"),
        ("index_reporter_aligned", "Aligning to Genome"), 
        ("align_genome", "Sorting reporter aligned reads"), 
        ("sort_genome_aligned", "Indexing Genome Aligned Reads"), 
        ("index_genome_aligned", "Done"),
    ])

steps = ["Start"]+[steps_mapper[step] for step in steps_mapper]
prefixes = list(thread_tracker.keys())

p = figure(height=200, width=800, y_range=prefixes, background_fill_color="lightgray", x_range=steps, title="Pipeline Progress", tools=[])
p.xgrid.grid_line_color = "gray"
p.xaxis.major_label_orientation = np.pi/8

p.ygrid.visible = False
source = ColumnDataSource(data={"y":prefixes, "right":["Start"]*len(prefixes), "height":[0.5]*len(prefixes), "left":[0]*len(prefixes)})
x = p.hbar(y="y", right="right", height="height", left="left", color=factor_cmap(field_name="y", palette=brewer["Spectral"][len(prefixes)], factors=prefixes), 
           line_color="black", source=source)

show(p, notebook_handle=True)

#Check every 1 sec whether the pipeline is finished and report wether it's been terminated.
try:
    while samples_done != len(fastq_files):
        for sample in thread_tracker:
            source.data = {"y":prefixes, "right":[steps_mapper[thread_tracker[sample]] for sample in thread_tracker], "height":[0.5]*len(prefixes), "left":[0]*len(prefixes)}
            x.view = CDSView(source=x.data_source)
            push_notebook()
        sleep(1)
    else:
        if kill_pipeline == False:
            source.data = {"y":prefixes, "right":[steps_mapper[thread_tracker[sample]] for sample in thread_tracker], "height":[0.5]*len(prefixes), "left":[0]*len(prefixes)}
            x.view = CDSView(source=x.data_source)
            push_notebook()
            print("Pipeline finished successfully!")
        else:
            print("Run terminated. Check for errors")
except:
    kill_pipeline = True  #allows KeyboardInterrupt to kill pipeline
    print("Run terminated. Check for errors")
    raise

    

runtime = time() - start_time
if runtime > 60:
    mins = round(runtime/60, 2)
    print("Run time:", mins, "minutes")
else:
    secs = round(runtime, 2)
    print("Run time:", secs, "seconds")


# In[59]:


def read_fate(sample_name):
    global output_dir
    with open(output_dir+"logs/"+sample_name+"_stdout.txt", "r") as log,     open(output_dir+"ncRNA_aligned/"+sample_name+"/"+sample_name+"_Log.final.out", "r") as ncRNA_log,     open(output_dir+"reporter_aligned/"+sample_name+"/"+sample_name+"_Log.final.out", "r") as reporter_log,     open(output_dir+"genome_aligned/"+sample_name+"/"+sample_name+"_Log.final.out", "r") as genome_log:
        file = log.readlines()[::-1]
        dedupe_end = file.index("Trim Reads\n")
        dedupe_start = file.index("--------------------\n", dedupe_end)
        for line in file[dedupe_end:dedupe_start]:
            if line.startswith("Input"):
                input_reads = int(line.split("\t")[1].split(" ")[0])
            elif line.startswith("Result"):
                deduplicated_reads = int(line.split("\t")[1].split(" ")[0])
        trim_end = file.index("Align to ncRNA\n")
        trim_start = file.index("--------------------\n", trim_end)
        for line in file[trim_end:trim_start]:
            if line.startswith("Result"):
                trimmed_reads = int(line.split("\t")[1].split(" ")[0])
        file = ncRNA_log.readlines()
        for line in file:
            if line.strip().startswith("Uniquely mapped reads number"):
                ncRNA_mapped_reads = int(line.split("|")[1].strip())
            elif line.strip().startswith("Number of reads mapped to multiple loci"):
                ncRNA_mapped_reads += int(line.split("|")[1].strip())
        file = reporter_log.readlines()
        for line in file:
            if line.strip().startswith("Uniquely mapped reads number"):
                reporter_mapped = int(line.split("|")[1].strip())
            elif line.strip().startswith("Number of reads mapped to multiple loci"):
                reporter_mapped += int(line.split("|")[1].strip())
        file = genome_log.readlines()
        for line in file:
            if line.strip().startswith("Number of input reads"):
                genome_input_reads = int(line.split("|")[1].strip())
            elif line.strip().startswith("Uniquely mapped reads number"):
                genome_mapped = int(line.split("|")[1].strip())
            elif line.strip().startswith("Number of reads mapped to multiple loci"):
                genome_mapped += int(line.split("|")[1].strip())
        unmapped_reads = genome_input_reads - genome_mapped
        print("Library was", str(round(deduplicated_reads/input_reads*100, 2))+"%", "unique.")
        print("Of those,", str(round(trimmed_reads/deduplicated_reads*100, 2))+"%", "survived trimming.")
        print(str(round(ncRNA_mapped_reads/trimmed_reads*100, 2))+"%", "mapped to ncRNA.")
        print(str(round(reporter_mapped/trimmed_reads*100, 2))+"%", "mapped to reporter.")
        print(str(round(genome_mapped/trimmed_reads*100, 2))+"%", "mapped to the genome.")
        print(str(round(unmapped_reads/trimmed_reads*100, 2))+"%", "remained unmapped.")
                
for barcode in barcode_name_mapping:
    print(barcode_name_mapping[barcode])
    read_fate(barcode_name_mapping[barcode])
    print("\n\n")

