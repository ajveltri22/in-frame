
# coding: utf-8

# Import Statements

# In[115]:


from threading import Thread
from subprocess import Popen
from time import sleep, time
import signal
import os
import re


# # Setup

# In[117]:


fastq_dir = "/home/anthony/cap_binding_pulldown/profiling/fastq/"
output_dir = "/home/anthony/profiling_pipeline/sample_data"
available_cores = 96
star_ncRNA_dir = "/home/anthony/reference_genome/STAR_index_ncRNA/"
star_genome_dir = "/home/anthony/reference_genome/STAR_index_ncRNA/"

fastq_dir += "/" if not fastq_dir.endswith("/") else ""
output_dir += "/" if not output_dir.endswith("/") else ""

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

# In[118]:


#defines process handler for running one sample through the pipeline
class run_one_sample():
    '''handles the running of one sample through the entire pipeline. 
    This class should be called by an individual thread. It will output stdout/stderr 
    readouts from the individual steps in the .../output/logs/ folder.
    
    It also uses a simple text log to pick up where the sample left off if the pipeline 
    fails. This log is stored in .../output/logs/pipeline_completion/
    
    The self.check method will terminate the run and kill all running process groups if
    one of the processes in this thread returns an exit code other than 0, or if any
    other thread changes the global kill_pipeline variable to True. Default behavior is 
    to set kill_pipeline to True whenever any thread fails, thus halting processing of 
    all samples.
    
    On completion, each thread will add one to the global samples_done counter, to 
    let the main thread know when all are finished.
    '''
    def __init__(self, filename):
        global kill_pipeline
        global samples_done
        global star_ncRNA_dir
        global star_genome_dir
        global fastq_dir
        global output_dir
        global numcores
        self.star_ncRNA_dir = star_ncRNA_dir
        self.star_genome_dir = star_genome_dir
        self.fastq_dir = fastq_dir
        self.filename = filename
        self.output_dir = output_dir
        self.numcores = numcores
        self.killflag = False
        match = re.search("(.+)(\.fastq|\.fastq\.gz)$", self.filename)
        self.prefix = match.group(1)
        self.suffix = match.group(2)
        self.steps = ["initiate", "deduplicate", "trim", "align_ncRNA", "align_genome", "sort", "index"] 
        if not os.access(output_dir+"logs/pipeline_completion/"+self.prefix+"_completion_log.txt", os.F_OK):
            self.current_step = 0
        else:
            self.completion = open(output_dir+"/logs/"+self.prefix+"_completion_log.txt", "r")
            self.current_step = self.steps.index(self.completion.readlines()[-1].strip())
            print("current_step", self.current_step)
            self.completion.close()
        
        with open(output_dir+"/logs/"+self.prefix+"_stdout.txt", "a") as self.stdout,             open(output_dir+"/logs/pipeline_completion/"+self.prefix+"_completion_log.txt", "a") as self.completion:
            try:
                self.proc = self.deduplicate(fastq_dir, output_dir, self.filename, numcores)
                self.check(self.proc)
                self.proc = self.trim_reads(output_dir, self.filename, numcores)
                self.check(self.proc)
                self.proc = self.align_ncRNA(output_dir, self.filename, numcores)
                self.check(self.proc)
                self.proc = self.align_genome(output_dir, self.filename, numcores)
                self.check(self.proc)
                self.proc = self.samtools_sort(output_dir, numcores)
                self.check(self.proc)
                self.proc = self.samtools_index(output_dir)
                self.check(self.proc)
                samples_done += 1
            except:
                self.killflag = True
                kill_pipeline = True
                samples_done += 1 
                try:
                    os.killpg(os.getpgid(self.proc.pid), signal.SIGTERM)
                except ProcessLookupError:
                    pass
                raise
    
    def deduplicate(self, fastq_dir, output_dir, filename, numcores):
        '''deduplicate ribosome profiling reads using dedupe.sh from 
        the BBTools suite.
        '''
        if self.killflag == False and self.current_step == 0:
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
            ], stdout=self.stdout, stderr=self.stdout, preexec_fn=os.setsid)
        else: 
            return 0
    
    def trim_reads(self, output_dir, filename, numcores):
        '''Trim adapters and low quality regions from reads using bbduk.sh
        from the BBTools suite.
        '''
        if self.killflag == False and self.current_step == 1:
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
                "out="+trimmed_dir+self.prefix+".trimmed"+self.suffix,
                "outm="+trimmed_dir+"failedQC/"+self.prefix+".failedQC"+self.suffix,
                "refstats="+trimmed_dir+"trimming_stats.txt",
                "literal=NNNNNNCACTCGGGCACCAAGGAC",
                "k=12", # this parameter sets the minimum kmer being trimmed. 
                                      #Longer = more specific, shorter = more sensitive
                "copyundefined=t",
                "ktrim=r",
                "forcetrimleft=4", #removes random barcode on left of reads.
                "minavgquality=10",
                "minlength=10",
                "threads="+str(numcores),
            ],
            stdout=self.stdout, stderr=self.stdout, preexec_fn=os.setsid)
        else:
            return 0

    def align_ncRNA(self, output_dir, filename, numcores):
        '''Align reads to ncRNA using STAR. ncRNA fasta downloaded from Ensembl.
        Output unaligned reads.
        '''
        if self.killflag == False and self.current_step == 2:
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
                "--readFilesIn", trimmed_dir+self.prefix+".trimmed"+self.suffix,
                "--outFileNamePrefix", ncRNA_aligned_dir+self.prefix+"_",
                "--outSAMtype", "BAM", "Unsorted",
                "--outReadsUnmapped", "Fastx",
            ]
            if self.suffix == ".fastq.gz":
                command.extend(["--readFilesCommand", "gunzip"])
            return Popen(command, stderr=self.stdout, stdout=self.stdout, 
                         preexec_fn=os.setsid)
        else:
            return 0
    
    def align_genome(self, output_dir, filename, numcores):
        '''Align remaining reads to genome.
        '''
        if self.killflag == False and self.current_step == 3:
            self.stdout.write("\n\nAlign to genome\n"+"".join(["-"]*20)+"\n\n")
            self.stdout.flush()
            ncRNA_aligned_dir = output_dir+"ncRNA_aligned/"+self.prefix+"/"
            tx_aligned_dir = output_dir+"genome_aligned/"+self.prefix+"/"
            if not os.access(tx_aligned_dir, os.F_OK):
                os.mkdir(tx_aligned_dir)
            command = [
                "STAR",
                "--runThreadN", str(self.numcores),
                "--genomeDir", self.star_genome_dir,
                "--readFilesIn", ncRNA_aligned_dir+self.prefix+"_Unmapped.out.mate1",
                "--outFileNamePrefix", tx_aligned_dir+self.prefix+"_",
                "--outSAMtype", "BAM", "Unsorted",
                "--outReadsUnmapped", "Fastx",
            ]
            return Popen(command, stderr=self.stdout, stdout=self.stdout,
                         preexec_fn=os.setsid)
        else:
            return 0
    
    def samtools_sort(self, output_dir, numcores):
        '''Sort BAM file from STAR ouput.
        '''
        if self.killflag == False and self.current_step == 4:
            self.stdout.write("\n\nSort BAM file\n"+"".join(["-"]*20)+"\n\n")
            self.stdout.flush()
            return Popen([
                "samtools",
                "sort",
                "-@", str(numcores),
                "-f", #use prefix as full file name
                output_dir+"genome_aligned/"+self.prefix+"/"+self.prefix+"_Aligned.out.bam",
                output_dir+"genome_aligned/"+self.prefix+"/"+self.prefix+"_Aligned.out.sorted.bam"
            ], stderr=self.stdout, stdout=self.stdout, preexec_fn=os.setsid)
        else:
            return 0
    
    def samtools_index(self, output_dir):
        '''Index BAM file from STAR output
        '''
        if self.killflag == False and self.current_step == 5:
            self.stdout.write("\n\nIndex BAM file\n"+"".join(["-"]*20)+"\n\n")
            self.stdout.flush()
            try:
                os.remove(output_dir+"genome_aligned/"+self.prefix+"/"+self.prefix+"_Aligned.out.bam",)
            except FileNotFoundError:
                pass
            return Popen([
                "samtools",
                "index",
                output_dir+"genome_aligned/"+self.prefix+"/"+self.prefix+"_Aligned.out.sorted.bam"
            ], stderr=self.stdout, stdout=self.stdout, preexec_fn=os.setsid)
        else:
            return 0
        
        
    def check(self, proc):
        '''Poll Popen processes returned by each method to determine if a nonzero
        error code was returned. If so, kill process group and set global kill_pipeline
        to True, signalling to other processing threads to shut down as well.
        Polling happens every 1 second.
        '''
        global kill_pipeline
        try:
            exit_code = proc.poll()
            while exit_code == None:
                exit_code = proc.poll()
                sleep(1)
                if kill_pipeline == True:
                    self.killflag = True
                if self.killflag == True:
                    os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
            else:
                if exit_code != 0:
                    self.killflag = True
                    kill_pipeline = True
                else:
                    self.current_step += 1
                    self.completion.write(self.steps[self.current_step]+"\n")
        except (AttributeError, ProcessLookupError): 
            #necessary to catch errors from methods returning 0 if self.killflag = True
            pass


# In[119]:


start_time = time()

kill_pipeline = False
samples_done = 0
sample_runs = {}

#start threads for running each sample through the pipeline
for filename in fastq_files:
    sample_runs[filename] = Thread(target=run_one_sample, args=(filename,))
    sample_runs[filename].start()
    
print("Pipeline Running...")

#Check every 1 sec whether the pipeline is finished and report wether it's been terminated.
try:
    while samples_done != len(fastq_files):
        sleep(1)
    else:
        if kill_pipeline == False:
            print("Pipeline finished successfully!")
        else:
            print("Run terminated. Check for errors")
except:
    kill_pipeline = True  #allows KeyboardInterrupt to kill pipeline
    raise

runtime = time() - start_time
if runtime > 60:
    mins = round(runtime/60, 2)
    print("Run time:", mins, "minutes")
else:
    secs = round(runtime, 2)
    print("Run time:", secs, "seconds")

