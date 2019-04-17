
# coding: utf-8

# # Profiling Metagene Generator  
# 
# Last updated 24 July 2018
# 
# ---
# **TODO:**  
# + Metagene Script
# 
#  - intron and exon reads
#  - how to deal with multimapping reads
#  - rethink all_[word]_df names and column names.
#  - normalization of codons ("codon centric" or "gene centric")
#  - allow genes without CDS - lincRNAs, etc.
#  - ~~Add individual codon analysis~~
#  - ~~normalization of genes by length (done for start_aligned)~~
#  - ~~add shifting to align to APE sites~~
#  - ~~read distribution plots~~
#  - ~~Just run pysam once~~
#  - ~~Density heatmaps~~
#  - ~~fix normalization of read sizes~~
# + pipeline
#  - dashboard allowing clicking on specific genes
#  - upstream steps (alignment/trimming)
#  - logging
#  - read length distribution plot for all reads
# ---
# 

# ## Import Statements

# In[1]:


import pandas as pd #used for constructing dataframes
import numpy as np
from multiprocessing import Pool #Multithreading pysam read lookup 
import multiprocessing as mp #Multithreading for metaplot summation
import pickle #used to save generated plots
import itertools

#plotting with Bokeh:
from bokeh.io import output_notebook, push_notebook, export, curdoc
from bokeh.plotting import figure, show
from bokeh.models import ColumnDataSource, CDSView, HoverTool, LinearColorMapper, FixedTicker, ColorBar, NumeralTickFormatter
from bokeh.layouts import column, gridplot
from bokeh.transform import linear_cmap
import colorcet as cc
output_notebook()

#adding iPython interaction to bokeh plots
from IPython.display import display
from ipywidgets import interact, IntSlider, ToggleButtons, SelectionSlider

#Reading bam files
import pysam

from time import strftime

prefix = "genome_aligned"


# ## Loading GTF annotations and filtering
# This cell loads a GTF file containing annotations and filters them for a subset of transcripts to use in metagene averaging. The hg38 Ensembl annotation GTF is used here and filtered down to those transcripts that are in the Gencode Basic set. 
# 
# In principle any GTF file can be used with some modification to this cell. The pandas column labels that are necessary downstream are: 
#     "feature": containing the values "transcript", "CDS", and "exon"
#     "transcript_id": containing an ID that associates "CDS" and "exon" features to a particular "transcript" feature
#     "start": containing the first position of the feature (inclusive)
#     "end": containing the last position of the feature (inclusive)
#     "strand": containing the value "+" or "-"
# 
# Here, the other columns are used for filtering down to Gencode Basic transcripts and one transcript per gene (with the highest Ensembl TSL number).

# In[9]:


pickled_gtf_df = "./parsed_gtf_20190405220439.pkl" #set to False to load from GTF file, set to filepath to restore pickled gtf_df
if not pickled_gtf_df:
    print("No pickled gtf_df set. Loading from GTF file")
    with open("/home/anthony/reference_genome/Homo_sapiens.GRCh38.92.gtf", "r") as inGTF:

        gtf_dict = {"chrom":[], "feature":[], "start":[], "end":[], 
                    "gene_id":[], "transcript_id":[], "exon_number":[], 
                    "strand":[], "annotation_type":[], "support_level":[],
                    "biotype":[]
                   }
        for line in inGTF:
            try:
                if not line.startswith("#"):
                    splitline = line.split("\t")
                    chrom = splitline[0]
                    feature = splitline[2]
                    if feature == "stop_codon": 
                        feature = "CDS"
                    if feature not in ["transcript", "CDS", "exon"] or chrom not in [str(num) for num in range(1,24)]+["X", "Y"]:
                        continue

                    if not 'tag "basic"' in line:
                        continue
                    start = int(splitline[3])
                    end = int(splitline[4])
                    gene_id = splitline[8].split("gene_id ")[1].split('"')[1]
                    transcript_id = splitline[8].split("transcript_id ")[1].split('"')[1]
                    biotype = splitline[8].split("transcript_biotype ")[1].split('"')[1]
                    try:
                        support_level = splitline[8].split("transcript_support_level")[1].split('"')[1]
                    except IndexError:
                        support_level = np.nan
                    annotation_type = splitline[1]
                    if feature == "exon":
                        exon_number = int(splitline[8].split("exon_number ")[1].split('"')[1])
                    else:
                        exon_number = np.nan
                    strand = splitline[6]
                    gtf_dict["chrom"].append(chrom)
                    gtf_dict["feature"].append(feature)
                    gtf_dict["start"].append(start)
                    gtf_dict["end"].append(end)
                    gtf_dict["gene_id"].append(gene_id)
                    gtf_dict["transcript_id"].append(transcript_id)
                    gtf_dict["exon_number"].append(exon_number)
                    gtf_dict["strand"].append(strand)
                    gtf_dict["annotation_type"].append(annotation_type)
                    gtf_dict["support_level"].append(support_level)
                    gtf_dict["biotype"].append(biotype)
            except:
                print(line)
                raise
    gtf_df = pd.DataFrame(gtf_dict)

    transcripts_deduped = gtf_df[(gtf_df.feature == "transcript") & (gtf_df.biotype == "protein_coding")
                                ].sort_values(by=["chrom", "gene_id", "support_level", "transcript_id", "start"]
                                             ).drop_duplicates(subset=["gene_id"], keep="first"
                                                              ).transcript_id
    gtf_df = gtf_df[gtf_df.transcript_id.isin(transcripts_deduped)
                   ].sort_values(by=["chrom", "gene_id", "support_level", "transcript_id", "start"])

    gtf_df["length"] = gtf_df.end.apply(int) - gtf_df.start.apply(int)
    max_gene = max(gtf_df.length)

    print("GTF loaded:", len(transcripts_deduped), "transcripts")

    #load dict of dicts containing spliced CDS coordinates for each gene
    # e.g. {"ENST00000321256":{"ATG":[1, 67, 160, ...], "ACG":[16, 19, ...], ...}, ...}
    codon_positions = pickle.load(open("/home/anthony/test_code/global_codon_positions.pkl", "rb"))
    gtf_df = gtf_df[gtf_df.transcript_id.isin(codon_positions.keys())]
    #bases = "TCAG"
    #all_codons = [a + b + c for a in bases for b in bases for c in bases] #create list of possible codons
    print("Loaded codon positions pkl")
    print("Filtered gtf_df down to", len(gtf_df[gtf_df.feature == "transcript"]), "matching transcripts.")


    filename = "./parsed_gtf_"+strftime("%Y%m%d%H%M%S")+".pkl"
    pickle.dump(gtf_df, open(filename, "wb"))
    print("Pickling gtf_df to",filename)

else:
    #load pickled gtf_df
    gtf_df = pickle.load(open(pickled_gtf_df, "rb"))
    print("Pickled gtf_df loaded:", len(gtf_df[gtf_df.feature == "transcript"]), "transcripts")
    #load dict of dicts containing spliced CDS coordinates for each gene
    # e.g. {"ENST00000321256":{"ATG":[1, 67, 160, ...], "ACG":[16, 19, ...], ...}, ...}
    #codon_positions = pickle.load(open("/home/anthony/test_code/global_codon_positions.pkl", "rb"))
    print("Loaded codon positions pkl")


# ### Downsample gene list for faster running/testing

# In[3]:


downsample_frac = 1 #input a value between 0 and 1 to downsample. Setting this to 1 will include all transcripts (no downsampling)
print(len(gtf_df[gtf_df.feature == "transcript"].sample(frac=downsample_frac, axis=0))) #print the number of genes in the downsampled dataframe
output_dict = {}


# ## Multithreaded density addition

# ### Define density calculation and summation functions
# ```count_density_2``` is the workhorse of the script that calculates normalized densities for each transcript  
# ```sum_density sums``` up all the results of ```count_density_2```

# In[10]:


pd.options.mode.chained_assignment = None  # default='warn'
def get_density(transcript):
    '''takes in a transcript line of gtf_df and creates a tx_density object. These will then be saved to a pkl file for fast parsing of downstream calculations.
    '''
    global prefix
    global gtf_df
    bamfile = pysam.AlignmentFile("/home/anthony/profiling_pipeline/sample_data/"+prefix+".sorted.bam", "rb")
    
    chrom = transcript[1].chrom
    strand = transcript[1].strand
    transcript_id = transcript[1].transcript_id
    tx_cds_exons_df = gtf_df[gtf_df.transcript_id == transcript_id] #forms a dataframe containing transcript, CDS, 
                                                                                  # and exon entries from gtf_df
    
    tx = tx_cds_exons_df[tx_cds_exons_df.feature == "transcript"] #transcript line of dataframe
    tx_left = tx.start.iloc[0] #Leftmost genome coordinate of transcript (could be 3' or 5' end depending on strand)
    tx_right = tx.end.iloc[0] #Rightmost genome coordinate of transcript

    exons = tx_cds_exons_df[tx_cds_exons_df.feature == "exon"] #Exon lines of dataframe
    tx_len = sum([(exon[1].end + 1) - exon[1].start for exon in exons.iterrows()]) #length of spliced transcript feature
    
    cds = tx_cds_exons_df[tx_cds_exons_df.feature == "CDS"] #CDS lines of dataframe
    cds_left_genomecoords = cds.iloc[0].start
    cds_right_genomecoords = cds.iloc[-1].end
    
    #creates a dataframe of all information on reads mapping anywhere in the transcript.
    read_iter = bamfile.fetch(chrom, int(tx_left), int(tx_right), multiple_iterators=True)
    read_dict = {"read_id":[], "length":[], "left_pos":[], "right_pos":[], "transcript_id":[]}
    for read in read_iter:
        read_dict["read_id"].append(read.query_name)
        read_dict["left_pos"].append(read.reference_start + 1) #pysam fetches coordinates 0-based, +1 to make 1-based
        read_dict["right_pos"].append(read.reference_end) #this gives nt 1 past end, so -1 to get end and +1 to make 1-based
        read_dict["length"].append(read.query_length)
        read_dict["transcript_id"].append(transcript_id)
    read_df = pd.DataFrame.from_dict(read_dict) 
    #read_df = all_read_df
    
    #creates two pd.Series to translate genome coordinates into spliced/unspliced transcript coords.
    unspliced_tx_len = tx_right - tx_left + 1
    exonic_positions = [] #pd.Series to translate genome coords to spliced transcript coords
    
    for exon in exons.iterrows():
        exonic_positions.extend(np.arange(exon[1].start, exon[1].end+1))
    exonic_positions = pd.Series(np.arange(1,len(exonic_positions)+1), index=exonic_positions) 
    cds_left = exonic_positions[cds_left_genomecoords]
    cds_right = exonic_positions[cds_right_genomecoords]
    global x
    x = (cds_left, cds_right, tx_len)
    if strand == "+":
        start_aligned = np.arange(-(cds_left-1), tx_len-cds_left+1)
        stop_aligned = np.arange(-(cds_right-1), tx_len-cds_right+1)
        transcript_positions = pd.Series(np.arange(1, unspliced_tx_len+1), index=np.arange(tx_left, tx_right+1))
    elif strand == "-":
        start_aligned = np.arange(cds_right-tx_len, cds_right)[::-1]
        stop_aligned = np.arange(cds_left-tx_len, cds_left)[::-1]
        #start_aligned = np.arange(-(cds_right-1), tx_len-cds_right+1)[::-1]
        #stop_aligned = np.arange(-(cds_left-1), tx_len-cds_left+1)[::-1]
        transcript_positions = pd.Series(np.arange(1, unspliced_tx_len+1)[::-1], index=np.arange(tx_left, tx_right+1))
    
    exonic_positions_start = pd.Series(start_aligned, index=exonic_positions.index)
    exonic_positions_stop = pd.Series(stop_aligned, index=exonic_positions.index)
    
    #only select reads contained entirely within annotated transcript and not contained entirely in exons
    unspliced_reads = read_df.loc[
        ((read_df.left_pos.isin(transcript_positions.index)) & 
         (read_df.right_pos.isin(transcript_positions.index)))] 
    tx_aligned_reads = read_df.loc[
        ((read_df.left_pos.isin(transcript_positions.index)) & 
         (read_df.right_pos.isin(transcript_positions.index)))].loc[:,["read_id", "length", "transcript_id"]]
    #all reads transformed to unspliced transcript coords
    if strand == "+":
        unspliced_reads["5p_pos_tx"] = transcript_positions[unspliced_reads.left_pos].values
        unspliced_reads["3p_pos_tx"] = transcript_positions[unspliced_reads.right_pos].values
    elif strand == "-":
        unspliced_reads["5p_pos_tx"] = transcript_positions[unspliced_reads.right_pos].values
        unspliced_reads["3p_pos_tx"] = transcript_positions[unspliced_reads.left_pos].values
    
    
    #only select reads that are contained entirely within introns
    exonic_reads = read_df.loc[
        (read_df.left_pos.isin(exonic_positions.index)) & 
        (read_df.right_pos.isin(exonic_positions.index))]
    cds_aligned_reads = read_df.loc[
        (read_df.left_pos.isin(exonic_positions.index)) & 
        (read_df.right_pos.isin(exonic_positions.index))].loc[:,["read_id"]]
    #change from genome coords to start aligned transcript coords
    if strand == "+":
        cds_aligned_reads["5p_pos_start_cds"] = exonic_positions_start[exonic_reads.left_pos].values
        cds_aligned_reads["3p_pos_start_cds"] = exonic_positions_start[exonic_reads.right_pos].values
        cds_aligned_reads["5p_pos_stop_cds"] = exonic_positions_stop[exonic_reads.left_pos].values
        cds_aligned_reads["3p_pos_stop_cds"] = exonic_positions_stop[exonic_reads.right_pos].values
    elif strand == "-":
        cds_aligned_reads["5p_pos_start_cds"] = exonic_positions_start[exonic_reads.right_pos].values
        cds_aligned_reads["3p_pos_start_cds"] = exonic_positions_start[exonic_reads.left_pos].values
        cds_aligned_reads["5p_pos_stop_cds"] = exonic_positions_stop[exonic_reads.right_pos].values
        cds_aligned_reads["3p_pos_stop_cds"] = exonic_positions_stop[exonic_reads.left_pos].values
    cds_aligned_reads["exonic"] = [True] * len(exonic_reads)
    
    #create dataframe containing all reads for transcript.
    read_df = pd.merge(cds_aligned_reads, unspliced_reads, on=["read_id"], suffixes=["_cds", "_tx"], how="outer")
    #read_df.transcript_id_tx.fillna(read_df.transcript_id_unspliced, inplace=True)

    
    #makes a list (called fivep or threep) containing all the read ends for transcript.
    # This has to take strandedness into account to determine 3' or 5' end.
    # numpy.histogram is used to turn list of ends into mapped-end depth at each position.
    tx_obj = pd.DataFrame()
    tx_obj["transcript_id"] = [transcript_id]
    tx_obj["cds_left_spliced"] = [cds_left]
    tx_obj["cds_right_spliced"] = [cds_right]
    tx_obj["cds_left_unspliced"] = [cds_left_genomecoords - tx_left + 1]
    tx_obj["cds_right_unspliced"] = [cds_right_genomecoords - tx_left + 1]
    tx_obj["strand"] = strand
    tx_obj["spliced_len"] = tx_len
    
    
    #create and return a tx_density object that contains all the read ends
    
    return((read_df, tx_obj))

def pickle_transcript_objects():
    global prefix
    global all_read_df
    global all_properties_df
    global codon_positions
    out_dict = {}
    pickle.dump({"tx_reads":all_read_df, "tx_properties":all_properties_df}, 
            open("./"+prefix+"_reads.pkl", "wb"))
    print("Pickled", prefix, "sample to ./"+prefix+"_reads.pkl")


# In[25]:


#%timeit -n1 -r1
file = False
if not file:
    p = Pool(96)
    map_out = p.map(get_density, gtf_df[gtf_df.feature == "transcript"].sample(frac=downsample_frac, axis=0).iterrows())
    #map_out = p.map(get_density, gtf_df[gtf_df.feature == "transcript"].iterrows())
    p.close()

    all_read_df = pd.concat([entry[0] for entry in map_out], axis=0, sort=False).reset_index(drop=True)
    all_properties_df = pd.concat([entry[1] for entry in map_out], axis=0, sort=False).reset_index(drop=True)
    #group = all_read_df.groupby(["transcript_id", "length"])
    #len_read_count = group.exonic.transform("count")
    group = all_read_df.groupby("transcript_id")
    read_count = group.exonic.transform("count")
    all_read_df["len_read_percent"] = 1/read_count
    
    pickle_transcript_objects()
else:
    in_df = pickle.load(open(file, "rb"))
    all_read_df = in_df["tx_reads"]
    all_properties_df = in_df["tx_properties"]
    print("loaded pickle", file)


# In[26]:


'''x = all_read_df.groupby(by="transcript_id_unspliced", as_index=False)


p = Pool(96)
new_list = p.starmap(myapply, zip(x[1], itertools.repeat(50), itertools.repeat(50), itertools.repeat("start")))
p.close
new_df = pd.concat(new_list, axis=0)'''

lengths = all_read_df.length.drop_duplicates()
nt_upstream = 50
nt_downstream = 50
mapped_end = "3p"
start_or_stop = "start"
align_to_site = "p_site"

column_name = mapped_end+"_pos_"+start_or_stop+"_cds"
plot_range = np.arange(-nt_upstream, nt_downstream + 1)

#normalize to lengths 
#rethink how this can be run just once
def get_lens(start_or_stop, up_or_downstream):
    global all_properties_df
    if start_or_stop == "start":
        plus_std = "cds_left_spliced"
        minus_std = "cds_right_spliced"
    elif start_or_stop == "stop":
        plus_std = "cds_right_spliced"
        minus_std = "cds_left_spliced"
    if up_or_downstream == "up":
        multiplier = -1
    elif up_or_downstream == "down":
        multiplier = 1
    plus = (all_properties_df[all_properties_df.strand == "+"][plus_std] - 1) * multiplier
    minus = (all_properties_df[all_properties_df.strand == "-"].spliced_len - all_properties_df[all_properties_df.strand == "-"][minus_std]) * multiplier
    lengths = pd.concat([plus, minus], axis=0)
    return(lengths)

upstream_lengths = get_lens(start_or_stop, "up")
downstream_lengths = get_lens(start_or_stop, "down")
upstream_lengths = np.cumsum(np.histogram(upstream_lengths, np.arange(min(upstream_lengths), max(upstream_lengths)+1))[0])
downstream_lengths = np.cumsum(np.histogram(downstream_lengths, np.arange(min(downstream_lengths), max(downstream_lengths)+1))[0][::-1])[::-1]
gene_lengths = np.concatenate([upstream_lengths[-nt_upstream:], downstream_lengths[1:nt_downstream+1]])

#histogramming
plotdf = all_read_df[all_read_df.length.isin(lengths)].dropna()
group_histo = all_read_df.dropna().groupby("length").apply(lambda x: np.histogram(x[column_name], plot_range, weights=x.len_read_percent)[0]/gene_lengths)
group_histo = pd.DataFrame(group_histo.to_dict())
group_histo.index = plot_range[:-1]
#histo = np.histogram(plotdf["5p_pos_start_cds"], np.arange(-50, 501), weights=plotdf.len_read_percent)[0]

#riboshifting
mode = group_histo.sum(axis=1).idxmax()
print("Highest peak is", mode, "from annotated start codon")
shift_dict = { #+1s to account for digestion heterogeneity
    "p_site":-mode + 1,
    "a_site":-mode + 3 + 1,
    "e_site":-mode - 3 + 1
}

#plots
mapper = linear_cmap(field_name='values', palette=cc.b_linear_blue_95_50_c20, low=0, high=group_histo.max().max())
y = [yval for xval in group_histo.index for yval in group_histo.columns]
x = [xval for xval in group_histo.index + shift_dict[align_to_site] for yval in group_histo.columns]
values = [group_histo.loc[xval,yval] for xval in group_histo.index for yval in group_histo.columns]
plot_df = pd.DataFrame({"x":x, "y":y, "values":values})
p = figure(width=1400)
p.add_tools(HoverTool(tooltips=[
    ("Read length", "@y"),
    ("Position (vs. start)", "@x"),
    ("Height", "@values")
]))
p.rect(x="x", y="y", color=mapper, line_color=mapper, width=1, height=1, source=plot_df, line_width=0.5)
show(p)

p2 = figure(title="Size distribution of mapped reads")
plot_df = pd.DataFrame(all_read_df[all_read_df.exonic == True].groupby("length").read_id.count().rename("number"))
plot_df.number = plot_df.number#/np.sum(plot_df.number)
p2.line(x="length", y="number", line_width=2, source=plot_df)
p2.add_tools(HoverTool(tooltips=[
    ("Read length", "@length"),
    ("Count", "@number"),
]))
show(p2)


# In[27]:


#annotate_codons_per_read
tx_codon_pos_dict = pickle.load(open("/home/anthony/test_code/tx_codon_pos_dict.pkl", "rb"))
all_read_df["5p_pos_start_cds_shifted"] = all_read_df["5p_pos_start_cds"] + shift_dict["p_site"]
all_read_df = pd.merge(all_read_df, tx_codon_pos_dict, how="left", 
                       left_on=["transcript_id", "5p_pos_start_cds_shifted"], 
                       right_on=["tx_id", "cds_pos"]
                      ).drop("tx_id", axis=1)
del tx_codon_pos_dict


# In[29]:


bases = "TCAG"
codons = [a + b + c for a in bases for b in bases for c in bases]
amino_acids = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
codon_table = dict(zip(codons, amino_acids))



'''
plot_list = []
three_ticker = FixedTicker(ticks=[1,2,3])
countsmax = codon_counts.max().max()
for codon in codon_counts.iterrows():
    title = codon[0]+" - "+codon_table[codon[0]]
    p = figure(title=title, y_range=(0, countsmax), tools=[])
    p.xaxis.ticker = three_ticker
    p.line(x=codon[1].index, y=codon[1].values, line_width=2)
    plot_list.append(p)
pl = plot_list


grid = gridplot([[pl[i], pl[i+1], pl[i+2], pl[i+3], pl[i+4], pl[i+5], pl[i+6], pl[i+7]] for i in range(0, len(plot_list), 8)], plot_width=170, plot_height=170)

show(grid)'''
codon_counts = all_read_df[all_read_df.length.isin([18,21,28])].groupby(["codon_nts", "codon_pos"])["5p_pos_start_cds_shifted"].count().unstack()
mapper = LinearColorMapper(palette=cc.b_linear_blue_95_50_c20, low=0, high=80000)#high=codon_counts.max().max())
mapper.high_color = "firebrick"

plotlist = []
for length in [[18],[21],[28]]:
    
    codon_counts = all_read_df[all_read_df.length.isin(length)].groupby(["codon_nts", "codon_pos"])["5p_pos_start_cds_shifted"].count().unstack()
    
    y = [yval for xval in codon_counts.index for yval in codon_counts.columns]
    x = [xval for xval in codon_counts.index for yval in codon_counts.columns]
    values = [codon_counts.loc[xval,yval] for xval in codon_counts.index for yval in codon_counts.columns]
    plot_df = pd.DataFrame({"x":x, "y":y, "values":values, "aa":[codon_table[nts] for nts in x]})
    p = figure(x_range=list(codon_counts.index.drop_duplicates()), width=1400, height=200, title="Read Length: "+str(length), y_axis_label="Position in Codon")
    p.xaxis.major_label_orientation = 1
    p.yaxis.ticker = FixedTicker(ticks=[1,2,3])
    p.add_tools(HoverTool(tooltips=[
        ("Codon", "@x"),
        ("Amino Acid", "@aa"),
        ("Position", "@y"),
        ("Height", "@values")
    ]))
    p.rect(x="x", y="y", color={'field': 'values', 'transform': mapper}, line_color={'field': 'values', 'transform': mapper}, width=1, height=1, source=plot_df, line_width=0.5)
    plotlist.append([p])
    color_bar = ColorBar(color_mapper=mapper, label_standoff=12, border_line_color=None, location=(0,0))
    p.add_layout(color_bar, 'right')
show(gridplot(plotlist))
#save(gridplot(plotlist), filename="./4E_Lys.html")


# In[30]:


all_read_df[
    (all_read_df.length == 21) &
    (all_read_df.exonic == True)
].groupby("transcript_id").read_id.count().sort_values(ascending=False)


# In[31]:


#IGV View
transcript = "ENST00000229239"
if transcript in gtf_df.transcript_id.values:
    tx_df = gtf_df[gtf_df.transcript_id == transcript]
    tx = tx_df[tx_df.feature == "transcript"]
    exons = tx_df[tx_df.feature == "exon"]
    CDS = tx_df[tx_df.feature == "CDS"]
    strand = tx.strand.values
    if strand == "+":
        start_codon = CDS.iloc[0].start
        stop_codon = CDS.iloc[-1].end
    elif strand == "-":
        start_codon = CDS.iloc[-1].end
        stop_codon = CDS.iloc[0].start
    

    hover = HoverTool(tooltips=[
        ("Reads", "@top")
    ], mode="vline")

    tx_reads = all_read_df[all_read_df.transcript_id == transcript].right_pos
    if len(tx_reads) > 0:
        histo = np.histogram(tx_reads, bins=np.arange(tx.start, tx.end+1))
        histo_df = pd.DataFrame({"top":histo[0], "left":histo[1][:-1], "right":histo[1][1:], "bottom":[0.1]*len(histo[0])})
        histo_df = histo_df[histo_df.top != 0]
        range_upper = histo_df.top.max()*10
    else:
        histo_df = pd.DataFrame()
        range_upper = 1

    read_p = figure(height=300, width=1500, y_axis_type="log", y_range = (0.1, range_upper), title=transcript, output_backend="webgl", tools=[hover, "xbox_zoom"])
    read_p.quad(left="left", right="right", top="top", bottom="bottom", color="gray", source=histo_df)
    read_p.xaxis.ticker = []
    read_p.xgrid.grid_line_color = None

    p = figure(height=100, width=1500, y_range=(-3,3), tools=["xbox_zoom", "reset"], output_backend="webgl")
    p.x_range = read_p.x_range
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.yaxis.ticker = []
    p.xaxis[0].formatter = NumeralTickFormatter(format="0,0")
    p.line(x=list(map(int, [tx.start, tx.end+1])), y=[0, 0], color="purple", line_width=5, level="glyph")
    p.quad(left=exons.start, right=exons.end+1, top=1, bottom=-1, color="white", line_color="purple", line_width=2, level="annotation")
    p.quad(left=CDS.start, right=CDS.end+1, top=1, bottom=-1, color="purple", line_color=None, level="annotation")
    if strand == "+":
        p.inverted_triangle(x=start_codon, y=2, color="green", size=10)
        p.inverted_triangle(x=stop_codon+1, y=2, color="red", size=10)
    elif strand == "-":
        p.inverted_triangle(x=start_codon+1, y=2, color="green", size=10)
        p.inverted_triangle(x=stop_codon, y=2, color="red", size=10)
        
    #p.quad(left=stop_codon, right=stop_codon+1, top=1, bottom=-1, color="red", level="annotation")


    grid = gridplot([[read_p],[p]])
    show(grid)
else:
    print(transcript, "not in gtf_df")


# # Old Code

# In[45]:



    
        
        
        
def count_density_2(transcript):
    '''This function is run by multiprocessing.Pool to count read densities along individual 
    genes according the the parameters set in this cell.
    Input is the transcript line from gtf_df
    This function returns a tuple of upstream density, downstream density, and the lengths of these regions (as lists).
    The output is normalized so that all read densities for a particular transcript add up to one (one transcript, one vote).
    '''
    global prefix
    bamfile = pysam.AlignmentFile("/home/anthony/cap_binding_pulldown/profiling/results/genome_mapped_reads/"+prefix+"Aligned.sortedByCoord.out.bam", "rb")
    global nt_upstream
    global nt_downstream
    global end_aligned
    global align_to
    global read_sizes
    global gtf_df
    global codon_positions
    strand = transcript[1].strand
    cds_codon_positions = codon_positions[transcript[1].transcript_id] #loads a dict of codon positions in spliced CDS coordinates
    tx_cds_exons_df = gtf_df[gtf_df.transcript_id == transcript[1].transcript_id] #forms a dataframe containing transcript, CDS, 
                                                                                  # and exon entries from gtf_df

    tx = tx_cds_exons_df[tx_cds_exons_df.feature == "transcript"] #transcript line of dataframe
    tx_left = tx.start.iloc[0] #Leftmost genome coordinate of transcript (could be 3' or 5' end depending on strand)
    tx_right = tx.end.iloc[0] #Rightmost genome coordinate of transcript

    exons = tx_cds_exons_df[tx_cds_exons_df.feature == "exon"] #Exon lines of dataframe of 
    tx_len = sum([(exon[1].end + 1) - exon[1].start for exon in exons.iterrows()]) #length of spliced transcript feature
    
    cds = tx_cds_exons_df[tx_cds_exons_df.feature == "CDS"] #CDS lines of dataframe
    cds_left_genomecoords = cds.iloc[0].start
    cds_right_genomecoords = cds.iloc[-1].end
    cds_len = sum([(cds[1].end + 1) - cds[1].start for cds in cds.iterrows()]) #spliced CDS length
    cds_started = False #variables to help in start and end detection of CDS during transcript splicing below:
    cds_ended = False
    
    histo = [] #initializing the histogram of mapped read ends along spliced transcript
    read_ends_dict = {} #dictionary containing ends of individual reads. Keys will be read IDs from BAM file to 
                        # ensure only single counting of reads. Value is list of [start, end] of mapped read in transcript coordinates
    offset = 0 #offset to convert genome coordinates into spliced transcript coordinates.
    
    #this for loop converts mapping read and CDS end coordinates from genome to transcript coordinates.
    for exon in exons.iterrows(): 
        offset += exon[1].start - 1 #the [1] is necessary beccause df.iterrows() returns a tuple of (index, line_of_df)
        #calculating CDS start and end in transcript coordinates 
        if exon[1].end >= cds.iloc[0].start and cds_started == False:
            #cds_left = exon[1].end - offset #test for aligning to ends of first exon
            cds_left = cds_left_genomecoords - offset
            cds_started = True
        if exon[1].end >= cds.iloc[-1].end and cds_ended == False:
            cds_right = cds_right_genomecoords - offset
            cds_len_calc = cds_right - cds_left + 1
            cds_ended = True
        read_iter = bamfile.fetch(exon[1].chrom, int(exon[1].start), int(exon[1].end), multiple_iterators=True)
        #finding exon-mapping reads:
        for read in read_iter:
            if "all" in read_sizes:
                read_ends_dict.setdefault(read.query_name, [read.reference_start-offset, read.reference_end-offset])
            elif read.reference_length in read_sizes:
                read_ends_dict.setdefault(read.query_name, [read.reference_start-offset, read.reference_end-offset])
        offset -= exon[1].end
    
    '''
    #sanity check to make sure cds_length calculated by transcript coordinates matches the 
    # value calculated by splicing CDS features together.
    try:
        assert cds_len_calc == cds_len
    except AssertionError:
        print(cds_len, cds_len_calc, offset, tx.iloc[0].transcript_id, exons)
        raise'''
    
    #makes a list (called fivep or threep) containing all the read ends for transcript.
    # This has to take strandedness into account to determine 3' or 5' end.
    # numpy.histogram is used to turn list of ends into mapped-end depth at each position.
    if end_aligned == 5:
        if strand == "+":
            fivep = [read_ends_dict[key][0] for key in read_ends_dict.keys()]
        elif strand == "-":
            fivep = [read_ends_dict[key][1] for key in read_ends_dict.keys()]
        histo = np.histogram(fivep, np.arange(0, tx_len+1, 1))[0]
    elif end_aligned == 3:
        if strand == "+":
            threep = [read_ends_dict[key][1] for key in read_ends_dict.keys()]
        elif strand == "-":
            threep = [read_ends_dict[key][0] for key in read_ends_dict.keys()]
        histo = np.histogram(threep, np.arange(0, tx_len+1, 1))[0]
    else:
        raise Exception("end_aligned variable should be 3 or 5...")

    #normalization (one transcript one vote):
    histosum = np.sum(histo)
    if histosum > 0:
        histo = list(map(lambda x: x/histosum, histo)) #divide each entry in histo by histosum
    
    #separates normalized counts upstream and downstream of cds_ends. 
    # Upstream will always be a reversed list to help downstream normalization to number of 
    # transcripts per position. Some transcripts might be shorter than given x_range of metaplot
    if align_to == "start":
        if strand == "+":
            histo_upstream = histo[:cds_left][::-1] #reversed to facilitate adding
            histo_downstream = histo[cds_left:]
        elif strand == "-":
            histo_upstream = histo[cds_right:] #reversed to facilitate adding
            histo_downstream = histo[:cds_right][::-1]
    elif align_to == "stop":
        if strand == "+":
            histo_upstream = histo[:cds_right][::-1] #reversed to facilitate adding
            histo_downstream = histo[cds_right:]
        elif strand == "-":
            histo_upstream = histo_upstream = histo[cds_left:] #reversed to facilitate adding
            histo_downstream = histo_downstream = histo[:cds_left][::-1]
    elif align_to in all_codons:
        raise Exception("you didn't code this part yet Anthony...")
    else:
        raise Exception("align_to variable must be 'start' or 'stop' or list of codons")

    return(histo_upstream[:nt_upstream], histo_downstream[:nt_downstream], [len(histo_upstream)], [len(histo_downstream)])

def sum_density(in_list, out_list, in_lock, out_lock):
    '''Adds up results of density calculation upstream. 
    This definition will be run by a worker process and manager processes handle
    in_list and out_list.
    '''
    mysum_upstream = [] #worker process's own sum of densities upstream of aligned feature
    mysum_downstream = [] #worker process's own sum of densities downstream of aligned feature
    lens_upstream = [] #list of lengths of worker process's added upstream transcript densities 
    lens_downstream = [] #list of lengths of worker process's added downstream transcript densities 
    newlist_upstream = [] #incoming individual transcript upstream density to be added to mysum_upstream
    newlist_downstream = [] #incoming individual transcript downstream density to be added to mysum_downstream
    while True:
        in_lock.acquire() #lock list so worker process is the only one that can add. Will wait for lock to become available
        if len(in_list) > 0:
            tup = in_list.pop()
            newlist_upstream = tup[0]
            newlist_downstream = tup[1]
            lens_upstream.extend(tup[2]) #keeps track of the lengths of all incoming newlists
            lens_downstream.extend(tup[3]) #keeps track of the lengths of all incoming newlists
        else: #break out of while loop if in_list is empty.
            in_lock.release()
            break
        in_lock.release()
        #gene_lens.append(len(newlist))

        if len(newlist_upstream) <= len(mysum_upstream): #adds newlist to mysum
            for n, val in enumerate(newlist_upstream):
                mysum_upstream[n] += newlist_upstream[n]
        else: #if newlist is longer than mysum, simply extend the list to include the extra positions.
            for n, val in enumerate(newlist_upstream[:len(mysum_upstream)]):
                mysum_upstream[n] += newlist_upstream[n]
            mysum_upstream.extend(newlist_upstream[len(mysum_upstream):])

        if len(newlist_downstream) <= len(mysum_downstream): #adds newlist to mysum
            for n, val in enumerate(newlist_downstream):
                mysum_downstream[n] += newlist_downstream[n]
        else:
            for n, val in enumerate(newlist_downstream[:len(mysum_downstream)]):
                mysum_downstream[n] += newlist_downstream[n]
            mysum_downstream.extend(newlist_downstream[len(mysum_downstream):])

    out_lock.acquire()
    out_list.append((mysum_upstream, mysum_downstream, lens_upstream, lens_downstream)) #gene_lens))
    out_lock.release()


# In[28]:


'''prefixes = ["CBP_elute", "4E_Lys", "4E_elute"] #prefix of input bamfile
align_to_iter = [["start", 50, 300], ["stop", 300, 50]] #list of lists. First position accepts "start" or "stop". 
                                                        # second position is upstream nucleotides for plot
                                                        # third position is downstream nucleotides for plot

#Accepts a list of lists containing integers or "all". Includes these read sizes in metagene:
# e.g. to calculate 21 and 28mer combination and each alone: [[21, 28], [21], [28]]
read_sizes_iter = [["all"], [20], [21], [22], [23], [24], [25], [26], [27], [28], [29], [30], [31], [32]]
end_aligned_iter = [5, 3] #Accepts a list integers 5 or 3. Maps the 5' or 3' end of the read.'''

prefixes = ["CBP_elute"] #prefix of input bamfile
align_to_iter = [["start", 50, 300], ["stop", 50, 300]] #list of lists. First position accepts "start" or "stop". 
                                                        # second position is upstream nucleotides for plot
                                                        # third position is downstream nucleotides for plot

#Accepts a list of lists containing integers or "all". Includes these read sizes in metagene:
# e.g. to calculate 21 and 28mer combination and each alone: [[21, 28], [21], [28]]
read_sizes_iter = [[21, 28], ["all"]]
end_aligned_iter = [5] #Accepts a list integers 5 or 3. Maps the 5' or 3' end of the read.


#put code inside this to generate multiple plots for interactivity. Uncomment code at bottom to save a pkl of densities.
for prefix in prefixes:
    for align_to, nt_upstream, nt_downstream in align_to_iter:
        for read_sizes in read_sizes_iter:
            for end_aligned in end_aligned_iter:
                print("\nSTARTING:", prefix, align_to, read_sizes, end_aligned)
                
                if ":".join([prefix, align_to, ", ".join([str(n) for n in read_sizes]), str(end_aligned)]) in output_dict.keys():
                    print("ALREADY DONE, SKIPPING: "+":".join([prefix, align_to, ", ".join([str(n) for n in read_sizes]), str(end_aligned)]))
                    continue
                

                print("\tStarted count_density...")
                p = Pool(96) #establishes a pool of 96 workers (should be <= to number of processor threads)
                try:
                    #creates list of results of count_density_2 for all transcripts, using "transcript" feature lines of gtf_df as input.
                    map_out = p.map(count_density_2,  gtf_df[gtf_df.feature == "transcript"].sample(frac=downsample_frac, axis=0).iterrows())
                except:
                    p.terminate()
                    raise
                p.close()

                print("\tStarted adding density...")
                adding_threads = 48 #number of threads used to add results of above density.


                #set up managers for in_list and out_list
                in_manager = mp.Manager()
                in_list = in_manager.list()
                in_lock = in_manager.Lock()

                out_manager = mp.Manager()
                out_list = out_manager.list()
                out_lock = out_manager.Lock()

                #add all outputs of count_density_2 to in_list
                in_lock.acquire()
                in_list.extend(map_out)
                in_lock.release()

                #start worker processes running sum_density.
                ps = [mp.Process(target=sum_density, args=(in_list, out_list, in_lock, out_lock)) for i in range(adding_threads)]
                [i.start() for i in ps]
                [i.join() for i in ps] #wait for all workers to finish.


                print("\tFinal sum...")
                in_lock.acquire()
                out_lock.acquire()

                assert(len(in_list) == 0) #make sure in_list is empty before reusing it.
                in_list.extend(out_list) #add the results of previous multithreaded adding to in_list for final sum.
                out_list = out_manager.list()

                in_lock.release()
                out_lock.release()

                #run sum_density one more time.
                p = mp.Process(target=sum_density, args=(in_list, out_list, in_lock, out_lock))
                p.start()
                p.join()

                print("\tNormalizing values to gene length and setting feature index to zero...")

                #make pandas series of upstream/downstream densities and lengths.
                H_s_upstream = pd.Series(list(out_list[0][0]))
                H_s_downstream = pd.Series(list(out_list[0][1]))
                gene_lens_upstream = pd.Series(list(out_list[0][2])).sort_values(ascending=False)
                gene_lens_downstream = pd.Series(list(out_list[0][3])).sort_values(ascending=False)

                #initialize pandas series of transcript "coverage" at each position of plot.
                count_s_upstream = pd.Series([0] * nt_upstream)
                count_s_downstream = pd.Series([0] * nt_downstream)

                #make transcript coverage lists upstream and downstream
                for length in gene_lens_upstream:
                    if length >= nt_upstream:
                        count_s_upstream = count_s_upstream + 1
                    else:
                        count_s_upstream[:length] = count_s_upstream[:length] + 1

                for length in gene_lens_downstream:
                    if length >= nt_downstream:
                        count_s_downstream = count_s_downstream + 1
                    else:
                        count_s_downstream[:length] = count_s_downstream[:length] + 1

                #normalize to number of transcripts covering each position of density.
                H_s_upstream = H_s_upstream / count_s_upstream
                H_s_downstream = H_s_downstream / count_s_downstream

                H_s_upstream = H_s_upstream.iloc[::-1] #reverse upstream pandas series to correct orientation.
                H_s_upstream.index = np.arange(-nt_upstream, 0)

                #concatenate upstream and downstream series.
                H_s = pd.concat([H_s_upstream, H_s_downstream])

                output_dict[":".join([prefix, align_to, ", ".join([str(n) for n in read_sizes]), str(end_aligned)])] = H_s
                print("\tDone.")

metagene_filename = "metagenes_"+strftime("%Y%m%d%H%M%S")+".pkl"
file = open(metagene_filename, "wb")
pickle.dump(output_dict, file)
file.close()
print("\nALL DONE.")
output_dict = {}


# ## Plots

# In[29]:


path_to_pkl = False #change this to open a specific pkl file, make it False to open the last processed pkl file.

if not path_to_pkl:
    pkl_to_open = open(metagene_filename, "rb")
else:
    pkl_to_open = open(path_to_pkl, "rb")
metagene_dict = pd.DataFrame.from_dict(pickle.load(pkl_to_open))

prefixes = set()
align_to = set()
read_lens = set()
end_aligned = set()
for key in metagene_dict.keys():
    options = key.split(":")
    prefixes.add(options[0])
    align_to.add(options[1])
    read_lens.add(options[2])
    end_aligned.add(options[3])

prefixes = list(prefixes)
prefixes.sort()
align_to = list(align_to)
align_to.sort()
end_aligned = list(end_aligned)
end_aligned.sort()
read_lens_ints = []
read_lens_combos = []
all_present = False
if "all" in read_lens:
    read_lens.remove("all")
    all_present = True
for entry in read_lens:
    if ", " in entry:
        read_lens_combos.append(entry)
    else:
        read_lens_ints.append(entry)
read_lens_ints.sort()
read_lens_combos.sort()
read_lens = read_lens_ints + read_lens_combos
        
        
plotdict = {}
for prefix in prefixes:
    plotdict[prefix] = figure(title = "Metagene "+prefix, y_axis_label = "Average Density per Nucleotide",
                              x_axis_label = "Position", plot_width = 1000, plot_height = 200)
linedict = {}
for plot in plotdict.keys():
    linedict[plot] = plotdict[plot].line(x=[], y=[], line_color="red", line_width=1)

cds_toggle = ToggleButtons(options=align_to, description="Feature Aligned: ")
alignment_toggle = ToggleButtons(options=[str(end)+"' Aligned" for end in end_aligned], description="End Aligned:")
if all_present:
    readlen_toggle = ToggleButtons(options=["All", "Select Below"], description="Read Length:")
else:
    readlen_toggle = ToggleButtons(options=["Select Below"], description="Read Length:")


read_len_slider = SelectionSlider(options=read_lens, selected=read_lens[0], continuous_update=True, 
                                  description="Select:", disabled=True)

@interact(cds_toggle=cds_toggle, alignment_toggle=alignment_toggle, 
          readlen_toggle=readlen_toggle, read_len=read_len_slider)
def Display(cds_toggle, alignment_toggle, readlen_toggle, read_len):
    if cds_toggle == "start":
        cds = "start"
    elif cds_toggle == "stop":
        cds = "stop"
    if alignment_toggle == "5' Aligned":
        end_map = "5"
    elif alignment_toggle == "3' Aligned":
        end_map = "3"
    if readlen_toggle == "All":
        read_len_slider.disabled = True
        rl = "all"
    elif readlen_toggle == "Select Below":
        read_len_slider.disabled = False
        rl = read_len
    colname = ":".join(["", cds, rl, end_map])
    
    if len(metagene_dict[prefix+colname]) == 0:
        print(metagene_dict[prefix+colname])
    for line in linedict.keys():
        data = {}
        data["x"] = list(metagene_dict[prefix+colname].index)
        data["y"] = list(metagene_dict[prefix+colname])
        linedict[line].data_source.data = data
        
    push_notebook()
    
    

show(column([plotdict[plot] for plot in plotdict.keys()]), notebook_handle=True)


# # Old Code

# In[34]:


names_dict = {
    "TTAGGC":"CBP Input",
    "GATCAG":"CBP IP",
    "AGTCAA":"4E Input",
    "CCGTCC":"4E IP"
}

#pickle.dump(H_df, open("./"+index+".pkl", "wb"))

p1 = figure(title = "Read 3' End Density Metagene - ", y_axis_label = "Average Density per Nucleotide",
            x_axis_label = "Position (Nucleotides from 5' End)", plot_width = 1000)


my_df = H_df

p1.line(y=my_df, x=my_df.index+1, 
        line_color="red",
        line_width=2)
show(p1)


# In[32]:


names_dict = {
    "TTAGGC":"CBP Input",
    "GATCAG":"CBP IP",
    "AGTCAA":"4E Input",
    "CCGTCC":"4E IP"
}

#pickle.dump(H_df, open("./"+index+".pkl", "wb"))

p1 = figure(title = "Read 3' End Density Metagene - ", y_axis_label = "Average Density per Nucleotide",
            x_axis_label = "Position (Nucleotides from 5' End)", plot_width = 1000)
p1.yaxis.axis_label_text_font_size="15pt"
p1.xaxis.axis_label_text_font_size="15pt"
p1.title.text_font_size="15pt"
p1.yaxis.major_label_text_font_size="12pt"
p1.xaxis.major_label_text_font_size="12pt"
p1.output_backend="svg"

my_df = H_df

p1.line(y=my_df, x=my_df.index+1, 
        line_color="red",
        line_width=2)
show(p1)


# In[7]:


#load pickled H_df



plot_index = "all" #'all' to show all plots

names_dict = {
    "TTAGGC":"CBP Input",
    "GATCAG":"CBP IP",
    "AGTCAA":"4E Input",
    "CCGTCC":"4E IP"
}

fig_dict = {}
H_fig_dict = {}
for idx in names_dict.keys():
    H_df_pkl = pickle.load(open("./"+idx+".pkl", "rb"))

    fig_dict[idx] = figure(title = "Read 3' End Density Metagene - "+names_dict[idx], y_axis_label = "Average Density per Nucleotide",
                x_axis_label = "Position (Nucleotides from 5' End)")
    fig_dict[idx].yaxis.axis_label_text_font_size="15pt"
    fig_dict[idx].xaxis.axis_label_text_font_size="15pt"
    fig_dict[idx].title.text_font_size="15pt"
    fig_dict[idx].yaxis.major_label_text_font_size="12pt"
    fig_dict[idx].xaxis.major_label_text_font_size="12pt"
    fig_dict[idx].output_backend="webgl"

    H_fig_dict[idx] = H_df_pkl
    fig_dict[idx].quad(top=H_fig_dict[idx], bottom=0, left=H_fig_dict[idx].index, right=H_fig_dict[idx].index+1,
            fill_color="darkgreen", 
            line_color="darkgreen",
            line_width=.5)

if plot_index == "all":
    for show_idx in fig_dict.keys():
        show(fig_dict[show_idx])
else:
    show(fig_dict[idx])


# ## Old code

# In[2]:


inq = Queue()
outq = Queue()


def sum_density(inq, outq):
    mysum = []
    while inq.empty() is False:
        newlist = inq.get()
        if len(newlist) <= len(mysum):
            for n, val in enumerate(newlist):
                mysum[n] += newlist[n]
        else:
            for n, val in enumerate(newlist[:len(mysum)]):
                mysum[n] += newlist[n]
            mysum.extend(newlist[len(mysum):])
    outq.put(mysum)
            
  
x = [[1,2], [1,2,3], [1,2,3], [1,2,3,4], [1,2], [1,2,3], [1,2,3], [1,2,1,1], [1,2,3], [1,2,3], [1,2,3,4]] * 10000 
[inq.put(i) for i in x]

ps = [mp.Process(target=sum_density, args=(inq, outq)) for i in range(96)]
[i.start() for i in ps]
[i.join() for i in ps]

while outq.empty() is False:
    inq.put(outq.get())

p = mp.Process(target=sum_density, args=(inq, outq))
p.start()
p.join()

while outq.empty() is False:
    print(outq.get())


# In[95]:


outq.qsize()


# In[39]:


while not q.empty():
    q.get()


# In[47]:


q2.empty()


# In[11]:


len([])

