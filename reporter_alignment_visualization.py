
# coding: utf-8

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
from bokeh.models import ColumnDataSource, CDSView, HoverTool, LinearColorMapper, FixedTicker, ColorBar, NumeralTickFormatter, ContinuousTicker
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

prefix = "mock_reporter"
graphs = []


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

# In[2]:


gtf_df = pd.DataFrame({
    "feature":["transcript", "CDS", "exon", "exon", "exon"],
    "chrom":["K20_reporter", "K20_reporter", "K20_reporter", "K20_reporter", "K20_reporter"],
    "start":[1, 1240, 1, 1240, 1300],
    "end":[2127, 1299, 1239, 1299, 2127],
    "strand":["+", "+", "+", "+", "+"],
    "transcript_id":["K20_reporter", "K20_reporter", "K20_reporter", "K20_reporter", "K20_reporter"],
})


# ## Multithreaded density addition

# ### Define density calculation and summation functions
# ```count_density_2``` is the workhorse of the script that calculates normalized densities for each transcript  
# ```sum_density sums``` up all the results of ```count_density_2```

# In[3]:


pd.options.mode.chained_assignment = None  # default='warn'
def get_density(transcript):
    '''takes in a transcript line of gtf_df and creates a tx_density object. These will then be saved to a pkl file for fast parsing of downstream calculations.
    '''
    global prefix
    global gtf_df
    bamfile = pysam.AlignmentFile("/home/anthony/siRNA_profiling/output/reporter_aligned/Mock/Mock_Aligned.out.plusstrand.sorted.bam")
    
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


# In[200]:


#%timeit -n1 -r1
file = False #"./genome_aligned_reads.pkl"
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


# In[7]:


'''x = all_read_df.groupby(by="transcript_id_unspliced", as_index=False)


p = Pool(96)
new_list = p.starmap(myapply, zip(x[1], itertools.repeat(50), itertools.repeat(50), itertools.repeat("start")))
p.close
new_df = pd.concat(new_list, axis=0)'''

lengths = all_read_df.length.drop_duplicates()
nt_upstream = 15
nt_downstream = 70
mapped_end = "5p"
start_or_stop = "start"
align_to_site = "p_site"

column_name = mapped_end+"_pos_"+start_or_stop+"_cds"
plot_range = np.arange(-nt_upstream, nt_downstream + 1)


gene_lengths = [1]*(nt_upstream+nt_downstream)

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


# In[4]:


## IGV View

file = "./mock_reporter_reads.pkl"
in_df = pickle.load(open(file, "rb"))
all_read_df = in_df["tx_reads"]
all_properties_df = in_df["tx_properties"]
print("loaded pickle", file)

minimum_read_length = 20

transcript = "K20_reporter"
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
        ("Reads", "@top"),
        ("Position", "@right")
    ], mode="vline", attachment="above", anchor="top_center")
    
    hover2 = HoverTool(tooltips=[
        ("Position", "@x")
    ], name="triangle")

    tx_reads = all_read_df[(all_read_df.length > 20) & (all_read_df.transcript_id == transcript)][["right_pos", "len_read_percent"]]
    if len(tx_reads) > 0:
        histo = np.histogram(tx_reads.right_pos, bins=np.arange(tx.start, tx.end+1), weights=tx_reads.len_read_percent)
        histo_df = pd.DataFrame({"top":histo[0]/max(histo[0]), "left":histo[1][:-1], "right":histo[1][1:], "bottom":[0]*len(histo[0])})
        histo_df = histo_df[histo_df.top != 0]
        range_upper = histo_df.top.max()*10
    else:
        histo_df = pd.DataFrame()
        range_upper = 1

    read_p = figure(height=300, width=1500, y_range=(0,1), title=transcript+" - "+file, output_backend="webgl", tools=[hover, "xbox_zoom"])#, y_axis_type="log", y_range = (0.1, range_upper))
    read_p.quad(left="left", right="right", top="top", bottom="bottom", color="gray", source=histo_df)
    read_p.xaxis.ticker = []
    read_p.xgrid.grid_line_color = None

    p = figure(height=100, width=1500, y_range=(-3,3), tools=["xbox_zoom", "reset", hover2], output_backend="webgl")
    p.x_range = read_p.x_range
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.yaxis.ticker = []
    p.xaxis[0].formatter = NumeralTickFormatter(format="0,0")
    p.line(x=list(map(int, [tx.start, tx.end+1])), y=[0, 0], color="purple", line_width=5, level="glyph")
    p.quad(left=exons.start, right=exons.end+1, top=1, bottom=-1, color="white", line_color="purple", line_width=2, level="annotation")
    p.quad(left=CDS.start, right=CDS.end+1, top=1, bottom=-1, color="purple", line_color=None, level="annotation")
    p.inverted_triangle(x="x", y="y", color="green", size=10, source=pd.DataFrame({"x":[start_codon], "y":[2]}))
    p.inverted_triangle(x="x", y="y",  color="red", size=10, source=pd.DataFrame({"x":[stop_codon+1], "y":[2]}))
        
    #p.quad(left=stop_codon, right=stop_codon+1, top=1, bottom=-1, color="red", level="annotation")


    grid = gridplot([[read_p],[p]])
    show(grid)
else:
    print(transcript, "not in gtf_df")
    
graphs.append(grid)


# In[5]:


## IGV View

file = "./ASCC3_reporter_reads.pkl"
in_df = pickle.load(open(file, "rb"))
all_read_df = in_df["tx_reads"]
all_properties_df = in_df["tx_properties"]
print("loaded pickle", file)

minimum_read_length = 20

transcript = "K20_reporter"
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
        ("Reads", "@top"),
        ("Position", "@right")
    ], mode="vline", attachment="above", anchor="top_center")
    
    hover2 = HoverTool(tooltips=[
        ("Position", "@x")
    ], name="triangle")

    tx_reads = all_read_df[(all_read_df.length > 20) & (all_read_df.transcript_id == transcript)][["right_pos", "len_read_percent"]]
    if len(tx_reads) > 0:
        histo = np.histogram(tx_reads.right_pos, bins=np.arange(tx.start, tx.end+1), weights=tx_reads.len_read_percent)
        histo_df = pd.DataFrame({"top":histo[0]/max(histo[0]), "left":histo[1][:-1], "right":histo[1][1:], "bottom":[0]*len(histo[0])})
        histo_df = histo_df[histo_df.top != 0]
        range_upper = histo_df.top.max()*10
    else:
        histo_df = pd.DataFrame()
        range_upper = 1

    read_p = figure(height=300, width=1500, y_range=(0,1), title=transcript+" - "+file, output_backend="webgl", tools=[hover, "xbox_zoom"])#, y_axis_type="log", y_range = (0.1, range_upper))
    read_p.quad(left="left", right="right", top="top", bottom="bottom", color="gray", source=histo_df)
    read_p.xaxis.ticker = []
    read_p.xgrid.grid_line_color = None

    p = figure(height=100, width=1500, y_range=(-3,3), tools=["xbox_zoom", "reset", hover2], output_backend="webgl")
    p.x_range = read_p.x_range
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None
    p.yaxis.ticker = []
    p.xaxis[0].formatter = NumeralTickFormatter(format="0,0")
    p.line(x=list(map(int, [tx.start, tx.end+1])), y=[0, 0], color="purple", line_width=5, level="glyph")
    p.quad(left=exons.start, right=exons.end+1, top=1, bottom=-1, color="white", line_color="purple", line_width=2, level="annotation")
    p.quad(left=CDS.start, right=CDS.end+1, top=1, bottom=-1, color="purple", line_color=None, level="annotation")
    p.inverted_triangle(x="x", y="y", color="green", size=10, source=pd.DataFrame({"x":[start_codon], "y":[2]}))
    p.inverted_triangle(x="x", y="y",  color="red", size=10, source=pd.DataFrame({"x":[stop_codon+1], "y":[2]}))
        
    #p.quad(left=stop_codon, right=stop_codon+1, top=1, bottom=-1, color="red", level="annotation")


    grid = gridplot([[read_p],[p]])
    show(grid)
else:
    print(transcript, "not in gtf_df")
    
graphs.append(grid)


# In[6]:


allgrid = gridplot([[grid] for grid in graphs])
show(allgrid)

