import numpy as np
from bokeh.plotting import figure, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import pdist, squareform
import scipy
import streamlit as st

import itertools
import pandas as pd
import matplotlib.pyplot as plt

def plot_dendrogram(uniq,model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
   
    dendrogram(linkage_matrix, **kwargs)

@st.experimental_memo
def getlinkage(uniqdf,distancemetric, optimal_ordering):
  return(linkage(uniqdf,method=distancemetric,optimal_ordering=optimal_ordering))

@st.cache()
def convert_df(fr_df):
  return fr_df.to_csv(index=False).encode('utf-8')

def view_alignment(aln, useconsensus,fr=1,fontsize="9pt", plot_width=1500,see_seq=False):
    """Bokeh sequence alignment view"""
    
    consensus_sequence=useconsensus
    
    #seqs = [[rec.seq] * int(np.ceil(int(rec.description))) for rec in (aln)]
    forfileseqs= [ rec.annotations['seq'] for rec in (aln)]
    forfileinsert = [rec.annotations['insert'] for rec in aln]
    seqs = [rec.seq for rec in (aln)]
    counts = [int(np.ceil(int(rec.description))) for rec in (aln)]
    alnscores = [np.float(rec.name) for rec in (aln)]

    
    ids = np.arange(0,len(seqs))
    ids=np.char.mod('%d', ids)
    ids=ids.tolist()

    colors=get_colors_consensus(seqs,consensus_sequence)
    
    N = len(seqs[0])
    S = len(seqs)    
    
    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()

    
    gg=np.array(colors).reshape((-1,np.max(gx)))
    ugg = gg
    countgg = np.array(counts)

    k=np.unique(gg)

    k=k[k != 'white'][::-1]

    v=np.arange(0,len(k))+1
        
    all_nnn=np.empty(shape=gg.shape)

    unique_nnn=np.empty(shape=k.shape)

    unique_nnn=np.empty(shape=ugg.shape)
    
    for key,val in zip(k,v):
      unique_nnn[ugg==key]=val
    
    
    unique_nnn[ugg=='white']=0
      
    col_score_kms=st.columns(11)
    st.write('Working on Framework '+str(fr))
    # setting distance_threshold=0 ensures we compute the full tree.
    dendo, axis =plt.subplots(1,1,figsize=(3,7))

    uniq,counts=unique_nnn,countgg

    if uniq.shape[0]>1:
      uniqdf=pd.DataFrame(uniq)


      uniqdf=np.array(uniqdf)
      countsforlink=counts.copy()
      countsforlink[np.argmax(countsforlink)]=1
      
      dfwcountnomain=np.repeat(uniqdf,countsforlink,axis=0)
      mismatchmutatedspotsuse=np.sum(dfwcountnomain != 0,axis=0)/dfwcountnomain.shape[0] > 0.01
      indelmutatedspotsuse=np.sum(dfwcountnomain == 5,axis=0)/dfwcountnomain.shape[0] > 0.001
      indexfirstletternotindel=np.where(~indelmutatedspotsuse)[0][0]
      if indexfirstletternotindel >0:
        indelmutatedspotsuse[0:indexfirstletternotindel]=False
      mutatedspotsuse=mismatchmutatedspotsuse+indelmutatedspotsuse
      with st.spinner('performing distance for '+str(uniqdf.shape[0])+' sequences with '+str(uniqdf.shape[1])+' positions and using'+str(np.sum(mutatedspotsuse))+' of those positions'):
        #uniqdftouseindendo=uniqdf[:,mutatedspotsuse]
        #uniqdftouseindendo,induniq = np.unique(uniqdftouseindendo,return_index=True)
        
        Z2 =getlinkage(uniqdf[:,mutatedspotsuse],'ward', optimal_ordering=True)
        #uniqdftouseindendo[scipy.cluster.hierarchy.leaves_list(Z2)]
        #induniq[scipy.cluster.hierarchy.leaves_list(Z2)]
      st.write('done with distance computation')
  
      countsforlink =countsforlink[scipy.cluster.hierarchy.leaves_list(Z2)] 
      counts=counts[scipy.cluster.hierarchy.leaves_list(Z2)] 
      
      seqs=np.array(seqs)[scipy.cluster.hierarchy.leaves_list(Z2)]
      forfileseqs=np.array(forfileseqs)[scipy.cluster.hierarchy.leaves_list(Z2)]
      forfileinsert=np.array(forfileinsert)[scipy.cluster.hierarchy.leaves_list(Z2)]
      strc=[str(hhh) if hhh>10 else '' for hhh in counts]
      ugg=ugg[scipy.cluster.hierarchy.leaves_list(Z2)] 

      R = dendrogram(Z2, orientation="left",leaf_font_size=6,labels=strc)

      uniqdf= np.array(uniqdf)[scipy.cluster.hierarchy.leaves_list(Z2)]
      
      countview = counts.copy()
      countview[np.argmax(counts)]=1
      countview_fordf = np.repeat(counts,countview)
      dfwcounts = np.repeat(uniqdf,countsforlink ,axis=0)
      

      uniqtocolors=np.repeat(ugg,countsforlink ,axis=0) 
      
      uniqtocolors[np.argmax(countview_fordf)]=['purple']*uniqtocolors.shape[1]


    col1,col2=st.columns([5,1])
    col1.write('Sequences')
    col1.write(seqs[np.argsort(-counts)])
    col2.write('Count of Sequences')
    col2.write(counts[np.argsort(-counts)])

    st.download_button(label='proper order Download',data=convert_df(pd.DataFrame({'sequence': seqs,'sequenceNotformatted': forfileseqs,'numberobserved':counts,'inserts':forfileinsert})),file_name='unique_data_'+str(fr)+'.csv',mime='text/csv') 

    seqs_for_view = np.repeat(seqs,countsforlink) 
    incrementforview = int(np.round(seqs_for_view.shape[0]/1000))
    if incrementforview > 1:
      indsview= np.arange(0,seqs_for_view.shape[0],incrementforview)
      indsview= np.append(indsview, np.argmax(countview_fordf))
      indsview=np.sort(np.unique(indsview))
      seqs_for_view =seqs_for_view[indsview]
      uniqtocolors = uniqtocolors[indsview]

    colors_for_view=uniqtocolors.ravel().tolist()
     
    N = len(seqs_for_view[0])
    S = len(seqs_for_view) 
    width = .4

    x = np.arange(1,N+1)
    y = np.arange(0,S,1)
    #creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    #flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
 
    #use recty for rect coords with an offset
    recty = gy+.5
    #recty2 = gy2+0.5
    #recty3 = gy3+0.5
    h= 1/S
    
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, colors=colors_for_view)) #text=text
    
    if see_seq:
      plot_height = len(seqs_for_view)*10+50
    x_range = Range1d(0,N+1, bounds='auto')
    view_range = (0,N)
    tools="xpan, xwheel_zoom, reset, save"

     #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=700,
                x_range=view_range, y_range=(0,S), tools=tools,#"xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                #text_font="monospace",
                text_font_size=fontsize)
               
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)
    p1.add_glyph(source, glyph)
    p1.add_glyph(source, rects)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    
    

    p = gridplot([[p1]], toolbar_location='below')
    cola,colb=st.columns([2,7])
      
    cola.pyplot(dendo)#, height=21)
    
    colb.bokeh_chart(p)


def make_seq(length=40):    
    return ''.join([random.choice(['A','C','T','G']) for i in range(length)])

def mutate_seq(seq):
    """mutate a sequence randomly"""
    seq = list(seq)
    pos = np.random.randint(1,len(seq),6)    
    for i in pos:
        seq[i] = random.choice(['A','C','T','G'])
    return ''.join(seq)

def get_colors_consensus(seqs,consensus_sequence):
    """make colors for bases in sequence"""
    #[i for i in range(len(s1)) if s1[i] != s2[i]]
    text = [[s[i],i] for s in list(seqs) for i in range(len(s))]
    def clrs(letter,pos):
      if consensus_sequence[pos]==letter and letter != '-':
        color='white'
      elif consensus_sequence[pos]==letter and letter =='-':
        color='gray'
      elif (letter != '-' and consensus_sequence[pos] =='-') or (letter == '-' and consensus_sequence[pos] !='-'):
        color='black'
      else:
        if letter=='A':  
          color='red' 
        elif letter=='T':
          color= 'blue'#'lightcoral'
        elif letter=='G':
          color='brown'#'crimson'
        elif letter=='C':
          color='turquoise'#firebrick'  
        elif letter=='N':
          color='#FFEFDB'   
        else:
          color='lightblue'
      return color
          
    colors = [clrs(s,i) for s,i in text]
    return colors




def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    clrs =  {'A':'red','T':'green','G':'orange','C':'blue','-':'white','N':'brown'}
    colors = [clrs[i] for i in text]
    return colors
