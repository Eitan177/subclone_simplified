from operator import sub
import numpy as np
from bokeh.plotting import figure, show
from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot

from scipy.cluster.hierarchy import dendrogram
from sklearn.datasets import load_iris
from sklearn.cluster import AgglomerativeClustering
#from scipy.cluster.hierarchy import linkage
from fastcluster import linkage
from scipy.spatial.distance import pdist, squareform
import scipy
import streamlit as st
from bokeh.palettes import Magma, Inferno, Plasma, Viridis, Cividis
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
import pdb

@st.experimental_memo
def getlinkage(uniqdf,distancemetric):#, optimal_ordering):
  return(linkage(uniqdf,method=distancemetric))#,optimal_ordering=optimal_ordering))

@st.cache()
def convert_df(fr_df):
  return fr_df.to_csv(index=False).encode('utf-8')

def make_alignmentplotwithcluster(aln_pd, useconsensus,fr=1,fontsize="9pt", plot_width=750,see_seq=False,reorder=False):
    """Bokeh sequence alignment view"""
    aln = list()
    for ind, row in aln_pd.iterrows():

        aln.append(SeqRecord(seq=row['sequence'],name=str(row['closest match']),description=str(row['numberobserved']),id=str(row['outliers']),annotations={'seq':row['sequence'],'insert':row['inserts']}))


    consensus_sequence=useconsensus
    
    #seqs = [[rec.seq] * int(np.ceil(int(rec.description))) for rec in (aln)]
    forfileseqs= [ rec.annotations['seq'] for rec in (aln)]
    forfileinsert = np.array([rec.annotations['insert'] for rec in aln])
    seqs = np.array([rec.seq for rec in (aln)])
    counts = [int(np.ceil(int(rec.description))) for rec in (aln)]
    alnscores = np.array([np.float(rec.name) for rec in (aln)])
    subclonecall = np.array([float(rec.id) for rec in aln])
    if len(set(alnscores)) <= 11 and len(set(alnscores)) >= 3:
        subclonecolors=Viridis[len(set(alnscores))]
    elif len(set(alnscores)) > 11 and len(set(alnscores))<14: 
        
        subclonecolors=Viridis[11]+Cividis[3][0:len(set(alnscores))-11]  
    elif len(set(alnscores)) >= 14 and len(set(alnscores)) <= 22:
        subclonecolors=Viridis[11]+Cividis[len(set(alnscores))-11]
    elif len(set(alnscores)) > 22 and len(set(alnscores))<= 24:  
        subclonecolors=Magma[11]+Viridis[11]+Cividis[3][0:len(set(alnscores))-22]
    elif len(set(alnscores)) >= 24 and len(set(alnscores)) <= 32:
        subclonecolors = Magma[11] + Viridis[11] + Cividis[len(set(alnscores)) - 22]
    elif len(set(alnscores)) >= 33:
        cc = Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + \
         Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11] + Cividis[10] + Magma[11] + Viridis[11]
        subclonecolors = cc[0:len(set(alnscores))]
    else:
        subclonecolors = Cividis[3][0:len(set(alnscores))]
     
    colors=get_colors_consensus(seqs,consensus_sequence)
    
    gg=np.array(colors).reshape((-1, len(seqs[0])))
    
    
    # pdb.set_trace()
    subclonecolormatch = np.repeat([subclonecolors[int(ii)] for ii in alnscores],5).reshape(-1,5)
    

    gg=np.hstack((gg,subclonecolormatch))


    for_ordering=subclonecolormatch[~pd.DataFrame(subclonecolormatch).duplicated()]

    fororder_list=[]
    
    if reorder:
        for ind, jj in pd.DataFrame(for_ordering).iterrows():
            
            fororder_list.append(np.where(np.all(subclonecolormatch==np.array(jj),axis=1))[0])
        contiguous_color_order=np.hstack(fororder_list)
        gg=gg[contiguous_color_order]
        forfileinsert=forfileinsert[contiguous_color_order]
        alnscores= alnscores[contiguous_color_order]
        counts=np.array(counts)[contiguous_color_order]
        seqs= seqs[contiguous_color_order]
        subclonecall=subclonecall[contiguous_color_order]
        subclonecolormatch=subclonecolormatch[contiguous_color_order]

    locofinsertion=[int(jj.split(',')[0].replace('[','').replace(']',''))  for jj in forfileinsert[np.where(forfileinsert != '[]')]]

    gg[np.where(forfileinsert != '[]'),locofinsertion]='teal'
    
    gg_subc=gg.copy()
    
    for ii in set(alnscores):
        
        gg_subc[alnscores==ii] =gg[np.logical_and(alnscores==ii ,subclonecall>0)]
    

    
    
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

    # setting distance_threshold=0 ensures we compute the full tree.
    dendo, axis =plt.subplots(1,1,figsize=(5,15))

    uniq,counts=unique_nnn,countgg

    if uniq.shape[0]>1:
      uniqdf=pd.DataFrame(uniq)


      uniqdf=np.array(uniqdf)
      countsforlink=counts.copy()
      countsforlink[np.argmax(countsforlink)]=np.ceil((sum(counts)-max(counts))/500)
      
      dfwcountnomain=np.repeat(uniqdf,countsforlink,axis=0)
      mismatchmutatedspotsuse=np.sum(dfwcountnomain != 0,axis=0)/dfwcountnomain.shape[0] > 0.01
      indelmutatedspotsuse=np.sum(dfwcountnomain == 5,axis=0)/dfwcountnomain.shape[0] > 0.001
      indexfirstletternotindel=np.where(~indelmutatedspotsuse)[0][0]
      if indexfirstletternotindel >0:
        indelmutatedspotsuse[0:indexfirstletternotindel]=False
      mutatedspotsuse=mismatchmutatedspotsuse+indelmutatedspotsuse
      
      # 03/26
      #with st.spinner('performing distance for '+str(uniqdf.shape[0])+' sequences with '+str(uniqdf.shape[1])+' positions and using'+str(np.sum(mutatedspotsuse))+' of those positions'):
        #uniqdftouseindendo=uniqdf[:,mutatedspotsuse]
        #uniqdftouseindendo,induniq = np.unique(uniqdftouseindendo,return_index=True)
        
        #Z2 =getlinkage(uniqdf[:,mutatedspotsuse],'ward')#, optimal_ordering=True)
        
        # 03/26
        #Z2 =getlinkage(uniqdf,'ward')
        

      #st.write('done with distance computation')
  
      #countsforlink =countsforlink[scipy.cluster.hierarchy.leaves_list(Z2)] 
      #counts=counts[scipy.cluster.hierarchy.leaves_list(Z2)] 
      
      #seqs=np.array(seqs)[scipy.cluster.hierarchy.leaves_list(Z2)]
      #forfileseqs=np.array(forfileseqs)[scipy.cluster.hierarchy.leaves_list(Z2)]
      #forfileinsert=np.array(forfileinsert)[scipy.cluster.hierarchy.leaves_list(Z2)]
      #strc=[str(hhh) if hhh>10 else '' for hhh in counts]
      #ugg=ugg[scipy.cluster.hierarchy.leaves_list(Z2)] 


       
      #subclonecall=subclonecall[scipy.cluster.hierarchy.leaves_list(Z2)]
      
      #gg_subc=gg_subc[scipy.cluster.hierarchy.leaves_list(Z2)]

      #R = dendrogram(Z2, orientation="left",leaf_font_size=6,labels=strc)

      #uniqdf= np.array(uniqdf)[scipy.cluster.hierarchy.leaves_list(Z2)]
      # 03/26        

      countview_fordf = np.repeat(counts,countsforlink)
      dfwcounts = np.repeat(uniqdf,countsforlink ,axis=0)
      
      #st.download_button(label='proper order Download',data=convert_df(pd.DataFrame({'sequence': seqs,'sequenceNotformatted': forfileseqs,'numberobserved':counts,'inserts':forfileinsert.tolist(),'subclone':subclonecall,'group':alnscores})),file_name='unique_data_'+str(fr)+'.csv',mime='text/csv') 

      uniqtocolors=np.repeat(ugg,countsforlink ,axis=0) 
      justsubc_tocolors=np.repeat(gg_subc,countsforlink ,axis=0) 
      forfileseqs = np.repeat(forfileseqs,countsforlink)
      forfileinsert=np.repeat(forfileinsert,countsforlink)
      
      alnscores=np.repeat(alnscores,countsforlink)
      uniqtocolors[np.where(countview_fordf==np.max(countview_fordf)),0:(uniqtocolors.shape[1]-5) ]=['purple']
      justsubc_tocolors[np.where(countview_fordf==np.max(countview_fordf)),0:(justsubc_tocolors.shape[1]-5) ]=['purple']


    seqs_for_view = np.repeat(seqs,countsforlink) 
    ## 03/26
    incrementforview = int(np.round(seqs_for_view.shape[0]/1000))
    ##incrementforview = 1
    ## 03/26
    
    if incrementforview > 1:
      indsview= np.arange(0,seqs_for_view.shape[0],incrementforview)
      indsview= np.append(indsview, np.argmax(countview_fordf))
      indsview=np.sort(np.unique(indsview))
      if len(np.where(justsubc_tocolors=='teal')[0]) > 0 and len(np.where(justsubc_tocolors=='teal')[0]) < 100:
          indsview=np.append(indsview,np.repeat(np.where(justsubc_tocolors=='teal')[0][0],5))
      indsview=np.sort(indsview) 
      seqs_for_view =seqs_for_view[indsview]
      uniqtocolors = uniqtocolors[indsview]
      alnscores=alnscores[indsview]
      justsubc_tocolors = justsubc_tocolors[indsview]
      
      forfileseqs = forfileseqs[indsview]
      forfileinsert=forfileinsert[indsview]
    
    st.download_button(label='sequences in alignment shown',data=convert_df(pd.DataFrame({'sequence': seqs_for_view,'sequenceNotformatted': forfileseqs,'inserts':forfileinsert.tolist(),'group':alnscores})),file_name='alignmenview_sequences_'+str(fr)+'.csv',mime='text/csv',key=np.random.rand()) 


    colors_for_view=np.flip(uniqtocolors,axis=1).ravel().tolist()
    just_subc_colors_for_view = np.flip(justsubc_tocolors,axis=1).ravel().tolist()
    
    N = len(seqs_for_view[0])
    S = len(seqs_for_view) 
    width = .4
   
    

    
    # 03/26
    #x = np.arange(1,N+1)
    x = np.arange(1,N+6)
    # 03/26

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
    
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, colors=colors_for_view[::-1])) #text=text
    subc_source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, colors=just_subc_colors_for_view[::-1])) #text=text
    
    if see_seq:
      plot_height = len(seqs_for_view)*10+50
    
    ## 03/26
    #x_range = Range1d(0,N+1, bounds='auto')
    x_range = Range1d(0,N+6, bounds='auto')

    view_range = (0,N)

    tools="xpan, xwheel_zoom, reset, save"
    
     #sequence text view with ability to scroll along x axis
    p1 = figure(title=None, plot_width=plot_width, plot_height=700,
                x_range=x_range, y_range=(0,S), tools=tools,#"xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                #text_font="monospace",
                text_font_size=fontsize)
               
    rects = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)



    p2 = figure(title=None, plot_width=plot_width, plot_height=700,
                x_range=x_range, y_range=(0,S), tools=tools,#"xpan,reset",
                min_border=0, toolbar_location='below')#, lod_factor=1)          
    glyph2 = Text(x="x", y="y", text="text", text_align='center',text_color="black",
                #text_font="monospace",
                text_font_size=fontsize)
               
    rects2 = Rect(x="x", y="recty",  width=1, height=1, fill_color="colors",
                line_color=None, fill_alpha=0.4)

    p2.add_glyph(source, glyph)
    p2.add_glyph(source, rects)
    
    p2.grid.visible = False
    p2.xaxis.major_label_text_font_style = "bold"
    p2.yaxis.minor_tick_line_width = 0
    p2.yaxis.major_tick_line_width = 0

    p1.add_glyph(subc_source, glyph2)
    p1.add_glyph(subc_source, rects2)

    p1.grid.visible = False
    p1.xaxis.major_label_text_font_style = "bold"
    p1.yaxis.minor_tick_line_width = 0
    p1.yaxis.major_tick_line_width = 0
    

    





    p = gridplot([[p1,p2]], toolbar_location='below')
    #p_just_subc = gridplot([[p2]], toolbar_location='below')
    #cola,colb=st.columns([3,7])
      
    #cola.pyplot(dendo)#, height=21)
    
    st.bokeh_chart(p)
    #colb.bokeh_chart(p_just_subc)
    return sum(counts)

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
