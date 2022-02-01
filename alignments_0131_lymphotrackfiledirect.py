from bdb import set_trace
import streamlit as st
#import neatbio.sequtils as utils
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from collections import Counter
import base64, zlib
# data pkgs
import matplotlib.pyplot as plt
import matplotlib
from Bio.SeqRecord import SeqRecord
import numpy as np
import pdb
import gzip
import shutil
from itertools import islice
import pandas as pd
st.set_page_config(layout="wide")
from streamlit import caching

import altair as alt  
import bokeh
from scipy.cluster.hierarchy import dendrogram, linkage
import pickle
import re
from Bio.Seq import Seq
from view_alignments_0131 import *
import time

# generator function 
def over_slice(test_str, K):
    itr = iter(test_str)
    res = tuple(islice(itr, K))
    if len(res) == K:
        yield res    
    for ele in itr:
        res = res[1:] + (ele,)
        yield res


def clone_kmer_show(fr3_clone,fr2_clone,fr1_clone):
    # initializing K
    
    clonelengthmax= np.max((len(fr3_clone),len(fr2_clone),len(fr1_clone)))
    K = st.slider('Kmer length',10,int(clonelengthmax),35)
    independentK = st.slider('Number of Independent Kmers',1,10,2)

    

    pd_res=pd.concat((pd.DataFrame(["".join(ele) for ele in over_slice(fr3_clone, K)]),\
    pd.DataFrame(["".join(ele) for ele in over_slice(fr2_clone, K)]),\
    pd.DataFrame(["".join(ele) for ele in over_slice(fr1_clone, K)])))    
    pd_res=pd_res.drop_duplicates()
   
    pd_res.rename(columns={0:'Seq_kmer'},inplace=True)

    # printing result
    st.write("Overlapping windows:")
    with st.expander('table of kmer with length '+str(K)):
        st.table(pd_res)
      
    return pd_res['Seq_kmer'].tolist(), K, independentK

def subextra(clone):
    clone=re.sub('[ ]',"",clone)
    clone=re.sub('\n','',clone) 
    return clone



with st.form(key='prep'):
    seq_file3 = st.file_uploader("Upload Framework 3 File",type = [".tsv",".txt"])
    seq_file2 = st.file_uploader("Upload Framework 2 File",type = [".tsv",".txt"])
    seq_file1 = st.file_uploader("Upload Framework 1 File",type = [".tsv",".txt"])

    fr3_clone = st.text_input("Framework3 Clone sequence ",\
        "ACACGGCTGTGTATTACTGTGCGAGAGATCTGAGGAACGAAATGGCGGGGGATGGTGCGGGGAGTTAGTCGACTACTACTACTACTACATGTACGTCTGGGGCAAAGGGACCAC")
  
    fr2_clone = st.text_input("Framework2 Clone sequence ",\
        "AGGGAAGGGGCTGGAGTGGGTGGCCAACATAAAGCAAGATGGATGTGAGAAATACTATGTGGACTCTGTGAAGGGCCGATTCCCCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGAGATCTGAGGAACGAAATGGCGGGGGATGGTGCGGGGAGTTAGTCGACTACTACTACTACTACATGTACGTCTGGGGCAAAGGGACCAC")    #fr2_clone = st.text_input("Framework2 Clone sequence 2",\

    fr1_clone = st.text_input("Framework1 Clone sequence ",\
        "CGCTGTCTATGGTGGGTCCTTCAGTGGTTACTATTGGACCTGGATCCGCCAGCTCCCAGGGGAGGGTCTGGAGTGGATTGGGGAAATCAGTCATAGTGGAAGCACCAACTACAATCCGTCCCTCAAGAGTCGAGTCACCATGTCAGTAGACACGTCCAAAAACCAGTTCTCCCTGAATTTGAGTTCTGTGACCGGCGCGGACACGGCTACATATTTTTGTGCGAGAGGCGGATTTTTCGGATGGGTGACTCGATCCAGGGCTCACTATCTTGACTACTGGGGCCAGGGAACCCT")

    fr3_clone=subextra(fr3_clone)
    fr2_clone=subextra(fr2_clone)
    fr1_clone=subextra(fr1_clone)


    to_rev_complement=st.radio('reverse complement input',['No','Yes'])
    if to_rev_complement == 'Yes':
        fr3_clone=str(Seq(fr3_clone).reverse_complement())
        fr2_clone=str(Seq(fr2_clone).reverse_complement())
        fr2_clone=str(Seq(fr1_clone).reverse_complement())
    
  
    res, K, independentK=clone_kmer_show(fr3_clone,fr2_clone,fr1_clone)
    restart=st.checkbox('start over')
    if restart:
        caching.clear_cache()
    submit_button = st.form_submit_button(label='Submit')


def replace_char_at_index(org_str, index, replacement):
    ''' Replace character at index in string org_str with the
    given replacement character.'''
    new_str = org_str
    if index < len(org_str):
        new_str = org_str[0:index] + replacement + org_str[index + 1:]
    return new_str
  
    
def remove_nonletter(objectalign):
    ## 01/10 changed seqA and seqB 
    seqA=objectalign[0][0]
    seqB=objectalign[0][1]
    seqB_forfile=seqB
    skips=re.finditer('-',seqA)
    counter=0
    for jj in skips:
       
        seqB=replace_char_at_index(seqB,jj.start()-counter,'')
        counter+=1
    return seqB, seqB_forfile    



def get_lymphotrackdirectfile(lrfile):
    lymphotrackfile=pd.read_table(lrfile,header=None)
    lrinfo = lymphotrackfile.iloc[np.arange(0,lymphotrackfile.shape[0],2)]
    lrsequence= lymphotrackfile.iloc[np.arange(1,lymphotrackfile.shape[0],2)]
    numseen=lrinfo[0].apply(lambda x: int(re.sub('_.+','',re.sub('.+Tcount','',str(x)))))
    lenseq = lrinfo[0].apply(lambda x: int(re.sub('.+length','',str(x))))
    lrdf = pd.DataFrame({'count':np.array(numseen),'seqlen':np.array(lenseq),'seq':lrsequence[0]})
    return lrdf
 
def process_seqs(seq_file):
    return get_lymphotrackdirectfile(seq_file)
    
class subclone:
    def __init__(self,fr3,fr2,fr1,resuse,independentK,fr3_clone,fr2_clone,fr1_clone):        
        self.res = resuse
        self.fr3_clone = fr3_clone
        self.fr2_clone = fr2_clone
        self.fr1_clone = fr1_clone
        self.fr3 = fr3
        self.fr2 = fr2
        self.fr1 = fr1
        self.independentK= independentK
        if fr3 is not None:
            self.clust0= process_seqs(fr3)
            st.write('framework 3 kmer match in progress')
            self.clust0_clone = self.process_match(self.clust0)
            if self.fr3_clone !='' and any(self.clust0_clone.seq.str.contains(self.fr3_clone)):
                self.abundant_seqfr3=self.fr3_clone
            else:        
                self.abundant_seqfr3=self.clust0_clone.seq.iloc[0]            
        if fr2 is not None:   
            self.clust1 = process_seqs(fr2)
            st.write('framework 2 kmer match in progress')
            self.clust1_clone = self.process_match(self.clust1)
            if self.fr2_clone !='' and any(self.clust1_clone.seq.str.contains(self.fr2_clone)):
                self.abundant_seqfr2=self.fr2_clone
            else:        
                self.abundant_seqfr2=self.clust1_clone.seq.iloc[0]              
        if fr1 is not None:   
            self.clust2 = process_seqs(fr1)
            st.write('framework 1 kmer match in progress')
            self.clust2_clone = self.process_match(self.clust2)
            if self.fr1_clone !='' and any(self.clust2_clone.seq.str.contains(self.fr1_clone)):
                self.abundant_seqfr1=self.fr1_clone
            else:        
                self.abundant_seqfr1=self.clust2_clone.seq.iloc[0]              

    def process_match(self,ses): 
           
        ses['foundsubstr']=0#np.array()   
        ses['foundindepKmer']=False#np.array()    
        pro=st.progress(0)
        countcol=0
        for indepKs in np.arange(0,len(self.res[0])):
            ses['foundsubstr']=0
            for indsingle_kmer in np.arange(indepKs, len(self.res),len(self.res[0])):
                single_kmer=res[indsingle_kmer]
                ses['foundsubstr']=ses['foundsubstr'] + ses.seq.str.contains(single_kmer)
                countcol+=1
               
            ses['foundindepKmer'] = ses['foundindepKmer']+ (self.independentK <= ses['foundsubstr'])  
                
            pro.progress(countcol/len(self.res))
                
        ses.reset_index(inplace=True)
        ses.rename(columns={0:'seq'},inplace=True)
        ses=ses[ses['foundindepKmer']]
        return(ses)
    
    def get_partclus_and_aln(self,partclus,clone,abundant_seq):
        if any(partclus['seq']==clone):
            abundant_seq_or_fr_clone=clone
        else:
            abundant_seq_or_fr_clone  = abundant_seq
        ii=0    
        ppp1=list()
        cols_show = st.columns(10)
        for x in partclus['seq']:  
            ppp1.append(pairwise2.align.globalms(abundant_seq_or_fr_clone,x,1,-1,-5,0))
            ii+=1
            if ii % 1000 == 0:
                cols_show[int(ii / 1000) % 10].write(str(ii) + 'aligned')   
        #ppp1=[pairwise2.align.globalms(abundant_seq_or_fr_clone,x,1,-1,-5,0) for x in partclus['seq']]
        
        New_seqBs=[remove_nonletter(aaa) for aaa in ppp1]
        
        partclus['score']=[ii[0].score for ii in ppp1]
        
        aln=[SeqRecord(seq=jjj[0][0],name=str(jjj[2]),description=str(jjj[1]),id=str(jjj[3]),annotations={'seq':jjj[0][1]}) for jjj in zip(New_seqBs,partclus['count'].tolist(),partclus['score'],partclus['seq'].tolist())]
        aln=[aln[jjj] for jjj in np.argsort(partclus['score'])]
        partclus=partclus.sort_values(['score'])
       
        score_touse=[pairwise2.align.localms(gg.id,clone,1,-1,-5,0,score_only=True) for gg in aln]
        
        scoremax=np.max(score_touse)
        argscores_touse=np.where(score_touse==scoremax)[0]
        
        ### 0108 this used to be pairwise2.align.localms(aln[argscore_touse].seq,clone,1,-1,-5,0)
        ids_backinto_sequences=[aln[jj].id for jj in argscores_touse]
        
        align_toclone_togettruncation=[pairwise2.align.localms(id,clone,1,-1,-5,0) for id in ids_backinto_sequences]
        startpoints =[jj[0].start for jj in align_toclone_togettruncation]
        argscore_touse=argscores_touse[np.argmax(startpoints)]
        align_toclone_togettruncation=align_toclone_togettruncation[np.argmax(startpoints)]
        clipright=0
       
        clipleft=align_toclone_togettruncation[0].start
        if clipleft == 1:
            clipleft = 0

        ## added 0108
        makesurewedonotneedtotrimseq=align_toclone_togettruncation=pairwise2.align.localms(aln[argscore_touse].seq,clone,1,-1,-5,0)
        #if makesurewedonotneedtotrimseq[0].start > 0:
            #st.write('SEQ MAY NEED TRIM')
            #pdb.set_trace()
        ###     
        aln_sortbyscore=[hh[0] for hh in zip(aln,score_touse) if not hh[1]<=len(clone)-15]
    
        if len(aln_sortbyscore)>1:
            aln=aln_sortbyscore
            partclus=partclus[np.array(score_touse)>len(clone)-30]
        
        return partclus, aln, clipleft, clipright

    def create_visual(self):

        if self.fr3 is not None:
            if self.fr3_clone != "":
                useclone=self.fr3_clone
            elif self.fr2_clone !="":
                useclone=self.fr2_clone  
            else:
                useclone=self.fr1_clone
            with st.spinner('formatting framework 3 sequences for '+ str(self.clust0_clone.shape[0])+' sequences and calculating alignment scores'):    
                part1, aln1, _cl,_cr = self.get_partclus_and_aln(self.clust0_clone,useclone,self.abundant_seqfr3)            
            if len(aln1)>0:
                view_alignment(aln1, self.abundant_seqfr3,3)
        if self.fr2 is not None:
            if self.fr2_clone != "":
                useclone=self.fr2_clone
            elif self.fr1_clone !="":
                useclone=self.fr1_clone  
            else:
                useclone=self.fr3_clone   
            with st.spinner('formatting framework 2 sequences for '+ str(self.clust1_clone.shape[0])+' sequences and calculating alignment scores'):                
                part2, aln2, _cl,_cr = self.get_partclus_and_aln(self.clust1_clone,useclone,self.abundant_seqfr2)
            if len(aln2)>0:    
                vv2,cluster_n2,fr2_ratioandreads=view_alignment(aln2, self.abundant_seqfr2,2)       
        if self.fr1 is not None:
            if fr1_clone != "":
                useclone=fr1_clone
            elif fr2_clone !="":
                useclone=fr2_clone  
            else:
                useclone=fr3_clone
            with st.spinner('formatting framework 1 sequences for '+ str(self.clust2_clone.shape[0])+' sequences and calculating alignment scores'):                     
                part3, aln3, _cl,_cr= self.get_partclus_and_aln(self.clust2_clone,useclone,self.abundant_seqfr1)
            if len(aln3)>0:
                view_alignment(aln3, self.abundant_seqfr1,1)


if seq_file1 is not None or seq_file2 is not None or seq_file3 is not None:
    
    if not(fr3_clone == '' and fr2_clone == '' and  fr1_clone == ''):
        st.write('input fr3 clone sequence is: '+fr3_clone)
        st.write('input fr2 clone sequence is: '+fr2_clone)
        st.write('input fr1 clone sequence is: '+fr1_clone)
                    
        subclone_obj=subclone(seq_file3,seq_file2,seq_file1,res,independentK, fr3_clone,fr2_clone,fr1_clone)
                
        subclone_obj.create_visual()

                
 