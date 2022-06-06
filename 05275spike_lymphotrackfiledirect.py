from bdb import set_trace
#from tkinter.messagebox import NO
import streamlit as st
#import neatBio.sequtils as utils
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
import sys
from view_alignments_0401 import *
from KevinsFunction import *
from Kevinsmerge import *
import time
sys.setrecursionlimit(5000)
# generator function 
def over_slice(test_str, K):
    itr = iter(test_str)
    res = tuple(islice(itr, K))
    if len(res) == K:
        yield res    
    for ele in itr:
        res = res[1:] + (ele,)
        yield res


def clone_kmer_show(K,fr_clone):
    # initializing K
    
    #clonelengthmax= np.max((len(fr3_clone),len(fr2_clone),len(fr1_clone)))
    #K = st.slider('Kmer length',10,int(clonelengthmax),35)
    #independentK = st.slider('Number of Independent Kmers',1,10,2)

    

    pd_res=pd.DataFrame(["".join(ele) for ele in over_slice(fr_clone, K)]) 
    pd_res=pd_res.drop_duplicates()
   
    pd_res.rename(columns={0:'Seq_kmer'},inplace=True)

    return pd_res['Seq_kmer'].tolist()

def subextra(clone):
    clone=re.sub('[ ]',"",clone)
    clone=re.sub('\n','',clone) 
    return clone

detectormerge=st.radio('What are we doing here?',('Detect and quantify subclones in one framework','merge the results of more than one framework'))


with st.form(key='prep'):
    seq_file = None
    files_to_merge = None
    if detectormerge == 'Detect and quantify subclones in one framework':
        #seq_file = st.file_uploader("Upload Framework File",type = [".tsv",".txt"])

        seq_file3a = st.file_uploader("Upload Framework 3 File",type = [".tsv",".txt"])
        seq_file3b = st.file_uploader("Upload 2nd Framework 3",type = [".tsv",".txt"])
        seq_file3c = st.file_uploader("Upload 3rd Framework 3",type = [".tsv",".txt"])        
        seq_file2a = st.file_uploader("Upload Framework 2 File",type = [".tsv",".txt"])
        seq_file2b = st.file_uploader("Upload 2nd Framework 2" ,type = [".tsv",".txt"])
        seq_file2c = st.file_uploader("Upload 3rd Framework 2",type = [".tsv",".txt"])         
        seq_file1a = st.file_uploader("Upload Framework  1 File",type = [".tsv",".txt"])
        seq_file1b = st.file_uploader("Upload 2nd Framework 1 File",type = [".tsv",".txt"])
        seq_file1c = st.file_uploader("Upload 3rd Framework 1",type = [".tsv",".txt"])         



        fr3_clone2 = st.text_input("Framework3 Clone sequence ",\
            "TCTGAGGACACGGCCGTGTATTACTGTGCGAGAGATAGGCGCGGGGAATGGCCTCCCTCGGATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCAC")
        fr3_clone = st.text_input("Framework3 Clone sequence 2",\
            "ACACGGCTGTGTATTACTGTGCGAGCCAGATATTGTAGTGGTGGTAGCTCCCTATCGGGGAGCTTTTGATATCTGGGGCCAAGGGACAAT")      
        fr2_clone2 = st.text_input("Framework2 Clone sequence ",\
            "TGGACAAGGGCTTGAGTGGATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATTACCGCGGACGAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGATAGGCGCGGGGAATGGCCTCCCTCGGATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCAC")
        fr2_clone = st.text_input("Framework2 Clone sequence 2",\
            "AGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGCCAGATATTGTAGTGGTGGTAGCTCCCTATCGGGGAGCTTTTGATATCTGGGGCCAAGGGACAAT")    
        fr1_clone2 = st.text_input("Framework1 Clone sequence ",\
            "CTTCTGGAGGCACCTTCAGCAGCTATGCTATCAGCTGGGTGCGACAGGCCCCTGGACAAGGGCTTGAGTGGATGGGAGGGATCATCCCTATCTTTGGTACAGCAAACTACGCACAGAAGTTCCAGGGCAGAGTCACGATTACCGCGGACGAATCCACGAGCACAGCCTACATGGAGCTGAGCAGCCTGAGATCTGAGGACACGGCCGTGTATTACTGTGCGAGAGATAGGCGCGGGGAATGGCCTCCCTCGGATTACTACTACTACTACTACATGGACGTCTGGGGCAAAGGGACCAC")
        fr1_clone = st.text_input("Framework1 Clone sequence 2",\
            "GCCTCTGGATTCACCTTCAGTAGCTATAGCATGAACTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTCTCATCCATTAGTAGTAGTAGTAGTTACATATACTACGCAGACTCAGTGAAGGGCCGATTCACCATCTCCAGAGACAACGCCAAGAACTCACTGTATCTGCAAATGAACAGCCTGAGAGCCGAGGACACGGCTGTGTATTACTGTGCGAGCCAGATATTGTAGTGGTGGTAGCTCCCTATCGGGGAGCTTTTGATATCTGGGGCCAAGGGACAAT")


        fr3_clone=subextra(fr3_clone)
        fr2_clone=subextra(fr2_clone)
        fr1_clone=subextra(fr1_clone)
        
        fr3_clone2=subextra(fr3_clone2)
        fr2_clone2=subextra(fr2_clone2)
        fr1_clone2=subextra(fr1_clone2)

        #fr_clone = st.text_input("Framework Clone sequence ",\
        #    "ACACGGCTGTGTATTACTGTGCGAGAGATCTGAGGAACGAAATGGCGGGGGATGGTGCGGGGAGTTAGTCGACTACTACTACTACTACATGTACGTCTGGGGCAAAGGGACCAC")
        #fr = st.radio('Framework',(1,2,3),2)
        #fr_clone=subextra(fr_clone)
        #to_rev_complement=st.radio('reverse complement input',['No','Yes'])
        #if to_rev_complement == 'Yes':
        #    fr_clone=str(Seq(fr_clone).reverse_complement())    
        K3 = st.slider('FR3 Kmer length',10,len(fr1_clone2),27)
        K2 = st.slider('FR3 Kmer length',10,len(fr1_clone2),57)
        K1 = st.slider('FR3 Kmer length',10,len(fr1_clone2),77)
        independentK = st.slider('Number of Independent Kmers',1,10,1)        
        res_list3a =clone_kmer_show(K3,fr3_clone)
        res_list3b =clone_kmer_show(K3,fr3_clone2)
        res_list2a =clone_kmer_show(K2,fr2_clone)
        res_list2b =clone_kmer_show(K2,fr2_clone2)
        res_list1a =clone_kmer_show(K1,fr1_clone)
        res_list1b =clone_kmer_show(K1,fr1_clone2)
    else:
        files_to_merge=st.file_uploader('Files to merge',type=[".csv"],accept_multiple_files=True)

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
    ind_pass=list()
    skips=re.finditer('-',seqA)
    counter=0
    for jj in skips:
       
        seqB=replace_char_at_index(seqB,jj.start()-counter,'')
        ind_pass.append(jj.start()-counter)
        counter+=1
    return seqB, seqB_forfile,ind_pass    



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



class merge_files:
    def __init__(self,seq_files):
        self.fr1_subclones = pd.read_csv(seq_files[0])
        self.fr2_subclones = pd.read_csv(seq_files[1])
        self.fr3_subclones = pd.read_csv(seq_files[2])
        self.merged_list = self.merge()
        self.merge_table= pd.concat(self.merged_list)
        st.download_button(label='merge subclone file',data=convert_df(self.merge_table),file_name='merge_data.csv',mime='text/csv') 

    def merge(self):
        outliers = [self.fr1_subclones, self.fr2_subclones, self.fr3_subclones]
        framework1_merge, framework2_merge, framework3_merge = merge_outliers(outliers)
        framework1_merge['framework'] = 1   
        framework2_merge['framework'] = 2
        framework3_merge['framework'] = 3
        
        return [framework1_merge, framework2_merge, framework3_merge]





class subclone:
    def __init__(self,seq_file,fr,resuse,independentK,fr_clone,must_be_identical_percent):      
        self.clust= process_seqs(seq_file)
        self.must_be_identical_percent = must_be_identical_percent
        self.first(fr,resuse,independentK,fr_clone)

    def first(self,fr,resuse,independentK,fr_clone):      
        self.res = resuse
        self.fr_clone = fr_clone
        self.fr = fr
        self.see_abundant=False
        self.fr_subcloneseqs=list()
        self.independentK= independentK

        st.write('framework kmer match in progress')
        self.fr_total, self.clust_clone = self.process_match(self.clust)
        st.write('formatting framework sequences for '+ str(self.clust_clone.shape[0]))
        if not self.see_abundant or self.fr_clone !='' and any(self.clust_clone.seq.str.contains(self.fr_clone)):
            self.abundant_seqfr=self.fr_clone
        else:        
            self.abundant_seqfr=self.clust_clone.seq.iloc[0]   
                        
        self.create_visual()
        cola,colb,colc = st.columns(3)
        cola.write('match to clone total '+str(self.fr_total)) 
        colb.write('total '+str(sum(self.clust['count'])))                 
        colc.write('clone fraction '+str(np.round(self.fr_total/sum(self.clust['count']),4)))

    def process_match(self,ses): 
           
        ses['foundsubstr']=0#np.array()   
        ses['foundindepKmer']=False#np.array()    
        pro=st.progress(0)
        countcol=0

        for indepKs in np.arange(0,len(self.res[0])):
            ses['foundsubstr']=0
            for indsingle_kmer in np.arange(indepKs, len(self.res),len(self.res[0])):
                single_kmer=self.res[indsingle_kmer]
                ses['foundsubstr']=ses['foundsubstr'] + ses.seq.str.contains(single_kmer)
                countcol+=1
               
            ses['foundindepKmer'] = ses['foundindepKmer']+ (self.independentK <= ses['foundsubstr'])  
                
            pro.progress(countcol/len(self.res))
        if ses.index[0] != 0:       
            ses.reset_index(inplace=True)
            ses.rename(columns={0:'seq'},inplace=True)
        ses_countsum=sum(ses[ses['foundindepKmer']]['count'])
        ses=ses[ses['foundindepKmer']].sort_values('count',ascending=False)#[0:5000]
        return(ses_countsum,ses)
    
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
        
        aln=[SeqRecord(seq=jjj[0][0],name=str(jjj[2]),description=str(jjj[1]),id=str(jjj[3]),annotations={'seq':jjj[0][1],'insert':jjj[0][2]}) for jjj in zip(New_seqBs,partclus['count'].tolist(),partclus['score'],partclus['seq'].tolist())]
        aln=[aln[jjj] for jjj in np.argsort(partclus['score'])]
        partclus=partclus.sort_values(['score'])
       
        score_touse=[pairwise2.align.localms(gg.id,clone,1,-1,-5,-1,score_only=True) for gg in aln]
        
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
        #makesurewedonotneedtotrimseq=align_toclone_togettruncation=pairwise2.align.localms(aln[argscore_touse].seq,clone,1,-1,-5,-1)
        #if makesurewedonotneedtotrimseq[0].start > 0:
            #st.write('SEQ MAY NEED TRIM')
            #pdb.set_trace()
        #### 
  
        aln_sortbyscore=[hh[0] for hh in zip(aln,score_touse) if not hh[1]<=(len(clone)*self.must_be_identical_percent )]
        
        st.write('keep ',str(len(aln_sortbyscore))+ 'kmer matches, filter '+str(np.sum([s <=(len(clone)* self.must_be_identical_percent ) for s in score_touse]))+' after pairwise alignment')
        

        if len(aln_sortbyscore)>1:
            aln=aln_sortbyscore
            #partclus=partclus[np.array(score_touse)>len(clone)-30]
        
        return partclus, aln, clipleft, clipright

    def create_visual(self):
        with st.spinner('formatting framework sequences for '+ str(self.clust_clone.shape[0])+' sequences and calculating alignment scores'):    
            part1, aln1, _cl,_cr = self.get_partclus_and_aln(self.clust_clone,self.fr_clone,self.abundant_seqfr)            
        
        if len(aln1)>1:
            tallysum, subcloneseqsfr=view_alignment(aln1, self.abundant_seqfr,self.fr) 
            self.fr_total = tallysum
            self.fr_subcloneseqs=subcloneseqsfr




if detectormerge == 'Detect and quantify subclones in one framework' and seq_file3a is not None and fr3_clone2 != '':
    
    for seq_filea, seqfileb,seqfilec,reslista,reslistb,fr_clone1,fr_clone2,must_be_identical_percent in zip([seq_file3a, seq_file2a, seq_file1a],[seq_file3b,seq_file2b,seq_file1b],[seq_file3c,seq_file2c,seq_file1c],[res_list3a,res_list2a,res_list1a],[res_list3b,res_list2b,res_list1b], [fr3_clone,fr2_clone,fr1_clone],[fr3_clone2,fr2_clone2,fr1_clone2],[0.75,0.85,0.89]):

        subclone_obj3a1=subclone(seq_filea,3,reslista,independentK, fr_clone1,must_be_identical_percent)
        subclone_obj3a1.first(3,reslistb,independentK, fr_clone2)
        subclone_obj3b1=subclone(seqfileb,3,reslista,independentK, fr_clone1,must_be_identical_percent)
        subclone_obj3b1.first(3,reslistb,independentK, fr_clone2)
        subclone_obj3c1=subclone(seqfilec,3,reslista,independentK, fr_clone1,must_be_identical_percent)
        subclone_obj3c1.first(3,reslistb,independentK, fr_clone2)      
    #subclone_obj3b1=subclone(seq_fileb ,3,res_list3a,independentK, fr3_clone)
    #subclone_obj3b1.first(3,res_list3b,independentK, fr3_clone2)

    #subclone_obj3c1=subclone(seq_file3c ,3,res_list3a,independentK, fr3_clone)
    #subclone_obj3c1.first(3,res_list3b,independentK, fr3_clone2)

    #subclone_obj2a1=subclone(seq_file2a,2,res_list2a,independentK, fr2_clone)
    #subclone_obj2a2=subclone(seq_file2a,2,res_list2b,independentK, fr2_clone2)   

    #subclone_obj2b1=subclone(seq_file2b,2,res_list2a,independentK, fr2_clone)
    #subclone_obj2b2=subclone(seq_file2b,2,res_list2b,independentK, fr2_clone2) 

    #subclone_obj2c1=subclone(seq_file2c,2,res_list2a,independentK, fr2_clone)
    #subclone_obj2c2=subclone(seq_file2c,2,res_list2b,independentK, fr2_clone2) 

    #subclone_obj1a1=subclone(seq_file1a,1,res_list3a,independentK, fr1_clone)
    #subclone_obj1a2=subclone(seq_file1a,1,res_list3b,independentK, fr1_clone2)  

    #subclone_obj1b1=subclone(seq_file1b,1,res_list3a,independentK, fr1_clone)
    #subclone_obj1b2=subclone(seq_file1b,1,res_list3b,independentK, fr1_clone2)  

    #subclone_obj1c1=subclone(seq_file1c,1,res_list3a,independentK, fr1_clone)
    #subclone_obj1c2=subclone(seq_file1c,1,res_list3b,independentK, fr1_clone2)  



elif files_to_merge is not None:
    resultsmerged_obj=merge_files(files_to_merge)

 
