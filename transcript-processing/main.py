from entropy import entropy
from bert import BERT

import pandas as pd
from fastGraph import fastGraph
from scipy.stats import ttest_ind as ttest
import glob
import random
import re
import csv

H = entropy()
wv = BERT()

def clean_text(text):
    # remove square brackets and their contents
    text = re.sub(r"\[.*?\]", "", text)
    # remove '+<' literal
    text = re.sub(r"\+<", "", text)
    # remove all commas and periods
    text = text.replace(",",'')
    text = text.replace(".",'')
    text = text.replace("?",'')
    # unintelligible
    text = text.replace("xxx",'')
    # remove anything after '+' (literal) to a space
    text = re.sub(r"\+.*?$", "", text)    
    # remove anythning that starts &= to space or end of line
    text = re.sub(r"&=.*?\s", "", text)
    text = re.sub(r"&=.*?$", "", text)
    text = re.sub(r"&-.*?\s", "", text)
    text = re.sub(r"&-.*?$", "", text)
    text = re.sub(r"&+.*?\s", "", text)
    text = re.sub(r"&+.*?$", "", text)
    # remove ‡
    text = text.replace("‡",'')
    # replace "_" with space
    text = text.replace("_",' ')
    # remove brackets but keep inside content
    text = re.sub(r"\((.*?)\)", r"\1", text)
    # remove chevron brackets but keep content
    text = re.sub(r"\<(.*?)\>", r"\1", text)
    # remove double spacings
    text = text.replace("  ", " ")
    text = text.replace("  ", " ")
    # trim white space
    text = text.strip()
    return text

def getTransript(fname, gram=False):

    fl = open(fname,'r')
    txt = fl.read()
    fl.close()

    txt = txt.replace("\n\t",' ')
    txt = txt.split("\n")

    tscpt = []
    gra = []
    ids = []

    for i in range(len(txt)-2):
        # let's get staggered * (lexical) content + grammatical tiers; only if aligned
        if txt[i].startswith('*') and (txt[i+2].startswith('%gra') or not gram):        
            ids.append(txt[i][txt[i].find('*'):txt[i].find(':')]) # get id for speaker        
            txt[i] = txt[i].replace(txt[i][txt[i].find('*'):txt[i].find(':')+1],'') 
            txt[i] = txt[i][0:txt[i].rfind(' ')]
            txt[i+2] = txt[i+2].replace(txt[i+2][txt[i+2].find('%'):txt[i+2].find(':')+1],'')
            tscpt.append(clean_text(txt[i]))
            gra.append(txt[i+2])
    
    # return trans, gra, ids
    return [tscpt,gra,ids]

def procChaFile(fname, outfile, baseline = '', mode = 'a', levels=[7,-1], k=10):

    print('Processing '+fname)
    [tscpt, gra, ids] = getTransript(fname)
    if baseline!='':
        [tscpt_bl, gra_bl, ids_bl] = getTransript(baseline)
        if len(tscpt_bl)<len(tscpt):
            # make tscpt the same size as tscpt_bl by cropping
            tscpt = tscpt[0:(len(tscpt_bl)-1)]

    vecs = []
    vecs_comp = []
    for i in range(len(tscpt)): # preload vectors; speeds things up
        if tscpt[i].strip()!="":
            vecs.append(wv(tscpt[i].strip(),level=levels)[0])
        else:
            vecs.append([])

    if baseline=='': # nb: baseline not implemented in tbi analysis
        vecs_comp = vecs

    # nb: raw text removed from data shared so as not to disseminate coelho corpus directly; see link below
    cols = ['fl','i','j','who_i','who_j','n_i','n_j','txt_i','txt_j', 'gra_i','gra_j','h_1','h_2','baseline', 'levels']
    df = pd.DataFrame(columns=cols) 

    for i in range(len(tscpt)):
        # get number of vectors in vecs_i set as n
        n_i = len(vecs[i])
        if n_i > 0:
            vecs_i = vecs[i]
            rg = range(max(i-k,0),min(i+k,len(tscpt)))
            for j in rg:                
                n_j = len(vecs_comp[j])
                if n_j>0:
                    vecs_j = vecs_comp[j]
                    h_val = H(vecs_i, vecs_j)
                    df = pd.concat([df, pd.DataFrame([[fname,i,j,ids[i],ids[j],
                        n_i,n_j,
                        tscpt[i],tscpt[j],gra[i],gra[j],
                        float(h_val[0]),float(h_val[1]), 
                        baseline, levels[0]]], columns=cols)])
    # df.to_csv(outfile, mode=mode, header=(mode=='w'), quoting=csv.QUOTE_NONE, quotechar='')
    df.to_csv(outfile, mode=mode, header=(mode=='w'), escapechar='\\')

# list *.cha files under the subfolder
# https://tbi.talkbank.org/access/English/Coelho.html
# nb: password protected, permission from talkbank.org easy to request
files_temp = glob.glob('n+tb/*.cha')

files = []
for i in range(0,len(files_temp)):
    [tscpt, gra, ids] = getTransript(files_temp[i])
    print(len(tscpt))
    longs = 0
    for j in range(len(tscpt)):
        # count the number of spaces in tscpt[j]
        if tscpt[j].count(' ')>=100:
            longs = longs + 1
    if longs == 0:
        files.append(files_temp[i])

print(files)
print(len(files)) # files for processing
lix = len(files)

for j in [1,2,3,4,5,6,7,8,9,10,11,12]:
    # \process the first file observed processes
    procChaFile(files[0],'processed.csv',levels=[j])
    for i in range(1,lix):
        print('Observed '+str(i))
        procChaFile(files[i],'processed.csv',levels=[j])
