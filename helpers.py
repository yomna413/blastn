"""
Created on Wed Dec 16 15:21:49 2020

"""
import numpy as np
""" Extracts every possible kmer from a sequence.

    Args:
        sequence(str): Sequence to extracts kmers from
        k(int): Word length
              
    Returns:
        kmers(list)(str):Every possible sequence kmer
 """
def extract_kmers(sequence,k):
    kmers=[]
    kmer = ""
    for i in range(len(sequence)-k+1):
        kmer = sequence[i]
        for j in range(1,k):
            kmer+=sequence[i+j]      
        kmers.append(kmer)
    return kmers



""" Preforms naive ungapped alignment between two kmers.

    Args:
        kmer1(str): First Sequence
        kmer2(str): Second Sequence
        match_score(float): Score added on match occurunce in ungapped alignment
        mismatch_score(float) : Score added on mismatch occurunce in ungapped alignment
              
    Returns:
        (float) Alignment Score
 """
def ungapped_alignment(kmer1,kmer2,match_score,mismatch_score):
    scores = np.zeros(len(kmer1))
    for i in range(len(kmer1)):
        if(kmer1[i-1]==kmer2[i-1]):
            scores[i] = scores[i-1]+match_score
        else:
            scores[i] = scores[i-1]+mismatch_score
    return scores[-1]

""" The core process of extention through Smith-Waterman algorithm.
     
    Args:
        hssp_table(data.frame) : A dataframe of hssps, the return of find_hssp function.
        querry_seq(str) : The querry sequence
        db_seq(str) : The database sequence
        k(int): Word length
        match_score(int): Socre given on match occurunce
        mismatch_score(int) : Socre given on mismatch occurunce
        gap_score(int) : Socre given on gap placment
        extention_threshold(int) : The threshold alignment score that if an extention falls
                                   behind it, the extention stops
    Returns : 
        raw_score(int): Raw alignment score between two extended kmers
        querry_alignment : querry kmer after extention
        db_alignment : db kmer after extention
"""


def compare(i,j, match_score, mismatch_score):
    if i==j :
        return match_score
    else:
        return mismatch_score




def SWM(querry_seq,db_seq,querry_index,db_index,k,
        match_score,mismatch_score,gap,
        extend_without_checking, extention_threshold,score):
    querry_alignment = ""
    db_alignment = ""
    raw_score = 0
    #dp = np.zeros((len(querry_seq),len(db_seq)))
   
    resume = True
    #row = len(querry_seq) - 1
    #col = len(db_seq) - 1

    n1 = len(querry_seq)+1
    n2 = len(db_seq)+1 
    dp = [[0]*n1 for i in range(n2)]
    maxI = 0
    maxJ = 0
    maxV = 0
    dp[db_index+k-1][querry_index+k-1] = score 
    for i in range(db_index+k,n2):
        for j in range(querry_index+k,n1):
             dp[i][j] = max(dp[i - 1][j - 1] + compare(db_seq[i - 1],querry_seq[j - 1], match_score, mismatch_score) ,  
                                 dp[i - 1][j] + gap ,  
                                 dp[i][j - 1] + gap , 
                                 0) 
             if dp[i][j] > maxV:
                 maxV = dp[i][j]
                 maxI = i
                 maxJ = j
   
    #backtrace
    
    output1 = ""
    output2 = ""
    i = maxI
    j = maxJ
   
    while (i>db_index+k-1) or (j>querry_index+k-1):
        if(i>db_index+k-1 and j> querry_index+k-1 and dp[i][j] == dp[i-1][j-1]+ compare(db_seq[i-1],querry_seq[j-1], match_score, mismatch_score)):
            output2 = db_seq[i-1] + output2
            output1 = querry_seq[j-1] + output1
            i -=1
            j -=1
        elif (i>0 and dp[i][j] ==dp[i-1][j]+gap):
            output2 = db_seq[i-1] + output2
            output1 = "_" + output1
            i-= 1
        elif (j>0 and dp[i][j] ==dp[i][j-1]+gap): 
            output2 = "_" + output2
            output1 = querry_seq[j-1] + output1
            j-= 1
        else:
            break
        
    
    return maxV,output1,output2
              
