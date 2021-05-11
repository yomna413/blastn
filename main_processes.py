"""
Created on Wed Dec 16 15:21:49 2020

"""
import pandas as pd
import math
import re
import sqlite3
exec(open('helpers.py').read())
""" Imports the database to search in.
    
    Args:
        filepath(str) : Database file path.
        
    Returns :
        listed_data(lst)(str) : A list of sequences found in database
"""

def import_data(filepath):
    if(re.search('.txt',filepath)):
       reader = open(filepath,'r')
       return re.split('\n',reader.read())
    else:
        conn = sqlite3.connect(filepath)
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM db_sequences")
        data = cursor.fetchall()
        listed_data = []
        for i in range(len(data)):
            listed_data.append(data[i][0])
        return listed_data

""" Finds the high scoring segment pairs.

    Args:
        querry_seq(str) : The querry sequence we search for
        db_seq(str) : One of the source database's sequences
        k(int) : Word length
        match_score(float) : Score added on match occurunce in ungapped alignment
        mismatch_score(float) : Score added on mismatch occurunce in ungapped alignment
        hssp_threshold(float) : The score threshold at which we decide whether the 
                         considered kmer is a hssp or not
       
    Returns:
        hssp_table(data.frame) : The high scoring segment pairs table, contains the kmer,
                     the starting index of querry kmer, the starting index of
                     db_kmer & the ungapped alignment score.
"""
def find_hssps(querry_seq,db_seq,k,match_score,mismatch_score,hssp_threshold):   
    querry_kmers = extract_kmers(querry_seq,k)
    for i in range(len(querry_kmers)-1):
        if(querry_kmers[i]==querry_kmers[i+1]):
            querry_kmers[i+1] = 'N'*k
    db_kmers = extract_kmers(db_seq,k)
    kmers =[]
    q_indicies=[]
    db_indicies=[]
    scores=[]
    q_kmer = []
    for i in range(len(querry_kmers)):
        for j in range(len(db_kmers)):
            score = ungapped_alignment(querry_kmers[i],db_kmers[j],match_score,mismatch_score)
            if(score>=hssp_threshold):
                kmers.append(db_kmers[j])
                q_kmer.append(querry_kmers[i])
                q_indicies.append(i)
                db_indicies.append(j)
                scores.append(score)
    data = {'db_kmer':kmers,
            'querry_index':q_indicies,
            'db_index':db_indicies,
            'score': scores,
            'q_kmer':q_kmer
            }
    hssp_table = pd.DataFrame(data,columns = ["db_kmer", "querry_index","db_index","score", "q_kmer"])        
    return hssp_table

""" Extends the seeded hssps

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
        (data.frame) : A data frame that contains the possible alignments & their raw alingment
                       scores ordered by that socre. 
                    
"""

def seed_extend(hssp_table,querry_seq,db_seq,k,
                match_score,mismatch_score,gap_score,
                extend_without_checking,extention_threshold):
    querry_alignments = []
    db_alignments = []
    raw_scores = []
    maxOutput1 = ""
    maxOutput2 = ""
    maxV =0
    best_i=0 
    for i in range(hssp_table.shape[0]):
       r,q,d = SWM(querry_seq,db_seq,hssp_table.iloc[i,1],hssp_table.iloc[i,2],k,
                    match_score,mismatch_score,gap_score,
                    extend_without_checking, extention_threshold,hssp_table.iloc[i,3])
       if( r>maxV ):
           maxOutput1 = q
           maxOutput2 = d
           maxV = r
           best_i = i
    final_db = hssp_table.iloc[best_i,0]+maxOutput2
    final_query = hssp_table.iloc[best_i,4]+maxOutput1
    #print(final_db)
   
    return final_query , final_db , maxV
"""    
    data = {'querry_alignment':querry_alignments,
            'db_alignment':db_alignments,
            'raw_score':raw_scores}
    df = pd.DataFrame(data)
    #return df.drop(df[df['raw_score'] < extention_threshold].index, inplace = True) 
    return df.sort_values(by=['raw_score'],ascending = 0).head(10)
"""
    
    
"""Assembles the algorithm for a database instance.
      Args:
        hssp_table(data.frame) : A dataframe of hssps, the return of find_hssp function.
        querry_seq(str) : The querry sequence
        db_seq(str) : The database sequence
        k(int): Word length
        match_score(int): Socre given on match occurunce
        mismatch_score(int) : Socre given on mismatch occurunce
        gap_score(int) : Socre given on gap placment
        hssp_threshold(float) : The score threshold at which we decide whether the 
                                considered kmer is a hssp or not
        extention_threshold(int) : The threshold alignment score that if an extention falls
                                   behind it, the extention stops
                                   
                                   
      Returns : 
         (data.frame) : A data frame that contains the possible alignments & their raw alingment
                        scores ordered by that socre. 

"""    
def blast_pipeline(querry_seq,db_seq,k,match_score,mismatch_score,gap_score,
                   hssp_threshold,extend_without_checking,extention_threshold):
    table = find_hssps(querry_seq,db_seq,k,match_score,mismatch_score,hssp_threshold)
    return seed_extend(table,querry_seq,db_seq,k,
                       match_score,mismatch_score,gap_score,
                       extend_without_checking,extention_threshold)
    
    

    
    
    
    
    
    
    