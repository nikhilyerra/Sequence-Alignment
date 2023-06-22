#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from random import choice
from Bio.Seq import Seq
import copy
import random
from operator import itemgetter
import sys


# In[2]:


inputFileLocation = sys.argv[1]
match_val = int(sys.argv[2])
replace_val = int(sys.argv[3])
gap_val = int(sys.argv[4])
outputFileLocation = sys.argv[5]


# In[3]:


filename = 'C:/Users/jeeta/OneDrive/Desktop/output2.fasta'


# In[4]:


sequences = [i for i in SeqIO.parse(inputFileLocation, 'fasta')]
fragment1 = list(sequences[0].seq)
fragment2 = list(sequences[1].seq)
#print(fragment1)
#print(fragment2)


# In[5]:


sequences_list=[]
for indexs in range(len(sequences)):
    fragment = list(sequences[indexs].seq)
    sequences_list.append(fragment)
sequences_list.sort(key=len)
#sequences_list = sequences_list[0:10]
len_list = len(sequences_list)
#print(len_list)


# In[6]:


#gap_val=-3
#match_val=+1
#replace_val=-1


# In[7]:


def dovetail_method(fragment1,fragment2, ind1,ind2, list_dovetail_alignment):
    
    len1 = len(fragment1)
    len2 = len(fragment2)
    dp = [[0 for i in range(len2+1)]for j in range(len1+1)]
    dp_pos=[[0 for i in range(len2+1)]for j in range(len1+1)]
    #dp[len1+1][len2+1]
    for x in range (len1+1):
        dp[x][0]=0
        dp_pos[x][0]=(x,0)
    for y in range (len2+1):
        dp[0][y]=0
        dp_pos[0][y]=(0,y)
    for x in range (1,len1+1):
        for y in range (1,len2+1):
            #print ('fragment1 val=',fragment1[x-1])
            #print ('fragment2 val=',fragment2[y-1])
            dp[x][y]=dp[x-1][y]+gap_val
            dp_pos[x][y]=(x-1,y)
            if(dp[x][y]<dp[x][y-1]+gap_val):
                dp_pos[x][y]=(x,y-1)
                dp[x][y]=dp[x][y-1]+gap_val
            if fragment1[x-1] == fragment2[y-1]:
                if(dp[x][y]<dp[x-1][y-1]+match_val):
                    dp[x][y]=dp[x-1][y-1]+match_val
                    dp_pos[x][y]=(x-1,y-1)
            else:
                if(dp[x][y]<dp[x-1][y-1]+replace_val):
                    dp[x][y]=dp[x-1][y-1]+replace_val
                    dp_pos[x][y]=(x-1,y-1)
        #print('\n')
    maxm = dp[len1][len2];
    dp_pos_x = len1
    dp_pos_y = len2
    last = 2
    
    ## 
    for x in range(1,len1):
        if(maxm<dp[x][len2]):
            maxm=dp[x][len2]
            dp_pos_x = x
            dp_pos_y = len2
    for y in range(1,len2):
        if(maxm<dp[len1][y]):
            dp_pos_x = len1
            dp_pos_y = y
            maxm=dp[len1][y]
            last = 1
            
    list_traverse=[]
    pos_x_last = dp_pos_x
    pos_y_last = dp_pos_y
    pos_x_first = dp_pos_x
    pos_y_first = dp_pos_y
    
    ## traversing for backtracking the dp
    while dp_pos_x>0 and dp_pos_y>0:
        pos_x_first = dp_pos_x
        pos_y_first = dp_pos_y
        list_traverse.append(tuple((dp_pos_x,dp_pos_y)))
        temp_x = dp_pos[dp_pos_x][dp_pos_y][0]
        temp_y = dp_pos[dp_pos_x][dp_pos_y][1]
        dp_pos_x = temp_x
        dp_pos_y = temp_y
        
    #print(dp)
    #print(list_traverse)
    
    ## calculating length of common alignment in fragment1 and fragment2
    len_align_first = pos_x_last - pos_x_first + 1
    len_align_second = pos_y_last - pos_y_first + 1
    
    ## preparing alignment values list
    list_single_pair =[]
    list_single_pair.append(maxm)
    list_single_pair.append(ind1)
    list_single_pair.append(ind2)
    list_single_pair.append(len_align_first)
    list_single_pair.append(len_align_second)
    list_single_pair.append(last)
    
    ## appending alignment values to the main list
    list_dovetail_alignment.append(list_single_pair)
    
    #print(list_single_pair)
    return maxm


# In[8]:


def generate_dovetail_scores_for_lists(sequences_list,list_dovetail_alignment):
    for index1 in range (len(sequences_list)):
        if index1+1 <= len(sequences_list):
            for index2 in range (index1+1,len(sequences_list)):
                fragment1 = sequences_list[index1]
                fragment2 = sequences_list[index2]
                dovetail_method(fragment1,fragment2,index1,index2,list_dovetail_alignment)


# In[9]:


def edit_list(sequences_list,list_dovetail_alignment):
    sorted(list_dovetail_alignment, key=itemgetter(0))
    len_alignment_lists = len(list_dovetail_alignment)
    maxm = list_dovetail_alignment[len_alignment_lists-1][0]
    if maxm <= 0:
        return False
    ind1 = list_dovetail_alignment[len_alignment_lists-1][1]
    ind2 = list_dovetail_alignment[len_alignment_lists-1][2]
    len_align_first = list_dovetail_alignment[len_alignment_lists-1][3]
    len_align_second = list_dovetail_alignment[len_alignment_lists-1][4]
    last = list_dovetail_alignment[len_alignment_lists-1][5]
    fragment1 = sequences_list[ind1]#list(sequences[ind1].seq)
    fragment2 = sequences_list[ind2]#list(sequences[index1].seq)
    new_fragment = copy.deepcopy(fragment1)
    if last == 1:
        new_fragment = fragment1 + fragment2[len_align_second:len(fragment2)]
    else:
        new_fragment = fragment2 + fragment1[len_align_first:len(fragment1)]
    del sequences_list[ind1]
    del sequences_list[ind2-1]
    sequences_list.append(new_fragment)
    #print(fragment1)
    #print(fragment2)
    #print(new_fragment)
    return True
        
        


# In[10]:


def fragment_assembler(sequences_list):
    list_dovetail_alignment = [[]]
    status = True
    while status == True and len(sequences_list) > 1:
        list_dovetail_alignment.clear()
        generate_dovetail_scores_for_lists(sequences_list,list_dovetail_alignment)
        status = edit_list(sequences_list,list_dovetail_alignment)
    


# In[11]:


#list_dovetail_alignmentx=[]
#val = dovetail_method(fragment1,fragment2,0,1,list_dovetail_alignmentx)
#print(val)


# In[12]:


fragment_assembler(sequences_list)
#print(sequences_list)
sequences_list.sort(key=len)
len_list = len(sequences_list)
#print(sequences_list[len_list-1])
final_fragment = sequences_list[len_list-1]


# In[13]:


final_sequence = copy.deepcopy(sequences)
new_sequence = copy.deepcopy(final_sequence[0])
final_fragment=''.join(final_fragment)
new_sequence.seq = Seq(final_fragment)
new_sequence.description = 'final fragment'
final_sequence.clear()
final_sequence.append(new_sequence)


# In[14]:


SeqIO.write(final_sequence, outputFileLocation , "fasta")

