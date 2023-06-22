#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from random import choice
from Bio.Seq import Seq
import copy
import sys


# In[2]:


import random


# In[ ]:


inputFileLocation = sys.argv[1]
k = int(sys.argv[2])
p = float(sys.argv[3])
outputFileLocation = sys.argv[4]


# In[3]:


#inputFileLocation = 'C:/Users/jeeta/OneDrive/Desktop/inputfile.fasta'


# In[4]:


sequences = [i for i in SeqIO.parse(inputFileLocation, 'fasta')]


# In[5]:


#print(sequences)


# In[6]:


sequences_input = sequences[0]


# In[7]:


sequences_input.name


# In[8]:


the_input_string = list(sequences_input.seq) 


# In[9]:


#print(the_input_string)


# In[10]:


#p = 0.005
#k =10


# In[11]:


the_input_string_length = len(the_input_string)
#print(the_input_string_length)


# In[12]:


length_list_mut_del = int (the_input_string_length*p)
#print(length_list_mut_del)
length_list_mut = int ((length_list_mut_del*4)/5)
length_list_del = length_list_mut_del - length_list_mut
#print(length_list_mut)
#print(length_list_del)


# In[13]:


list_mut_del = random.sample(range(0, the_input_string_length), length_list_mut_del)
list_mut_del.sort()
#print(list_mut_del)


# In[14]:


#print(len(list_mut_del))


# In[15]:


listA = ['G','C','T']
listG = ['A','C','T']
listC = ['A','G','T']
listT = ['A','G','C']


# In[16]:


def mutation_process(input_string):
    
    cnt = 0
    #print(len(input_string))
    for x in list_mut_del:
        if cnt < 4 :
            cnt += 1
            temp_str =""
            temp_str = input_string[x]
            mut_rand= random.randint(0, 2)
            if temp_str == 'A':
                temp_str = listA[mut_rand]
            elif temp_str == 'G':
                temp_str = listG[mut_rand]
            elif temp_str == 'C':
                temp_str = listC[mut_rand]
            else:
                temp_str = listT[mut_rand]
        
        #print(input_string[x])    
            input_string[x] = temp_str
       # print(input_string[x]) 
       # print("\n")
        else:
            cnt = 0
    #print(len(input_string))
 


# In[17]:


def del_process(input_string):
    
    #print(len(input_string))
   # print(length_list_del)
   # print(len(input_string) - length_list_del)
    cnt = 0
    cnt_del = 0
    for x in list_mut_del:
  
        if cnt < 4 :
            cnt += 1
            #print(cnt)
        else :
            del_pos = x - (cnt_del)
         #   #print(del_pos)
            cnt_del +=1
            del input_string[del_pos]
            cnt = 0

    #print(len(input_string))


# In[18]:


for seq_num in range(k-1):
    new_sequence = copy.deepcopy(sequences_input)
    new_input_string = list(new_sequence.seq)
   
    mutation_process(new_input_string)
    del_process(new_input_string)
    
   
    new_input_string=''.join(new_input_string)
    new_sequence.seq = Seq(new_input_string)
  
    seq_num = str(seq_num+1)
    new_sequence.description = 'Nucleotide mutation ' + seq_num
    #print(new_sequence.description)
    
    #print(sequences_input.seq)
    #print(new_sequence.seq)
    #print('end')
    sequences.append(new_sequence)
    
    


# In[19]:


SeqIO.write(sequences,outputFileLocation, "fasta")


# In[ ]:




