#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio import SeqIO
from random import choice
from Bio.Seq import Seq
import copy
import random
import sys


# In[ ]:


inputFileLocation = sys.argv[1]
x = int(sys.argv[2])
z = int(sys.argv[3])
y = int(sys.argv[4])
outputFileLocation = sys.argv[5]


# In[2]:


#filename = 'C:/Users/jeeta/OneDrive/Desktop/output1.fasta'


# In[3]:


sequences = [i for i in SeqIO.parse(inputFileLocation, 'fasta')]
chopped_sequences = copy.deepcopy(sequences)
chopped_sequences.clear()


# In[4]:


length_sequences = len(sequences)
#print(length_sequences)


# In[5]:


#x = 100
#y = 200
#z = 220


# In[6]:


sequences_input = sequences[0]


# In[7]:


the_input_string = list(sequences_input.seq)
#print(the_input_string)


# In[8]:


def chopping_function(new_input_string):
    #print(new_input_string)
    while len(new_input_string) > x:
        len_chopping = random.randint(x, y)
        #print('x')
        #print(x)
        #print(len_chopping)
        chopped_list = new_input_string[0:len_chopping]
        #print(chopped_list)
        #print(new_input_string)
        del new_input_string[0:len_chopping]
        if len_chopping <= z:
            chopped_seq_list = copy.deepcopy(sequences_input)
            chopped_list = ''.join(chopped_list)
            chopped_seq_list.seq = Seq(chopped_list)
            chopped_sequences.append(chopped_seq_list)


# In[9]:


for ind in range(length_sequences-1) :
    new_sequence = copy.deepcopy(sequences[ind])
    new_input_string = list(new_sequence.seq)
    chopping_function(new_input_string)


# In[10]:


SeqIO.write(chopped_sequences, outputFileLocation , "fasta")

