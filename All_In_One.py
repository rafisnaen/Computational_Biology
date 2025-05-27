#!/usr/bin/env python
# coding: utf-8

# In[21]:


# Sequence A : ATGCGTACGTTAGCCTAGGCTAACGTAGGCTTACGGTAGCTAGCTGATCGTACGTAGCTAG
# Sequence B : GCTAGCATCGGATACGTTAGGCCATGCGTACCTGGAATCGTACCGTGGATCGATCAGGTAC
"""
•	Find the length of both sequences.
•	Determine the number of times the codon (CGC) appears in both sequences.
•	Locate the first occurrence index position of the sub-sequence (CAGTC) in both sequences.
•	Combine the first 15 nucleotides from sequence A and the last 10 nucleotides from sequence B to create a new DNA sequence called sequence C.
•	Reverse the sequence of sequence C.
•	Plot the nucleotide base (A, C, G, T) frequency of each of the sequence A, sequence B, and sequence C using the Matplotlib library.
"""


# In[22]:


from Bio.Seq import Seq
Seq_A = Seq("ATGCGTACGTTAGCCTAGGCTAACGTAGGCTTACGGTAGCTAGCTGATCGTACGTAGCTAG")
Seq_B = Seq("GCTAGCATCGGATACGTTAGGCCATGCGTACCTGGAATCGTACCGTGGATCGATCAGGTAC")


# In[23]:


# 1
len_A = len(Seq_A)
len_B = len(Seq_B)
print(len_A, len_B)


# In[24]:


#2
print(Seq_A.count("CGC"))
print(Seq_B.count("CGC"))


# In[25]:


#3
print(Seq_A.find("CAGTC"))
print(Seq_B.find("CAGTC"))


# In[26]:


#4
Seq_C = (Seq_A[:15] + Seq_B[-10:])
print(Seq_C)


# In[27]:


#5
Reversed_Seq_C = Seq_C[::-1]
print(Reversed_Seq_C)


# In[28]:


#6 
import matplotlib.pyplot as plt
from collections import Counter

freq_A = Counter(Seq_A)
freq_B = Counter(Seq_B)
freq_C = Counter(Seq_C)

freq_keys_A = freq_A.keys()
freq_keys_B = freq_B.keys()
freq_keys_C = freq_C.keys()

freq_values_A = freq_A.values()
freq_values_B = freq_B.values()
freq_values_C = freq_C.values()

# _,axs = plt.subplots(2,2)
# axs[0,0].bar(freq_keys_A, freq_values_A)
# axs[0,1].bar(freq_keys_B, freq_values_B)
# axs[1,0].bar(freq_keys_C, freq_values_C)
# plt.show()

plt.bar(freq_keys_A, freq_values_A)
plt.title("Sequence A")
plt.show()

plt.bar(freq_keys_B, freq_values_B)
plt.title("Sequence B")
plt.show()

plt.bar(freq_keys_C, freq_values_C)
plt.title("Sequence C")
plt.show()


# •	Determine the molecular weight of both sequence.
# •	Refer to the Wallace rule to calculate the melting point of sequence A and sequence B and plot them using Matplotlib library.
# •	Calculate the GC content of sequence A and sequence B.
# •	Calculate the AT content of sequence A and sequence B.
# •	Plot the GC and AT content of sequence A and sequence B using Matplotlib library.

# In[29]:


from Bio.SeqUtils import molecular_weight
weight_A = molecular_weight(Seq_A)
weight_B = molecular_weight(Seq_B)

print(weight_A, weight_B)


# In[30]:


from Bio.SeqUtils import MeltingTemp
MeltPoint_A = MeltingTemp.Tm_Wallace(Seq_A)
MeltPoint_B = MeltingTemp.Tm_Wallace(Seq_B)
print(MeltPoint_A, MeltPoint_B)

plt.bar(["A","B"], [MeltPoint_A, MeltPoint_B])
plt.xlabel("Sequence")
plt.ylabel("Melting Point")
plt.show()


# In[31]:


from Bio.SeqUtils import GC
GC_A = GC(Seq_A)
GC_B = GC(Seq_B)

AT_A = 100 - GC_A
AT_B = 100 - GC_B

_,axs = plt.subplots(1,2)
axs[0].bar(["A", "B"], [GC_A, GC_B])
axs[0].set_title("GC Content")
axs[1].bar(["A", "B"], [AT_A, AT_B])
axs[1].set_title("AT Content")
plt.show()


# •	Transcribe the DNA sequence A and sequence B into mRNA.
# •	Combine 24 starting nucleotides of mRNA sequence A and 21 ending nucleotides of mRNA sequence B and set it as mRNA sequence C.
# •	Translate the mRNA sequence A, sequence B, and sequence C into amino acids sequences.
# •	Transcribe back the mRNA sequence C into DNA.
# •	Calculate the complementary sequences of DNA sequence A, sequence B, and sequence C.

# In[32]:


mRNA_A = Seq_A.transcribe()
mRNA_B = Seq_B.transcribe()

mRNA_C = Seq_A[:24] + Seq_B[-21:]
print(mRNA_C)


# In[33]:


from Bio.SeqUtils import seq3

translate_A = mRNA_A.translate()
translate_B = mRNA_B.translate()
translate_C = mRNA_C.translate()

Amino_Acid_A = seq3(translate_A)
Amino_Acid_B = seq3(translate_B)
Amino_Acid_C = seq3(translate_C)

print(Amino_Acid_A)
print(Amino_Acid_B)
print(Amino_Acid_C)


# In[34]:


Seq_C = mRNA_C.back_transcribe()
print(Seq_C)


# In[35]:


complement_A = Seq_A.complement()
complement_B = Seq_B.complement()
complement_C = Seq_C.complement()

print(complement_A)
print(complement_B)
print(complement_C)


# •	Determine score of best local alignment.
# •	Calculate similarity of DNA sequence A and sequence B using hamming distance method.
# •	Calculate similarity of DNA sequence A and sequence B using levenshtein distance method
# •	Map DNA sequence A and sequence B similarity using dot plot technique.

# In[36]:


from Bio import pairwise2
from Bio.pairwise2 import format_alignment

local_score = pairwise2.align.localxx(Seq_A, Seq_B, one_alignment_only = True, score_only = True)
print(local_score)

# for i in local_score:
#     print(i)
#     print(format_alignment(*i))


# In[37]:


def hamming(Seq_A, Seq_B):
    result = 0
    for x,y in zip(Seq_A, Seq_B):
        if(x != y):
            result += 1
    return result

print(hamming(Seq_A, Seq_B))


# In[38]:


from Levenshtein import distance
print(distance(Seq_A, Seq_B))


# In[39]:


print(" |"+Seq_A)
for i in Seq_B:
    result = ""
    for j in Seq_A:
        if(j == i):
            result += "X"
        else:
            result += " "
    print(i + "|" + result)


# Load fasta/genbank file
# - Take the first half of sequence, and the last half of sequence
# - Calculate melting point using NN
# - Search the first index of "GCAT", and how many is "GCAT" Appears
# - Calculate GC Content
# - Calculate global alignment with score with match = 2 mismatch = -1
# - plot the melt point and GC content

# In[46]:


from Bio import SeqIO
extracted = SeqIO.read("./session_5/sequence.fasta", "fasta")
print(extracted.seq)

Seq_X = extracted.seq
length = round(len(Seq_X) / 2)
Seq_D = Seq_X[:length]
Seq_E = Seq_X[length:]


# In[ ]:


MeltPoint_D = MeltingTemp.Tm_NN(Seq_D)
MeltPoint_E = MeltingTemp.Tm_NN(Seq_E)
print(MeltPoint_D, MeltPoint_E)


# In[56]:


print(Seq_D.find("GCAT"))
print(Seq_E.find("GCAT"))

print(Seq_D.count("GCAT"))
print(Seq_E.count("GCAT"))


# In[51]:


GC_D = GC(Seq_D)
GC_E = GC(Seq_E)
print(GC_D, GC_E)


# In[53]:


global_score = pairwise2.align.globalmx(Seq_D, Seq_E, 2, -1, one_alignment_only = True, score_only = True)
print(global_score)

# for i in global_score:
#     print(i)
#     print(format_alignment(*i))


# In[55]:


plt.bar(["D", "E"], [GC_D, GC_E])
plt.xlabel("Sequence")
plt.ylabel("GC Composition")
plt.title("GC Content")
plt.show()

plt.bar(["D", "E"], [MeltPoint_D, MeltPoint_E])
plt.title("Melting Temperature")
plt.xlabel("Sequence")
plt.ylabel("Melt Temp")
plt.show()


# In[ ]:




