#!/usr/bin/env python
# coding: utf-8

# In[1]:


# データを総合して解析を行う
import os
from Neptune_Function_d56Fe import read, Ratio, Delta #呼び出す関数は.py形式で保存しておく。
import pandas as pd
import numpy as np
get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
os.chdir("自分のフォルダのパス")
print(os.getcwd())


# In[2]:


#blk + standard1
blk1_name=input('blk file.exp:')
std1_name=input('standard file.exp:')

blk1,std1 = read(blk1_name,std1_name)


# In[3]:


#blk + sample
blk2_name=input('blk file.exp:')
sample_name=input('sample file.exp:')

blk2,sample = read(blk2_name,sample_name)


# In[4]:


#blk + standard2
blk3_name=input('blk file.exp:')
std2_name=input('standard file.exp:')

blk3,std2 = read(blk3_name,std2_name)


# In[5]:


IRMM1 = Ratio(blk1,std1)
Sample = Ratio(blk2,sample)
IRMM2 = Ratio(blk3,std2)


# In[9]:


delta_values = Delta(IRMM1,Sample,IRMM2)
delta_values


# In[ ]:




