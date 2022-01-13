#!/usr/bin/env python
# coding: utf-8

# In[1]:


# データを総合して解析を行う
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def XANESPlot(file, name, color):
    spectrum=[[],[]]
    n,m = 1, 1
    for k in range(len(file)):
        p = pd.read_csv(file[k])
        for i in range(len(p)):
            if p[p.columns[0]][i]=="[BG_BEGIN]":
                n = i
            elif p[p.columns[0]][i]=="[BG_END]":
                m = i
            else: continue
        data = pd.read_csv(file[k],encoding="SHIFT_JIS",sep=('\t'),usecols=[1,2,3,4],names=[1,2,3,4],skiprows=n+3)[1:m-n]
        spectrum[0].append(data[1])
        spectrum[1].append(data[2])
    
    fig = plt.figure(figsize=(8,10))
    ax = fig.add_subplot(111)
    ax.set_xlabel("energy [eV]",fontsize=16)
    ax.set_ylabel(" ")
    ax.axes.yaxis.set_ticklabels([])
    ax.grid()
    ax.set_xlim(7090, 7200)

    for n in range(len(spectrum[1])):
        energy = [float(i) for i in spectrum[0][n]]
        line = spectrum[1][n] + 5 - n/2 
        ax.plot(energy,line,color = color[n],linewidth = 3,label = name[n])
    ax.legend(loc="lower right")






