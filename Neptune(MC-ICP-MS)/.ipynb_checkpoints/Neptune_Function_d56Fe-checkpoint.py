import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Read(blk_name,sample_name):
    #blk, sampleファイルの読み込み
    blk=pd.read_csv(blk_name+'.exp',encoding="SHIFT_JIS",sep=('\t'),skiprows=16,usecols=[2,3,4,5,6,7,8,9], 
                names=['52Cr','54Fe','56Fe','57Fe','58Fe','60Ni','63Cu','65Cu'])[1:21]
    sample=pd.read_csv(sample_name+'.exp',encoding="SHIFT_JIS",sep=('\t'),skiprows=16,usecols=[2,3,4,5,6,7,8,9], 
                   names=['52Cr','54Fe','56Fe','57Fe','58Fe','60Ni','63Cu','65Cu'])[1:61]
    return(blk, sample)


def Ratio(blk,sample):
    #原子量
    Fe54_mass=53.939615
    Fe56_mass=55.9349375
    Fe57_mass=56.935394
    Fe58_mass=57.9332756
    Cu62_mass=61.92835
    Cu63_mass=62.929601
    Cu64_mass=63.92797
    Cu65_mass=64.927794
    #標準試料IRMMの同位体比(56Fe54Feは56Fe/54Feという意味)
    IRMM56Fe54Fe=15.6985871271586
    IRMM57Fe54Fe=0.362574568288854
    IRMM58Fe54Fe=0.0482103610675039
    IRMM65Cu63Cu=0.44513
    IRMM54Cr52Cr=0.0282256620797479
    
    #外れ値検定用の下限と上限を求める関数: Low_lim, High_lim
    def Low_lim(X):
        min = np.average(X)-2*np.std(X)
        return min
    def High_lim(X):
        max = np.average(X)+2*np.std(X)
        return max
    
    # このままではDataframe型なので行列(numpy配列)へ変換: blk1
    blk1=blk.to_numpy()
    blk1=blk1.astype('float64')
    
    # blk1 の外れ値検定: reject_judge1
    #平均±2sd外のデータを外し、除去後に平均±2sdをもう一度求めるという操作を複数回行う。
    a=20
    global reject_judge12 #外から呼び出せるようグローバル変数にしておく
    reject_judge1=[blk1[:,0],blk1[:,1],blk1[:,2],blk1[:,3],blk1[:,4],blk1[:,5],blk1[:,6],blk1[:,7]]
    MIN, MAX=np.zeros(8), np.zeros(8)
    reject_judge12=[blk1[:,0],blk1[:,0],blk1[:,1],blk1[:,1],blk1[:,2],blk1[:,2],
                blk1[:,3],blk1[:,3],blk1[:,4],blk1[:,4],blk1[:,5],blk1[:,5],blk1[:,6],blk1[:,6],blk1[:,7],blk1[:,7]]
    reject_judge12=pd.DataFrame(reject_judge12).T
    
    #blkの外れ値を計算
    for m in range(8):
        for n in range(5):
            MIN[m]=Low_lim(reject_judge1[m])
            MAX[m]=High_lim(reject_judge1[m])
            for i in range(a):
                if (reject_judge1[m][i]<=MIN[m])or(reject_judge1[m][i]>=MAX[m]):
                    reject_judge1[m][i]=float('nan')
                else: continue
            reject_judge1[m]=[x for x in reject_judge1[m] if np.isnan(x)==False]
            a=len(reject_judge1[m])
            
    #blk の外れ値検定結果の表示: reject_judge12
    for m in range(8):
        for i in range(20):
            if (reject_judge12[2*m][i]<=MIN[m])or(reject_judge12[2*m][i]>=MAX[m]):
                reject_judge12[2*m+1][i]='reject'
            else: reject_judge12[2*m+1][i]='ok'
    reject_judge12.columns=['52Cr','Judge52Cr','54Fe','Judge54Fe','56Fe','Judge56Fe','57Fe','Judge57Fe','58Fe','Judge58Fe','60Ni','Judge60Ni','63Cu','Judge63Cu','65Cu','Judge65Cu']
    
    # blk1の平均と標準偏差: blk1_ave,blk1_2sd
    blk1_ave, blk1_2sd=np.zeros(8), np.zeros(8)
    for i in range(8):
        blk1_ave[i]=np.average(reject_judge1[i])
        blk1_2sd[i]=2*np.std(reject_judge1[i])
    sample1=sample.to_numpy()
    sample1=sample1.astype('float64')
    
    # sample1 = sample1 -(blk1_ave)
    for i in range(60):
        for j in range(8):
            sample1[i,j]=sample1[i,j]-blk1_ave[j]
            
    # 比率の計算: Ratio1
    Ratio1=np.zeros((60,5))
    Ratio1[:,0]=sample1[:,2]/sample1[:,1]
    Ratio1[:,1]=sample1[:,3]/sample1[:,1]
    Ratio1[:,2]=sample1[:,4]/sample1[:,1]
    Ratio1[:,3]=sample1[:,7]/sample1[:,6]
    Ratio1[:,4]=sample1[:,0]/sample1[:,1]
    
    # β値の計算: Beta1
    Beta1=np.zeros((60,4))
    Beta1[:,0]=np.log(Ratio1[:,0]/IRMM56Fe54Fe)/np.log(Fe56_mass/Fe54_mass)
    Beta1[:,1]=np.log(Ratio1[:,1]/IRMM57Fe54Fe)/np.log(Fe57_mass/Fe54_mass)
    Beta1[:,2]=np.log(Ratio1[:,2]/IRMM58Fe54Fe)/np.log(Fe58_mass/Fe54_mass)
    Beta1[:,3]=np.log(Ratio1[:,3]/IRMM65Cu63Cu)/np.log(Cu65_mass/Cu63_mass)
    
    # interference isotopes: Iso1
    Iso1=np.zeros((60,5))
    for i in range(60):
        if sample1[i,0]>0:
            Iso1[i,0]=sample1[i,0]*IRMM54Cr52Cr*(Cu64_mass/Cu62_mass)**Beta1[i,3]
        else: Iso1[i,0]=0
            
    Iso1[:,1]=sample1[:,1]-Iso1[:,0]
    Iso1[:,2]=sample1[:,2]/Iso1[:,1]
    Iso1[:,3]=sample1[:,3]/Iso1[:,1]
    Iso1[:,4]=sample1[:,4]/Iso1[:,1]
    
    # Cu補正: Corre1
    Corre1=np.zeros((60,3))
    Corre1[:,0]=Iso1[:,2]/((Fe56_mass/Fe54_mass)**Beta1[:,3])
    Corre1[:,1]=Iso1[:,3]/((Fe57_mass/Fe54_mass)**Beta1[:,3])
    Corre1[:,2]=Iso1[:,4]/((Fe58_mass/Fe54_mass)**Beta1[:,3])
    
    # 外れ値の処理
    z=60
    global reject_judge22 #外から呼び出せるようグローバル変数にしておく
    reject_judge2=[Corre1[:,0],Corre1[:,1],Corre1[:,2]]
    reject_judge22=[Corre1[:,0],Corre1[:,0],Corre1[:,1],Corre1[:,1],Corre1[:,2],Corre1[:,2]]
    reject_judge22=pd.DataFrame(reject_judge22).T
    MIN, MAX = np.zeros(3), np.zeros(3)
    
   # blkと同様にして外れ値除去を 5回行う(nが試行回数)
    for m in range(3):
        for n in range(5):
            MIN[m]=Low_lim(reject_judge2[m])
            MAX[m]=High_lim(reject_judge2[m])
            for i in range(z):
                if (reject_judge2[m][i]<=MIN[m])or(reject_judge2[m][i]>=MAX[m]):
                    reject_judge2[m][i]=float('nan')
                else: continue
            reject_judge2[m]=[x for x in reject_judge2[m] if np.isnan(x)==False]
            z=len(reject_judge2[m])
        
    #Corre1 の外れ値検定結果の表示: reject_judge22
    for m in range(3):
        for i in range(60):
            if (reject_judge22[2*m][i]<=MIN[m])or(reject_judge22[2*m][i]>=MAX[m]):
                reject_judge22[2*m+1][i]='reject'
            else: reject_judge22[2*m+1][i]='ok'
    reject_judge22.columns=['56Fe/54Fe','Judge56Fe','57Fe/54Fe','Judge57Fe','58Fe/54Fe','Judge58Fe'] 
    
    
    def PlotResult(x,y,z,t):
        plt.subplot(2,2,z)
        plt.subplots_adjust(wspace=0.4, hspace=0.4)
        plt.plot(np.arange(0,t,1),x,color='gray')
        for i in range (t):
            plt.scatter(i,x[i], color='royalblue' if y[i]=='ok' else 'red')
            plt.xlabel('cycle')
            plt.ylabel(x.name)
    
    fig = plt.subplots(figsize=(10,8))
    blk=PlotResult(reject_judge12['56Fe'],reject_judge12['Judge56Fe'],1,20)
    test=PlotResult(reject_judge22['56Fe/54Fe'],reject_judge22['Judge56Fe'],2,60)
    test2=PlotResult(reject_judge22['57Fe/54Fe'],reject_judge22['Judge57Fe'],3,60)
    test3=PlotResult(reject_judge22['58Fe/54Fe'],reject_judge22['Judge58Fe'],4,60)
    
    
    #平均,2SD,2SEを算出
    ave1, SD1, SE1, SE_rate1 = np.zeros(3),np.zeros(3),np.zeros(3),np.zeros(3)
    for i in range(3):
        ave1[i]=np.average(reject_judge2[i])
        SD1[i]=2*np.std(reject_judge2[i])
        SE1[i]=SD1[i]/(len(reject_judge2[i]))**0.5
        SE_rate1[i]=SE1[i]/ave1[i]
    result = np.c_[ave1,SD1,SE1,SE_rate1]
    result = pd.DataFrame(result).T
    result.columns=['56Fe/54Fe','57Fe/54Fe','58Fe/54Fe']
    result.index=['average','2sd','2se','2se%']
    return result


def Delta(X,Y,Z):
    X,Y,Z=X.to_numpy(),Y.to_numpy(),Z.to_numpy()
    IRMM=(X[0]+Z[0])/2
    SD=2*((X[0]-IRMM)**2+(Z[0]-IRMM)**2)**0.5
    delta=(Y[0]/IRMM-1)*1000
    sigma_2=(Y[3]**2+(X[3]**2+Z[3]**2)/4)**0.5*1000
    
    result=np.c_[delta,sigma_2]
    result = pd.DataFrame(result)
    result.index=['δ56Fe','δ57Fe','δ58Fe']
    result.columns=['‰','2sd']
    return(result)