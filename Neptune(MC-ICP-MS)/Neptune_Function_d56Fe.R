# Ratio: 56Fe/54Feを計算する
Ratio<-function(x1,x2){ 

library(ggpubr)
  
#原子量
Fe54_mass<-53.939615
Fe56_mass<-55.9349375
Fe57_mass<-56.935394
Fe58_mass<-57.9332756
Cu62_mass<-61.92835
Cu63_mass<-62.929601
Cu64_mass<-63.92797
Cu65_mass<-64.927794

#標準試料IRMMの同位体比(56Fe54Feは56Fe/54Feという意味)
IRMM56Fe54Fe<-15.6985871271586
IRMM57Fe54Fe<-0.362574568288854
IRMM58Fe54Fe<-0.0482103610675039
IRMM65Cu63Cu<-0.44513
IRMM54Cr52Cr<-0.0282256620797479

#外れ値検定用の下限と上限を求める関数: Low_lim, High_lim
Low_lim <- function(X){mean(X)-2*sd(X)}
High_lim <- function(X){mean(X)+2*sd(X)}

# 直前のblkデータ: blk1
blk1<-matrix(nrow=20,ncol=8)
for(i in 1:20){
  for(j in 1:8){
    blk1[i,j]<-as.numeric(x1[i,j+1])}}
colnames(blk1)<-c("52Cr","54Fe","56Fe","57Fe","58Fe","60Ni","63Cu","65Cu")
blk56Fe<-blk1[,3]

# blk1 の外れ値検定: judge1
#平均±2sd外のデータを外し、除去後に平均±2sdをもう一度求める操作を複数回行う。
a<-20
judge1<-list(blk1[,1],blk1[,2],blk1[,3],blk1[,4],blk1[,5],blk1[,6],blk1[,7],blk1[,8])
judge1.2<-data.frame(number=rep(1:20),Cr52=blk1[,1],Judge52Cr=1,Fe54=blk1[,2],Judge54Fe=1,
                     Fe56=blk1[,3],Judge56Fe=1,Fe57=blk1[,4],Judge57Fe=1,Fe58=blk1[,5],
                     Judge58Fe=1,Ni60=blk1[,6],Judge60Ni=1,Cu63=blk1[,7],Judge63Cu=1,Cu65=blk1[,8],Judge65Cu=1)
MIN<-list(1,1,1,1,1,1,1,1)
MAX<-list(1,1,1,1,1,1,1,1)

# 外れ値除去を 5回行う(nが試行回数)
for (m in 1:8){
  for (n in 1:5){
    MIN[[m]]<-Low_lim(judge1[[m]])
    MAX[[m]]<-High_lim(judge1[[m]])
    for (i in 1:a){
      if ((judge1[[m]][i] <= MIN[[m]])||(judge1[[m]][i] >= MAX[[m]])){
        judge1[[m]][i]<-NA}
      else{next()}}
    judge1[[m]] <- na.omit(judge1[[m]])
    a<-length(judge1[[m]])}}

for (m in 1:8){  
  for(i in 1:20){
    if((judge1.2[i,2*m] <= MIN[[m]])||(judge1.2[i,2*m] >= MAX[[m]])){
      judge1.2[i,2*m+1]<-"reject"}
    else{judge1.2[i,2*m+1]<-"ok"}}}

# blk1の平均と標準偏差: blk1_ave,blk1_2sd
blk1_ave<-numeric()
blk1_2sd<-numeric()
for(i in 1:8){
  blk1_ave[i]<-mean(judge1[[i]])
  blk1_2sd[i]<-2*sd(judge1[[i]])}



# 試料データ:sample1
sample1<-matrix(nrow=60,ncol=8)
for(i in 1:60){
  for(j in 1:8){sample1[i,j]<-as.numeric(x2[i,j+1])}}

# sample1 <- sample1 -(blk1_ave)
for(i in 1:60){
  for(j in 1:8){
    sample1[i,j]<-sample1[i,j]-blk1_ave[j]}}
colnames(sample1)<-c("52Cr","54Fe","56Fe","57Fe","58Fe","60Ni","63Cu","65Cu")

# 比率の計算: Ratio1
Ratio1<-matrix(nrow=60,ncol=5)
Ratio1[,1]<-sample1[,3]/sample1[,2]
Ratio1[,2]<-sample1[,4]/sample1[,2]
Ratio1[,3]<-sample1[,5]/sample1[,2]
Ratio1[,4]<-sample1[,8]/sample1[,7]
Ratio1[,5]<-sample1[,1]/sample1[,2]
colnames(Ratio1)<-c("56Fe/54Fe","57Fe/54Fe","58Fe/54Fe","65Cu/63Cu","52Cr/54Fe")

# β値の計算: Beta1
Beta1<-matrix(nrow=60,ncol=4)
Beta1[,1]<-log(Ratio1[,1]/IRMM56Fe54Fe)/log(Fe56_mass/Fe54_mass)
Beta1[,2]<-log(Ratio1[,2]/IRMM57Fe54Fe)/log(Fe57_mass/Fe54_mass)
Beta1[,3]<-log(Ratio1[,3]/IRMM58Fe54Fe)/log(Fe58_mass/Fe54_mass)
Beta1[,4]<-log(Ratio1[,4]/IRMM65Cu63Cu)/log(Cu65_mass/Cu63_mass)
colnames(Beta1)<-c("β(56Fe/54Fe)","β(57Fe/54Fe)","β(58Fe/54Fe)","β(65Cu/63Cu)")

# interference isotopes: Iso1
Iso1<-matrix(nrow=60,ncol=5)
for(i in 1:60){
  if (sample1[i,1]>0){
    Iso1[i,1]<-sample1[i,1]*IRMM54Cr52Cr*(Cu64_mass/Cu62_mass)^Beta1[i,4]}
  else {Iso1[i,1]<-0}}
Iso1[,2]<-sample1[,2]-Iso1[,1]
Iso1[,3]<-sample1[,3]/Iso1[,2]
Iso1[,4]<-sample1[,4]/Iso1[,2]
Iso1[,5]<-sample1[,5]/Iso1[,2]
colnames(Iso1)<-c("54Cr","54Fe","56Fe/54Fe","57Fe/54Fe","58Fe/54Fe")

# Cu補正: Corre1
Corre1<-matrix(nrow=60,ncol=3)
Corre1[,1]<-Iso1[,3]/((Fe56_mass/Fe54_mass)^Beta1[,4])
Corre1[,2]<-Iso1[,4]/((Fe57_mass/Fe54_mass)^Beta1[,4])
Corre1[,3]<-Iso1[,5]/((Fe58_mass/Fe54_mass)^Beta1[,4])
colnames(Corre1)<-c("56Fe/54Fe","57Fe/54Fe","58Fe/54Fe")

# 外れ値の処理
z<-60
judge2<-list(Corre1[,1],Corre1[,2],Corre1[,3])
judge2.2<-data.frame(number=rep(1:60),
                            Fe56Fe54=Corre1[,1],Judge56Fe=1,
                            Fe57Fe54=Corre1[,2],Judge57Fe=1,
                            Fe58Fe54=Corre1[,3],Judge58Fe=1)
MIN<-list(1,1,1)
MAX<-list(1,1,1)

# blkと同様にして外れ値除去を 5回行う(nが試行回数)
for (m in 1:3){
  for (n in 1:5){
    MIN[[m]]<-Low_lim(judge2[[m]])
    MAX[[m]]<-High_lim(judge2[[m]])
    for (i in 1:z){
      if ((judge2[[m]][i] <= MIN[[m]])||(judge2[[m]][i] >= MAX[[m]])){
        judge2[[m]][i]<-NA}
      else{next()}}
    judge2[[m]] <- na.omit(judge2[[m]])
    z<-length(judge2[[m]])}}

for (m in 1:3){  
  for(i in 1:60){
    if((judge2.2[i,2*m] <= MIN[[m]])||(judge2.2[i,2*m] >= MAX[[m]])){
      judge2.2[i,2*m+1]<-"reject"}
    else{judge2.2[i,2*m+1]<-"ok"}}}

#結果の描画
PlotResult<-function(data,x,y,name){
  ggplot()+
    geom_line(data=data,aes(x=number,y=data[,x],group=1),color="gray")+
    geom_point(data=data,aes(x=number,y=data[,x],color=data[,y]),size=4)+
    scale_color_manual(values = c("Royalblue", "red"))+
    scale_x_continuous("cycle")+
    scale_y_continuous(name)+
    theme_bw()+
    theme(legend.title = element_blank())
}

Judge1<-PlotResult(judge1.2,6,7,"blk 56Fe CPS")
Judge2<-PlotResult(judge2.2,2,3,"56Fe/54Fe")
Judge3<-PlotResult(judge2.2,4,5,"57Fe/54Fe")
Judge4<-PlotResult(judge2.2,6,7,"58Fe/54Fe")

#判定結果を出力
gridExtra::grid.arrange(Judge1,Judge2,Judge3,Judge4,ncol=2)

#平均,2SD,2SEを算出
ave1<-c(1,1,1)
SD1<-c(1,1,1)
SE1<-c(1,1,1)
for (i in 1:3){
  ave1[i]<-mean(judge2[[i]])
  SD1[i]<-2*sd(judge2[[i]])
  SE1[i]<-SD1[i]/length(judge2[[i]])^0.5}

SE_rate1<-SE1/ave1
result<-rbind(ave1,SD1,SE1,SE_rate1)
rownames(result)<-c("average","2sd","2se","2se%")
colnames(result)<-c("56Fe/54Fe","57Fe/54Fe","58Fe/54Fe")

return(result)}


# Delta: デルタ値を計算する(X: standard1, Y: Sample, Z: standard2)
Delta<-function(X,Y,Z){
  IRMM<-(X[1,]+Z[1,])/2
  SD<-2*sqrt(((X[1,]-IRMM)^2+(Z[1,]-IRMM)^2))
  delta <- (Y[1,]/IRMM-1)*1000
  sigma_2<-(Y[4,]^2+(X[4,]^2+Z[4,]^2)/4)^0.5*1000
  
  result<-rbind(delta,sigma_2)
  colnames(result)<-c("d56Fe","d57Fe","d58Fe")
  rownames(result)<-c("‰","2σ")
  result<-data.frame(result)
  return(result)}




