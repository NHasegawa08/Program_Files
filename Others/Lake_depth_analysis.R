setwd("フォルダのパス")
rawdata<-read.csv("saiko-dem-2016_10m.csv")#ファイル名

data<-matrix(nrow=length(rawdata[,1]))
for (i in 1:length(rawdata[,1])){
  if(rawdata[i,3] > 0){data[i]<-rawdata[i,3]}
  else{data[i]<-NA}}
data<-na.omit(data)
data<-data[,1]

len<-length(data)

data2<-matrix(nrow=len,ncol=2)
data2[,2]<-data

for(j in 1:floor(max(data))+1){
  for (i in 1:len){
    if((data[i] >= j-1) && (data[i] < j)){data2[i,1] <- j}
    else {next()}}
  for(i in 1:len){#なぜかここを付け足さないと下限の値がNAになる。
    if(is.na(data2[i,1])==TRUE){data2[i,1]<-1}
    else{next()}
  }}

summary<-matrix()

for(j in 1:max(data2[,1])){
    summary[j]<-(sum(data2[,1] == j))}
names(summary)<-rep(1:(floor(max(data))+1))
summary

