setwd("フォルダのパス")

x<-read.csv("ファイル名.dat",skip=16,sep="")
name1<-"Energy"
name2<-"Sample"
Energy<-12398.52/2/3.13551/sin(x[,2]/180*pi)
Trans<- (-1)*log(x[,5]/x[,4])
ut<-x[,6]/x[,4]
plot(Energy,ut)
Map<- data.frame(Energy=Energy, ut=ut)
names(Map)<-c(name1,name2)

#1回目の時
result<-data.frame(Map)

#2回目以降ほかのデータとmergeするとき
result<-data.frame(result,Map)

#結果をエクセルにまとめて保存
write.table(result,"result.csv",append=T,col.names=NA,row.names=T,sep=",",quote=F)
