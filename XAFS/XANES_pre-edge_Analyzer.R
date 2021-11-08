library(ggpubr)

#ここにディレクトリ名、ファイル名を入力
setwd("フォルダのパス")
filename<-"ファイル名.xan"#現状REXのXANファイルを読み込む形になっています。

#xanファイルの入力
#pre-edge解析の時はEX_BEGINの列（normalize前のデータ）を抽出する(以下のコードで自動でやってくれる)
p<-read.csv(filename)
n<-1
m<-1
for (i in 1:length(p[,1])){
  if (p[i,1]=="[EX_BEGIN]"){n<-i}
  else if (p[i,1]=="[EX_END]"){m<-i}
  else{next()}}
x1<-read.csv(filename, skip = n+2 ,nrows= m-n-1,sep="",header=F)

par(mfrow=c(1,1))
plot(x1,type="l")
print("拡大範囲を選択")
c1<-identify(x1[,1],x1[,2])

c1<-x1[c1[1]:c1[2],1:2]#pre_edgeの範囲
plot(c1, type="l")

print("pre-edgeの範囲を指定")
c2<-identify(c1[,1],c1[,2])
c2<-c1[c(1:c2[1], c2[2]:length(c1[,1])),1:2]

#xのデータを各々個別のデータフレームに入れてタグ付け
df1 <- data.frame(x1=c1[,1],y1=c1[,2],group="raw data")
df2 <- data.frame(x2=c2[,1],y2=c2[,2])
Energy<-df1$x1/1000

#regression: スプラインにより求めた回帰値
reg<-smooth.spline(df2$x2,df2$y2,spar = 0.7)#とりあえず平滑化パラメータは自動設定にしておく
regy<-predict(reg,df1$x1)
regression_data<-data.frame(regy, group="spline")

#Subtract: 生データ-回帰値
Subtract<-df1$y1-regy$y
Subtract_data<-data.frame(x=df1$x1 , y=Subtract , group="subtraction")

#Gauss関数でフィッティングを行う
Fit<-function(v){
  Slope<-v[1]
  FWHM<-v[2]
  centroid<-v[3]
  Gauss<-numeric(length(df1$x1))
  
  for (i in 1:length(df1$x1)){
    Gauss[i]<-Slope/FWHM/(sqrt(pi/2))*exp(-2*((Energy[i]-7.105)*1000-centroid)^2/(FWHM^2))}
  
  Sum<-sum((Subtract - Gauss)^2)
  return(Sum)}

#最適化関数optim:デフォルトで目的関数の最小化をするようになっている
#Sumの値が最も小さくなるSlope（尖度）,FWHM（半値幅）,centroid（中心）の値を決定する
Opt<-optim(par=c(0.001, 2 ,5),fn=Fit)#parの中にはパラメータの初期値を入れる。的外れすぎると収束しなくなる

#最適化されたGauss関数
result<-Opt[["par"]]
names(result)<-c("Slope","FWHM","Centroid")
Gauss<-data.frame(x=df1$x1,y=rep(1:length(df1$x1)),group="Gauss")
for (i in 1:length(df1$x1)){
  Gauss[i,2]<-result[1]/result[2]/(sqrt(pi/2))*exp(-2*((Energy[i]-7.105)*1000-result[3])^2/(result[2]^2))}

#結果の図示
ggplot()+
  geom_point(data=df1,aes(x = x1, y = y1, color=group),alpha = 1, size=3.3)+
  geom_line(data=regression_data,aes(x = x, y = y,color=group),alpha = 0.4, size=4)+
  geom_point(data=Subtract_data,aes(x = x, y = y,color=group),alpha = 1, size=3.3)+
  geom_line(data=Gauss,aes(x = x, y = y,color=group),alpha = 0.4, size=4)+
  scale_x_continuous("Energy [eV]") + 
  scale_y_continuous(" ")+
  scale_color_manual(labels=c("raw data","spline","subtraction","Gaussian"),
                     breaks=c("raw data","spline","subtraction","Gauss"),
                     values=c("tomato","#ffb6c1","#6495ed","#40e0d0"))+
  theme_bw() +
  theme(axis.text.y=element_text(size=15),
        axis.text.x=element_text(size=15),
        axis.title = element_text(size=15))+
  theme(legend.position =c(0,1),
        legend.justification = c(0,1),
        legend.text=element_text(size=15),
        legend.background=element_rect(fill="white",colour = "gray"),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank())

par(mfrow=c(1,1))
print(result)
