library(ggpubr)

setwd("フォルダのパス")
data<-read.csv("result.csv")#自分で纏めたファイルです。
data<-data[,-1]
data_num<-length(data)/2

#dataを行列からデータフレームに変換 
sp2<-data.frame(data)
#各列に名前を付ける(念のため)
names(sp2)<-c("energy1","試料１","energy2","試料２","energy3","試料３",
              "energy4","試料４","energy5","試料５","energy6","試料６",
              "energy7","試料７","energy8","試料８","energy9","試料９",
              "energy10","試料１０","energy11","試料１１","energy12","試料１２",
              "ebergy13","試料１３","energy14","試料１４","energy15","試料１５")

#dataのデータを各々個別のデータフレームに入れてタグ付け
label <- c()
df<-mapply(rep,1:data_num,0)
for (i in 1:data_num){
  label[i] <- names(sp2[2*i])
  df[[i]] <- data.frame(x=sp2[,2*i-1],y=sp2[,2*i],group=label[i])}

#異常点除去
abnormal<-mapply(rep,1:data_num,0)
for (i in 1:data_num){
  x<-df[[i]][,1]
  y<-df[[i]][,2]
  plot(x,y, main = df[[i]][1,3], xlab = "energy [eV]",type="l")
  print(c("異常点をクリックで入力:図",i))
  abnormal[[i]]<-identify(x,y)}#グラフの点の配列番号を示してくれる(ESCキーで操作を終了する)

for (i in 1:data_num){
  for (j in 1: length(df[[i]][,1])){
      df[[i]][abnormal[[i]],] <- NA}
  df[[i]] <- na.omit(df[[i]])}

#異常点除去後のデータを表示
par(mfrow=c(3,5))
for (i in 1:data_num){
  x<-df[[i]][,1]
  y<-df[[i]][,2]
  plot(x,y, main = df[[i]][1,3], xlab = "energy [eV]",type="l")}
print("異常点除去後のデータ")
readline("Push Enter:")#待機
par(mfrow=c(1,1))

#BGを指定する
range<-rep(1:data_num)
for (i in 1:data_num){
  x<-df[[i]][,1]
  y<-df[[i]][,2]
  plot(x,y, main = df[[i]][1,3], xlab = "energy [eV]",type="l")
  print(c("BG終点をクリックで入力:図",i))
  range[i]<-identify(x,y)}#グラフの点の配列番号を示してくれる(ESCキーで操作を終了する)

#BG補正をかける
par(mfrow=c(3,5))
BG<-mapply(rep,1:data_num,0)
for (i in 1:data_num){
  plot(df[[i]][,1],df[[i]][,2],main = df[[i]][1,3],xlab = "energy [eV]",type="l")
  x<-df[[i]][1:range[i],1]
  y<-df[[i]][1:range[i],2]
  model<-lm(y~x)
  BG[[i]]<-predict(model,df[[i]])
  par(new=T)
  plot(df[[i]][,1],BG[[i]],ylim=c(min(df[[i]][,2]),max(df[[i]][,2])),ann = F,
       main = df[[i]][1,3],xlab = "energy [eV]",type="l",col="red")
  par(new=F)}
print("BGの表示")
readline("Push Enter:")#待機

#規格化位置の選択
par(mfrow=c(1,1))
top<-rep(1:data_num)
for (i in 1:data_num){
  x<-df[[i]][,1]
  y<-df[[i]][,2]
  plot(x,y, main = df[[i]][1,3], xlab = "energy [eV]",type="l")
  print(c("Normalize位置をクリックで入力:図",i))
  top[i]<-identify(x,y)}#グラフの点の配列番号を示してくれる(ESCキーで操作を終了する)


#規格化
norm_df<-df
par(mfrow=c(3,5))
for (i in 1:data_num){
  norm_df[[i]][,1] <- df[[i]][,1]
  norm_df[[i]][,3] <- df[[i]][,3]
  norm_df[[i]][,2] <- (df[[i]][,2]-BG[[i]])*(1/(df[[i]][top[i],2]-min(df[[i]][,2]))) #最大値が１になるようにする
  plot(norm_df[[i]][,1],norm_df[[i]][,2],main = norm_df[[i]][1,3],xlab = "energy [eV]",type = "l")
}
print("規格化完了")
readline("Push Enter:")#待機


#修正したデータを保存
max <- c()
for(i in 1:data_num){
  max<-max(length(norm_df[[i]][,1]))}
summary <- matrix(nrow=max,ncol=2*data_num)

for (k in 1:data_num){
for(i in 1:max){
  if (is.na(norm_df[[k]][i,]) == FALSE){
    summary[i,2*k-1]<-norm_df[[k]][i,1]
    summary[i,2*k]<-norm_df[[k]][i,2]}
  else {summary[i,] <- NA}}}

colnames(summary) <- names(sp2)

write.table(summary,"normalized_spectrum.csv",append=F,col.names=NA,row.names=T,sep=",",quote=F)


#スペクトルを5本ずつ描画
data<-norm_df
for (i in 1:length(data)){
  for (j in 1:length(norm_df[[i]][,1])){
    data[[i]][j,2]<-norm_df[[i]][j,2]+8-i/2}}


ResultPlot<-function(a,b,c,d,e,title){
  ggplot()+
    geom_line(data=data[[a]],aes(x = x, y = y, color=group),alpha = 1, size=2)+
    geom_line(data=data[[b]],aes(x = x, y = y, color=group),alpha = 1, size=2)+
    geom_line(data=data[[c]],aes(x = x, y = y, color=group),alpha = 1, size=2)+
    geom_line(data=data[[d]],aes(x = x, y = y, color=group),alpha = 1, size=2)+
    geom_line(data=data[[e]],aes(x = x, y = y, color=group),alpha = 1, size=2)+
    scale_x_continuous("Energy [eV]",limits = c(7090,7250)) + 
    scale_y_continuous(" ")+
    labs(title=title)+
    scale_color_manual(labels=c("試料名1","試料名2","試料名3","試料名4","試料名5"),
                       breaks=c(data[[a]][1,3],data[[b]][1,3],
                                data[[c]][1,3],data[[d]][1,3],data[[e]][1,3]),
                       values=c("#ff6347",
                                "orange",
                                "limegreen",
                                "#66cdaa",
                                "dodgerblue"))+
    theme_bw() +
    theme(axis.text.x=element_text(size=15),
          axis.text.y=element_blank(),
          axis.title = element_text(size=15),
          plot.margin= unit(c(1, 5, 1, 1), "lines"))+
    theme(legend.position = c(1.03,0.02),
          legend.justification = c(0,0),
          legend.text=element_text(size=15),
          legend.background=element_rect(fill="white",colour = "gray"),
          legend.title = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
}

p1 <- ResultPlot(1,2,3,4,5,"title1")
p2 <- ResultPlot(6,7,8,9,10,"title2")
p3 <- ResultPlot(11,12,13,14,15,"title3")

gridExtra::grid.arrange(p1,p2,p3,nrow=1)

