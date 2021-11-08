#Neptune_functionの読み込み
.myfunc.env = new.env()
sys.source('~/自分のフォルダのパス/Neptune_function(d56Fe).R', envir=.myfunc.env)
attach(.myfunc.env)

#フォルダ選択
setwd("生データのフォルダのパス")

#生データの入力（基本的にはskipは固定でOK）
x1 <- read.csv("01-blk.exp",skip=16,nrows=20,sep="",row.names=1) # blk file
x2 <- read.csv("02-IRMM.exp",skip=16,nrows=60,sep="",row.names=1) # IRMM file
x3 <- read.csv("03-blk.exp",skip=16,nrows=20,sep="",row.names=1) # blk file
x4 <- read.csv("04-sample.exp",skip=16,nrows=60,sep="",row.names=1) # sample file
x5 <- read.csv("05-blk.exp",skip=16,nrows=20,sep="",row.names=1) # blk file
x6 <- read.csv("06-IRMM.exp",skip=16,nrows=60,sep="",row.names=1) # IRMM file

#計算
IRMM1  <- Ratio(x1,x2)
Sample <- Ratio(x3,x4)
IRMM2  <- Ratio(x5,x6)
delta_values <- Delta(IRMM1,Sample,IRMM2)
print(delta_values)

#結果の保存
write.table(delta_values,"result.csv",append=T,col.names=NA,row.names=T,sep=",",quote=F)
