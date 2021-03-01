#required packages-----
library(randomForest) 
library(ROCR)

#functions---------
BS_percent<-function(X){
  X[,ncol(X)+1]<-X[,3]/X[,4];X[,ncol(X)+1]<-X[,5]/X[,6];X[,ncol(X)+1]<-X[,7]/X[,8];XX<-X[,c(1,11,12,13)];colnames(XX)[1]<-"ID"
  return(XX)
}

getChr5<-function(DF){
  C1to4<-DF[c(grep("AT5G",DF$ID,invert=T)),]
  C5<-DF[c(grep("AT5G",DF$ID)),]
  C<-list(C1to4,C5)
  names(C)<-c("C1to4","C5")
  return(C)
}


labelDF<-function(DATA,cutoff){
  DATA[,1][DATA[,1]>=cutoff]<-0 ;DATA[,1][DATA[,1]< cutoff]<-1
  print(table(DATA[,1]))
  DATA[,1]<-as.factor(DATA[,1])
  return(DATA)
}

omitID<-function(DF){
  DF2<-DF[,c(grep("ID",colnames(DF),invert=T))]
  return(DF2)
}

balance<-function(DATA){
  k<-nrow(DATA[DATA[,1]==1,])
  N_DATA<-DATA[DATA[,1]==0,];NK_DATA<-N_DATA[sample(1:nrow(N_DATA),k),];DATA.sc<-rbind(NK_DATA,DATA[DATA[,1]==1,])
  DATA[,1]<-as.factor(DATA[,1])
  return(DATA.sc)
}


plot_ROC<-function(predicted,groundtruth){
  pred<-prediction(as.numeric(predicted),as.numeric(groundtruth))
  perf <- performance(pred, "tpr", "fpr")
  auct<- performance(pred, "auc");auc<-auct@y.values[[1]]
  plot(perf);text(0.8,0.1,paste("AUC = ",round(auc,digits=3)))
  return(invisible(list(perf=perf,auc=auc)))
}

F_score<-function(predicted,groundtruth){ # input predicted class and ground truth class. output F.score,Precision and Recall
  TBL<-table(predicted,groundtruth) #true_false cross table
  Precision<-TBL[1,1]/sum(TBL[1,])
  Recall<-TBL[1,1]/sum(TBL[,1])
  F.score<-2*Precision*Recall/(Precision+Recall)
  FPC<-data.frame(F.score,Precision,Recall)
  return(FPC)
}

#Train Random Forest-----
#read data
DATA<-read.table('path_to_folder/RF_data_fld.txt',header=T)


#prep the data; divide the data into Chr1-4 and Chr5, and label the data with FLD enrichment > 20 (first column of the 'DATA'.RPM normalized control-FLD) 
C<-getChr5(labelDF(DATA,-20)) 
#randam sampling for negatively labelled data
train=balance(omitID(C$C1to4));CV=balance(omitID(C$C5))
#training
model<-randomForest(-train[,2:ncol(train)],y=train[,1],ntree=1000)

#see accuracy
prediction<-predict(model,-CV[,2:ncol(train)]);table(prediction)
F_score(prediction,CV[,1])


#visualise 
IMP<-data.frame(model$importance)
sortbythis<-c(49,50,51,45,46,47,72,25,32,26,28,30,29,27,35,36,71,48,31,33,34,62,63,64,56,57,58,52,53,54,74,13,20,14,16,18,17,15,23,24,73,55,19,21,22,65,66,67,42,43,44,38,39,40,70,1,8,2,4,6,5,3,11,12,69,41,7,9,10,59,60,61,37,68)
IMP[,ncol(IMP)+1]<-sortbythis
IMPs<-IMP[sortbythis,]
I<-IMPs$MeanDecreaseGini

#Fig.2C
#plot Cumulative density 
bin=max(I)/99
plot(NULL,xlim=c(0,max(I)+bin),ylim=c(0,1),cex=0,type='n',xlab="Importance",ylab="Cumulative frequency") # custamize x
cols = colorRamp(c("#FFFFFF","#6D97BF","#19202E"))
bin=max(I)/99;rect(seq(0,max(I),length=100), 0, seq(0,max(I),length=100)+bin, 1, col=rgb(cols(0:99/99)/255), border=NA)
plot(ecdf(IMPs$MeanDecreaseGini),add=T)


#plot heatmap of features except mRNA,length
seps<-c(0,40,150,250,350,500); 
i<-NROW(I)-2
x<-matrix(I[1:i],nrow=i/3)
y <- t(x[nrow(x):1, ncol(x):1])[ncol(x):1, ]
par(mar=c(0,0,0,0));image(y, col=rgb(cols(0:99/99)/255), xaxt = "n", yaxt = "n")
#plot heatmap of mRNA and length(3 rd and 4 th row). 1st adn second rows are max/min color
par(mar=c(0,0,0,0));image(t(as.matrix(c(I[(NROW(I)-1):NROW(I)],min(I),max(I)))),col =rgb(cols(0:99/99)/255), xaxt = "n", yaxt = "n")


#Extended Data Fig. 5a
#plot ROC using RNAP2 around TTS as the only predictor.
groundtruth<-CV$TTS_fld
par(mar=c(1,1,2,1))
RC<-plot_ROC(CV$TTS_WT_CMA601,groundtruth);mtext("CMA601_TTS",side=3,line=1)
plot(RC$perf@x.values[[1]],RC$perf@y.values[[1]],lwd=3,xlab="",ylab="",type="l");text(0.75,0.1,paste("AUC = ",round(RC$auc,digits=2)),cex=1.2);mtext("CMA601_TTS",side=3,line=0.5)
RC<-plot_ROC(CV$TTS_WT_CMA602,groundtruth);mtext("CMA602_TTS",side=3,line=1)
plot(RC$perf@x.values[[1]],RC$perf@y.values[[1]],lwd=3,xlab="",ylab="",type="l");text(0.75,0.1,paste("AUC = ",round(RC$auc,digits=2)),cex=1.2);mtext("CMA602_TTS",side=3,line=0.5)
RC<-plot_ROC(CV$TTS_WT_CMA603,groundtruth);mtext("CMA603_TTS",side=3,line=1)
plot(RC$perf@x.values[[1]],RC$perf@y.values[[1]],lwd=3,xlab="",ylab="",type="l");text(0.75,0.1,paste("AUC = ",round(RC$auc,digits=2)),cex=1.2);mtext("CMA603_TTS",side=3,line=0.5)

