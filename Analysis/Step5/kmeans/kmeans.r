#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
S=read.table(file=paste("sum_top_bottom_",args[1],"_cluster",sep=""),header=TRUE,sep="\t")
prcomp(S[,-1], center = TRUE,scale. = TRUE)->PCAS
data.frame(PC1=PCAS$x[,1],PC2=PCAS$x[,2],sam=S$sample)->hun
rownames(S)<-S$sample
kmeans(S[,2:ncol(S)], 2, nstart = 20)->tSclusters
as.data.frame(tSclusters$cluster)->Q
rownames(Q)->Q$sample
colnames(Q)<-c("Cluster","sample")
merge(hun,Q,by.x="sam",by.y="sample")->hunQ
tiff(paste("kmeans_",args[1],".tiff",sep=""),units="in",height=5,width=5,res=300)
plot(hunQ$PC1,hunQ$PC2,col=ifelse(hunQ$Cluster==2,"red","blue"),pch=16,xlab="PC1",ylab="PC2")
dev.off()
nrow(S)/2->n
data.frame(Sample=S$sample,Truth=c(rep(0,n),rep(1,n)))->ROC
merge(ROC,hunQ,by.x="Sample",by.y="sam")->ROCp
write.table(file=paste("truth_predicted_ROC_",args[1],".txt",sep=""),ROCp,col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")
