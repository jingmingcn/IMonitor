args <- commandArgs(TRUE)
infile<-args[1]
out<-args[2]
name<-args[3]

rt<-read.table(infile,header=T)
pdf(width=8, height=6,file=out)
y<-rt[1:nrow(rt),1:26]
y<-as.matrix(y)

par(mar=c(6,5,5.5,2.5))
if(nrow(rt)==4) {	
	barplot(y, xlab="Length(bp)", cex.lab=1.3,cex.main=1.4,  ylab="Percentage(%)", main=paste("V(D)J Gene Deletion Length(",name,")"),cex.axis=1.2,beside=TRUE,col=c("blue","green","yellow","grey"),legend=rownames(rt),names=0:25)
}else {
	barplot(y, xlab="Length(bp)", cex.lab=1.3,cex.main=1.4,  ylab="Percentage(%)", main=paste("V(D)J Gene Deletion Length(",name,")"),cex.axis=1.2,beside=TRUE,col=c("blue","green"),legend=rownames(rt),names=0:25)
}

dev.off()
