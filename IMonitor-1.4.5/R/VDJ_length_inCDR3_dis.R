args <- commandArgs(TRUE)
infile<-args[1]
out<-args[2]
name<-args[3]

rt<-read.table(infile,header=T)
pdf(width=8, height=6,file=out)
y<-rt[1:nrow(rt),]
y<-as.matrix(y)

par(mar=c(6,5,5.5,2.5))
if(nrow(rt)==3) {	
	barplot(y, xlab="Length(bp)", cex.lab=1.3,cex.main=1.4,  ylab="Percentage(%)", main=paste("V(D)J Gene Length In CDR3(",name,")"),cex.axis=1.2,beside=TRUE,col=c("blue","green","yellow"),legend=rownames(rt),names=0:30)
}else {
	barplot(y, xlab="Length(bp)", cex.lab=1.3,cex.main=1.4,  ylab="Percentage(%)", main=paste("V(D)J Gene Length In CDR3(",name,")"),cex.axis=1.2,beside=TRUE,col=c("blue","green"),legend=rownames(rt),names=0:30)
}

dev.off()
