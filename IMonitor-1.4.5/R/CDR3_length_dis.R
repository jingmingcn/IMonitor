args <- commandArgs(TRUE)
infile<-args[1]
out<-args[2]
name<-args[3]

rt<-read.table(infile,header=F)
pdf(width=8, height=6,file=out)
y<-rt[,2]
names(y)<-rt[,1]
par(mar=c(6,5,5.5,2.5))
barplot(y,col="light blue", xlab="CDR3(nucleotide) Length(bp)", cex.lab=1.3,cex.main=1.4, ylab="Percentage(%)", main=paste("CDR3 Length(",name,")"),cex.axis=1.2)
dev.off()
