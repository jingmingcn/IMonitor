args <- commandArgs(TRUE)
infile<-args[1]
out<-args[2]
name<-args[3]

rt <- read.table(infile,header=F)
pdf(width=8, height=6,file=out)
par(mar=c(6,5,5.5,2.5))

barplot(rt[,3], names.arg=rt[,1], xlab="Length (bp)", ylab="Frequency(%)", main=paste("Insert-size length(",name,")",sep=""), col="light blue", cex.lab=1.3,cex.main=1.4,)
dev.off()

