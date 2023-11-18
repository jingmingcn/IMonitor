args <- commandArgs(TRUE)
infile<-args[1]
out<-args[2]
name<-args[3]

data<-read.table(infile,header = F,fill = T)
pdf(out,height = 6 , width=8)

yrange<-range(c(data[,2],data[,4]))
yrange[1]<-yrange[1]*0.8
yrange[2]<-yrange[2]*1.2
xrange<-range(data[,1])
xl<-xrange[1]+(xrange[2]-xrange[1])/2
yl<-yrange[1]+(yrange[2]-yrange[1])/2


par(mar=c(7,7,5,2)+0.1)
	plot(data[,1],data[,2],xlab="",ylab="",axes=F,pch=2,cex=2,main=name,cex.main=2,ylim=yrange)
	lines(spline(data[,1],data[,2]),col="black",cex=1.8,lwd=2)
	par(new=T)
	plot(data[,1],data[,4],xlab="Effective sequence",ylab="The number of Unique CDR3",col="red",cex.lab=1.7,pch=11,cex=2,cex.axis=1.4,ylim=yrange,bty="n")
	lines(spline(data[,1],data[,4]),col="red",cex=1.8,lwd=2)
	
	legend(xl,yl,c("Observed","Chao1_correction"),col=c("black","red"),pch=c(2,11),cex=1.5)
dev.off()
