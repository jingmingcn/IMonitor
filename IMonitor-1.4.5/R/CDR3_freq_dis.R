args <- commandArgs(TRUE)
infile1<-args[1]
infile2<-args[2]
infile3<-args[3]
out<-args[4]
name<-args[5]

d1<-read.table(infile1,header = F,fill = T)
d2<-read.table(infile2,header = F,fill = T)
d3<-read.table(infile3,header = F,fill = T)
#x1_range<-range(d1[,1])
y12_range<-range(log10(d1[,2]))
y<-trunc(y12_range[2]+1)


pdf(out,height = 8 , width=6)
layout(matrix(c(1,1,2,3),2,2,byrow=TRUE))
	par(mar=c(4,5,4,4))
	D<-matrix(c(1:nrow(d1),rep(0,nrow(d1))),nrow=nrow(d1),ncol=2)
	D[,2]<-sort(log10(d1[,2]))
	plot(D[,1],D[,2],bty="u",axes=F,ylim=c(0,y),xlab="CDR3 Sequence",ylab="",main=paste("Abundance Distribution(",name,")"),cex.lab=1.4,cex.main=1.5,cex.axis=1.3,xlim=c(1,nrow(d1)),type="l")
#	barplot(sort(log10(d1[,2])),bty="u",axes=F,ylim=c(0,y),xlab="CDR3 Sequence",ylab="",main=paste("Abundance Distribution(",name,")"),cex.lab=1.3,cex.main=1.4,cex.axis=1.2)
	x=range(D[,1])
	ymin=min(D[,2])
	polygon(c(x[1],D[,1],x[2]),c(ymin,D[,2],ymin),col="black")
	grid(col="grey")
	abline(h=0)
	abline(h=y)
	y_new<-0;
	for (i in 1:y) {y_new<-c(y_new,10^i)}

	axis(2,at=seq(0,y,1),labels=y_new,las=2)
	mtext(side=2,"CDR3 Abundance",cex.lab=1.3,line=3.5)
	axis(4,at=seq(0,y,1),labels=y_new,las=2)

	rate<-paste(round(d2[,3]),"%",sep="")
	lb<-paste(d2[,1],rate,sep=":")
	pie(d2[,3],labels=lb,main="CDR3-AA",cex=1.3,col=rainbow(4))

	par(mar=c(5,5,5,1.5))
	barplot(d3[1:10,3],col="light blue",xlab="Clone",ylab="Frequency(%)",main="Top10 of CDR3-AA",cex.lab=1.4,cex.axis=1.3,cex.main=1.4,names=1:10)


dev.off()
