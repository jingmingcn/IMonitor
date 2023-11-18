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
	barplot(sort(log10(d1[,2])),bty="u",axes=F,ylim=c(0,y),xlab="CDR3 Sequence",ylab="",main=paste("Abundance Distribution(",name,")"),cex.lab=1.3,cex.main=1.4,cex.axis=1.2)
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
	pie(d2[,3],labels=lb,main="CDR3-AA",cex=1.2)

	par(mar=c(5,5,5,1.5))
	barplot(d3[1:10,3],col="light blue",xlab="Clone",ylab="Frequency(%)",main="Top10 of CDR3-AA",cex.lab=1.3,cex.axis=1.2,cex.main=1.4,names=1:10)


dev.off()
