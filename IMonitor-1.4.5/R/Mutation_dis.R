args <- commandArgs(TRUE)
infile<-args[1]
out<-args[2]
name<-args[3]

rt<-read.table(infile,header=F)
pdf(width=8, height=6,file=out)
#y<-as.matrix(rt[,2:3])
#names(y)<-rt[,1]
ymax=0
par(mar=c(6,5,5.5,5))

	for (i in 1:nrow(rt))
	{
		rt[i,3] = rt[i,3]/10;
	}
	ymax = max(c(rt[,2],rt[,3]))
	ymax = trunc(ymax)+2
	y<-as.matrix(rt[,2:3])
	barplot(y,col=c("light blue","grey","blue","black"), xlab="", cex.lab=1.3,cex.main=1.4, ylab="Base Mutation Rate(%)", main=paste("Hyper-Mutation(",name,")"),cex.axis=1.2,beside=TRUE,ylim=c(0,ymax),names=c("Base Mutation","Sequence Mutation"),legend=c("V-gene","D-gene","J-gene","Overall"),space=c(0,1),cex.names=1.3)
	
	if(ymax>10) ymax=10
	axis(side=4,at=seq(0,ymax,1),seq(0,ymax*10,10))
	mtext(side=4,"Sequence Mutation Rate(%)",line=3,cex=1.3,cex.axis=1.2)

	for (j in 1:nrow(rt))
	{
		text(j+0.5,rt[j,2]+0.1,round(rt[j,2],digits=2),cex=1.1)
	}
	for (j in 1:nrow(rt))
	{
		text(j+5.5,rt[j,3]+0.2,round(rt[j,3],digits=2),cex=1.1)
	}
dev.off()
