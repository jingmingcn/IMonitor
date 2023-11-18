args <- commandArgs(TRUE)
out<-args[1]
name<-args[2]
#	d1<-read.table(infile1,header = F,fill = T)

pdf(out,height = 21 , width=12, bg="gray92")
#grid(nx=3,ny=6,col = "lightgray")
layout(matrix(c(1,2,3,4,5,6,7,8,9,10,11,12,13,13,13,14,14,14,14,14,14),7,3,byrow=TRUE))


#--------------------------------------------	1 row	----------------------------
# 1.1 Insert size Distribution
rt<-read.table(args[3])
if(nrow(rt)==1 & ncol(rt)==1 & rt[1,1]==0){
	plot(NA,NA,xlab="",ylab="",axes=F,main="Insert-size length",xlim=c(0,4),ylim=c(0,3),cex.main=2)
	text(2,1.5,"NULL",cex=2)
}else{
	par(mar=c(6,5,5.5,2.5))
	barplot(rt[,3], names.arg=rt[,1], xlab="Length (bp)", ylab="Frequency(%)", main="Insert-size length", col="light blue", cex.lab=1.7,cex.main=2,cex.axis=1.4,cex.names=1.4)
}

# 1.2 Rarefraction curve 
data<-read.table(args[4],header = F,fill = T)
if(nrow(data)==1 & ncol(data)==1 & data[1,1]==0){
	plot(NA,NA,xlab="",ylab="",axes=F,main="Saturation Curve",xlim=c(0,4),ylim=c(0,3),cex.main=2)
	text(2,1.5,"NULL",cex=2)
}else{

yrange<-range(c(data[,2],data[,4]))
yrange[1]<-yrange[1]*0.8
yrange[2]<-yrange[2]*1.2
xrange<-range(data[,1])
xl<-xrange[1]+(xrange[2]-xrange[1])/2
yl<-yrange[1]+(yrange[2]-yrange[1])/2


par(mar=c(7,7,5,2)+0.1)
        plot(data[,1],data[,2],xlab="",ylab="",axes=F,pch=2,cex=2,main="Saturation Curve",cex.main=2,ylim=yrange,cex.axis=1.4)
        lines(spline(data[,1],data[,2]),col="black",cex=1.8,lwd=2)
        par(new=T)
        plot(data[,1],data[,4],xlab="Effective sequence",ylab="The number of Unique CDR3",col="red",cex.lab=1.7,pch=11,cex=2,cex.axis=1.4,ylim=yrange,bty="n")
        lines(spline(data[,1],data[,4]),col="red",cex=1.8,lwd=2)

        legend(xl,yl,c("Observed","Chao1_correction"),col=c("black","red"),pch=c(2,11),cex=1.5)
}

# 1.3 CDR3 NT length Distribution
rt<-read.table(args[5],header=F)
	y<-rt[,2]
	names(y)<-rt[,1]
	par(mar=c(6,5,5.5,2.5))
	barplot(y,col="light blue", xlab="CDR3(nucleotide) Length(bp)", cex.lab=1.7,cex.main=2, ylab="Percentage(%)", main="CDR3 Length",cex.axis=1.4)

#------------------------------------------- 2 row	-----------------------
# 2.1 CDR3 AA frequency Distribution
d1<-read.table(args[6],header = F,fill = T)
y12_range<-range(log10(d1[,2]))
y<-trunc(y12_range[2]+1)
        par(mar=c(4,5,4,5))
#        barplot(sort(log10(d1[,2])),bty="u",axes=F,ylim=c(0,y),xlab="CDR3 Sequence",ylab="",main="Abundance Distribution",cex.lab=1.7,cex.main=2,cex.axis=1.4)
	D<-matrix(c(1:nrow(d1),rep(0,nrow(d1))),nrow=nrow(d1),ncol=2)
	D[,2]<-sort(log10(d1[,2]))
	plot(D[,1],D[,2],bty="u",axes=F,ylim=c(0,y),xlab="CDR3 Sequence",ylab="",main="Abundance Distribution",cex.lab=1.4,cex.main=1.5,cex.axis=1.3,xlim=c(1,nrow(d1)),type="l")
	x=range(D[,1])
	ymin=min(D[,2])
	polygon(c(x[1],D[,1],x[2]),c(ymin,D[,2],ymin),col="black")
        grid(col="grey")
        abline(h=0)
        abline(h=y)
        y_new<-0;
        for (i in 1:y) {y_new<-c(y_new,10^i)}

        axis(2,at=seq(0,y,1),labels=y_new,las=2,cex.axis=1.4)
        mtext(side=2,"CDR3 Abundance",cex.lab=1.7,line=3.5)
        axis(4,at=seq(0,y,1),labels=y_new,las=2,cex.axis=1.4)

# 2.2 CDR3 AA fraction Distribution
d2<-read.table(args[7],header = F,fill = T)
        rate<-paste(round(d2[,3]),"%",sep="")
	lb<-paste(d2[,1],rate,sep=":")
	pie(d2[,3],labels=lb,main="CDR3-AA",cex=1.5,cex.main=2,cex.axis=1.5,cex.lab=1.7,col=rainbow(4),border="NA")

# 2.3 CDR3 AA top10 clones
d3<-read.table(args[8],header = F,fill = T)
        par(mar=c(5,5,5,1.5))
	barplot(d3[1:10,3],col="light blue",xlab="Clone",ylab="Frequency(%)",main="Top10 of CDR3-AA",cex.lab=1.7,cex.axis=1.4,cex.main=2,names=1:10,cex.names=1.4)

#---------------------------------------   3 row	-------------------------------
# 3.1 VDJ length in CDR3
rt<-read.table(args[9],header=T)
y<-rt[1:nrow(rt),]
y<-as.matrix(y)

par(mar=c(6,5,5.5,2.5))
if(nrow(rt)==3) {
        barplot(y, xlab="Length(bp)", cex.lab=1.7,cex.main=2,  ylab="Percentage(%)", main="V(D)J Gene Length In CDR3",cex.axis=1.4,beside=TRUE,col=c("blue","green","yellow"),legend=rownames(rt),names=0:30,cex.names=1.4,cex=1.5)
}else {
        barplot(y, xlab="Length(bp)", cex.lab=1.7,cex.main=2,  ylab="Percentage(%)", main="V(D)J Gene Length In CDR3",cex.axis=1.4,beside=TRUE,col=c("blue","green"),legend=rownames(rt),names=0:30,cex.names=1.4,cex=1.5)
}

# 3.2 Deletion length Distribution
rt<-read.table(args[10],header=T)
y<-rt[1:nrow(rt),1:26]
y<-as.matrix(y)

par(mar=c(6,5,5.5,2.5))
if(nrow(rt)==4) {
        barplot(y, xlab="Length(bp)", cex.lab=1.7,cex.main=2,  ylab="Percentage(%)", main="V(D)J Gene Deletion Length",cex.axis=1.4,beside=TRUE,col=c("blue","green","yellow","grey"),legend=rownames(rt),names=0:25,cex.names=1.4,cex=1.5)
}else {
        barplot(y, xlab="Length(bp)", cex.lab=1.7,cex.main=2,  ylab="Percentage(%)", main="V(D)J Gene Deletion Length",cex.axis=1.4,beside=TRUE,col=c("blue","green"),legend=rownames(rt),names=0:25,cex.names=1.4,cex=1.5)
}

# 3.3 Insertion length Distribution
rt<-read.table(args[11],header=T)
	y<-rt[1:nrow(rt),1:26]
	y<-as.matrix(y)

	par(mar=c(6,5,5.5,2.5))
	if(nrow(rt)==2) {
		        barplot(y, xlab="Length(bp)", cex.lab=1.7,cex.main=2,  ylab="Percentage(%)", main="V(D)J Gene Insertion Length",cex.axis=1.4,beside=TRUE,col=c("blue","green"),legend=rownames(rt),names=0:25,cex.names=1.4,cex=1.5)
	}else {
		        barplot(y, xlab="Length(bp)", cex.lab=1.7,cex.main=2,  ylab="Percentage(%)", main="V(D)J Gene Insertion Length",cex.axis=1.4,beside=TRUE,col=c("light blue"),legend=rownames(rt),names=0:25,cex.names=1.4,cex=1.5)
	}

#-------------------------------------------	4 row	----------------------
# 4.1 Hyper - mutation
rt<-read.table(args[12],header=F)
if(nrow(rt)==1 & ncol(rt)==1 & rt[1,1]==0){
	plot(NA,NA,xlab="",ylab="",axes=F,main="Hyper-Mutation",xlim=c(0,4),ylim=c(0,3),cex.main=2)
	text(2,1.5,"NULL",cex=2)
}else{
	ymax=0
	par(mar=c(6,5,5.5,5))

        for (i in 1:nrow(rt))
        {
                rt[i,3] = rt[i,3]/10;
        }
        ymax = max(c(rt[,2],rt[,3]))
        ymax = trunc(ymax)+2
        y<-as.matrix(rt[,2:3])
        barplot(y,col=c("light blue","grey","blue","black"), xlab="", cex.lab=1.7,cex.main=2, ylab="Base Mutation Rate(%)", main="Hyper-Mutation",cex.axis=1.4,beside=TRUE,ylim=c(0,ymax),names=c("Base Mutation","Sequence Mutation"),legend=c("V-gene","D-gene","J-gene","Overall"),space=c(0,1),cex.names=1.6,cex=1.5)

        if(ymax>10) ymax=10
        axis(side=4,at=seq(0,ymax,1),seq(0,ymax*10,10),cex.axis=1.4)
        mtext(side=4,"Sequence Mutation Rate(%)",line=3,cex.lab=1.7)

        for (j in 1:nrow(rt))
        {
                text(j+0.5,rt[j,2]+0.1,round(rt[j,2],digits=2),cex=1.3)
        }
        for (j in 1:nrow(rt))
        {
                text(j+5.5,rt[j,3]+0.2,round(rt[j,3],digits=2),cex=1.3)
        }

}

# 4.2 J gene Usage 
rt<-read.table(args[13],header=F)
	y<-rt[,3]
	names(y)<-rt[,1]
	par(mar=c(6,5,5.5,2.5))
	barplot(y,col="light blue", xlab="", cex.lab=1.7,cex.main=2, ylab="Percentage(%)", main="J Gene Usage",cex.axis=1.4,las=2,cex.names=1.5)

# 4.3 NULL
plot(NA,NA,xlab="",ylab="",axes=F,xlim=c(0,4),ylim=c(0,3))

#--------------------------	5 row  ------------------------------
# V gene Usage
rt<-read.table(args[14],header=F)
	y<-rt[,3]
	names(y)<-rt[,1]
	par(mar=c(6,5,5.5,2.5))
	barplot(y,col="light blue", xlab="", cex.lab=1.7,cex.main=2, ylab="Percentage(%)", main="V Gene Usage",cex.axis=1.4,las=2,cex.names=1.5)

#----------------------------	6 row	--------------------------------
# V-J Pairing

fre<-read.table(args[15])
fre<-fre[,c(1:2,4)]
sort(unique(fre[,2]))->fre.l2
sort(unique(fre[,1]))->fre.l1
fre.m<-matrix(0,nrow=length(fre.l1),ncol=length(fre.l2))
rownames(fre.m)<-fre.l1
colnames(fre.m)<-fre.l2
for(i in 1:nrow(fre)){
 fre.m[fre[i,1],fre[i,2]]=fre.m[fre[i,1],fre[i,2]]+fre[i,3]
}
sort(fre.m,decreasing=T)->sort.l
plot.red=0;


#=======================================================

#pdf("R14_A_barplot_3d.pdf",width=25,height=10)
#out.pdf<-paste(args[2],".pdf", sep="")
#pdf(out.pdf,width=25,height=10)
bar.size=4
space.size=12
bs.s<-bar.size+space.size
anglexy=45
sin.a=sin(anglexy/180*3.1415926)
cos.a=cos(anglexy/180*3.1415926)
col.line="grey"
lwd.line=0.2
col.grid="black";
#col.bg=grey(5/8)
#col.bg="lightcyan"
#col.lab="white"
col.lab="black"
cex.lab=1.7
font.lab=2 # 1 or 2 or 3 or 4
xlab.pos=par("usr")[3]-2*bs.s
ylab.pos=0.73


if(plot.red){
l.f1=as.numeric(substr(sort.l[2],1,1))
l.h=nchar(trunc(sort.l[2]))
l.min=(l.f1+1)*10^(l.h-1)
}else{
l.f1=as.numeric(substr(sort.l[1],1,1))
l.h=nchar(trunc(sort.l[1]))
l.min=(l.f1+1)*10^(l.h-1)
}



par(mar=c(8.5,2,2,0))

title <- "V-J pairing Usage"
#plot(0,xaxs="i",yaxs="i",xlab=NA,ylab=NA,axes=F,xlim=c(0,bs.s*nrow(fre.m)+bs.s*cos.a*ncol(fre.m))*1.2,ylim=c(0,bs.s*(ncol(fre.m)*sin.a+10))*1.1,cex=0,main="R14_A_Barplot3d")
plot(0,xaxs="i",yaxs="i",xlab=NA,ylab=NA,axes=F,xlim=c(0,bs.s*nrow(fre.m)+bs.s*cos.a*ncol(fre.m))*1.2,ylim=c(0,bs.s*(ncol(fre.m)*sin.a+10))*1.1,cex=0)
title(title,line=-1.3,cex.main=2.5)

h.left=(bar.size+space.size)*10
#}
lwd.grid=0.5

lines(c(0,0),c(0,h.left),col=col.grid,lwd=lwd.grid)
lines(c(ncol(fre.m)*(bar.size+space.size)*cos.a,ncol(fre.m)*(bar.size+space.size)*cos.a),c(ncol(fre.m)*(bar.size+space.size)*sin.a,ncol(fre.m)*(bar.size+space.size)*cos.a+h.left),col=col.grid,lwd=lwd.grid)

x0=ncol(fre.m)*(bar.size+space.size)*cos.a+nrow(fre.m)*(bar.size+space.size)
y0=ncol(fre.m)*(bar.size+space.size)*sin.a
lines(c(x0,x0),c(y0,y0+h.left),col=col.grid,lwd=lwd.grid)

polygon(c(0,ncol(fre.m)*(bar.size+space.size)*cos.a,x0,nrow(fre.m)*(bar.size+space.size)),c(0,ncol(fre.m)*(bar.size+space.size)*sin.a,y0,0),col="grey",border=NA)

for(j in 0:10){
x1=0
y1=j*(bar.size+space.size)
x2=x1+ncol(fre.m)*(bar.size+space.size)*cos.a
y2=y1+ncol(fre.m)*(bar.size+space.size)*sin.a
lines(c(x1,x2),c(y1,y2),col="grey",lwd=lwd.grid)
}
for(j in 0:10){
x1=ncol(fre.m)*(bar.size+space.size)*cos.a
y1=j*(bar.size+space.size)+ncol(fre.m)*(bar.size+space.size)*sin.a
x2=x1+(bar.size+space.size)*nrow(fre.m)
y2=y1
lines(c(x1,x2),c(y1,y2),col="grey",lwd=lwd.grid)
}


for(j in ncol(fre.m):1){
#col.bar=rgb(0,0,1,alpha=round((j/ncol(fre.m))*4/10+0.6,2))
col.ceiling="grey"

#col.bar="purple"
col.bar="black"
#col.ceiling="blue"
for(i in 1:nrow(fre.m)){

z_r=(bs.s)/(l.min/10)
h=z_r*fre.m[i,j]
#print(h)
if(plot.red){
if(fre.m[i,j]==sort.l[1]){
h=1.1*h.left
}
}
x1=(i-1)*(bar.size+space.size)+space.size/2+space.size/2*cos.a+(j-1)*(bar.size+space.size)*cos.a
x2=x1+bar.size
x3=x1+cos.a*bar.size
x4=x3+bar.size
y1=((j-1)*(bar.size+space.size)+space.size/2)*sin.a
y2=y1
y3=y1+sin.a*bar.size
y4=y1+sin.a*bar.size


x5=x1-space.size/2-space.size/2*cos.a
x6=x5+space.size+bar.size
x7=x5+(space.size+bar.size)*cos.a
x8=x6+(space.size+bar.size)*cos.a
y5=y1-space.size/2*sin.a
y6=y5
y7=y5+(space.size+bar.size)*sin.a
y8=y7


lines(c(x5,x6),c(y5,y6),col=col.grid,lwd=lwd.line)
lines(c(x5,x7),c(y5,y7),col=col.grid,lwd=lwd.line)
lines(c(x6,x8),c(y6,y8),col=col.grid,lwd=lwd.line)
lines(c(x7,x8),c(y7,y8),col=col.grid,lwd=lwd.line)

if(fre.m[i,j] == 0){
next
}

lines(c(x1,x2),c(y1,y2),col=col.line,lwd=lwd.line)
lines(c(x1,x3),c(y1,y3),col=col.line,lwd=lwd.line)
lines(c(x2,x4),c(y2,y4),col=col.line,lwd=lwd.line)
lines(c(x3,x4),c(y3,y4),col=col.line,lwd=lwd.line)
lines(c(x1,x1),c(y1,y1+h),col=col.line,lwd=lwd.line)
lines(c(x2,x2),c(y2,y2+h),col=col.line,lwd=lwd.line)
lines(c(x3,x3),c(y3,y3+h),col=col.line,lwd=lwd.line)
lines(c(x4,x4),c(y4,y4+h),col=col.line,lwd=lwd.line)
lines(c(x1,x2),c(y1+h,y2+h),col=col.line,lwd=lwd.line)
lines(c(x1,x3),c(y1+h,y3+h),col=col.line,lwd=lwd.line)
lines(c(x2,x4),c(y2+h,y4+h),col=col.line,lwd=lwd.line)
lines(c(x3,x4),c(y3+h,y4+h),col=col.line,lwd=lwd.line)



if(plot.red==1){
if(fre.m[i,j]==sort.l[1]){
text(x1,y1+h+0.1*h.left,fre.m[i,j],col=2,cex=0.6,font=2,pos=4,offset=-0.3)
}
}



rect(x1,y1,x2,y2+h,col=col.bar,border=NA)
polygon(c(x2,x4,x4,x2),c(y2,y4,y4+h,y2+h),col=col.bar,border=NA)
polygon(c(x1,x2,x4,x3),c(y1+h,y2+h,y4+h,y3+h),col="ivory",border=NA)

}}

text(seq(0.5*(bar.size+space.size),(nrow(fre.m))*(bar.size+space.size),by=(bar.size+space.size)),rep(-0.2*(bar.size+space.size),nrow(fre.m)),labels=rownames(fre.m),font=font.lab,col=col.lab,cex=cex.lab,xpd=TRUE,srt=90,adj=1)

x=seq(nrow(fre.m)*(bar.size+space.size)+0.5*(bar.size+space.size)*cos.a,nrow(fre.m)*(bar.size+space.size)+(ncol(fre.m)-0.5)*(bar.size+space.size)*cos.a,by=(bar.size+space.size)*cos.a)
y=seq(0.5*(bar.size+space.size)*sin.a,(ncol(fre.m)-0.5)*(bar.size+space.size)*sin.a,by=(bar.size+space.size)*sin.a)
text(x,y,pos=4,colnames(fre.m),font=font.lab,col=col.lab,offset=ylab.pos,cex=cex.lab)
x=nrow(fre.m)*(bar.size+space.size)+ncol(fre.m)*(bar.size+space.size)*cos.a
y=seq(ncol(fre.m)*(bar.size+space.size)*sin.a,ncol(fre.m)*(bar.size+space.size)*sin.a+10*(bar.size+space.size),by=(bar.size+space.size))
text(x,y+0.1*bs.s,pos=4,paste(seq(0,l.min,by=l.min/10),"%"),font=font.lab,col=col.lab,cex=cex.lab)


