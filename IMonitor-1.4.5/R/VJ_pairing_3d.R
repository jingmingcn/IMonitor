

#===========================================
#filename="R14_A.plot"
args <- commandArgs(TRUE)
filename<-args[1]
fre<-read.table(filename)
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
#if((sort.l[1]/sort.l[2])>1.5){
#plot.red=1
#}



#=======================================================

#pdf("R14_A_barplot_3d.pdf",width=25,height=10)
out.pdf<-paste(args[2],".pdf", sep="")
pdf(out.pdf,width=25,height=10)
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

title <- paste("V-J pairing usage in", args[3])
#plot(0,xaxs="i",yaxs="i",xlab=NA,ylab=NA,axes=F,xlim=c(0,bs.s*nrow(fre.m)+bs.s*cos.a*ncol(fre.m))*1.2,ylim=c(0,bs.s*(ncol(fre.m)*sin.a+10))*1.1,cex=0,main="R14_A_Barplot3d")
plot(0,xaxs="i",yaxs="i",xlab=NA,ylab=NA,axes=F,xlim=c(0,bs.s*nrow(fre.m)+bs.s*cos.a*ncol(fre.m))*1.2,ylim=c(0,bs.s*(ncol(fre.m)*sin.a+10))*1.1,cex=0)
title(title,line=-1.3,cex.main=2.5)

h.left=(bar.size+space.size)*10
#for(j in 0:ncol(fre.m)){
#x1=j*(bar.size+space.size)*cos.a
#y1=j*(bar.size+space.size)*sin.a
#x2=x1
#y2=y1+h.left
#lines(c(x1,x2),c(y1,y2),col=col.grid)
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


#for(j in 0:nrow(fre.m)){
#x1=ncol(fre.m)*(bar.size+space.size)*cos.a+j*(bar.size+space.size)
#y1=ncol(fre.m)*(bar.size+space.size)*sin.a
#x2=x1
#y2=y1+h.left
#lines(c(x1,x2),c(y1,y2),col=col.grid)
#}



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

#lines(c(x1,x2),c(y1,y2),col=col.line,lwd=lwd.line)
#lines(c(x2,x4),c(y2,y4),col=col.line,lwd=lwd.line)
#lines(c(x1,x1),c(y1,y1+h),col=col.line,lwd=lwd.line)
#lines(c(x2,x2),c(y2,y2+h),col=col.line,lwd=lwd.line)
#lines(c(x4,x4),c(y4,y4+h),col=col.line,lwd=lwd.line)
#lines(c(x1,x2),c(y1+h,y2+h),col=col.line,lwd=lwd.line)
#lines(c(x1,x3),c(y1+h,y3+h),col=col.line,lwd=lwd.line)
#lines(c(x2,x4),c(y2+h,y4+h),col=col.line,lwd=lwd.line)
#lines(c(x3,x4),c(y3+h,y4+h),col=col.line,lwd=lwd.line)



}}


#text(seq((bar.size+space.size),(nrow(fre.m))*(bar.size+space.size),by=(bar.size+space.size)),xlab.pos,labels=rownames(fre.m),xpd=TRUE,font=font.lab,pos=4,offset=0,col=col.lab,cex=cex.lab)

text(seq(0.5*(bar.size+space.size),(nrow(fre.m))*(bar.size+space.size),by=(bar.size+space.size)),rep(-0.2*(bar.size+space.size),nrow(fre.m)),labels=rownames(fre.m),font=font.lab,col=col.lab,cex=cex.lab,xpd=TRUE,srt=90,adj=1)

x=seq(nrow(fre.m)*(bar.size+space.size)+0.5*(bar.size+space.size)*cos.a,nrow(fre.m)*(bar.size+space.size)+(ncol(fre.m)-0.5)*(bar.size+space.size)*cos.a,by=(bar.size+space.size)*cos.a)
y=seq(0.5*(bar.size+space.size)*sin.a,(ncol(fre.m)-0.5)*(bar.size+space.size)*sin.a,by=(bar.size+space.size)*sin.a)
text(x,y,pos=4,colnames(fre.m),font=font.lab,col=col.lab,offset=ylab.pos,cex=cex.lab)
x=nrow(fre.m)*(bar.size+space.size)+ncol(fre.m)*(bar.size+space.size)*cos.a
y=seq(ncol(fre.m)*(bar.size+space.size)*sin.a,ncol(fre.m)*(bar.size+space.size)*sin.a+10*(bar.size+space.size),by=(bar.size+space.size))
text(x,y+0.1*bs.s,pos=4,paste(seq(0,l.min,by=l.min/10),"%"),font=font.lab,col=col.lab,cex=cex.lab)

dev.off()



