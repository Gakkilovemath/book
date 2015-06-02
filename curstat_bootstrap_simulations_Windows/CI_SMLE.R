# R-code for confidence intervals for hazards.
# Competing risk model with current status data
# Piet Groeneboom
# 3-26-2014

    C<-read.table("MLE.txt")
   	D<-read.table("CI_SMLE.txt")
   	E<-read.table("SMLE.txt")
   	f <- function(x) {(1-exp(-x))/(1-exp(-2))}
   	x0<-seq(0,2,by=0.01)
   	y0<-f(x0)
   	x1<-C[,1]
   	y1<-C[,2]
   	x2<-D[,1]
   	y2<-D[,2]
   	t1<-D[,2]
   	u1<-D[,3]
   	z1<-E[,1]
   	z2<-E[,2]
   	
   	
   	plot(c(-10000,-10000),xlim=c(0,2), ylim=c(0.0,1.0), main= "", ylab="",xlab="",bty="n",las=1)
   	lines(x0, y0,lty=2,lwd=2,col="red")
	lines(x1, y1,lty=1,lwd=2,type="s")
	lines(z1, z2,lty=1,lwd=2,col="blue")
	segments(x2,t1,x2,u1)
   	segments(x2-0.005,t1,x2+0.005,t1)
   	segments(x2-0.005,u1,x2+0.005,u1)
 

  