# R-code for confidence intervals for hazards.
# Competing risk model with current status data
# Piet Groeneboom
# 3-26-2014

   A<-read.table("MLE.txt")
   B<-read.table("SMLE.txt")
   C<-read.table("CI_SMLE.txt")
   x1<-A[,1]
   y1<-A[,2]
   y2<-A[,3]
   y3<-A[,4]
   x2<-B[,1]
   z1<-B[,2]
   z2<-B[,3]
   z3<-B[,4]
   x3<-C[,1]
   t1<-C[,2]
   u1<-C[,3]
   t2<-C[,4]
   u2<-C[,5]
   t3<-C[,6]
   u3<-C[,7]
   
   
   plot(c(-100,-100),xlim=c(min(x2),max(x2)), ylim=c(0,max(u1)), main= "", ylab="",xlab="",bty="n",las=1)
   lines(x1, y1,lwd=2,lty=2,type ="s")
   lines(x2, z1,lwd=2,col="red")
   segments(x3,t1,x3,u1)
   segments(x3-0.05,t1,x3+0.05,t1)
   segments(x3-0.05,u1,x3+0.05,u1)

   plot(c(-100,-100),xlim=c(min(x2),max(x2)), ylim=c(0,max(u2)), main= "", ylab="",xlab="",bty="n",las=1)
   lines(x1, y2,lwd=2,lty=2,type ="s")
   lines(x2, z2,lwd=2,col="red")
   segments(x3,t2,x3,u2)
   segments(x3-0.05,t2,x3+0.05,t2)
   segments(x3-0.05,u2,x3+0.05,u2)
   
   
   plot(c(-100,-100),xlim=c(min(x2),max(x2)), ylim=c(0,0.1), main= "", ylab="",xlab="",bty="n",las=1)
   lines(x1, y3,lwd=2,lty=2,type ="s")
   lines(x2, z3,lwd=2,col="red")
   segments(x3,t3,x3,u3)
   segments(x3-0.05,t3,x3+0.05,t3)
   segments(x3-0.05,u3,x3+0.05,u3)


  