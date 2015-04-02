# The present script expects two columns: the first column
# contains the (ordered) observation times and the second column the risk 1,...,K
# for which failure is observed or 0, if there is no failure.


	library(Rcpp)
	A<-read.table("dataThai.txt") 
	icm<-sourceCpp("icm.cpp")
	output <- ComputeMLE(A)
	
	x1 <- output$MLE[,1:4]
	x2 <- output$SMLE[,1:4]
  
    output$MLE
		
   	plot(c(-100,-100),xlim=c(min(x2[,1]),max(x2[,1])), ylim=c(0,max(x1[,2:4],x2[,2:4])),main= "",ylab="",xlab="",bty="n",las=1)
   	lines(x1[,1], x1[,2],type="s",col="blue")
   	lines(x1[,1], x1[,3],type="s",lty=2,col="red")
   	lines(x1[,1], x1[,4],type="s",lty=3)
   	lines(x2[,1], x2[,2])
   	lines(x2[,1], x2[,3])
   	lines(x2[,1], x2[,4])
   	 
   save.image("dataThai.RData")
		