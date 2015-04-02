# The present script has as input the (ordered) observation times and
# the frequencies of the observations for the different risks. 
# It expects K+2 columns: the first column
# contains the observation times, the next column the frequencies
# for the observations where no failure is observed, followed by the
# columns with the frequencies for the observed failures of type k, k=1,...,K.
# One has to specify the number of risks K.


	library(Rcpp) 
	A<-as.matrix(read.table("dataThai2.txt"))
	T<-A[,1]
	K<-3
	frequencies <- A[,2:(K+2)]
	icm<-sourceCpp("icm2.cpp")
	output <- ComputeMLE(T,frequencies,K)
	
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
   	 
   save.image("dataThai2.RData")
		
