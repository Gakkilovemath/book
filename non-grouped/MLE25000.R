# The present script expects two columns: the first column
# contains the (ordered) observation times and the second column the risk 1,...,K
# for which failure is observed or 0, if there is no failure.

	library(Rcpp) 
	A<-read.table("data25000.txt")
	icm<-sourceCpp("icm.cpp")
	
	output <- ComputeMLE(A)
	output
	x1 <- output$MLE[,1:4]
	x2 <- output$SMLE[,1:4]
   	
   	 plot(c(-100,-100),xlim=c(min(x2[,1]),max(x2[,1])), ylim=c(0,max(x1[,2:4],x2[,2:4])),main= "",ylab="",xlab="",bty="n",las=1)
   	lines(x1[,1], x1[,2],type="s",col="blue")
   	lines(x1[,1], x1[,3],type="s",lty=2,col="red")
   	lines(x1[,1], x1[,4],type="s",lty=3)
   	f1 <- function(x) {(1/3)*(1-exp(-x))}
   	f2 <- function(x) {(2/3)*(1-exp(-2*x))}
   	f <- function(x) {f1(x)+f2(x)}
   	z <- x2[,1]
   	lines(z, f1(z))
   	lines(z, f2(z))
   	lines(z, f(z))
   	
  
    save.image("data25000.RData")
		
