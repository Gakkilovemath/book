#####
# Figure 2.1 in Nonparametric Estimation under Shape Constraints
# Piet Groeneboom and Geurt Jongbloed

# The function "convexmin" uses the classical pool adjacent violators
# algorithm for computing the isotonic estimate and has been around
# since the appearance of the Book by Barlow et al. (the 4 B's or B^4),
# see Barlow et al. (1972) Statistical Inference under Order Restrictions.
# Wiley Series in Probability and Mathematical Statistics. 
#####



convexmin<-function(m,cumw,cs,y)
{
 	y[1] <- cs[1]/cumw[1]
 		
 	for (i in 2:m)
 	{
 		y[i] <- (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1])
 		if (y[i-1]>y[i])
		{
			j<-i
			while (y[j-1] > y[i] && j>1)
			{
				j<-(j-1)
				if (j>1)
					y[i] <- (cs[i]-cs[j-1])/(cumw[i]-cumw[j-1])
				else
					y[i] <- cs[i]/cumw[i]
				for (n in j:(i-1))	y[n] <- y[i]
			}
		}
 	}
    return(y)
}

m<-7
x<-c(1.59,1.62,1.59,1.64,1.65,1.68,1.70)
cumw <- c(1:7)
x<-x-mean(x)
cs<-cumsum(x)
y<-vector(length=m)
y<-convexmin(m,cumw,cs,y)
cumy<-cumsum(y)

plot(0:m,c(0,cs),main="",xlab="",ylab="",bty="n",type="n",las=1)
points(0:m,c(0,cs),pch=19)
lines(0:m,c(0,cumy),col="blue")
