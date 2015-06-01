#####
# Figure 2.8 in Nonparametric Estimation under Shape Constraints
# Piet Groeneboom and Geurt Jongbloed 
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


#####
# Cumulative sum diagram current status data (a)
####

#n<-200
#t<-10*rexp(n) # simulated data for illustration
#d<-(5*rexp(n)<=t) # simulated data for illustration
#to<-order(t)
#ts<-t[to]
#ds<-d[to]

t<-scan("RubellaTime.txt")
d<-scan("RubellaDelta.txt")
to<-order(t)
ts<-t[to]
ds<-d[to]
n<-length(ds)

cumw<-1:n
cs<-cumsum(ds)
y<-vector(length=n)
y<-convexmin(n,cumw,cs,y)
cumy<-cumsum(y)

plot(0:n,c(0,cs),main="",xlab="",ylab="",bty="n",type="n",las=1)
plot(0:n,c(0,cs),main="",xlab="",ylab="",bty="n",pch=20,las=1)
lines(0:n,c(0,cumy),lwd=2,col="blue")

#####
# MLE for current status (b) 
####

 j<-0
 if (y[1]>0)
	j<-1	
	
 for (i in 2:n)
 {
	if (y[i]>y[i-1])
		j<-j+1
 }
 m<-j

 u<-vector(length=m)
 v<-vector(length=m)

 j<-0
 if (y[1]>0)
 {
	j<-1
	u[j]<-ts[1]
	v[j]<-y[1]
 }


 for (i in 2:n)
 {
	if (y[i]>y[i-1])
	{
		j<-j+1
		u[j]<-ts[i]
		v[j]<-y[i]
	}
 }


plot(ts,y,main="",xlab="",ylab="",bty="n",type="n",las=1)
points(c(0,u),c(0,v),pch=19)
lines(c(0,u,ts[n]),c(0,v,v[m]),lwd=2,type = "s",col="blue")








