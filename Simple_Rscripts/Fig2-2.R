#####
# Figure 2.2 in Nonparametric Estimation under Shape Constraints
# Piet Groeneboom and Geurt Jongbloed
# The input file LakeMonona.txt is used. It was taken from
# http://www.aos.wisc.edu/~sco/lakes/Monona-ice.html


# The function "concavemaj" uses the classical pool adjacent violators
# algorithm, just like the function "convexmin" for Figure 2.1 en is
# derived from that algorithm by reversing inequality signs
#####


concavemaj<-function(m,cumw,cs,y)
{
 	y[1] <- cs[1]/cumw[1]
 		
 	for (i in 2:m)
 	{
 		y[i] <- (cs[i]-cs[i-1])/(cumw[i]-cumw[i-1])
 		if (y[i-1]<y[i])
		{
			j<-i
			while (y[j-1] < y[i] && j>1)
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

x1<-scan("LakeMonona.txt")
x<-x1-mean(x1)
m<-length(x)
cs<-cumsum(x)
cumw<-c(1:m)

y<-vector(length=m)
y<-concavemaj(m,cumw,cs,y)
cumy<-cumsum(y)

z<-1855+c(1:m)
y<-y+mean(x1)
plot(z,x1,main="",xlab="",ylab="",bty="n",type="n",las=1)
points(z,x1)
lines(z,y,lwd=2,type="s",col="blue")

plot(0:m,c(0,cs),main="",xlab="",ylab="",bty="n",type="n",las=1)
points(0:m,c(0,cs),pch=19)
lines(0:m,c(0,cumy),lwd=2,col="blue")
