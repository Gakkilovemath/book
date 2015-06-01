#####
# Figure 2.4 in Nonparametric Estimation under Shape Constraints
# Piet Groeneboom and Geurt Jongbloed 
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

#n<-15
#x<-sort(rexp(n)), can be used to generate own data
x<-c(0.1327366, 0.1387036, 0.1841652, 0.2732703, 0.4060988, 0.5234132, 0.5563788,
0.6316339, 0.6461327, 1.0647460, 1.5012775, 1.6169203, 1.6999337, 4.2833256,
6.0971440) # the data used in the book
n<-length(x)

######
# The empirical cdf with its least concave majorant (a)
######

plot(ecdf(x),verticals=TRUE,do.p=FALSE,main="",xlab="",ylab="",bty="n",las=1,lty=2)
lines(ecdf(x),verticals=FALSE,do.p=FALSE)
y<-vector(length=n+1)
cumw<-x
cs<-seq(from=1/n, to=1,by=1/n)
y<-concavemaj(n,cumw,cs,y)
y[n+1]=0
cumy<-vector(length=n+1)
cumy[1]<-y[1]*x[1]
for (i in 2:(n+1))
 	cumy[i] <- cumy[i-1]+y[i]*(x[i]-x[i-1])
lines(c(0,x,7),c(0,cumy),lwd=2,col="blue")

######
# The resulting left continuous Grenander estimate (b)
######

j<-0
for (i in 1:n)
{
	if (y[i+1]<y[i])
		j<-j+1
}
m<-j
t<-vector(length=m)
u<-vector(length=m)
j<-0
for (i in 1:n)
{
	if (y[i+1]<y[i])
	{
		j<-j+1
		t[j]<-x[i]
		u[j]<-y[i]
	}
}


plot(c(x,7),y,main="",xlab="",ylab="",bty="n",type="n",las=1)
points(t,u,pch=19)
lines(c(0,x,7),c(y[1],y),lwd=2,type = "S",col="blue")




