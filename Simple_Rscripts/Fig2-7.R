#####
# Figure 2.7 in Nonparametric Estimation under Shape Constraints
# Piet Groeneboom and Geurt Jongbloed 
#####

#####
# Visualization of current status data
####

# t<-scan("RubellaTime.txt")
# d<-scan("RubellaDelta.txt")
m<-200
t<-10*rexp(m) # simulated data for illustration
d<-(5*rexp(m)<=t) # simulated data for illustration
tzero<-t[d==0]
tone<-t[d==1]
t0<-sort(tzero)
t1<-sort(tone)
nt0<-length(t0)
nt1<-length(t1)
t1<-t1[nt1:1]
t0<-t0[nt0:1]

plot(c(0,100),c(1,m),type="n",main="",xlab="",ylab="",bty="n",las=1)
for (i in (1:nt0)){
  lines(c(t0[i],100),c(i,i))}
for (i in (1:nt1)){
  lines(c(0,t1[i]),c(i+nt0,i+nt0))}









