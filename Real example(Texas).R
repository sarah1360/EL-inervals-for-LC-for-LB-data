library(MASS)
library(summarytools)
library(tseries)
library(boot)
library(invgamma);
#--------------------------------drawing lorenz curve
shape=1; rate=1;
t0=seq(0.05,0.95,0.05)
rho=function(t) pinvgamma(qinvgamma(t, shape+1, rate), shape, rate)
#plot(c(0,1),c(0,1), type="n", xlab = "Proportion of population", ylab = "Income population", xlim=c(0,1), ylim=c(0,1))
#abline(a=0, b=1, col="red")
#lines(t0, rho(t0), type = "l", col="blue")
#text(0.6,0.2, "Lorenz Curve", col = "blue", adj = c(0, -.1))
#text(0.4, 0.4, "Perfect equality", col = "red", adj = c(-0.2, -.15))
#--------------------------------data
data<-read.table("c:\\Users\\Nikam\\Documents\\SARAH\\My works\\Vejdani\\Data\\Texas_County.txt",header=F)
y=data[[1]]
n=length(y)
#--------------------------------descriptive statistics
#descr(y,style== st_options("grid"))
#jarque.bera.test(y)
#truehist(y,nbins=50,col=27, xlab="Total personal income", ylab="Density", sub="(a)")
#lines(density(y,bw=8,n=520))

#qqnorm(y,sub="(b)")
#qqline(y,col=2)

#boxplot(y, xlab="whole data", ylab = "Time", col="violet")
#--------------------------------Entries
B=400
yb=matrix(nrow=B, ncol=n)
psib=rho_nb=c(); psi=rho_n=c();
alpha=0.05; 
z=qnorm(1-alpha/2);
q=qchisq(p=alpha, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE);

CI=function(t){
#--------------------------------Computing 
y01<-sort(y)
yy_na=1/y01
sum=sum(yy_na)

if(t<=yy_na[1]/sum){psi=y01[1]}  else{kn=2
                                 while(sum(yy_na[1:kn])/sum <t){kn=kn+1}
                                 psi=y01[kn]}

rho_n=mean(1*(y <=psi))
rho_tilde=mean(psi*(t-1*(y<=psi))/y) + rho_n
sigma2=rho_n*(1-rho_n)
sigma1=max(0,(psi^2)*(((1-t)^2)*mean((1/(y^2))*(y<=psi))+(t^2)*mean((1/(y^2))*(y>psi)))+sigma2+2*psi*t*(t-1)*mean(1/y))
k=sigma1/sigma2

#--------------------------------Bootstrap
for(i in 1:B){yb[i,]=sample(y,n,replace=TRUE)}
y0<-apply(yb,1,sort);
yy=1/y0
sum=apply(yy,2,sum)
for(i in 1:B){
           if(t<=sum(yy[1,i])/sum[i]) {psib[i]=y0[1,i]}  else{kn=2 
                                                         while(sum(yy[1:kn,i])/sum[i] < t) {kn=kn+1}
                                                         psib[i]=y0[kn,i]
                                                            }
                                                        
           rho_nb[i]=mean(1*(yb[i,]<=psib[i]))}
        
M=mean(rho_nb)
V=sum((rho_nb-M)^2)/(B-1)
BVL=M-z*sqrt(V); BVU=M+z*sqrt(V)
#---------------------------------NA1
NA1L=rho_n-z*sqrt(sigma1/n)
NA1U=rho_n+z*sqrt(sigma1/n)
#---------------------------------NA2
NA2L=rho_tilde-z*sqrt(sigma1/n)
NA2U=rho_tilde+z*sqrt(sigma1/n)
#--------------------------------EL1
f3<-function(gamma){a=(1*(y<=psi)-gamma)^(-1)
                      -2*sum(log(n*a/sum(a)))-k*q}

g1=uniroot(f3, interval=c(-5,-0.01), extendInt="yes")$root
g2=uniroot(f3, interval=c(1.01,5), extendInt="yes")$root

w1=((1*(y<=psi)-g1)^(-1))/sum((1*(y<=psi)-g1)^(-1))
w2=((1*(y<=psi)-g2)^(-1))/sum((1*(y<=psi)-g2)^(-1))

ELL1=sum(w1*(y<=psi))
ELU1=sum(w2*(y<=psi))
#--------------------------------EL2
x=psi*(t-1*(y<=psi))/y + 1*(y<=psi)
x0=sort(x)

f2<-function(gamma){a=(x-gamma)^(-1)
                      -2*sum(log(n*a/sum(a)))-q}

d1=uniroot(f2, interval=c(-5,x0[1]-0.0000000001))$root
d2=uniroot(f2, interval=c(x0[n]+0.01,5))$root

w1=((x-d1)^(-1))/sum((x-d1)^(-1))
w2=((x-d2)^(-1))/sum((x-d2)^(-1))

ELL2=sum(w1*x)
ELU2=sum(w2*x)
#--------------------------------print
cat("t0:",t,"\n")
cat("rho_hat(t0)=",rho_n,"\n") 
cat("EL1 Confidence interval= (", ELL1, "," , ELU1, "),", "Length=", ELU1-ELL1, "\n")
cat("EL2 Confidence interval= (", ELL2, "," , ELU2, "),", "Length=", ELU2-ELL2, "\n")
cat("NA1 Confidence interval= (", NA1L, "," , NA1U, "),", "Length=", NA1U-NA1L, "\n")
cat("NA2 Confidence interval= (", NA2L, "," , NA2U, "),", "Length=", NA2U-NA2L, "\n")
cat("Bootstrap Confidence interval= (", BVL, "," , BVU, "),", "Length=", BVU-BVL, "\n")
100*c(ELL1-rho_n, ELU1-rho_n, ELL2-rho_n,ELU2-rho_n,NA1L-rho_n,NA1U-rho_n,NA2L-rho_n,NA2U-rho_n,BVL-rho_n,BVU-rho_n)
}

band=sapply(t0, CI)
#-------------------------------Plots
#------Length
Length<-read.csv("c:\\Users\\Nikam\\Documents\\SARAH\\My works\\Vejdani\\Data\\Length.csv",header=T)
t<-seq(0.05,0.95,0.05)
plot(t,Length[[1]], col="orange", type='b',  lwd='2', xlab='t', ylab='Length', ylim=c(0.01,0.08))
lines(t,Length[[2]], col="forestgreen", type='b', pch=2, lwd='2')
lines(t,Length[[3]], col="blue", lty=2, type='b', pch=1, cex=3, lwd='2')
lines(t,Length[[4]], col="black", lty=2, type='b', pch=3, cex=2, lwd='2')
lines(t,Length[[5]], col="dark violet", type='b',  pch=3, lwd='2')
legend(0.7,0.04, c("EL1", "EL2","NA1","NA2","Boot"), col=c("orange","forestgreen","blue", "black", "dark violet"), pch=c(3,2,1,3,3), lty=c(1,1,2,2,1), lwd=c(2,2,2,2,2), merge=FALSE)
#------Confidence intervals
plot(t[-1], band[1,][-1], col="orange", type='b', lwd='2', xlab='t', ylab='Confidence Limits', ylim=c(-5,5))
lines(t[-1],band[2,][-1],col="orange", type='b',  lwd='2')
lines(t[-1],band[3,][-1],col="forestgreen", type='b', pch=2, lwd='2')
lines(t[-1],band[4,][-1],col="forestgreen", type='b', pch=2, lwd='2')
lines(t[-1],band[5,][-1], col="blue", lty=2, type='b', cex=3, lwd='2')
lines(t[-1],band[6,][-1], col="blue", lty=2, type='b', cex=3, lwd='2')
lines(t[-1],band[7,][-1], col="black", lty=2, type='b', pch=3, cex=2, lwd='2')
lines(t[-1],band[8,][-1], col="black", lty=2, type='b', pch=3, cex=2, lwd='2')
lines(t[-1],band[9,][-1], col="dark violet", type='b',  pch=3, lwd='2')
lines(t[-1],band[10,][-1],col="dark violet", type='b',  pch=3, lwd='2')
lines(t[-1],rep(0,18), col="dark red", type='b', pch=4, lwd='2')
legend(0.7,-0.1, bty='n', c("EL1", "EL2","NA1","NA2", "Boot"), col=c("orange","forestgreen","blue", "black", "dark violet"), lty=c(1,1,2,2,1), lwd=c(2,2,2,2,2), pch=c(1,2,1,3,3), merge=FALSE)

