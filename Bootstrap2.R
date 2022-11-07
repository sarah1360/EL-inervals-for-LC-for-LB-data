rm(list=ls(all=TRUE))
library(tictoc)
library(invgamma);
tic(boot)
out=function(t){
set.seed(2021);
N=10000; B=400; n=50; 
shape=3; rate=1;
yb=matrix(nrow=B, ncol=n)
psi=rho_n=c(); BVL=BVU=c();
alpha=0.05; z=qnorm(1-alpha/2); 
#-------------------------------
rho=pinvgamma(qinvgamma(t, shape+1, rate), shape, rate)
for(j in 1:N){
y=rinvgamma(n, shape, rate)

y00<-sort(y)
yyy=1/y00
sum1=sum(yyy)

if(t<=yyy[1]/sum1){psi=y00[1]} else{kn=2
                                 while(sum(yyy[1:kn])/sum1 <t){kn=kn+1}
                                 psi0=y00[kn]}
rho_n0=mean(1*(y<=psi0))
#---------------------------------bootstrap resampling
for(i in 1:B){yb[i,]=sample(y,n,replace=TRUE)}
y0<-apply(yb,1,sort);
yy=1/y0
sum=apply(yy,2,sum)
for(i in 1:B){
           if(t<=sum(yy[1,i])/sum[i]) {psi[i]=y0[1,i]}  else{kn=1  
                                                        while(sum(yy[1:kn,i])/sum[i] < t) {kn=kn+1}
                                                        psi[i]=y0[kn,i]
                                                            }
                                                        
           rho_n[i]=mean(1*(yb[i,]<=psi[i]))}
        
#-------------------------------
M=mean(rho_n)
V=mean((rho_n-M)^2)
T=(rho_n-rho_n0)/sqrt(V/n)

#-------------------------------
BVL[j]=M-z*sqrt(V); BVU[j]=M+z*sqrt(V)
}
L=mean(BVU-BVL); C=mean(1*(rho>=BVL)*(BVU>=rho))

cat("Rho:",rho,"\n")

cat("Bootstrap-based Confidence Interval = (",  mean(BVL)," , ", mean(BVU),"); Average Length =",L," and Coverage Probability=",C )                               
cat("\n")

}
out(0.1)
toc()
t=seq(0.1,0.9,0.1)
Boot=sapply(t,out)
write.table(Boot,"D:\\bootstrap.txt",col.names=F,row.names=F)
                   
                  