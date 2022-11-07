rm(list=ls(all=TRUE))
library(invgamma);library(rmutil); library(rootSolve);library(tictoc)
tic(EL1)
out=function(t){
set.seed(2021);
N=10000; n=50; 
shape=3; rate=1;
psi=rho_n=c();
alpha=0.05; q=qchisq(p=alpha, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE);
cb=c()
#-------------------------------
rho=pinvgamma(qinvgamma(t, shape+1, rate), shape, rate)
for(j in 1:N){
y=rinvgamma(n, shape, rate)
#-------------------------------
y0<-sort(y)
yy=1/y0
sum=sum(yy)
                                                    
#--------------------------------estimator
if(t<=yy[1]/sum){psi=y0[1]} else{kn=2
                                 while(sum(yy[1:kn])/sum <t){kn=kn+1}
                                 psi=y0[kn]}
#--------------------------------
rho_n=mean(1*(y <=psi))
sigma2=rho_n*(1-rho_n)
sigma1=max(0,(psi^2)*(((1-t)^2)*mean((1/(y^2))*(y<=psi))+(t^2)*mean((1/(y^2))*(y>psi)))+sigma2+2*psi*t*(t-1)*mean(1/y))
k=sigma1/sigma2
#--------------------------------
 f3<-function(gamma){a=(1*(y<=psi)-gamma)^(-1)
                      -2*sum(log(n*a/sum(a)))-k*q}
#--------------------------------
g1=uniroot(f3, interval=c(-5,-0.01), extendInt="yes")$root
g2=uniroot(f3, interval=c(1.01,5), extendInt="yes")$root

w1=((1*(y<=psi)-g1)^(-1))/sum((1*(y<=psi)-g1)^(-1))
w2=((1*(y<=psi)-g2)^(-1))/sum((1*(y<=psi)-g2)^(-1))

e1=sum(w1*(y<=psi))
e2=sum(w2*(y<=psi))

cb=rbind(cb, c(e1,e2, e2-e1,sign((rho-e1)*(e2-rho  )   )  +1  )  )
                                      
  }

cat("Rho:",rho,"\n")

cat("EL1-based Confidence Band = (",  mean(cb[,1])," , ", mean(cb[,2]),"); Average Length =",  mean(cb[,3])," and Coverage Probability=",mean(cb[,4]) /2 )                               
cat("\n")

                                         
     }
out(0.1)
toc()
t=seq(0.1,0.9,0.1)
EL1=sapply(t,out)
write.table(EL1,"D:\\EL1.txt",col.names=F,row.names=F)





