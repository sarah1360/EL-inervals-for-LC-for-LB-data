rm(list=ls(all=TRUE))
library(invgamma);library(rmutil); library(rootSolve); library(tictoc);

out=function(t){
set.seed(2021);
N=10000; n=400; 
shape=3; rate=1;
psi=rho_n=c();
alpha=0.05; 
z=qnorm(1-alpha/2)
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
if(t<=yy[1]/sum){psi=y0[1]} else{kn=1
                                 while(sum(yy[1:kn])/sum <t){kn=kn+1}
                                 psi=y0[kn]}

#--------------------------------
rho_n=mean(1*(y <=psi))
sigma2=rho_n*(1-rho_n)
sigma1=max(0,(psi^2)*(((1-t)^2)*mean((1/(y^2))*(y<=psi))+(t^2)*mean((1/(y^2))*(y>psi)))+sigma2+2*psi*t*(t-1)*mean(1/y))
#--------------------------------
rho_n=mean(1*(y <=psi))
x=psi*(t-1*(y<=psi))/y + 1*(y<=psi)
#--------------------------------
e1=mean(x)-z*sqrt(sigma1/n)
e2=mean(x)+z*sqrt(sigma1/n)
cb=rbind(cb, c(e1,e2, e2-e1,sign((rho-e1)*(e2-rho  )   )  +1  )  )
                                      
  }

cat("Rho:",rho,"\n")
cat("NA-based Confidence Band = (",  mean(cb[,1])," , ", mean(cb[,2]),"); Average Length =",  mean(cb[,3])," and Coverage Probability=",mean(cb[,4]) /2 )                               
cat("\n")

                                        
     }

t=seq(0.1,0.9,0.1)
NAI=sapply(t,out)
write.table(NAI2,"D:\\NA2.txt",col.names=F,row.names=F)
