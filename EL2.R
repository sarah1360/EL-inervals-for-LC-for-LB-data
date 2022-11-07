rm(list=ls(all=TRUE))
library(invgamma);library(rmutil); library(rootSolve);library(tictoc)
tic(EL2)
out=function(t){
set.seed(2021);
N=10000; n=50; 
shape=5; rate=1;
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
#--------------------------------Akbari estimator
if(t<=yy[1]/sum){psi=y0[1]} else{kn=1
                                 while(sum(yy[1:kn])/sum <t){kn=kn+1}
                                 psi=y0[kn]}
#--------------------------------
rho_n=mean(1*(y <=psi))
x=psi*(t-1*(y<=psi))/y + 1*(y<=psi)
x0=sort(x)
#--------------------------------
 f3<-function(gamma){a=(x-gamma)^(-1)
                      -2*sum(log(n*a/sum(a)))-q}
#--------------------------------
g1=uniroot(f3, interval=c(-5,x0[1]-0.01), extendInt="yes")$root
g2=uniroot(f3, interval=c(x0[n]+0.01,5), extendInt="yes")$root

w1=((x-g1)^(-1))/sum((x-g1)^(-1))
w2=((x-g2)^(-1))/sum((x-g2)^(-1))

e1=sum(w1*x)
e2=sum(w2*x)

cb=rbind(cb, c(e1,e2, e2-e1,sign((rho-e1)*(e2-rho  )   )  +1  )  )
                                      
  }

cat("Rho:",rho,"\n")

cat("EL2-based Confidence Band = (",  mean(cb[,1])," , ", mean(cb[,2]),"); Average Length =",  mean(cb[,3])," and Coverage Probability=",mean(cb[,4]) /2 )                               
cat("\n")

                                         
     }
out(0.1)
toc()
t=seq(0.1,0.9,0.1)
EL2=sapply(t,out)
write.table(EL2,"D:\\EL2.txt",col.names=F,row.names=F)





