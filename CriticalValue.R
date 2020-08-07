### computing the critical value of the test, for given parameters
rm(list=ls())
### initialization
  sigma=sqrt(3); #under H_1
  sigmahat=6;
  lambda=0.25; #under H_1
  n=20;
  Delta=0.2;
  K=1.1;
  nr=10; #within VW
  mystep=0.1;
  ### generating the random data also depends on the user
### functions
fd=function(z)
  {out=-z*dnorm(z);
   list(out=out)}#
f3=function(z)
  {out=(-z^3+3*z)*dnorm(z);
   list(out=out)}#
gd=function(z)
  {out=(lambda-1)*z*dnorm(z) - (lambda*z/sigma^2)*dnorm(z,sd=sigma);
   list(out=out)}#
g3=function(z)
  {out= (3*z+(1-lambda)*z^3-3*z*lambda)*dnorm(z) + (3*z*lambda/sigma^4 - z^3*lambda/sigma^6)*dnorm(z,sd=sigma);
   list(out=out)}#
F=function(z,delta,kappa)
  {out=pnorm(z) + delta *fd(z)$out/2 + kappa * delta^2 * f3(z)$out/24;
   list(out=out)}#
G0=function(z)
  {out = (1-lambda)*pnorm(z) + lambda*pnorm(z, mean=0, sd=sigma);
   list(out=out)}# 
G=function(z,delta,kappa)
  {out = G0(z)$out + delta *gd(z)$out/2 + kappa * delta^2 * g3(z)$out/24;
   if (out>1) print(c("attention",out));
   list(out=out)}# 

VW=function(z) ##### finding the supremum under H0 and H1
  {#print("VW");
   pv=pw=rep(0,nr+4); #random choices + boundaries
   for (i in 1:nr)
     {de=runif(1, min=0, max=Delta); 
      ka=runif(1, min=1, max=K);
      pv[i]=F(z,de,ka)$out;
      pw[i]=G(z,de,ka)$out;
     }#for
    pv[nr+1]=F(z,de=0,ka=1)$out;     pv[nr+2]=F(z,de=0,ka=K)$out; 
    pv[nr+3]=F(z,de=Delta,ka=1)$out; pv[nr+4]=F(z,de=Delta,ka=K)$out;
    pw[nr+1]=G(z,de=0,ka=1)$out;     pw[nr+2]=G(z,de=0,ka=K)$out; 
    pw[nr+3]=G(z,de=Delta,ka=1)$out; pw[nr+4]=G(z,de=Delta,ka=K)$out;
    vout=max(pv); wout=max(pw); #print(c(z,vout,wout));
    list(vout=vout, wout=wout);   
  }#VW   

derivdobre=function(a,b,myint) #derivative of vectors a,b
  {p=length(a);
   derv=derw=rep(0.5,p);
   nula=p/2+0.5; #print(c("nula", nula));
   for (i in 2:(p-1))
      derv[i]=(a[i+1]-a[i-1])/(2*mystep);
   for (i in 2:(nula-1))
     {#derv[i]=(a[i+1]-a[i-1])/(2*mystep);
      derw[i]=(b[i+1]-b[i-1])/(2*mystep);
     }
   for (i in (nula+1):(p-1))
     {#derv[i]=(a[i+1]-a[i-1])/(2*mystep);
      derw[i]=(b[i+1]-b[i-1])/(2*mystep);
    }
  po=derw/derv; 
  # print(cbind(a,b,derv,derw,derw/derv));
  # x11();
  # jpeg("U:/obr11.jpg", height=5,width=5,units="in", res=600); 
  # plot(myint[2:(p-1)],derv[2:(p-1)], xlab=" ", ylab=" ", ylim=c(-0.1,0.5), xlim=c(-5,5)); lines(myint[2:(p-1)],derv[2:(p-1)]);
  # plot(myint[2:(nula-1)],derv[2:(nula-1)], xlab=" ", ylab=" ", ylim=c(-0.05,2), xlim=c(-5,5)); lines(myint[2:(nula-1)],derv[2:(nula-1)]);
  # points(myint[(nula+1):(p-1)],derv[(nula+1):(p-1)], xlab=" ", ylab=" "); lines(myint[(nula+1):(p-1)],derv[(nula+1):(p-1)]);
  # points(myint[2:(nula-1)],derw[2:(nula-1)], col="red"); lines(myint[2:(nula-1)],derw[2:(nula-1)], col="red");
  # points(myint[(nula+1):(p-1)],derw[(nula+1):(p-1)], col="red"); lines(myint[(nula+1):(p-1)],derw[(nula+1):(p-1)], col="red");
  # points(myint[2:(nula-1)],po[2:(nula-1)], col="brown"); lines(myint[2:(nula-1)], po[2:(nula-1)], col="brown");
  # points(myint[(nula+1):(p-1)], po[(nula+1):(p-1)], col="brown"); lines(myint[(nula+1):(p-1)], po[(nula+1):(p-1)], col="brown");
  # dev.off() 
 list(podil=po)
}#derivdobre    

testh=function(myint, po, zz)
  {if (myint[1]<zz) {lo=which.max(myint[myint<zz]);} 
    else {lo=1;}
   out=po[lo];
   list(out=out)
}#testh

mytest=function(myint, po) ### for real data
  {mydata=rnorm(n); #H0 ############################### generating the data
   for (i in 1:18) mydata[i]=rnorm(1,mean=0,sd=sqrt(sigmahat)); #possible H1
   z=myval=rep(0,length(mydata));
   for (i in 1:length(z))
     myval[i]=testh(myint, po, mydata[i])$out;#the closest value, which is smaller
  krith=prod(myval);
  #print(c("critical value", krith));
  list(krith=krith)
}#mytest  

main=function(iter=1000)
  {kh=rep(-1,n);
  myint=seq(-5,5,mystep); #=-5:5; #or enumeration 
   #index of zero: length(myint)/2+0.5;
  vv=ww=rep(0,length(myint)); 
  for (i in 1:length(myint))
    {rm=VW(myint[i]);
     vv[i]=rm$vout;
     ww[i]=rm$wout;}
    #ko=ww[2:length(myint)]-ww[1:(length(myint)-1)]; #ko>0
    #print(cbind(myint,vv,ww,rm));
    #jpeg("U:/obr11.jpg", height=5,width=5,units="in", res=600);   
    #plot(myint, vv, ylim=c(0,1.2), xlab=" ", ylab=" "); 
    #  points(myint, ww, col="red", pch=3);
  sm=derivdobre(vv,ww,myint);
  for (i in 1:iter)
    kh[i]=mytest(myint, sm$podil)$krith;
  list(out=kh)
}#main   

### critical value
#   my=main(10000); plot(sort(my$out)); sort(my$out)[9500]
### power of the test
#   my=main(10000); plot(sort(my$out)); sum(my$out>0.37)/10000;