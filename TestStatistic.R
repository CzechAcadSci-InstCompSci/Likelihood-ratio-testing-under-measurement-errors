### we consider a given random sample
### the first step is the computation of the test statistic 
### the test considers a simple (i.e. not composite) null hypothesis and a simple alternative hypothesis
### the implementation is technically difficult but conceptually straightforward, following the paper Broniatowski et al. (2018)
rm(list=ls())
### init.
  Delta=0.1;
  K=0.1; #parameter corresponding to the variability of the measurement errors
  sigma=3;
  alpha=0.25;
  nr=10; #within the function VW
  mystep=0.1;
### functions
fd=function(z)
  {out=-z*dnorm(z);
   list(out=out)}#
f3=function(z)
  {out=(-z^3+3*z)*dnorm(z);
   list(out=out)}#
gd=function(z)
  {out=(alpha-1)*z*dnorm(z) - (alpha*z/sigma^2)*dnorm(z,sd=sigma);
   list(out=out)}#
g3=function(z)
  {out= (3*z+(1-alpha)*z^3-3*z*alpha)*dnorm(z) + (3*z*alpha/sigma^4 - z^3*alpha/sigma^6)*dnorm(z,sd=sigma);
   list(out=out)}#
F=function(z,delta,kappa)
  {out=pnorm(z) + delta *fd(z)$out/2 + kappa * delta^2 * f3(z)$out/24;
   list(out=out)}#
G0=function(z)
  {out = (1-alpha)*pnorm(z) + alpha*pnorm(z, mean=0, sd=sigma);
   list(out=out)}# 
G=function(z,delta,kappa)
  {out = G0(z)$out + delta *gd(z)$out/2 + kappa * delta^2 * g3(z)$out/24;
   if (out>1) print(c("attention",out));
   list(out=out)}# 
###### finding the supremum, required by the test statistic
VW=function(z)
  {#print("VW");
   pv=pw=rep(0,nr+4); #random choices + boundaries
   for (i in 1:nr)
     {delta=runif(1, min=0, max=Delta); 
      kappa=runif(1, min=0, max=K);
      pv[i]=F(z,delta,kappa)$out;
      pw[i]=G(z,delta,kappa)$out;
     }#for
    pv[nr+1]=F(z,delta=0,kappa=0)$out; pv[nr+2]=F(z,delta=0,kappa=K)$out; 
    pv[nr+3]=F(z,delta=Delta,kappa=0)$out; pv[nr+4]=F(z,delta=Delta,kappa=K)$out;
    pw[nr+1]=G(z,delta=0,kappa=0)$out; pw[nr+2]=G(z,delta=0,kappa=K)$out; 
    pw[nr+3]=G(z,delta=Delta,kappa=0)$out; pw[nr+4]=G(z,delta=Delta,kappa=K)$out;
    vout=max(pv); wout=max(pw); #print(c(z,vout,wout));
    list(vout=vout, wout=wout);   
  }#VW   

main=function()
  {myint=seq(-5,5,mystep); #=-5:5; #or enumeration 
   #index of zero: length(myint)/2+0.5;
  vv=ww=rep(0,length(myint)); 
  for (i in 1:length(myint))
    {rm=VW(myint[i]);
     vv[i]=rm$vout;
     ww[i]=rm$wout;}
  rm=ww[2:length(myint)]-ww[1:(length(myint)-1)]; #nesmi byt zaporne
 # print(cbind(myint,vv,ww,rm));
 # x11();
 # jpeg("U:/obr11.jpg", height=5,width=5,units="in", res=600);   
  plot(myint, vv, ylim=c(0,1.2), xlab=" ", ylab=" "); 
    points(myint, ww, col="red", pch=3);
  sm=derivdobre(vv,ww,myint);
  podil=sm$podil; x=sm$x;
  sm=mytest(x,podil,length(myint));
 # dev.off()  
}#main   

mytest=function(x,po,p) ### for real data
  {#plot(x,po, ylim=c(-0.05,2)); #simplistic, not needed
   mydata=c(-0.77, 0.43, 0.86, 0.58, -1.29);
   z=myval=rep(0,length(mydata));
   print("in mytest"); print(cbind(x[2:(p-1)],po[2:(p-1)]));
   for (i in 1:length(z))
     myval[i]=testh(mydata[i],x,po)$out;#the closest value, which is smaller
  krith=prod(myval);
  print(c("test statistic for fixed data", krith));
  list(krith=krith)
}#mytest  
#nejbl=function(poz,x,po)
#  {print(cbind(x,po));
#  }#nejbl
deriv=function(a,b,myint)
  {#derivative of vectors a,b
   p=length(a);
   derv=derw=rep(0.5,p);
   nula=p/2+0.5; print(c("nula", nula));
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
   print(cbind(a,b,derv,derw,derw/derv));
# x11();
# jpeg("U:/obr13.jpg", height=5,width=5,units="in", res=600); 
 #  plot(myint[2:(p-1)],derv[2:(p-1)], xlab=" ", ylab=" ", ylim=c(-0.1,0.5), xlim=c(-5,5)); lines(myint[2:(p-1)],derv[2:(p-1)]);
   plot(myint[2:(nula-1)],derv[2:(nula-1)], xlab=" ", ylab=" ", ylim=c(-0.05,5), xlim=c(-5,5)); lines(myint[2:(nula-1)],derv[2:(nula-1)]);
   points(myint[(nula+1):(p-1)],derv[(nula+1):(p-1)], xlab=" ", ylab=" "); lines(myint[(nula+1):(p-1)],derv[(nula+1):(p-1)]);
   points(myint[2:(nula-1)],derw[2:(nula-1)], col="red"); lines(myint[2:(nula-1)],derw[2:(nula-1)], col="red");
   points(myint[(nula+1):(p-1)],derw[(nula+1):(p-1)], col="red"); lines(myint[(nula+1):(p-1)],derw[(nula+1):(p-1)], col="red");
   po=derw/derv;
   points(myint[2:(nula-1)],po[2:(nula-1)], col="brown"); lines(myint[2:(nula-1)], po[2:(nula-1)], col="brown");
   points(myint[(nula+1):(p-1)], po[(nula+1):(p-1)], col="brown"); lines(myint[(nula+1):(p-1)], po[(nula+1):(p-1)], col="brown");
 # dev.off() 
 list(x=myint,podil=po)
}#deriv
    
testh=function(zz,x,po)
  {if (x[1]<zz) {lo=which.max(x[x<zz]);} 
    else {lo=1;}
  #print(c(zz,lo,lo+1)); 
  out=po[lo];
  #if (abs(lo-zz)>abs((lo+1)-zz))
  #  out=po[up];#nedef.
  #print("testh"); print(c(zz,out));
  list(out=out)
}#testh
### main()