### computation of the power under H_1, for an already known critical value of the test
### i.e. the power does not use the data any more
### here, the functions, which were introduced in CriticalValue.R, are used in a different context
### the power is evaluated for a specific H_1, while the critical value is obtained based on the distribution of the test statistic under H_0
rm(list=ls())
### init.
  Delta=0.2;
  K=1.1;
  sigma=2;
  lambda=0.2;
  nr=10; #within VW
  mystep=0.1;
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

VW=function(z) ###### finding the supremum
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

deriv=function(a,b,myint)
  {#derivative of vectors a,b
   p=length(a);
   derv=derw=rep(0.5,p);
   myzero=p/2+0.5; #print(c("zero", myzero));
   for (i in 2:(p-1))
      derv[i]=(a[i+1]-a[i-1])/(2*mystep);
   for (i in 2:(myzero-1))
     {#derv[i]=(a[i+1]-a[i-1])/(2*mystep);
      derw[i]=(b[i+1]-b[i-1])/(2*mystep);
     }
   for (i in (myzero+1):(p-1))
     {#derv[i]=(a[i+1]-a[i-1])/(2*mystep);
      derw[i]=(b[i+1]-b[i-1])/(2*mystep);
     }
 # print(cbind(a,b,derv,derw,derw/derv));#values of the derivatives
 # x11();
 # jpeg("U:/obr11.jpg", height=5,width=5,units="in", res=600); 
 #  plot(myint[2:(p-1)],derv[2:(p-1)], xlab=" ", ylab=" ", ylim=c(-0.1,0.5), xlim=c(-5,5)); lines(myint[2:(p-1)],derv[2:(p-1)]);
 #  plot(myint[2:(myzero-1)],derv[2:(myzero-1)], xlab=" ", ylab=" ", ylim=c(-0.05,1.2), xlim=c(-5,5)); lines(myint[2:(myzero-1)],derv[2:(myzero-1)]);
 #  points(myint[(myzero+1):(p-1)],derv[(myzero+1):(p-1)], xlab=" ", ylab=" "); lines(myint[(myzero+1):(p-1)],derv[(myzero+1):(p-1)]);
 #  points(myint[2:(myzero-1)],derw[2:(myzero-1)], col="red"); lines(myint[2:(myzero-1)],derw[2:(myzero-1)], col="red");
 #  points(myint[(myzero+1):(p-1)],derw[(myzero+1):(p-1)], col="red"); lines(myint[(myzero+1):(p-1)],derw[(myzero+1):(p-1)], col="red");
   po=derw/derv;
 #  points(myint[2:(myzero-1)],po[2:(myzero-1)], col="brown"); lines(myint[2:(myzero-1)], po[2:(myzero-1)], col="brown");
 #  points(myint[(myzero+1):(p-1)], po[(myzero+1):(p-1)], col="brown"); lines(myint[(myzero+1):(p-1)], po[(myzero+1):(p-1)], col="brown");
 #  dev.off() 
 list(x=myint,podil=po)
}#deriv    

testh=function(zz,x,po)#technical
  {if (x[1]<zz) {lo=which.max(x[x<zz]);} 
    else {lo=1;}
  #print(c(zz,lo,lo+1)); 
  out=po[lo];
  #if (abs(lo-zz)>abs((lo+1)-zz))
  #  out=po[up];#nedef.
  #print("testh"); print(c(zz,out));
  list(out=out)
}#testh

mytest=function(x,po,p) ### for real data
  {#plot(x,po, ylim=c(-0.05,2)); #simplistic, not needed
   mydata=myval=rep(0,20); mydata[1:15]=rnorm(15); mydata[16:20]=rnorm(5,mean=0,sd=2.5);
   for (i in 1:length(mydata))
     myval[i]=testh(mydata[i],x,po)$out;
  krith=prod(myval);
  #print("in mytest, generated data"); print(mydata);  
  #print("in mytest, individual likelihood values");
  #print(myval);
  #print(c("in mytest, critical value", krith));
  list(krith=krith)
}#mytest  

main=function()
  {myint=seq(-5,5,mystep); #=-5:5; #or enumeration 
   #index of zero: length(myint)/2+0.5;
  vv=ww=rep(0,length(myint)); 
  for (i in 1:length(myint))
    {rm=VW(myint[i]);
     vv[i]=rm$vout;
     ww[i]=rm$wout;}
  rm=ww[2:length(myint)]-ww[1:(length(myint)-1)]; #must be non-negative
 # print(cbind(myint,vv,ww,rm));
 # x11();
 # jpeg("U:/obr11.jpg", height=5,width=5,units="in", res=600);   
 # plot(myint, vv, ylim=c(0,1.2), xlab=" ", ylab=" "); 
 #   points(myint, ww, col="red", pch=3);
  sm=deriv(vv,ww,myint);
  podil=sm$podil; x=sm$x;
  sm=mytest(x,podil,length(myint));
 # dev.off()  
}#main   

myrepeat=function(nc)
  {out=rep(0,nc);
  print(c("i, critical value"));
  for (i in 1:nc)
    {out[i]=main()$krith;
    print(c(i,out[i]));#critical value
    }
  #plot(out);
  list(out=out)
}#myrepeat

### critical value
myrun=function()
{nc=100;
my=opakovat(nc)$out;
plot(sort(my));
print(sum(my>1.29)/nc);
}# myrun()