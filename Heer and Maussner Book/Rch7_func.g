@ ------------------------------ Rch7_func.g ----------------------------

date: March 20, 2007

author: Burkhard Heer

Chapter 7.2., distribution function

computation of the policy functions and the distribution function
in the heterogenous-agent neoclassical growth model
with value function iteration

FUNCTION Approximation

monotonocity of a' in a

concavity of v

Linear interpolation between grid points

Maximization: golden section search method

input procedures:
golden -- routine for golden section search
lininter -- linear interpolation

-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph;
Macheps=1e-20;

#include ch8_toolbox.src; 

graphset;
/*
save path="C:\\Dokumente und Einstellungen\\ba6wp1\\Eigene Dateien\\buch\\prog\\ch521";
*/

h0=hsec;


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

/* numerical computation */
tol=0.001;          /* stopping criterion for final solution of K */
tol1=1e-7;          /* stopping criterion for golden section search */
neg=-1e10;
eps=0.05;
nexp=2;         /* number of coefficients in functional approximation */
nint=20;        /* nint*nexp - number of nodes for Gauss-Chebyshev quadrature */
psi=0.5;        /* update of abare, avare, afold, .. */

/* calibration of parameters */
alpha=0.36;
beta=0.995;
delta=0.005;
sigma=2;
r=0.005;
w=5;
tau=0.02;
rep=0.25;
pp=(0.9565~0.0435|0.5~0.5);
/* stationary employment /unemployment */
pp1=equivec1(pp);

nn0=pp1[1];
kk0=(alpha/(1/beta-1+delta))^(1/(1-alpha))*nn0;
w0=(1-alpha)*kk0^alpha*nn0^(-alpha);
b=w0*rep;

amin1=-2;                 /* asset grid */
amax1=1000;
na=201;
astep=(amax1-amin1)/(na-1);
a=seqa(amin1,astep,na);
/* initialization of the value function: */
v=u((1-tau)*r*a+(1-tau)*w)~u((1-tau)*r*a+b);
v=v/(1-beta);                /* agents consume their income */
copt=zeros(na,2);           /* optimal consumption */
aopt=zeros(na,2);           /* optimal next-period assets */

nit=100;                /* number of maximum iterations over value function */
ngk=0;               /* initial number of iterations over distribution function */
crit1=1;
kritg=1;
nq=100;

psik=0.95;

kritw=1;
nn0=pp1[1];

kk0=(alpha/(1/beta-1+delta))^(1/(1-alpha))*nn0;
kkbar=kk0;
kk1=kk0-1;
kbarq=zeros(nq,1);      /* convergence of kk0 */
k1barq=zeros(nq,1);      /* convergence of kk1 */
kritwq=zeros(nq,2);     /* convergence of value function v */
crit=1;


q=0;
do until q==nq or (q>20 and abs((kk1-kk0)/kk0)<tol); /* iteration over K */
    q=q+1;
    w=(1-alpha)*kk0^alpha*nn0^(-alpha);
    r=alpha*kk0^(alpha-1)*nn0^(1-alpha)-delta;
    kbarq[q]=kk0;
    if ngk<25000;
        ngk=ngk+500;
    endif;
    nit=nit+2;

    crit=1;
    j=0;
    do until (j>50 and crit<tol) or (j==nit);
        j=j+1;
        vold=v;
        "iteration j~q: " j~q;
	    "kk0~kk1" kk0~kk1;
        "time elapsed: " etstr(hsec-h0);
        "error value function: " crit;
        "error k: " crit1;
        "error distribution function: " kritg;
        "kk0: " kk0;

@ ----------------------------------------------------------------

Step 2: Iteration of the value function

----------------------------------------------------------------- @

        e=0;    /* iteration over the employment status */
        do until e==2;  /* e=1 employed, e=2 unemployed */
            e=e+1;
            i=0;        /* iteration over asset grid a in period t */
            l0=0;
            do until i==na;
                i=i+1;
                l=l0;
                v0=neg;
                ax=0; bx=-1; cx=a[na];
                do until l==na; /* iteration over a' in period t*1 */
                    l=l+1;
                    if e==1;
                        c=(1+(1-tau)*r)*a[i]+(1-tau)*w-a[l];
                     else;
                        c=(1+(1-tau)*r)*a[i]+b-a[l];
                    endif;
                    if c>0;
                        v1=bellman(a[i],a[l],e);
                        if v1>v0;
                            v[i,e]=v1;
                            if l==1;
                               ax=a[1]; bx=a[1]; cx=a[2];
                            elseif l==na;
                                ax=a[na-1]; bx=a[na]; cx=a[na];
                            else;
                                ax=a[l-1]; bx=a[l]; cx=a[l+1];
                            endif;
                            v0=v1;
                            l0=l-1;
                        else;
                            l=na;   /* concavity of value function */
                        endif;
                    else;
                        l=na;
                    endif;
                endo;   /* l=1,..,na */
                if ax==bx;  /* boundary optimum, ax=bx=a[1]  */
                    bx=ax+eps*(a[2]-a[1]);
                    if value(bx,e)<value(ax,e);
                        aopt[i,e]=a[1];
                        v[i,e]=bellman(a[i],a[1],e);
                    else;
                        aopt[i,e]=golden(&value1,ax,bx,cx,tol1);
                        v[i,e]=bellman(a[i],aopt[i,e],e);
                    endif;
                elseif bx==cx;  /* boundary optimum, bx=cx=a[n] */
                    bx=cx-eps*(a[na]-a[na-1]);
                    if value(bx,e)<value(cx,e);
                        aopt[i,e]=a[na];
                    else;
                        aopt[i,e]=golden(&value1,ax,bx,cx,tol1);
                    endif;
                else;
                    aopt[i,e]=golden(&value1,ax,bx,cx,tol1);
                endif;
                v[i,e]=bellman(a[i],aopt[i,e],e);
            endo;   /* i=1,..na */
        endo;   /* e=1,2 */
        crit=meanc(abs(vold-v));
    endo;   /* j=1,..nit */

    copt[.,1]=(1+(1-tau)*r)*a+(1-tau)*w-aopt[.,1];
    copt[.,2]=(1+(1-tau)*r)*a+b-aopt[.,2];

@ --------------------------------------------------------------------------------

computation of the stationar distribution

------------------------------------------------------------------------------ @


af=zeros(nexp+1,2);
afold=af;
afold[1,1]=pp1[1]/(amax1-amin1);
afold[1,2]=pp1[2]/(amax1-amin1);

/* computation of the initial mean and variances */
abareini=amean1(afold[1:3,1])/pp1[1];
abaruini=amean1(afold[1:3,2])/pp1[2];
avareini=avar1(afold[1:3,1],abareini)/pp1[1];
avaruini=avar1(afold[1:3,2],abaruini)/pp1[2];

/* solve for the coefficients */
abare=abareini;
abaru=abaruini;
avare=avareini;
avaru=avaruini;

    x1=afold[.,1]|afold[.,2];
    {xf,jcode}=FixVMN1(x1,&mom2);
af[.,1]=xf[1:nexp+1];
af[.,2]=xf[nexp+2:2*nexp+2];


am=zeros(ngk,2);
kkit4=zeros(ngk,1);
crit=tol+1;

nround=ngk/100;
j=0;
do until j==ngk;
    j=j+1;
    afe=afold[.,1];
    afu=afold[.,2];
    abareneu=amean2(afe,afu,1)/pp1[1];
    abaruneu=amean2(afe,afu,2)/pp1[2];
    avareneu=avar(afe,afu,1,abare)/pp1[1];
    avaruneu=avar(afe,afu,2,abaru)/pp1[2];
    abare=psi*abare+(1-psi)*abareneu;
    abaru=psi*abaru+(1-psi)*abaruneu;
    avare=psi*avare+(1-psi)*avareneu;
    avaru=psi*avaru+(1-psi)*avaruneu;
     /* solve for the coefficients */
    x1=afe|afu;
    {xf,jcode}=FixVMN1(x1,&mom2);
     af[.,1]=xf[1:nexp+1];
     af[.,2]=xf[nexp+2:2*nexp+2];
     am[j,1]=amean1(af[.,1])/pp1[1]; am[j,2]=amean1(af[.,2])/pp1[2]; 
     kkit4[j]=am[j,1]*pp1[1]+am[j,2]*pp1[2];

        if round(j/nround)==j/nround;
            "iteration q~j~capital stock: " q~j~kkit4[j];
            "k0~k1: ";
            kbarq[1:q]~k1barq[1:q];
        endif;     
     crit=meanc(abs(afold-af));
     afold=psi*af+(1-psi)*afold;    
endo;   /* j=1,..nit */

    save disfunc=af;


    kk1=kkit4[ngk];
    k1barq[q]=kk1;


    tau0=b*pp1[2]/(w*nn0+r*kk0);
    tau=psik*tau+(1-psik)*tau0;

    crit1=abs((kk1-kk0)/kk0);
    kk0=psik*kk0+(1-psik)*kk1;

    
    "iteration: " q;
    "time elapsed: " etstr(hsec-h0);
    "error value function: " crit;
    "error distribution: " kritg;
    "error capital stock: " crit1;
    "kbarq~k1barq'";
    kbarq[1:q]~k1barq[1:q];
    save afunc=a;
    save k1func=k1barq;
    if round(q/2)==q/2; cls; endif;
endo;   /* q=1,..,nq */

"iteration: " q;
"time elapsed: " etstr(hsec-h0);
"error value function: " crit;
"error distribution: " kritg;
"error capital stock: " crit1;
wait;

    "press any key to plot solution of projection method... "; wait;
    "mean capital stock employed/unemployed: ";
    amean1(af[.,1])/pp1[1]~amean1(af[.,2])/pp1[2];
    "mean of the economoy: ";
    kkk4=amean1(af[.,1])+amean1(af[.,2]);
    kkk4; wait;
    title("convergence of amean");
    xy(seqa(1,1,ngk),am[1:ngk,1]*pp1[1]+am[1:ngk,2]*pp1[2]); wait;
    /* plotting the solution */
    title("wealth distribution of the employed agent");
    xlabel("asset level a");
    xy(a,af[1,1]*exp(af[2,1].*a+af[3,1].*a^2)); wait;
    title("wealth distribution of the unemployed agent");
    xy(a,af[1,2]*exp(af[2,2].*a+af[3,2].*a^2)); wait;
    wait;
    title("wealth distribution");
    xy(a,af[1,1]*exp(af[2,1].*a+af[3,1].*a^2)~af[1,2]*exp(af[2,2].*a+af[3,2].*a^2));
    wait;
    

/* plotting the solution */
_plegstr="employed\000unemployed";
title("value function");
xy(a,v);
wait;
xy(a,v[.,1]);
wait;
xy(a,v[.,2]);
wait;

title("consumption");
xy(a,copt);
wait;
xy(a,copt[.,1]);
wait;
xy(a,copt[.,2]);
wait;

title("change in asset level");
xy(a,aopt-a);
wait;

xy(a,aopt[.,1]-a);
wait;
xy(a,aopt[.,2]-a);
wait;

title("next-period assets");
_plegstr="employed\000a";
xy(a,aopt[.,1]~a);
wait;

_plegstr="unemployed\000a";
xy(a,aopt[.,2]~a);
wait;


@  ----------------------------  procedures -----------


u(x) -- utility function

value(a,e) -- returns the value of the value function for asset
                a and employment status e

value1(x) -- given a=a[i] and epsilon=e, returns the
             value of the bellman equation for a'=x

bellman -- value for the right-hand side of the Bellman equation

------------------------------------------------------- @

proc u(x);
   retp(x^(1-sigma)/(1-sigma));
endp;

proc value(x,y);
    retp(lininter(a,vold[.,y],x));
endp;

proc value1(x);
    retp(bellman(a[i],x,e));
endp;


proc bellman(a0,a1,y);
   local c;
   if y==1;
      c=(1+(1-tau)*r)*a0+(1-tau)*w-a1;
   else;
      c=(1+(1-tau)*r)*a0+b-a1;
   endif;
   if c<0;
      retp(neg);
   endif;
   if a1>=a[na];
      retp(a1^2*neg);
   endif;
   retp(u(c)+beta*(pp[y,1]*value(a1,1)+pp[y,2]*value(a1,2)));
endp;

proc lininter(xd,yd,x);
  local j;
  j=sumc(xd.<=x');
  retp(yd[j]+(yd[j+1]-yd[j]).*(x-xd[j])./(xd[j+1]-xd[j]));
endp;

proc golden(&f,ay,by,cy,tol);
    local f:proc,x0,x1,x2,x3,xmin,r1,r2,f1,f2;
    r1=0.61803399; r2=1-r1;
    x0=ay;
    x3=cy;
    if abs(cy-by)<=abs(by-ay);
        x1=by; x2=by+r2*(cy-by);
    else;
        x2=by; x1=by-r2*(by-ay);
    endif;
    f1=-f(x1);
    f2=-f(x2);
    do until abs(x3-x0)<=tol*(abs(x1)+abs(x2));
        if f2<f1;
            x0=x1;
            x1=x2;
            x2=r1*x1+r2*x3;
            f1=f2;
            f2=-f(x2);
        else;
            x3=x2;
            x2=x1;
            x1=r1*x2+r2*x0;
            f2=f1;
            f1=-f(x1);
        endif;
    endo;
    if f1<=f2;
        xmin=x1;
        else;
        xmin=x2;
    endif;
    retp(xmin);
endp;

proc equivec1(p);
    local n,x;
    n=rows(p);
    p=diagrv(p,diag(p)-ones(n,1));
    p=p[.,1:n-1]~ones(n,1);
    x=zeros(n-1,1)|1;
    retp((x'*inv(p))');
endp;


proc amean2(afe,afu,x);
    local n,m,z,y,sum1,y1e,y1u,inte,intu,i,l;
    n=rows(afe);
    m=nint*n;  /* m>n quadrature points */
    z=-cos((2.*seqa(1,1,m)-1)./(2*m).*pi);
    y=(z+1).*(amax1-amin1)/2+amin1;
    sum1=0;
    y1e=lininter(a,aopt[.,1],y);  
    y1u=lininter(a,aopt[.,2],y);
    if n==2;  
        inte=afe[1]*exp(afe[2].*y);
        intu=afu[1]*exp(afu[2].*y);
    endif;
    if n==3;
        inte=afe[1]*exp(afe[2].*y+afe[3].*y^2);
        intu=afu[1]*exp(afu[2].*y+afu[3].*y^2);
    endif;
    l=0;
    do until l==m;
        l=l+1;
        sum1=sum1+pp[1,x]*maxc(amin1|y1e[l])*inte[l]*sqrt(1-z[l]^2);
        sum1=sum1+pp[2,x]*maxc(amin1|y1u[l])*intu[l]*sqrt(1-z[l]^2);
    endo;     
    sum1=sum1*pi/2*(amax1-amin1)/m;
    retp(sum1);
endp;

proc amean1(afx);
    local n,m,z,y,sum1,intx,i,l,j;
    n=rows(afx);
    m=nint*n;  /* m>n quadrature points */
    z=-cos((2.*seqa(1,1,m)-1)./(2*m).*pi);
    y=(z+1).*(amax1-amin1)/2+amin1;
    sum1=0;
    if n==2;
        intx=afx[1]*exp(afx[2].*y);
    endif;
    if n==3;
        intx=afx[1]*exp(afx[2].*y+afx[3].*y^2);
    endif;
    l=0;
    do until l==m;
       l=l+1;
       sum1=sum1+maxc(amin1|y[l])*intx[l]*sqrt(1-z[l]^2);
    endo;     
    sum1=sum1*pi/2*(amax1-amin1)/m;
    retp(sum1);
endp;


proc avar1(afx,x);
    local n,m,z,y,sum1,intx,i,l,j;
    n=rows(afx);
    m=nint*n;  /* m>n quadrature points */
    z=-cos((2.*seqa(1,1,m)-1)./(2*m).*pi);
    y=(z+1).*(amax1-amin1)/2+amin1;
    sum1=0;
    if n==2;
        intx=afx[1]*exp(afx[2].*y);
    endif;
    if n==3;
        intx=afx[1]*exp(afx[2].*y+afx[3].*y^2);
    endif;
    l=0;
    do until l==m;
        l=l+1;
        sum1=sum1+(maxc(amin1|y[l])-x)^2*intx[l]*sqrt(1-z[l]^2);
    endo;     
    sum1=sum1*pi/2*(amax1-amin1)/m;
    retp(sum1);
endp;

proc avar(afe,afu,x,abar);
    local n,m,z,y,sum1,y1e,y1u,inte,intu,i,l;
    n=rows(af);
    m=nint*n;  /* m>n quadrature points */
    z=-cos((2.*seqa(1,1,m)-1)./(2*m).*pi);
    y=(z+1).*(amax1-amin1)/2+amin1;
    sum1=0;
    y1e=lininter(a,aopt[.,1],y);  
    y1u=lininter(a,aopt[.,2],y);
    if n==2;  
        inte=afe[1]*exp(afe[2].*y);
        intu=afu[1]*exp(afu[2].*y);
    endif;
    if n==3;
        inte=afe[1]*exp(afe[2].*y+afe[3].*y^2);
        intu=afu[1]*exp(afu[2].*y+afu[3].*y^2);
    endif;

    l=0;
    do until l==m;
        l=l+1;
        sum1=sum1+pp[1,x]*(maxc(amin1|y1e[l])-abar)^2*inte[l]*sqrt(1-z[l]^2);
        sum1=sum1+pp[2,x]*(maxc(amin1|y1u[l])-abar)^2*intu[l]*sqrt(1-z[l]^2);
    endo;     
    sum1=sum1*pi/2*(amax1-amin1)/m;
    retp(sum1);
endp;


proc mom2(x);   /* moments computation for n=2 */
    local  n,m,z,y,sum1,sum2,afe,afu,inte,intu,l,res;
    n=rows(x)/2;
    m=nint*n;  /* m>n quadrature points */
    z=-cos((2.*seqa(1,1,m)-1)./(2*m).*pi);
    y=(z+1).*(amax1-amin1)/2+amin1;
    sum1=zeros(n,1);
    sum2=zeros(n,1);
    afe=x[1:n];
    afu=x[n+1:2*n];
    inte=afe[1]*exp(afe[2].*y+afe[3].*y^2);
    intu=afu[1]*exp(afu[2].*y+afu[3].*y^2);
    l=0;
    do until l==m;
         l=l+1;
         sum1[1]=sum1[1]+inte[l]*sqrt(1-z[l]^2);
         sum2[1]=sum2[1]+intu[l]*sqrt(1-z[l]^2);
         sum1[2]=sum1[2]+inte[l]*sqrt(1-z[l]^2)*maxc(y[l]|amin1);
         sum2[2]=sum2[2]+intu[l]*sqrt(1-z[l]^2)*maxc(y[l]|amin1);
         sum1[3]=sum1[3]+inte[l]*sqrt(1-z[l]^2)*(maxc(y[l]|amin1)-abare)^2;
         sum2[3]=sum2[3]+intu[l]*sqrt(1-z[l]^2)*(maxc(y[l]|amin1)-abaru)^2;
     endo;     
     sum1=sum1*pi/2*(amax1-amin1)/m;
     sum2=sum2*pi/2*(amax1-amin1)/m;
     res=x;
     res[1]=(pp1[1]-sum1[1])/pp1[1];
     res[2]=(pp1[2]-sum2[1])/pp1[2];
     res[3]=(abare*pp1[1]-sum1[2])/abare;
     res[4]=(abaru*pp1[2]-sum2[2])/abaru;
     res[5]=(avare*pp1[1]-sum1[3])/avare;
     res[6]=(avaru*pp1[2]-sum2[3])/avaru;
     retp(res);
endp;

