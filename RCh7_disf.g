@ ------------------------------ rch7_disf.g ----------------------------

Date: March 19, 2007

author: Burkhard Heer

Chapter 7.2., distribution function

computation of the policy functions and the distribution function
in the heterogenous-agent neoclassical growth model
with value function iteration

DISTRIBUTION FUNCTION

computation of the inidividual policy functions assumes:
- monotonocity of a' in a
- concavity of v

Linear interpolation between grid points

Maximization: golden section search method

input procedures:
golden -- routine for golden section search
lininter -- linear interpolation


-------------------------------------------------------------------------@


new;
clear all;
cls;
Macheps=1e-20;
library pgraph;
graphset;


/*
load path="C:\\Dokumente und Einstellungen\\ba6wp1\\Eigene Dateien\\buch\\prog\\ch521";
save path="C:\\Dokumente und Einstellungen\\ba6wp1\\Eigene Dateien\\buch\\prog\\ch521";
*/
h0=hsec;


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

/* numerical computation */
tol=0.001;          /* stopping criterion for final solution of K */
tol1=1e-7;          /* stopping criterion for golden section search */
tolg=0.0000001;       /* distribution function */
neg=-1e10;
eps=0.05;


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
amax1=3000;
na=201;
astep=(amax1-amin1)/(na-1);
a=seqa(amin1,astep,na);
nk=3*na;            /* asset grid for distribution */
agstep=(amax1-amin1)/(nk-1);
ag=seqa(amin1,agstep,nk);

nk0=sumc(ag.<=kk0);
gk=zeros(nk,2);     /* distribution */

/* initialization of the value function: */
v=u((1-tau)*r*a+(1-tau)*w)~u((1-tau)*r*a+b);
v=v/(1-beta);                /* agents consume their income */
copt=zeros(na,2);           /* optimal consumption */
aopt=zeros(na,2);           /* optimal next-period assets */

nit=100;                /* number of maximum iterations over value function */
ngk=0;               /* initial number of iterations over distribution function */
crit1=1;
kritg=1;
nq=70;

psi=0.95;

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
do until q==nq or (abs((kk1-kk0)/kk0)<tol and q>30 and ngk>=25000); /* iteration over K */
    q=q+1;
    w=(1-alpha)*kk0^alpha*nn0^(-alpha);
    r=alpha*kk0^(alpha-1)*nn0^(1-alpha)-delta;
    kbarq[q]=kk0;
    if ngk<25000;
        ngk=ngk+500;
    endif;

    
    kt=zeros(ngk,1);

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


    /* iteration to find invariant distribution */
    q1=0;
    kritg=1;

    /* Step 2: initialization of the distribution functions */
        /*
        gk=(ag-amin1)/(amax1-amin1);
        gk=gk.*pp1';
        */
        gk=zeros(nk,2);
        gk[nk0+1:nk,.]=ones(nk-nk0,2);
        gk=gk.*pp1';

    /* Step 3: computation of the inverse of a'(a,e) and a'(a,u) */

    "computation of a'^-1";
    /* elimination of the rows with the same entry as the subsequent row */
    aopte=delif(a~aopt[.,1],aopt[.,1].==aopt[2:na,1]|0);
    aoptu=delif(a~aopt[.,2],aopt[.,2].==aopt[2:na,2]|0);

    a1=zeros(nk,2); /* a'^(-1) at ag[i], i=1,..,nk and epsilon=e,u */
    l=0;
    do until l==2;  /* l=1 employed, l=2 unemployed in period t */
        l=l+1;
        i=0;
        do until i==nk; /* ag[i] -- wealth in t+1 */
            i=i+1;
            if l==1;
                if ag[i]<=aopte[1,2]; /* a'> ag[i] for all a */
                    a1[i,l]=a[1]-0.1;
                elseif ag[i]>=maxc(aopte[.,2]);     /* ag[i]>a' for all a */
                    a1[i,l]=maxc(aopte[.,1]);
                else;
                    a1[i,l]=lininter(aopte[.,2],aopte[.,1],ag[i]);
                endif;
            else;
                if ag[i]<=aoptu[1,2]; 
                    a1[i,l]=a[1]-0.1;
                elseif ag[i]>=maxc(aoptu[.,2]); 
                    a1[i,l]=maxc(aoptu[.,1]);
                else;
                    a1[i,l]=lininter(aoptu[.,2],aoptu[.,1],ag[i]);
                endif;
            endif;
        endo;
    endo;

    /* step 4: iteration of the distribution function */

    "computation of invariant distribution of wealth..";
	"iteration q: " q;
	"kbarq~k1barq: ";
	kbarq[1:q]~k1barq[1:q];

    q1=0;
    do until (q1>ngk);
        q1=q1+1; 

        gk0=gk;
        gk=zeros(nk,2);
        l=0;
        do until l==2;
            l=l+1;
            i=0;
            do until i==nk;
                i=i+1;
                k0=a1[i,l];
                l1=0;
                do until l1==2;
                    l1=l1+1;
                    if k0<=ag[1]; 
                        gk[i,l1]=gk[i,l1]+0;
                    elseif k0>=ag[nk];
                        gk[i,l1]=gk[i,l1]+pp[l,l1]*pp1[l];
                    else;
                        gk[i,l1]=gk[i,l1]+lininter(ag,gk0[.,l],k0)*pp[l,l1];
                    endif;
                endo;
            endo;
        endo;


        gk=gk/(gk[nk,1]+gk[nk,2]);
        kritg=sumc(abs(gk0-gk));
        nround=ngk/100;
        if round(q1/nround)==q1/nround;

            kkt=(gk[1,1]+gk[1,2])*ag[1];
            gk1=(gk[2:nk,1]+gk[2:nk,2])-(gk[1:nk-1,1]+gk[1:nk-1,2]);
            ag1=(ag[2:nk]+ag[1:nk-1])/2;
            kkt=kkt+gk1'*ag1;
            "iteration~capital stock: " q1~kkt;
            qt=q1/nround;
            kt[qt]=kkt;
        endif;

    endo;   /* q1=1,.., invariant distribution */

    
    /* computation of the mean capital stock */
    kk1=(gk[1,1]+gk[1,2])*ag[1];
    gk1=(gk[2:nk,1]+gk[2:nk,2])-(gk[1:nk-1,1]+gk[1:nk-1,2]);
    ag1=(ag[2:nk]+ag[1:nk-1])/2;
    kk1=kk1+gk1'*ag1;

    ag1=ag[1]|ag1;
    gk1=gk[1,.]|(gk[2:nk,.]-gk[1:nk-1,.]);

    k1barq[q]=kk1;


    tau0=b*pp1[2]/(w*nn0+r*kk0);
    tau=psi*tau+(1-psi)*tau0;

    crit1=abs((kk1-kk0)/kk0);
    kk0=psi*kk0+(1-psi)*kk1;

    
    "iteration: " q;
    "time elapsed: " etstr(hsec-h0);
    "error value function: " crit;
    "error distribution: " kritg;
    "error capital stock: " crit1;
    "kbarq~k1barq'";
    kbarq[1:q]~k1barq[1:q];
    save adis=a,aopt,a1,ag,r,kk0,w,b,q,nn0,tau;
    save aoptdis=aopt,vdis=v;
    save gkdis=gk,gk1dis=gk1,agdis=ag,ag1dis=ag1;
    save k1disf=k1barq;
    save kt;
    if q/2==round(q/2); cls; endif;


endo;   /* q=1,..,nq */

"iteration: " q;
"time elapsed: " etstr(hsec-h0);
"error value function: " crit;
"error distribution: " kritg;
"error capital stock: " crit1;
wait;



title("Mean of the distribution function");
xlabel("Number of iterations over distribution function F(.)");
xy(seqa(1,1,ngk),kt);
wait;

/* plotting the solution */
xlabel("asset level a");
title("invariant distribution");
_plegstr="employed\000unemployed";
xy(ag,gk);
wait;
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
