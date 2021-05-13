@ ------------------------------ RCh7_hug.g ----------------------------

Chapter 7.3, risk-free rate of return
Author: Burkhard Heer

Date: August 2004

computes Huggett (1993), calibration from Huggett

computation of the policy functions and the distribution function
in the heterogenous-agent neoclassical growth model
with value function iteration

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
graphset;
h0=hsec;


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

/* computational parameters */
tolw=1e-6;              /* stopping criterion for value function */
tol=0.001;  /* stopping criterion for absolute divergence in capital stock */
krita=1;
tol1=1e-7;            /* stopping criterion for golden section search */
tolg=1e-10;         /* distribution function */
neg=-1e10;
delr=0.01;          /* if excess demand is 1, the interest rate is decreased by zeta */
eps=0.05;
nit=250;                /* number of maximum iterations over value function */
ngk=5000;              /* number of iterations over distribution */
crit=1;
nq=100;             /* number of iterations over r */                   
psi=0.5;
kritg=1;

/* parameter values */
beta=0.99322;
sigma=1.5;

/* initial guess for q,r,equilibrium assets */
q0=1;
r=1/q0-1;
r01=0.01;
a0=1;   
/* endowment */
eh=1;
el=0.1;
endow=eh|el;
pp=(0.925~0.075|0.5~0.5);

/* ergodic distribution */
pp1=equivec1(pp);

amin1=-2;                 /* asset grid */
amax1=4;
na=101;
astep=(amax1-amin1)/(na-1);
a=seqa(amin1,astep,na);

nk=3*na;            /* asset grid for distribution */
agstep=(amax1-amin1)/(nk-1);
ag=seqa(amin1,agstep,nk);

/* initialization of the value function: */
v=u(r*maxc(a'|zeros(1,na))+eh)~u(r*maxc(a'|zeros(1,na))+el);
v=v/(1-beta);                /* agents consume their income */
copt=zeros(na,2);           /* optimal consumption */
aopt=zeros(na,2);           /* optimal next-period assets */


nn0=pp1[1];

abarq=zeros(nq,1);      /* convergence of kk0 */
rbarq=zeros(nq,1);      /* convergence of kk1 */
delrq=zeros(nq,1);
kritwq=zeros(nq,2);     /* convergence of value function v */
crita=1;

q=0;
do until q==nq or (q>30 and kritg<tolg and crit<tolw and crita<tol); 
    q=q+1;
    cls;
    if q==2;
        r=r01;
    endif;

    "iteration: " q;
    "time elapsed: " etstr(hsec-h0);
    "error value function: " crit;
    "error a: " crita;
    "error distribution function: " kritg;
    "a0: " a0;
    "r: " r;

    crit=1;
    j=0;
    do until (j>100 and crit<tol) or (j==nit);
        j=j+1;
        q~j;
        abarq[1:q]';
        rbarq[1:q]';
        vold=v;


@ ----------------------------------------------------------------

Step 2: Iteration of the value function

----------------------------------------------------------------- @

        e=0;    /* iteration over the endowment */
        do until e==2;  /* endow[1]=eh, endow[2]=el */
            e=e+1;
            i=0;        /* iteration over asset grid a in period t */
            l0=0;
            do until i==na;
                i=i+1; 
                l=l0;
                v0=neg;
                ax=a[1]; bx=a[1]; cx=a[na];
                do until l==na; /* iteration over a' in period t*1 */
                    l=l+1; 
                    c=a[i]+endow[e]-a[l]/(1+r);
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


    copt[.,1]=a+endow[1]-aopt[.,1]/(1+r);
    copt[.,2]=a+endow[2]-aopt[.,2]/(1+r);


    /* iteration to find invariant distribution */
    q1=0;
    kritg=1;
    /* initialization of the distribution functions */
    if q<=10; 
            gk=ones(nk,2)/nk;
            gk=gk.*pp1';
    endif; 
    if q==10; ngk=5*ngk; endif;
    "computation of invariant distribution of wealth..";
    do until (q1>ngk);
        q1=q1+1; q~q1~kritg';
        if q>1;
             abarq[q-1];
        endif;
        gk0=gk;
        gk=zeros(nk,2);
        l=0;
        do until l==2;
            l=l+1;
            i=0;
            do until i==nk;
                i=i+1;
                k0=ag[i];
                if k0<=amin1; 
                    k1=aopt[1,l]; 
                elseif k0>=amax1;
                    k1=aopt[na,l];
                else;    
                    k1=lininter(a,aopt[.,l],k0);
                endif;
                if k1<=amin1;
                    gk[1,1]=gk[1,1]+gk0[i,l]*pp[l,1];
                    gk[1,2]=gk[1,2]+gk0[i,l]*pp[l,2];
                elseif k1>=amax1;
                    gk[nk,1]=gk[nk,1]+gk0[i,l]*pp[l,1];
                    gk[nk,2]=gk[nk,2]+gk0[i,l]*pp[l,2];
                elseif (k1>amin1) and (k1<amax1);
                    j=sumc(ag.<=k1)+1;
                    n0=(k1-ag[j-1])/(ag[j]-ag[j-1]);
                    gk[j,1]=gk[j,1]+n0*gk0[i,l]*pp[l,1];
                    gk[j,2]=gk[j,2]+n0*gk0[i,l]*pp[l,2];
                    gk[j-1,1]=gk[j-1,1]+(1-n0)*gk0[i,l]*pp[l,1];
                    gk[j-1,2]=gk[j-1,2]+(1-n0)*gk0[i,l]*pp[l,2];
                endif;
            endo;
        endo;

        gk=gk/sumc(sumc(gk));
        kritg=sumc(abs(gk0-gk));
    endo;   /* q1=1,.., invariant distribution */


a0=(gk[.,1]+gk[.,2])'*ag;

crita=abs(a0);

    abarq[q]=a0;
    rbarq[q]=r;
    delrq[q]=delr;
    if q>1;
	    r1=rbarq[q];
  	    r0=rbarq[q-1];
  	    delr=(r1-r0)/(abarq[q]-abarq[q-1]);
 	    r0=r-delr*a0;
	    r=psi*r+(1-psi)*r0;
    endif;

    "iteration: " q;
    "time elapsed: " etstr(hsec-h0);
    "error value function: " crit;
    "error distribution: " kritg;
    "error assets: " crita;
endo;   /* q=1,..,nq */



"iteration: " q;
"time elapsed: " etstr(hsec-h0);
"error value function: " crit;
"error distribution: " kritg;
"error assets: " crita;
wait;

/* plotting the solution */
xlabel("asset level");
_plegstr="eh\000el";
title("value function");
xy(a,v);
xy(a,v[.,1]);
xy(a,v[.,2]);

title("consumption");
xy(a,copt);

xy(a,copt[.,1]);
xy(a,copt[.,2]);

title("change in asset level");
xy(a,aopt-a);
wait;

xy(a,aopt[.,1]-a);
xy(a,aopt[.,2]-a);
"Figures 7.5 and 7.6";
wait;

save ahug=a,aopthug=aopt;
_plctrl=0;
_plegctl=1;
title("");
xlabel("Asset level");
ylabel("Next-period assets a'");
_plegstr="a'=a'(1.0,a)\000a'=a";
xy(a,aopt[.,1]~a);
wait;

_plegstr="a'=a'(0.1,a)\000a'=a";
xy(a,aopt[.,2]~a);
wait;

save hugag=ag, huggk=gk;

title("stationary distribution function");
xy(ag,cumsumc(gk));
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
   c=a0+endow[y]-a1/(1+r);
   if c<0;
      retp(neg);
   endif;
   if a1>=a[na];
      retp(a1^2*neg);
   endif;
   if a1<a[1];
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
