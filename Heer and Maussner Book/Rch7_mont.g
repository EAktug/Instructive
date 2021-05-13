@ ------------------------------ Rch7_mont.g ----------------------------


Date: March 19, 2007

Chapter 7.2., Monte Carlo Simulation

computation of the policy functions and the distribution function
in the heterogenous-agent neoclassical growth model
with value function iteration


MONTE CARLO SIMULATION

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

/*
load path="d:\\buch\\prog";
save path="d:\\buch\\prog";
*/

h0=hsec;


@ ----------------------------------------------------------------

Step 1: Parameterization and Initialization

----------------------------------------------------------------- @

/* numerical computation */
tol=0.001;          /* stopping criterion for final solution of K */
tol1=1e-7;          /* stopping criterion for golden section search */
tolg=0.00001;       /* distribution function */
neg=-1e10;
eps=0.05;

ni=10000;       /* number of individuals */

/* calibration of parameters */
alpha=0.36;
beta=0.995;
delta=0.005;
sigma=2;
r=0.005;
w=5;
tau=0.02;
rep=0.25;
b=rep*(1-tau)*w;
pp=(0.9565~0.0435|0.5~0.5);
/* stationary employment /unemployment */
pp1=equivec1(pp);


amin1=-2;                 /* asset grid */
amax1=3000;
na=201;
astep=(amax1-amin1)/(na-1);
a=seqa(amin1,astep,na);
nk=3*na;            /* asset grid for distribution */
agstep=(amax1-amin1)/(nk-1);
ag=seqa(amin1,agstep,nk);
gk=zeros(nk,2);     /* distribution */

/* initialization of the value function: */
v=u((1-tau)*r*a+(1-tau)*w)~u((1-tau)*r*a+b);
v=v/(1-beta);                /* agents consume their income */
copt=zeros(na,2);           /* optimal consumption */
aopt=zeros(na,2);           /* optimal next-period assets */

nit=1000;                /* number of maximum iterations over value function */
ngk=0;               /* initial number of iterations over distribution function */
crit1=1;
kritg=1;
nq=50;

psi=0.95;

kritw=1;
nn0=pp1[1];
kk0=(alpha/(1/beta-1+delta))^(1/(1-alpha))*nn0;
w0=(1-alpha)*kk0^alpha*nn0^(-alpha);
b=w0*rep;
kkbar=kk0;
kk1=kk0-1;
kbarq=zeros(nq,1);      /* convergence of kk0 */
k1barq=zeros(nq,1);      /* convergence of kk1 */
kritwq=zeros(nq,2);     /* convergence of value function v */
crit=1;


q=0;
do until q==nq or (abs((kk1-kk0)/kk0)<tol and q>50); /* iteration over K */
    q=q+1; 

    if q/2==round(q/2); cls; endif;
    w=(1-alpha)*kk0^alpha*nn0^(-alpha);
    r=alpha*kk0^(alpha-1)*nn0^(1-alpha)-delta;
    kbarq[q]=kk0;

    if ngk<25000;  
        ngk=ngk+500;
    endif;
    
    crit=1;
    j=0;
    do until (j>50 and crit<tol) or (j==nit);
        j=j+1;
        vold=v;
        "iteration j~q: " j~q;
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
    /*          and the moments nu1 (mean) und nu2 (standard deviation) */
    /* if q<5; */
    ae=zeros(ni,2);
    nie=round(ni*pp1[1]);
    ae[.,1]=(amax1-amin1)*seqa(1,1,ni)/ni;    /* uniform distribution */
    ae[.,1]=ones(ni,1)*kk0;
    aeneu=zeros(ni,2);
    i=0;
    do until i==ni;
        i=i+1;
        if rndu(1,1)<pp1[1];
            ae[i,2]=1;      
        else;
            ae[i,2]=2;        
        endif;
    endo;
    /* endif; */
    
    mu1=kk0;
    mu2=0;
    
    
    "computation of invariant distribution of wealth..";
    do until (q1>ngk);    
        q1=q1+1;
        p=rndu(ni,1);
        i=0;
        do until i==ni;  
            i=i+1;
            k0=ae[i,1];     /* wealth in period t */
            empl=ae[i,2];
            if k0<=amin1;
                 k1=aopt[1,empl];
            elseif k0>=amax1;
                 k1=aopt[na,empl];
            else;
                 k1=lininter(a,aopt[.,empl],k0);
            endif;
            aeneu[i,1]=k1;


            /* step 4 */ 
            if ae[i,2]==1;  /* employed in t */
                if p[i]<=pp[1,1];
                    aeneu[i,2]=1;                    
                else;
                    aeneu[i,2]=2;
                endif;
            else;
                if p[i]<=pp[2,1];
                    aeneu[i,2]=1;                    
                else;
                    aeneu[i,2]=2;
                endif;
            endif;
        endo;

        /* step 5: set of statistics */
        mu1neu=meanc(aeneu[.,1]);
        mu2neu=stdc(aeneu[.,1]);


        nround=ngk/100;
        if round(q1/nround)==q1/nround;
            "computation of invariant distribution of wealth..";
	        "iteration q1~q: " q1~q;
	        "kbarq~k1barq: ";
	        kbarq[1:q]~k1barq[1:q];
            "mu1~mu2: " mu1~mu2;
            "mu1neu~mu2neu: " mu1neu~mu2neu;
        endif;

        kritg=meanc(abs(mu1-mu1neu|mu2-mu2neu));
        mu1=mu1neu;
        mu2=mu2neu;
        
        nuempl=sumc(aeneu[.,2])-ni;    /* number of unemployed people */
        xuempl=nuempl-round(pp1[2]*ni);  
        if xuempl<0;
            iuempl=0;
            do until iuempl==-xuempl;
                ii=round(rndu(1,1)*ni);
                if aeneu[ii,2]==1;
                    aeneu[ii,2]=2;
                    iuempl=iuempl+1;
                endif;
            endo;       
        elseif xuempl>0;

            iempl=0;
            do until iempl==xuempl;
                ii=round(rndu(1,1)*ni);
                if aeneu[ii,2]==2;
                    aeneu[ii,2]=1;
                    iempl=iempl+1;
                endif;
            endo;       
        endif;            
        ae=aeneu;
        
    endo;   /* q1=1,.., invariant distribution */
    save fmont=ae;

        
    /* computation of the mean capital stock */
    kk1=mu1;
    k1barq[q]=kk1;

    save k1mon=k1barq;

    /* invariant distribution */
    ag=sortc(ae,1);
    age=delif(ag,ag[.,2].==2);  /* employed agents */
    agu=delif(ag,ag[.,2].==1);

    numbere=rows(age);
    numberu=rows(agu);
    age=zeros(ni-numbere-1,2)|age[.,1]~seqa(1/ni,1/ni,numbere)|amax1~pp1[1];    
    agu=zeros(ni-numberu-1,2)|agu[.,1]~seqa(1/ni,1/ni,numberu)|amax1~pp1[2];
    
    save agemont=age,agumont=agu;

    if q==1 or q==nq;
        xlabel("assets a");
        title("relative frequency");
        xy(age[.,1]~agu[.,1],age[.,2]~agu[.,2]);
    endif;

    /* update of the tax rate */
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
endo;   /* q=1,..,nq */

"iteration: " q;
"time elapsed: " etstr(hsec-h0);
"error value function: " crit;
"error distribution: " kritg;
"error capital stock: " crit1;
wait;

/* plotting the solution */
        xlabel("assets a");
        title("relative frequency");
        xy(age[.,1]~agu[.,1],age[.,2]~agu[.,2]);
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
