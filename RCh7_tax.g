@ ------------------------------ RCh7_tax.G ----------------------------

    date: March 21, 2007

    author: Burkhard Heer

    solves the tax reform model in chapter 7.3, heer/maussner, 2nd ed.

    value function iteration, two-dimensional grid over k,n
        
    you need to adjust the save path before running the program

-------------------------------------------------------------------------@

library pgraph;
graphset;

/*
save path="d:\\buch\\prog\\";
*/

cls;
clear all;

@ ---------------------------

 parameter

---------------------------- @

h0=hsec;

beta=0.96;      /* discount factor */
r=0.07;         /* initial value of the interest rate */
sigma=2;        /* coefficient of relative risk aversion */
alpha=0.36;     /* production elasticity of capital */
delta=0.04;     /* rate of depreciation */
gam1=10;         /* disutility from working */
gam0=0.13;
tauc=0.15;      /* consumption tax, initial guess */
tauy=0.174;       /* income tax rate */

repl=0.52;       /* replacement ratio of unemployment insurance, lowest
           productivity */
gexpy=0.196;     /* government expenditure share in GDP */

eps=0|0.447559|0.785109|1.054437|1.712894;  /* productivities */

pp={0.35000  0.65000  0.00000  0.00000  0.00000,
    0.08000  0.67506  0.17023  0.03638  0.03833,
    0.08000  0.16514  0.51624  0.20025  0.03836,
    0.08000  0.04223  0.19946  0.52237  0.15593,
    0.08000  0.03706  0.03446  0.16059  0.68789};

nk=1000;
kmin=0;
kmax=12;
kstep=(kmax-kmin)/(nk-1);
k=seqa(kmin,kstep,nk);
nn=100;
nmin=0;
nmax=1;
nstep=(nmax-nmin)/(nn-1);
ngrid=seqa(nmin,nstep,nn);

nep=5;          /* number of relative productivity types */
nit=10;          /* number of iteration over value function */
nq=50;          /* number of iterations over k */

ngk=500;        /* minimum number of iterations for invariant distribution */
tol=0.0001;      /* percentage deviation of final solution */
tol1=1e-5;      /* percentage deviaton of value function */
tolg=0.00001;     /* tolerance for distribution function */

kritw=1; kritg=1;   
neg=-1e20;      /* initial value for value function */
phi1=0.9;       /* update: percentage of last iteration kbar,nbar */
phi=0.5;        /* update of v */
psi=0.001;      /* constant in utility function */



@ ---------------------------------

initialization

--------------------------------- @

/* initial guess for aggregate capital and employment */
kk0=3;
nn0=0.18;
kk1=kk0;
nn1=nn0;
lsuppl=0.3;
lstd=0;

wgini=0;
lgini=0;

w0=(1-alpha)*(kk0/nn0)^alpha;
y0=kk0^alpha*nn0^(1-alpha);
b0=repl*0.3*w0*eps[2];     /* unemployment compensation */

copt=zeros(nk,nep);     /* optimal consumption */
nopt=ones(nk,nep)*0.2;  /* optimal labor supply */
aopt=zeros(nk,nep);     /* optimal next period capital stock */
aoptn=ones(nk,nep);    /* grid point index of aopt */
pep1=equivec(pp);       /* invariant productivity distribution */

@ -----------------------------------------------

iteration until capital converges to steady state

------------------------------------------------- @

kbarq=zeros(nq,1);      /* convergence of kk0 */
k1barq=zeros(nq,1);      /* convergence of kk1 */

nbarq=zeros(nq,1);
kritwq=zeros(nq,1);     /* convergence of value function v */

ls1=zeros(nq,1);        /* average labor supply */
ls2=zeros(nq,1);        /* variational coefficient of labor */
wg1=zeros(nq,1);        /* wealth gini */
lg1=zeros(nq,1);        /* labor gini */

q=0;
do until q==nq or (abs((kk1-kk0)/kk0)<tol and kritw<tol1 and kritg<tolg); 
  q=q+1;
  w=(1-alpha)*kk0^alpha*nn0^(-alpha);
  r=alpha*kk0^(alpha-1)*nn0^(1-alpha);
  kbarq[q]=kk0;
  nbarq[q]=nn0;
 
  kk1=kk0;
  nn1=nn0;

  /* initialization of value function */
  if q==1;
    v=u((r-delta)*k+b0,0)/(1-beta);     /* value function of unemployed */
    j=1;
    do until j==nep;
      j=j+1;
      v=v~u((r-delta)*k+w*eps[j]*nn0,nn0)/(1-beta);     /* value function */
    endo;
  endif;

  /* value function iteration */

  it=0;
  kritw=100;
  if q==10; nit=nit*2; endif;
  do until it==nit or abs(kritw)<tol1;
    it=it+1;
    vnew=zeros(nk,nep);
    /* compuation of the decsion rules for the worker */
    l=0;
    do until l==nep;     /* productivity types */
      l=l+1;
      "q~it~l~kk0~nn0~kritw: " q~it~l~kk0~nn0~kritw;

    i=0;
    m0=0; j0=nn+1;
    do until i==nk;     /* asset holding in t */
    i=i+1;
    k0=k[i]; 
    v0=neg;
    m=m0;
    do until m==nk;  /* asset holding in t+1 */
      m=m+1;
      k1=k[m];
      j=j0;
      do until j==1;
        j=j-1;
        if l==1; j=1; endif;
        empl=ngrid[j];
        y=k0*(r-delta)+w*empl*eps[l];
        c=((1-tauy)*y+k0+b0*(l==1)-k1)/(1+tauc);
        if c>0 and empl<1;
            v1=u(c,empl)+beta*v[m,.]*pp[l,.]';
            if v1>v0;
               v0=v1;
               copt[i,l]=c;
               aopt[i,l]=k1;
               aoptn[i,l]=m;
               nopt[i,l]=empl;
               m0=maxc(0|m-2);  /* concavity of value function */
               j0=minc(nn|j+2);
            endif;
        endif;
        if i>1;         /* monotonocity of c in k */
            if c<copt[i-1,l];
            j=1;
            endif;
        endif;
      endo;   /* j=1,..,nn */
      endo; /* m=1,..,na */
      vnew[i,l]=v0;
      endo;     /* i=1,..,nk, asset holdings in t */
    endo;       /* l=1,..,nep, productivity type in t */

    kritw=meanc(meanc(abs(1-vnew.*(1/v))));
    kritwq[q]=kritw;
    v=phi*v+(1-phi)*vnew;
  endo;       /* it=1,.., */
  save ngrid,k,copt,v,nopt,aopt;

  /* iteration to find invariant distribution */
  q1=0;
  kritg=1;
  /* initialization of the distribution functions */
  if q<=10; 
      gk=ones(nk,nep)/nk;
      gk=gk.*pep1';
  endif; 
    if q==10; ngk=ngk*10; endif;
  "computation of invariant distribution of wealth..";
  do until q1>ngk or (ngk>50 and kritg<tolg);
    q1=q1+1;
    gk0=gk;
    gk=zeros(nk,nep);
    l=0;
    do until l==nep;
      l=l+1;
      i=0;
      do until i==nk;
        i=i+1;
        k1=aopt[i,l];
        if k1<=kmin;
            gk[1,1]=gk[1,1]+gk0[i,l]*pp[l,1];
            gk[1,2]=gk[1,2]+gk0[i,l]*pp[l,2];
            gk[1,3]=gk[1,3]+gk0[i,l]*pp[l,3];
            gk[1,4]=gk[1,4]+gk0[i,l]*pp[l,4];
            gk[1,5]=gk[1,5]+gk0[i,l]*pp[l,5];
        elseif k1>=kmax;
            gk[nk,1]=gk[nk,1]+gk0[i,l]*pp[l,1];
            gk[nk,2]=gk[nk,2]+gk0[i,l]*pp[l,2];
            gk[nk,3]=gk[nk,3]+gk0[i,l]*pp[l,3];
            gk[nk,4]=gk[nk,4]+gk0[i,l]*pp[l,4];
            gk[nk,5]=gk[nk,5]+gk0[i,l]*pp[l,5];
        elseif (k1>kmin) and (k1<kmax);
            j=sumc(k.<k1)+1;
            n0=(k1-k[j-1])/(k[j]-k[j-1]);
            gk[j,1]=gk[j,1]+n0*gk0[i,l]*pp[l,1];
            gk[j,2]=gk[j,2]+n0*gk0[i,l]*pp[l,2];
        gk[j,3]=gk[j,3]+n0*gk0[i,l]*pp[l,3];
        gk[j,4]=gk[j,4]+n0*gk0[i,l]*pp[l,4];
        gk[j,5]=gk[j,5]+n0*gk0[i,l]*pp[l,5];
        gk[j-1,1]=gk[j-1,1]+(1-n0)*gk0[i,l]*pp[l,1];
        gk[j-1,2]=gk[j-1,2]+(1-n0)*gk0[i,l]*pp[l,2];
        gk[j-1,3]=gk[j-1,3]+(1-n0)*gk0[i,l]*pp[l,3];
        gk[j-1,4]=gk[j-1,4]+(1-n0)*gk0[i,l]*pp[l,4];
        gk[j-1,5]=gk[j-1,5]+(1-n0)*gk0[i,l]*pp[l,5];
        endif;
      endo;
    endo;
    gk=gk/sumc(sumc(gk));
    kritg=sumc(abs(gk0-gk));
  endo;   /* q1=1,.., invariant distribution */

@
  if q==1 or q==nq;      /* test of the two first-order conditions */
    res1=ones(nk,nep);
    res2=zeros(nk,nep-1);
    j=0;
    do until j==nep;
      j=j+1;
      i=0;
      do until i==nk;
        i=i+1;
        k0=k[i];
        c0=copt[i,j];
        n0=nopt[i,j];
        y0=k0*(r-delta)+w*n0*eps[j]+b0*(j==1);
        w0=(1-tauy)*w;
        if j>1;
            res2[i,j-1]=(1+tauc)*dul(c0,n0)/(duc(c0,n0)*eps[j]*w0)-1;
        endif;
        kind=aoptn[i,j];
        k1=aopt[i,j];
        l=0;
        do until l==nep;
            l=l+1;
            c1=copt[kind,l];
            n1=nopt[kind,l];
            y1=k1*(r-delta)+w*n1*eps[l]+b0*(l==1);
            r1=1+(r-delta)*(1-tauy);
            res1[i,j]=res1[i,j]-duc(c1,n1)/duc(c0,n0)*pp[j,l]*beta*r1;
        endo; /* l=1,..,nep */
      endo; /* i=1,..,na */
    endo; /* j=1,..,nep */
    xlabel("capital stock");
    ylabel("");
    title("residual function for intertemporal foc");
    xy(k[10:nk-10],res1[10:nk-10,.]);
    title("residual function for intratemporal foc");
    xy(k[10:nk-10],res2[10:nk-10,.]);
  endif;
@
    /* average labor supply of all working agents */
    lsuppl=sumc(sumc(gk.*nopt))/(1-pep1[1]);
    "average labor supply: " lsuppl;
    /* average labor supply of the eps[2] for unemployment compensation */
    lsuppl2=sumc(gk[.,2].*nopt[.,2])/pep1[2];
    ls1[q]=lsuppl;
    /* computation of the standard deviation of labor supply */
    lstd=sumc(sumc(gk.*(nopt-lsuppl)^2))/(1-pep1[1]);
    "coefficient of variation of labor supply" sqrt(lstd)/lsuppl;
    ls2[q]=sqrt(lstd)/lsuppl;
    
    /* computation of the standard deviation of effective labor supply */
    elstd=sumc(sumc(gk[.,2:5].*((nopt[.,2:5]'.*eps[2:5]-nn0)')^2))/(1-pep1[1]);
    "coefficient of variation of effective labor supply" sqrt(elstd)/nn0;
    vel=sqrt(elstd)/nn0;

  /* computation of aggregate consumption and taxes, the mean capital stock and employment */
  taxes=0;
  ytaxes=0;
  cc0=0;
  kk0=0;
  nn0=0;
  mhour=0;  /* mean hours */
  i=0;
  do until i==nk;
    i=i+1;
    l=0;
    do until l==nep;
      l=l+1;
      kk0=kk0+gk[i,l]*k[i];
      nn0=nn0+gk[i,l]*nopt[i,l]*eps[l];
      mhour=mhour+gk[i,l]*nopt[i,l];
       y1=nopt[i,l]*eps[l]*w+(r-delta)*k[i];
      taxes=taxes+tauc*copt[i,l]*gk[i,l]+tauy*y1*gk[i,l];
      ytaxes=ytaxes+tauy*y1*gk[i,l];
      cc0=cc0+copt[i,l]*gk[i,l];
    endo;
  endo;
   mhour=mhour/(1-pep1[1]);
   k1barq[q]=kk0;
   tauc0=(gexpy*kk0^alpha*nn0^(1-alpha)+b0*pep1[2]-ytaxes)/cc0;
  b01=repl*lsuppl2*w*eps[2];
  "new capital stock~employment~mean hours~tauc0~ui: ";
  kk0~nn0~mhour~tauc0~b01;
  kk0=phi1*kk1+(1-phi1)*kk0;
  nn0=phi1*nn1+(1-phi1)*nn0;
  b0=phi1*b0+(1-phi1)*b01;
  tauc=phi1*tauc+(1-phi1)*tauc0;
  gexpy0=gexpy*kk0^alpha*nn0^(1-alpha);
/*
  save taxes0=ytaxes,gexpy0;
  save b0,nn0,kk0,aopt,copt,nopt,v;
  save gk,tauc,w,mhour;
*/  

/* computation of the gini coefficient of capital distribution */
gk1=sumc(gk');
j=1;
fj=zeros(nk,1);
fj[1]=gk1[1];
hj=zeros(nk,1);
hj[1]=fj[1]*k[1]/kk0;
gini=1-hj[1]*gk1[1];
do until j==nk;
  j=j+1;
  fj[j]=gk1[j];
  hj[j]=hj[j-1]+fj[j]*k[j]/kk0;
  gini=gini-(hj[j]+hj[j-1])*fj[j];
endo;


wgini=gini;
wg1[q]=wgini;
save wg1;

naly=500;
ylmax=w*eps[5];
yl=seqa(0,ylmax/(naly-1),naly);
gyl=zeros(naly,1);

l=1;
do until l==nep;
  l=l+1;
  i=0;
  do until i==nk;
    i=i+1;
    l1=nopt[i,l]*eps[l]*w;
    j=sumc(yl.<=l1)+1;
    n0=(l1-yl[j-1])/(yl[j]-yl[j-1]);
    gyl[j,1]=gyl[j,1]+n0*gk[i,l];
    gyl[j-1,1]=gyl[j-1,1]+(1-n0)*gk[i,l];
  endo;
endo;

/* nopt=0 */
gyl[1]=0;
gyl=gyl/sumc(gyl);

/* computation of the labor gini coefficient */
j=1;
ylbar=sumc(yl'*gyl);
fj=zeros(naly,1);
fj[1]=gyl[1];
hj=zeros(naly,1);
hj[1]=fj[1]*yl[1]/ylbar;
gini=1-hj[1]*gyl[1];
do until j==naly;
  j=j+1;
  fj[j]=gyl[j];
  hj[j]=hj[j-1]+fj[j]*yl[j]/ylbar;
  gini=gini-(hj[j]+hj[j-1])*fj[j];
endo;

lgini=gini;
lg1[q]=gini;
save lg1;
save blgini=lgini;

endo;   /* q=1,.. over kk0 */



"Parameters: ";
"gam0~gam1:";
gam0~gam1;
"w~b0~r:";
w~b0~r;
"kbar: ";
kbarq';
"kritw: ";
kritwq';
"nbarq: ";
nbarq';
wait;
    res1=ones(nk,nep);
    res2=zeros(nk,nep-1);
    j=0;
    do until j==nep;
      j=j+1;
      i=0;
      do until i==nk;
     i=i+1;
     k0=k[i];
     c0=copt[i,j];
     n0=nopt[i,j];
     y0=k0*(r-delta)+w*n0*eps[j]+b0*(j==1);
     w0=(1-tauy)*w;
     if j>1;
         res2[i,j-1]=(1+tauc)*dul(c0,n0)/(duc(c0,n0)*eps[j]*w0)-1;
     endif;
     kind=aoptn[i,j];
     k1=aopt[i,j];
     l=0;
     do until l==nep;
       l=l+1;
       c1=copt[kind,l];
       n1=nopt[kind,l];
       y1=k1*(r-delta)+w*n1*eps[l]+b0*(l==1);
       r1=1+(r-delta)*(1-tauy);
       res1[i,j]=res1[i,j]-duc(c1,n1)/duc(c0,n0)*pp[j,l]*beta*r1;
     endo; /* l=1,..,nep */
    endo; /* i=1,..,na */
   endo; /* j=1,..,nep */
    xlabel("capital stock");
    ylabel("");
    title("residual function for intertemporal foc");
    xy(k[10:nk-10],res1[10:nk-10,.]);
    title("residual function for intratemporal foc");
    xy(k[10:nk-10],res2[10:nk-10,.]);

"labor supply: ";
ls1';
"variational coefficient of labor supply: ";
ls2';
"wealth gini: ";
wg1';
"labor income gini: ";
lg1';

"aggregate capital stock: kbarq"; kbarq[1:q]';
"value function approximation: "; kritwq[1:q]';
"computational time:      "; etstr(hsec-h0);




@ ------------------------------------------------------

  Procedures:

  lininter - linear interpolation
  u(.) - utility function
  duc - marginal utility of consumption
  dul - marginal utility of leisure

  ------------------------------------------------------  @



proc lininter(xd,yd,x);
  local j;
  j=sumc(xd.<=x');
  retp(yd[j]+(yd[j+1]-yd[j]).*(x-xd[j])./(xd[j+1]-xd[j]));
endp;


proc equivec(P);
  local c,va,ve;
  {va,ve}=eigv(P');
  c=abs(ve[.,maxindc(va)]);
  retp(c./sumc(c));
endp;

proc u(c,l);
  local u;
  u=(c+psi)^(1-sigma)/(1-sigma)+gam0*(1-l)^(1-gam1)/(1-gam1);
  retp(u);
endp;

proc duc(c,l);  /* marginal utility of consumption */
  local duc;
  duc=(c+psi)^(-sigma);
  retp(duc);
endp;

proc dul(c,l);  /* marginal utility of labor */
  local dul;
  dul=gam0*(1-l)^(-gam1);
  retp(dul);
endp;
