@ ------------------------------ rch7_fig.g ----------------------------

Date: May 1, 2008

author: Burkhard Heer

computes the figures in Chapter 7, 2nd edition, Heer/Maussner


-------------------------------------------------------------------------@


new;
clear all;
cls;
library pgraph,user;
GraphSettings;
_plwidth=7;
_pxsci=5;

@
/* input from RCh7_dis.g */
graphsettings;
"Figure F7_1:"; wait;
load kt;
nt=sumc(kt.>10);
title("");
xtics(0,25000,5000,2000);
ylabel("K");
xlabel("Number of iterations over distribution function F(.)");
xy(seqa(250,250,100),kt[1:100]);
wait;
@

/* input from RCH7_hug.g */
graphsettings;
"Figure 7.5: " wait;
load ahug,aopthug;
_plctrl=0;
_plegctl={1 7 2 -1.8};
_pltype={6,1};
title("");
xlabel("Asset level");
ylabel("Next-period assets a'");
_plegstr="a'=a'(1.0,a)\000a'=a";
xy(ahug,aopthug[.,1]~ahug);
wait;

"Figure 7.6: " wait;
_plegstr="a'=a'(0.1,a)\000a'=a";
xy(ahug,aopthug[.,2]~ahug);
wait;

"Figure 7.7: " wait;
ylabel("Change in asset level a'-a");
_plegctl={1 7 -1.8 -1};
_pltype={6, 1};
_plegstr="e=1.0\000e=0.1";
xy(ahug,aopthug-ahug);
wait;

load hugag, huggk;

"Figure 7.8: " wait;
_plegctl={1 7 1 0.3};
_pltype={6, 1};
_plegstr="e=1.0\000e=0.1";
title("");
ylabel("Invariant distribution");
xy(hugag,cumsumc(huggk));
wait;


/* Monte-Carlo Simulation: RCh7_mont.g */
load age=agemont,agu=agumont;
nage=rows(age); age=age[1:nage-1,.];
nagu=rows(agu); agu=agu[1:nagu-1,.];

/* distribution function: RCh7_dis.g */
load adis,gkdis=gk1dis,agdis;

n=rows(agdis);
gmont=zeros(n,2);
age1=delif(age,age[.,1].==0);
agu1=delif(agu,agu[.,1].==0);
nmont=rows(age1)+rows(agu1);
da=(agdis[2]-agdis[1])/2;
i=0;
i0=0;
do until i==n;
    i=i+1;
    if sumc(age1[.,1].>=agdis[i]-da)>0;
        x=selif(age1,age1[.,1].>=agdis[i]-da); 
        if sumc(x.<=agdis[i]+da)>0;
            x=selif(x,x[.,1].<=agdis[i]+da);
            gmont[i,1]=rows(x);
        endif;
    endif;
    if sumc(agu1[.,1].>=agdis[i]-da)>0;
        x=selif(agu1,agu1[.,1].>=agdis[i]-da); 
        if sumc(x.<=agdis[i]+da)>0;
            x=selif(x,x[.,1].<=agdis[i]+da);
            gmont[i,2]=rows(x);
        endif;
    endif;
endo;

gmont=gmont/nmont;
@
title("frequency distribution: Monte Carlo");
xy(agdis,gmont);
wait;
@
mf=1/(agdis[2]-agdis[1]);
fkdis=gkdis*mf;
n1=sumc(agdis.<=1000);

/* Density Function from RCh7_den.g */
load aden,gkden,agden;
mf=1/(agden[2]-agden[1]);
fkden=gkden*mf;

@
Graphsettings;
"Figure F7_3:" wait;
_plegctl={1 7 600 0.002};
_plegctl=1; 
_pltype={6, 2};
_plegstr="employed\000unemployed";
title("");
xlabel("Individual wealth a");
ylabel("Density");
xy(agden[1:n1],fkden[1:n1,.]);
wait;
cls;
@

/* Input from RCh7_func.g */
load af=disfunc;
fkfunc=af[1,1]*exp(af[2,1].*agdis[1:n1]+af[3,1].*agdis[1:n1]^2)+af[1,2]*exp(af[2,2].*agdis[1:n1]+af[3,2].*agdis[1:n1]^2);

"Figure F7_4: "; wait;
_plegctl={1 5 500 0.0015};
_plegstr="Distribution\000Density\000Monte-Carlo\000Function Approx.";
_plctrl={10,15,0,0};
_pltype={6,6,6,3};
_pstype={2,8,0,0};
_psymsiz={3,3,0,0};
title("");
xlabel("Individual wealth a");
ylabel("Density");
xy(agdis[1:n1],(fkdis[1:n1,1]+fkdis[1:n1,2])~(fkden[1:n1,1]+fkden[1:n1,2])~(gmont[1:n1,1]+gmont[1:n1,2])*mf~fkfunc);
wait;


/*
title("wealth distribution");
xy(afunc,af[1,1]*exp(af[2,1].*afunc+af[3,1].*afunc^2)+af[1,2]*exp(af[2,2].*afunc+af[3,2].*afunc^2));
wait;
*/
Graphsettings;
"Figure F7_3:" wait;
_plegctl={1 9 600 0.002};
_pltype={6, 1};
_plegstr="employed\000unemployed";
title("");
xlabel("Individual wealth a");
ylabel("Density");
xy(agdis[1:n1],fkdis[1:n1,.]);
wait;

    graphsettings;
    load k1disf;
    nk=sumc(k1disf.>10);
    "Figure F7_2:"; wait;
    ylabel("K");
    xlabel("Number of iterations over K");
    xy(seqa(1,1,nk),k1disf[1:nk]);
    wait;





end;



/* GraphSettings: Changes Gauss' default initialization of graphic routines
**
** Usage: GraphSettings;
*/
  
proc(0)=GraphSettings;

 external matrix  _pltype, _ptitlht, _pnumht,
         _paxht, _pmcolor, _plwidth, _pcsel;
 external proc GraphSet;
 GraphSet;
   _pdate="";        @ do not plot date               @
  _pltype=6;         @ solid lines only               @
 _plwidth=4;         @ line width always 4            @
   _paxht=0.15;      @ size of axis labels            @        
  _pnumht=0.15;      @ size of axis numbering         @
 _ptitlht=0.20;      @ size of titles                 @
 _pmcolor={0,        @ color of axis: black           @
           0,        @ color of axis' numbers         @
           0,        @ color of x-axis label          @
           0,        @ color of y-axis label          @
           0,        @ color of z-axis label          @
           0,        @ color of title                 @
           0,        @ color of boxes                 @
           0,        @ color of date and author       @
           15};      @ background color:        white @
   _pcsel={0,        @ color of first line: black     @
           1,        @ color of second line: dark blue@
           2,        @ color of third line: green     @
           4,        @ color of fourth line: red      @
           3,        @ color of fifth line: light blue@
          10,        @ color of sixth line: light green@
          8};        @ color of seventh line: dark grey@
          
retp;
endp;           

