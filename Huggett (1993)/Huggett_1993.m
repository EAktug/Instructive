%% %Emrehan AKTUG 30.10.2016
clear
clc
format shortG

tic        %to calculate duration of process

amin=-2;    %min of grid
ainc=0.01;  %increment value=> if change, change gini coef part!!!
amax=5;     %max of grid
a=(amin:ainc:amax);     %asset grid
nk=size(a,2);  % size of asset grid

Y=[1 0.1];      %income in employment and unemployment state
nz=size(Y,2);   %size of states
Mark=[0.80 0.20; 0.50 0.50]; %Markov-Chain matrix
alpha=1.5;
beta=0.96;
A0=5;           %just to enter the loop, arbitrary
q0=0.99;        %just selected value between 1 and beta.

bigiter=0;      %to keep track of loop
incr=0.0001;    %increment for q, bond price

while (A0>0.0001 || A0<-0.0001) %0.0001 is tolerance value for asset holdings
    
    %% consumption
    c=zeros(nk,nk,nz);
    for j=1:nk
        for i=1:nk
            for m=1:nz %for every state
                c(i,j,m)=a(j)+Y(m)-q0*a(i); %from budget constraint
                %if a(i)>((a(j)+Y(m))/q0)
                %   c(i,j,m)=NaN;
                %end
            end
        end
    end
    
    %% utility
    c(c<=0)=NaN;        %negative consmption punished, for borrowing constraint
    u=(c.^(1-alpha))/(1-alpha);
    u(isnan(u))=-99999;
    
    %% Value Function
    v=zeros(nk,nz);     %initial guess
    vnext=zeros(nk,nz); %value function obtained as a result of our guess
    arule=zeros(nk,nz);   %to keep location->decision rule will be obtained
    viter=zeros(nk,nk,nz); %intermediate step, vnext will be formed by this matrix
    diff=5; %just to enter the loop
    iter=0; %keeping track of the process
    
    while diff>0.00001  %tolerance value for Value Function
        
        for t=1:nk
            for b=1:nk
                for z=1:nz
                    viter(t,b,z)=u(t,b,z); %for z^th value of state
                    for m=1:nz             %expectation operation
                        viter(t,b,z)=viter(t,b,z)+beta*Mark(z,m)*v(t,m);
                    end
                end
            end
        end
        
        for x=1:nk
            for m=1:nz
                [vnext(x,m),arule(x,m)]=max(viter(:,x,m),[],1);
            end
        end
        diff=max(abs(vnext-v));
        
        %% McQueen-Porteus Algorithm
        bhigh=(beta/(1-beta))*max(max(vnext-v));
        blow=(beta/(1-beta))*min(min(vnext-v));
        v=vnext+((blow+bhigh)/2);
        
        iter=iter+1;
        %if mod(iter,50)==0
        %    disp(['Iter is at ',num2str(iter),' difference is ',num2str(max(diff))])
        %end
    end
    
    %% Probability Matrix
    probs= (1/(2*nk))*ones(nk,2);   %uniform distribution is created
    new=zeros(nk,nz);
    new2=zeros(nk,nz);
    tol=1;              %just to enter the loop
    %{
for i=1:nk
   for z=1:nz
      new(arule(i,z),z)=new(arule(i,z),z)+probs(i,z)*Mark(z,z);
      new(arule(i,z),nz-z+1)=new(arule(i,z),nz-z+1)+probs(i,z)*Mark(z,nz-z+1);
   end
end
    %}
    while tol>0.00001   %tolerance for time-invariant distribution
        new2=zeros(nk,nz);
        for i=1:nk
            for z=1:nz
                for m=1:nz
                    new2(arule(i,z),m)=new2(arule(i,z),m)+probs(i,z)*Mark(z,m);
                end  %distribution is updated with movement of agents according to decision rule
            end
        end
        tol = max(max(abs(new2-probs)));    %until it is less than 0.0001
        probs = new2;
    end
    %% A0 is obtained
    rulemat=amin+arule.*((amax-amin)/nk);
    A0mat=probs.*rulemat;       %asset matrix is the probability times decision rule
    A0=sum(sum(A0mat));         %sum of asset matrix gives the total asset holdings
    
    %% update q0
    if A0<0         %not enough demand for bond, reduce price
        q0=q0-incr;
    elseif A0>0     %too much demand, increase price
        q0=q0+incr;
    end
    
    newA(bigiter+1,1)=A0;   %for displaying the result
    newA(bigiter+1,2)=q0;
    newA(bigiter+1,3)=incr;
    
    %changing increment size if stuck around 0
    if bigiter>2
        if newA(bigiter+1,1)*newA(bigiter,1)<0  %if oscillates around 0, reduce increment
            incr=incr/2;
        end
    end
    %% track iteration
    bigiter=bigiter+1;
    if mod(bigiter,10)==0
        disp(['Iter is at ',num2str(bigiter),' and q0 is ',num2str(q0), ' and A0 is ',num2str(A0)])
    end
    
end
toc        %end of main loop, so duration is calculated
disp(['Final:Iter is at ',num2str(bigiter),' and q0 is ',num2str(q0), ' and A0 is ',num2str(A0)])

%% Steady-States for each state
for z=1:nz
    s=1;
    for i=1:nk
        if a(arule(i,z))==a(i)
            disp(['SS is at ', num2str(a(i)), ' with z=', num2str(z)]);
            steady(s,z)=a(i); %I formed a matrix, which shows the points where a-grid and decision rule intersects
            s=s+1;
        end
    end
end

for i=1:nk
    if probs(i,1)<probs(i,2)
        disp(['Prob. of being employed is lower at asset', num2str(a(i))])
    end     %there is no such point, that is why it is always higher than or equal to prob. of being unemployed
end
%% total resource (Question-4)
for i=1:nk
    for z=1:nz
        wealth(i,z)=probs(i,z)*(Y(z)+a(i));
    end
end
Lorenz=[wealth(:) probs(:)];    %created for LorenzCurve
LorenzOrd=sortrows(Lorenz,1);   %sorted
LorenzFinal=cumsum(LorenzOrd);  %cumulative sum for graph
totwealth=sum(sum(wealth))      %total wealth is calculated

%% Gini Coeff (Question-4)
SharePeop=LorenzFinal(:,2)*100;
ShareWealth=LorenzFinal(:,1)*(100/totwealth);

% negative integral
area1=trapz(SharePeop(1:1180),ShareWealth(1:1180));
area1_1=trapz(SharePeop(1:1180),SharePeop(1:1180));
area1half=area1_1-area1; %area between 45d and Lorenz

%positive integral
area2=trapz(SharePeop(1181:end),ShareWealth(1181:end));
area2_1=trapz(SharePeop(1181:end),SharePeop(1181:end));
area2half=area2_1-area2;  %area between 45d and Lorenz

areatotal=area1half+area2half;
gini=areatotal*2/10000;

%% Lorenz Curve (Question-4)
hold all
plot(SharePeop,ShareWealth)
plot(SharePeop,SharePeop)
xlabel('Cumulative Share of People')
ylabel('Cumulative Share of Wealth')
legend('Lorenz Curve','45d line')

%% total welfare (Question-5)
for i=1:nk
    for z=1:nz
        welfare(i,z)=probs(i,z)*(vnext(i,z));
    end
end
welfareIns=sum(sum(welfare));
FullMat=Mark^100000;
c_aver=FullMat(1,1)*Y(1)+FullMat(1,2)*Y(2);
welfareFull=(c_aver^(1-alpha))/((1-beta)*(1-alpha));
welfareDiff=welfareFull-welfareIns;
disp(['Welfare difference is ',num2str(welfareDiff)]);
%% All Figures

% For Q-1, the Value Functions drawn for each state
%{
hold all
plot(a,vnext(:,1),'-.')
plot(a,vnext(:,2),'--')
legend('Employed','Unemployed')
title('Value Function over asset grid')
xlabel('Asset Grid')
ylabel('Value Function')
%}

% For Q-1, the decision rules drawn for each state

%{
hold all
plot(a,a(arule(:,1)),'-.')
plot(a,a(arule(:,2)),'--')
plot(a,a)
legend('Employed','Unemployed','45d line')
title('Decision Rule over asset grid')
xlabel('Asset Grid')
ylabel('Decision Rule')
%}

% For Q-2, the stationary distribution

%{
surf(probs)
xlabel('Employment Status (Binary)')
ylabel('Asset Grid')
zlabel('Probability')
%}

