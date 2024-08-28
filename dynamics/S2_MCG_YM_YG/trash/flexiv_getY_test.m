clear;clc;
load flexiv_get_tau.mat
load flexiv_getY.mat
load T_flexiv.mat

%% 
betas = 0.5*ones(70,1);
betas(61) = 0;   %xx7  % XXR7=I7xx - I7yy;
betas(64) = 0;   %yy7
betas(62) = 0;   %xy7
betas(63) = 0;   %xz7
betas(65) = 0;   %yz7
betas(67) = 0;   %mx7
betas(68) = 0;   %my7



beta =  T_37_by_70*betas;

L1=0.365;
L2=0.395;
L3=0.385;
L4=0.11;
L5=0.02;
g=-9.81;

cnt = 1;
beta1  = beta(cnt);cnt=cnt+1;
beta2  = beta(cnt);cnt=cnt+1;
beta3  = beta(cnt);cnt=cnt+1;
beta4  = beta(cnt);cnt=cnt+1;
beta5  = beta(cnt);cnt=cnt+1;
%beta6  = beta(cnt);cnt=cnt+1;
beta7  = beta(cnt);cnt=cnt+1;
beta8  = beta(cnt);cnt=cnt+1;
beta9  = beta(cnt);cnt=cnt+1;
beta10 = beta(cnt);cnt=cnt+1;
beta11 = beta(cnt);cnt=cnt+1;
%beta12 = beta(cnt);cnt=cnt+1;
beta13 = beta(cnt);cnt=cnt+1;
beta14 = beta(cnt);cnt=cnt+1;
beta15 = beta(cnt);cnt=cnt+1;
beta16 = beta(cnt);cnt=cnt+1;
beta17 = beta(cnt);cnt=cnt+1;
%beta18 = beta(cnt);cnt=cnt+1;
beta19 = beta(cnt);cnt=cnt+1;
beta20 = beta(cnt);cnt=cnt+1;
beta21 = beta(cnt);cnt=cnt+1;
beta22 = beta(cnt);cnt=cnt+1;
beta23 = beta(cnt);cnt=cnt+1;
%beta24 = beta(cnt);cnt=cnt+1;
beta25 = beta(cnt);cnt=cnt+1;
beta26 = beta(cnt);cnt=cnt+1;
beta27 = beta(cnt);cnt=cnt+1;
beta28 = beta(cnt);cnt=cnt+1;
beta29 = beta(cnt);cnt=cnt+1;
beta30 = beta(cnt);cnt=cnt+1;
beta31 = beta(cnt);cnt=cnt+1;
beta32 = beta(cnt);cnt=cnt+1;
beta33 = beta(cnt);cnt=cnt+1;
beta34 = beta(cnt);cnt=cnt+1;
beta35 = beta(cnt);cnt=cnt+1;
beta36 = beta(cnt);cnt=cnt+1;
%beta37 = beta(cnt);cnt=cnt+1;
beta38 = beta(cnt);cnt=cnt+1;
beta39 = beta(cnt);cnt=cnt+1;
beta40 = beta(cnt);cnt=cnt+1;
beta41 = beta(cnt);cnt=cnt+1;
beta42 = beta(cnt);%cnt=cnt+1;
%beta43 = beta(cnt);cnt=cnt+1;


%% get_tau
q = ones(7,1);
dq = zeros(7,1);
ddq = zeros(7,1);

q1 = q(1);
dq1 = dq(1);
ddq1 = ddq(1);
q2 = q(2);
dq2 = dq(2);
ddq2 = ddq(2);
q3 = q(3);
dq3 = dq(3);
ddq3 = ddq(3);
q4 = q(4);
dq4 = dq(4);
ddq4 = ddq(4);
q5 = q(5);
dq5 = dq(5);
ddq5 = ddq(5);
q6 = q(6);
dq6 = dq(6);
ddq6 = ddq(6);
q7 = q(7);
dq7 = dq(7);
ddq7 = ddq(7);

tau = flexiv_get_tau(L2,L3,L5,beta1,beta2,beta3,beta4,beta5,beta7,beta8,beta9,beta10,beta11,...
    beta13,beta14,beta15,beta16,beta17,beta19,beta20,beta21,beta22,beta23,beta25,beta26,...
    beta27,beta28,beta29,beta30,beta31,beta32,beta33,beta34,beta35,beta36,beta38,beta39,...
    beta40,beta41,beta42,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,g,q2,q3,q4,q5,q6);
Y = flexiv_getY(L2,L3,L5,ddq1,ddq2,ddq3,ddq4,ddq5,ddq6,ddq7,dq1,dq2,dq3,dq4,dq5,dq6,dq7,g,q2,q3,q4,q5,q6);

[M,cg,tau2,Ys]=flexiv_get_M_cg_tau_Ys(q,dq,ddq,betas);

%%
norm(Y-Ys*T_pinv_70_by_37)


%% cal differences
Y*beta-Ys*betas
tau-Y*beta
tau-tau2
tau2-(M*ddq+cg)









