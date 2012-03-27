% VRFT method using a non-linear model
% real system: y(k)=0.9y(k-1)+tanh(u(t-1))
%========================================================================
close all; clear all;
clc;
P=path;
path(P,'../functions')
P=path;
path(P,'../functions/signal')
P=path;
path(P,'../functions/plots')
P=path;
path(P,'../functions/rational')

% simulation parameters
rho_size=2;
cut=1;
m=5;
a=0.7; b=0.9; c=0.3; d=0.2; e=0.5; 

Ts=1;
exper = 10;

%model.a = [1 -(b+e) b*e]; 
%model.b = [d -a*d];
%model.c = [1 -c];
%model.d = [1 0];
% Model M=1/z
model.mn = [1];
model.md = [1 0];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.0005;
C_den=[1 -1];

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
beta=[tf([1 0], C_den, model.TS); tf([1],C_den , model.TS)];

%Instrumental Variables
IV=beta*L*(inv(M)-1);

theta = zeros(exper, rho_size);
[u1 N]=f_get_prbs(m);
%u = ones(N, 1)*0.5;
u=u1*.2;
%========================================================================
% real system
%========================================================================
y=zeros(N, 1);
y(1)=0;
for k=2:N
    y(k)=0.2*y(k-1)^3+tanh(u(k-1)); 
end

if max(y) > 10000
    error('y divergiu');
end
    
for i = 1: exper
    [el y1] = f_get_vrft_e_nl(model, u, y);
    [el2 y2] = f_get_vrft_e_nl(model, u, y);
    ul=lsim(L,u);
    phy2=lsim(IV, y1);
    instr2=lsim(IV, y2);
    phy=phy2(cut:max(size(phy2)),:,1);
    instr=instr2(cut:max(size(instr2)),:,1);
    theta(i,:)=inv(instr'*phy)*instr'*ul(cut:max(size(ul)));
end
mtheta=mean(theta);
C=tf(mtheta,C_den, model.TS)
Jvr=f_get_vrft_Jvr(C, el, u(1:max(size(u))-1))
zero(C)
%Jmr=f_get_vrft_Jmr(C, model)
khsgdjg


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model parameter definition
m_rat.n_dim   = 6;
m_rat.dim     = 6;
% to indo do
m_rat.texp    = [1 1 1 1 1 1];
m_rat.yu      = [2 2 2 1 1 1];
m_rat.regr    = [0 1 2 1 2 3];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r = 3 u = 2 y=1 none =0
m_rat.yplus_yur =[0 0 0 0 0 0];
% tels the d param
m_rat.yplus_exp =[0 0 0 0 0 0];
% tels the C param
m_rat.yplus_regr = [0 0 0 0 0 0];

m_rat.err_m_rat   = 0;
m_rat.err_enable = true
m_rat.err_model = 0;
%% Simulation parameters
simul=struct('N', N, 'nEstimates', exper, 'np', model.noise_std, 'maxError', 0.001, 'l', 100, 'diffConv', .1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[e y] = f_get_vrft_el(model, u);
theta = f_rational_model(simul, m_rat, u, [u(1)], e);

variance =var(theta);
cda=0.8; cdb=0.9; cdc=0.5; cdd=1; cde=0.36; cdf=0.7;
expect=[cda -cda*(cdb+cdc) cda*cdb*cdc (cdd+cde+cdf) -(cdd*cde+cdd*cdf+cde*cdf) cdd*cde*cdf];
ret=mean(theta);
% controller with no filter
C=tf([ret(1) ret(2) ret(3) 0], [1 -ret(4) -ret(5) -ret(6)], model.TS);
% optimal controller 
Cd=zpk([0 cdb cdc],[cdd cde cdf],cda, 1);

% costs
Jvr=f_get_vrft_Jvr(C, e, u)
Jmr=f_get_vrft_Jmr(C, model)

G=tf(model.b,model.a, model.TS);
T=feedback(C*G, 1);
step(T, M)
legend;
grid;

figure
bode(C, Cd)
legend;

