% VRFT method using betas as a transfer function array
% example of luciola thesys page 92 
% where the controller function isn't inside controller class

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
rho_size=3;
cut=1;
m=10;
a=0.7; b=0.9; c=0.3; d=0.2; e=0.5; 

Ts=1;
exper = 100;
%% BL system: 
% G_0(z)=d(z-a)/((z-b)(z-e))
% H_0(z)=z/(z-c)
% M(z)=0.36z/(z-0.6)^2

%% Ideal Controler
% C_0(z)=0.8z(z-0.9)(z-0.5)
%       (z-1)(z-0.36)(z-0.7)
% y(t)=G_0(z)*u(t)+H_0(z)*e(t)

model.a = [1 -(b+e) b*e]; 
model.b = [d -a*d];
model.c = [1 -c];
model.d = [1 0];
model.mn = [0.16 0];
model.md = [1 -1.2 0.36];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.0005;
C_den=[1 -1 0];

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
beta=[tf([1 0 0], C_den, model.TS); tf([1 0],C_den , model.TS);tf([1],C_den , model.TS)];

%Instrumental Variables
IV=beta*L*(inv(M)-1);

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);

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
% r =3 u = 2 y=1 none =0
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
theta = f_rational_model(simul, m_rat, [u(1)], u, e, 0 , 0);

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

