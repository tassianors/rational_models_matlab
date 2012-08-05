% VRFT method example extracet and adapted from campestrini Doc. Thesys page 92 
% where the ideal controller function isn't inside controller class
%===================================

close all; clear all;
clc;
P=path;
path(P,'../functions')
P=path;
path(P,'../functions/signal')
P=path;
path(P,'../functions/plots')
P=path;
path(P,'../functions/vrft')

%===================================
% simulation parameters
rho_size=3;
m=20;
exper = 100;
cut=1;
%===================================
% variables initialization
theta   = zeros(exper, rho_size);
thetaL  = zeros(exper, rho_size);
[u N]   = f_get_prbs(m);

% Real system parameters
a=0.7; b=0.9; c=0.3; d=0.2; e=0.5;
%===================================
model.a          = [1 -(b+e) b*e]; 
model.b          = [d -a*d];
model.c          = [1 -c];
model.d          = [1 0];
model.mn         = [0.16 0];
model.md         = [1 -1.2 0.36];
model.TS         = 1;
model.delay      = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std  = 0.1;

C_den = [1 -1 0];
M     = tf(model.mn, model.md, model.TS);
L     = (1-M)*M;
beta  = [tf([1 0 0], C_den, model.TS); tf([1 0],C_den , model.TS); tf([1],C_den , model.TS)];
%Instrumental Variables
IV=beta*L*(inv(M)-1);

for i = 1: exper
    [el y]   = f_get_vrft_el(model, u);
    [el2 y2] = f_get_vrft_el(model, u);
    ul       = lsim(L,u);
    phi2     = lsim(IV, y);
    instr2   = lsim(IV, y2);
    phi      = phi2(cut: max(size(phi2)),:,1);
    instr    = instr2(cut: max(size(instr2)),:,1);
    thetaL(i,:) = inv(instr'*phi)*instr'* ul(cut:max(size(ul)));
end

for i = 1: exper
    [el y] = f_get_vrft_el(model, u);
    phi2   = lsim(beta, el);
    phi    = phi2(cut:max(size(phi2)),:,1);
    theta(i,:) = inv(phi'*phi)*phi'* u(cut:max(size(u)));
end

%===================================
variance  = var(theta);
varianceL = var(thetaL);

% controller with no filter
C=tf(mean(theta),C_den, model.TS)
% controller with L filter and IV 
CL=tf(mean(thetaL),C_den, model.TS)
% ideal controller 
Cd=zpk([0 0.9 0.5],[1 0.36 0.7],0.8, model.TS);

% costs
Jvr  = f_get_vrft_Jvr(C, el, u)
Jmr  = f_get_vrft_Jmr(C, model)
JvrL = f_get_vrft_Jvr(CL, el, u)
JmrL = f_get_vrft_Jmr(CL, model)

%===================================
% graphic ploting
G  = tf(model.b,model.a, model.TS);
T  = feedback(C*G, 1);
TL = feedback(CL*G, 1);

step(T, TL, M)
legend; grid; figure;
bode(C, CL, Cd)
legend;

