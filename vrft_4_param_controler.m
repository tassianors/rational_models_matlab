%% Aguirre 10.4
close all; clear all;
clc;
P=path;
path(P,'./functions')
%% Real system - variables
a1=.1; a2=-.5; a3=1; b1=.1;alpha=0.7;

%% model parameter definition
m_rat.n_dim   = 4;
m_rat.dim     = 4;
m_rat.texp    = [2 1 1 1];
m_rat.yu      = [1 1 1 1];
m_rat.regr    = [1 2 0 0];

% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r = 3 u = 2 y=1 none =0
m_rat.yplus_yur = [0 0 0 1];
% tels the d param
m_rat.yplus_exp = [0 0 0 3];
% tels the C param
m_rat.yplus_regr = [0 0 0 2];
m_rat.err_m_rat   = 0;
m_rat.err_enable = true
%% Simulation parameters
simul=struct('N', 300, 'nEstimates', 10, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);

% initial conditions
y=zeros(simul.N, 1);
%u=ones(simul.N, 1);
u=f_get_square_signal(simul.N);

%% Simulation of real system
for k=3:simul.N
    y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k))/(1+b1*y(k-2)^3);
end
stem(y)

%% VRFT parameters
model.Ts =1;
model.Tf = simul.N-1;
t=[0:model.Ts:model.Tf];
% M is the desired transfer function in Closed Loop
M=zpk([],[-alpha], 1+alpha, model.Ts);
atraso=tf([1],[1 0], model.Ts);

W=1/M*atraso;
rl_=lsim(W, y, t);
rl=rl_(2:size(rl_,1));
rl(size(rl_,1))=0;
% control input
el=rl-y;
% yt=y(1:size(y,1));
% yt(size(y,1))=0;
yt=zeros(size(y));
yt(2:size(y,1))=y(1:size(y,1)-1);

%% Rational model - get the rational m_rat estimative
ret = f_rational_model(simul, m_rat, yt, [y(1)], rl)
c1=-a1/a3;
c2=-a2/a3;
c3=1/a3;
c4=b1/a3;
CC=[c1 c2 c3 c4 ]
f_draw_elipse(ret(:,1), ret(:,2), c1, c2);
f_draw_elipse(ret(:,3), ret(:,4), c3, c4);