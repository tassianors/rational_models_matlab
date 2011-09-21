%==========================================================================
% Projeto - VRFT
% Tassiano Neuhaus
% tassianors@gmail.com
%==========================================================================
clear all; close all;clc;

% Sample time
alpha = 0.8;
beta = 0.9;

% Final time [s]
model.Ts=1;
model.N=100;
model.dim=2;
model.regr = [0 1];
model.eul= [1 0];

M=tf([1+alpha+beta+alpha*beta],[1+alpha+beta+2*alpha*beta (alpha+beta) 1],model.Ts,'Variable','z^-1');
G=tf([1+alpha],[alpha 1],model.Ts,'Variable','z^-1');

% Time vector
t = 0:1:model.N-1;
u = (square(0.02*pi*t)');

% response of unknown plant to u input signal
% Controller output signal
y=lsim(G, u, t);
yl = f_get_wnoise(y, 0);
% get the signal rl whose generate the same yl, but considering M TF.
W=1/M;
rl=lsim(W, yl, t);

% Controller input signal
el=rl-yl;

teta=f_calc_mmq_theta(model, u, el)
% to be used in graphic plotting
c1=(1+beta)/beta;
c2=-1/beta;
C=[c1 c2]
