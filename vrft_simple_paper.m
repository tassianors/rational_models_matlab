%==========================================================================
% Projeto - VRFT
% Tassiano Neuhaus
% tassianors@gmail.com
%==========================================================================
clear all; close all;clc;

% Sample time
a=.6;
b=.8;
c=.9;
k=10;

% Final time [s]
model.Ts=1;
model.N=500;
model.dim=4;
model.regr = [1 0 1 2];
model.eul= [0 1 1 1];

M=tf([1-c],[1 -c],model.Ts,'Variable','z');
G=tf([k 0],[1 -(a+b) a*b],model.Ts,'Variable','z');

% Time vector
t = 0:1:model.N-1;
u = (square(0.02*pi*t)');

% response of unknown plant to u input signal
% Controller output signal
yl=lsim(G, u, t);
%yl = f_get_wnoise(y, 0);
% get the signal rl whose generate the same yl, but considering M TF.
atraso=tf([1],[1 0], model.Ts);
W=1/M*atraso;
rl=lsim(W, yl, t);
rl=rl(2:max(size(rl)));
rl(max(size(yl)))=0;

% Controller input signal
el=rl-yl;

teta=f_calc_mmq_theta(model, u, el)
% to be used in graphic plotting
d=(1-c)/k;
c1=1;
c2=d;
c3=-d*(a+b);
c4=d*a*b;
C=[c1 c2 c3 c4]
