%% Aguirre 10.4
close all; clear all;
clc;
%% Real system - variables
a1=.3; a2=-2; a3=1.5; b1=.2;

%% model parameter definition
model.n_dim   = 3;
model.dim     = 4;
model.texp    = [2 1 1 3];
model.yu      = [1 1 0 1];
model.regr    = [1 2 1 2];
model.err_model   = 0;
model.err_enable = true;
%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 30, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);
y=zeros(simul.N, 1);
% initial conditions
y(1)=0;
y(2)=0;
u=ones(1, simul.N);

%% Simulation of real system
for k=max(abs(model.regr))+1:simul.N
    y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(1+b1*y(k-2)^3);
end
ret = f_rational_model(simul, model, y, [y(1) y(2)], u);
f_draw_elipse(ret(:,1), ret(:,2), a1, a2);
f_draw_elipse(ret(:,3), ret(:,4), a3, b1);