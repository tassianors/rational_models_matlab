%% Aguirre 10.4
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
%% Real system - variables
a1=.3; a2=-2; a3=1.5; b1=.2;

%% model parameter definition
m_rat.n_dim   = 3;
m_rat.dim     = 4;
m_rat.a_exp    = [2 1 1 3];
m_rat.a_signal_type      = [1 1 2 1];
m_rat.a_regress    = [1 2 1 2];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.b_signal_type = [0 0 0 0];
% tels the d param
m_rat.b_exp = [0 0 0 0];
% tels the C param
m_rat.b_regress = [0 0 0 0];

m_rat.error_in_account = true;
%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 10, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);

% initial conditions
y=zeros(simul.N, 1);
u=ones(simul.N, 1);

%% Simulation of real system
for k=max(abs(m_rat.a_regress))+1:simul.N
    y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(1+b1*y(k-2)^3);
end
%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;
expected=[8.658 0.001223 -0.0441 -0.08381 0.001766];

for k=max(abs(m_rat.a_regress))+1:simul.N
    y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1));
end

%% Rational model - get the rational m_rat estimative
ret = f_rational_model(simul, m_rat, [y(1)], y, u, 0, 0)

std(ret)
f_draw_elipse(ret(:,1), ret(:,2), a1, a2);
% f_draw_elipse(ret(:,3), ret(:,4), a3, b1);