%% Aguirre 10.4
close all; clear all;
clc;
%% Real system - variables
a1=.3; a2=-2; a3=1.5; b1=.2;

%% model parameter definition
m_rat.n_dim   = 3;
m_rat.dim     = 4;
m_rat.texp    = [2 1 1 3];
m_rat.yu      = [1 1 2 1];
m_rat.regr    = [1 2 1 2];
m_rat.err_m_rat   = 0;
m_rat.err_enable = true;
%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 2, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);

% initial conditions
y=zeros(simul.N, 1);
u=ones(simul.N, 1);

%% Simulation of real system
for k=max(abs(m_rat.regr))+1:simul.N
    y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(1+b1*y(k-2)^3);
end
%% Rational model - get the rational m_rat estimative
ret = f_rational_model(simul, m_rat, y, [y(1) y(2)], u)
%yf_draw_elipse(ret(:,1), ret(:,2), a1, a2);
% f_draw_elipse(ret(:,3), ret(:,4), a3, b1);