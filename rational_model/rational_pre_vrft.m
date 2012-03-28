%% Model expansion
close all; clear all;
clc;
P=path;
path(P,'./functions')

%% Real system - variables
a=0.8; b=0.9; c=0.5; d=1; e=0.36; f=0.7;

Cd=zpk([0 b c],[d e f],a, 1);

m_rat.n_dim   = 6;
m_rat.dim     = 6;
% to indo do
m_rat.texp    = [1 1 1 1 1 1];
m_rat.yu      = [3 3 3 1 1 1];
m_rat.regr    = [0 1 2 1 2 3];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r =3 u = 2 y=1 none =0
m_rat.yplus_yur =[0 0 0 0 0 0];
% tels the d param
m_rat.yplus_exp =[0 0 0 0 0 0];
% tels the C param
m_rat.yplus_regr = [0 0 0 0 0 0];

m_rat.err_enable = true
m_rat.err_model = 0;

% initial conditions
[u N]=f_get_prbs(3);
y=lsim(Cd, u);


%% Simulation parameters
simul=struct('N', N, 'nEstimates', 5, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);

expected=[a -a*(b+c) a*b*c (d+e+f) -(d*e+d*f+e*f) d*e*f];

%% Rational model - get the rational m_rat estimative
ret = f_rational_model(simul, m_rat, [y(1)], y, zeros(size(y)), u', 0,0)
f_draw_elipse(ret(:,1), ret(:,2), a, a2);
f_draw_elipse(ret(:,3), ret(:,4), a3, b1);