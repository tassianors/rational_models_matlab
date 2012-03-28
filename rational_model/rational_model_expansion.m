%% Model expansion
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
a1=.3; a2=-.5; a3=1.5; b1=.2;

%% model parameter definition
% model example 
%y(k) = (y(k-a1)^b1)*(y(k-c1)^d1)+...+(y(k-an)^bn)*(y(k-cn)^dn)+(y(k-ua1)^ub1)*(u(k-uc1)^ud1)+...+(y(k-uan)^ubn)*(u(k-ucn)^udn)
%       1+(y(k-an1)^bn1)*(y(k-cn1)^dn1)+...+(y(k-am)^bm)*(y(k-cm)^dm)+(y(k-ua1)^ub1)*(u(k-uc1)^ud1)+...+(y(k-uan)^ubn)*(u(k-ucn)^udn)
m_rat.n_dim   = 3;
m_rat.dim     = 4;
% to indo do
m_rat.texp    = [2 1 1 3];
m_rat.yu      = [1 1 2 1];
m_rat.regr    = [1 2 1 2];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r = 3 u = 2 y=1 none =0
m_rat.yplus_yur = [0 1 2 0];
% tels the d param
m_rat.yplus_exp = [0 1 1 0];
% tels the C param
m_rat.yplus_regr = [0 1 1 0];

m_rat.err_enable = true
%% Simulation parameters
simul=struct('N', 300, 'nEstimates', 50, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);

% initial conditions
y=zeros(simul.N, 1);
u=ones(simul.N, 1);
%u=f_get_square_signal(simul.N);

%% Simulation of real system
for k=3:simul.N
    y(k)=(a1*y(k-1)^2+a2*y(k-2)*y(k-1)+a3*u(k-2)*u(k-1))/(1+b1*y(k-2)^3);
end
stem(y)
%% Rational model - get the rational m_rat estimative
ret = f_rational_model(simul, m_rat, [y(1) y(2)], y, u, 0, 0)
f_draw_elipse(ret(:,1), ret(:,2), a1, a2);
f_draw_elipse(ret(:,3), ret(:,4), a3, b1);