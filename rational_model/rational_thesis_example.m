% Rational Model example to be used inside my thesis
% example extracted from : rational model identification - zhu and billings 1991
%========================================================================
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
P=path;
path(P,'../functions/vrft')

% simulation parameters
rho_size=5;
m=5;

% static non linearity
a1=.2;a2=0.1;a3=1;
b1=1;b2=1

Ts=1;
exper = 10;
theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);

%========================================================================
% plant simul
%========================================================================
y=zeros(N, 1);
for k=3:N
    y(k)=(a1*y(k-1)+a2*y(k-1)*u(k-1)+a3*u(k-1))/(1+b1*y(k-1)^2+b2*y(k-2)^2);
end
if max(y) > 10000
    error('y divergiu');
end
%========================================================================
% Controller model definition
%========================================================================
m_rat.n_dim   = 3;
m_rat.dim     = 5;
m_rat.err_model =0;
m_rat.texp    = [1 1 1 2 2];

m_rat.yu      = [1 1 2 1 1];
m_rat.regr    = [1 1 1 1 2];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.yplus_yur = [0 2 0 0 0];
% tels the d param
m_rat.yplus_exp = [0 1 0 0 0];
% tels the C param
m_rat.yplus_regr = [0 1 0 0 0];

m_rat.err_enable = true
m_rat.err_size = 1;
%% Simulation parameters
simul=struct('N', N, 'nEstimates', exper, 'np', 0.005, 'maxError', 0.01, 'l', 100, 'diffConv', .1);

theta = f_rational_model(simul, m_rat, [y(1)], y, u, 0, 0);

m_theta = mean(theta)
v_theta = var(theta)
c_theta = cov(theta)

y1=zeros(N, 1);
for k=3:N
    y1(k)=(m_theta(1)*y1(k-1)+m_theta(2)*y1(k-1)*u(k-1)+m_theta(3)*u(k-1))/(1+m_theta(4)*y1(k-1)^2+m_theta(5)*y1(k-2)^2);
end

err=0;
for k=1:N
    err = err+abs(y(k)-y1(k))^2;
end

f_draw_elipse(theta(:,1), theta(:,2), a1, a2);
