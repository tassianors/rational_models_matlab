% VRFT method using a non-linear model
% real system: y(k)=0.9y(k-1)+tanh(u(t-1))
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

% simulation parameters
rho_size=7;
cut=1;
m=3;

% static non linearity
a1=.5;a2=1;b1=.25;
mn=0.4;

Ts=1;
exper = 1;

model.mn = [mn];
model.md = [1 -(1-mn)];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0;

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);
%u=f_get_square_signal(N);
%========================================================================
% plant simul
%========================================================================
y=zeros(N, 1);
for k=3:N
    y(k)=(a1*u(k-1)*y(k-1)+a2*u(k-1))/(1+b1*y(k-2)^2);
end

plot(y)

if max(y) > 10000
    error('y divergiu');
end
%========================================================================
% Controller model definition
%========================================================================

%% model parameter definition
% model example 
%y(k) = (y(k-a1)^b1)*(y(k-c1)^d1)+...+(y(k-an)^bn)*(y(k-cn)^dn)+(y(k-ua1)^ub1)*(u(k-uc1)^ud1)+...+(y(k-uan)^ubn)*(u(k-ucn)^udn)
%       1+(y(k-an1)^bn1)*(y(k-cn1)^dn1)+...+(y(k-am)^bm)*(y(k-cm)^dm)+(y(k-ua1)^ub1)*(u(k-uc1)^ud1)+...+(y(k-uan)^ubn)*(u(k-ucn)^udn)
m_rat.n_dim   = 6;
m_rat.dim     = 7;
m_rat.err_model =0;
% to indo do
m_rat.texp    = [1 1 1 1 3 2 1];
m_rat.yu      = [2 3 3 2 3 3 3];
m_rat.regr    = [0 1 0 0 1 1 0];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.yplus_yur = [0 0 0 3 0 3 0];
% tels the d param
m_rat.yplus_exp = [0 0 0 2 0 1 0];
% tels the C param
m_rat.yplus_regr = [0 0 0 1 0 0 0];

m_rat.err_m_rat   = 0;
m_rat.err_enable = true
%% Simulation parameters
simul=struct('N', N-1, 'nEstimates', exper, 'np', model.noise_std, 'maxError', 0.001, 'l', 100, 'diffConv', .1);

%========================================================================
% vrft
%========================================================================
for i = 1: exper
    [e yl rl] = f_get_vrft_e_nl(model, u, y);
    ul=lsim(L,u);
    theta = f_rational_model(simul, m_rat, u(1:max(size(u))-1), [u(1)], e, y(1:max(size(y))-1))
end
ee=e;
yy=yl;
for k=2:N-1
    uu(k)=(mn/a2*ee(k)+mn/a2*yy(k-1)+(1-mn)/a2*yy(k)+mn*b1/a2*ee(k)*yy(k-1)^2+mn*b1/a2*yy(k-1)^3+(1-mn)*b1/a2*yy(k-1)^2*yy(k))/(1+a1/a2*yy(k));
end

theta0=[mn/a2 mn/a2 (1-mn)/a2 mn*b1/a2 mn*b1/a2 (1-mn)*b1/a2 a1/a2]
theta
u2=f_y_model([0], e, y, theta0, m_rat);
u3=f_y_model([0], e, y, theta, m_rat);

plot(u3-u2)

