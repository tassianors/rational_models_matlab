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
rho_size=11;
cut=1;
m=1;

% static non linearity
b1=1.5;b3=0.5;
ga1=0.5;
gb1=0.9;

Ts=1;
exper = 1;

model.a = [1 -gb1]; 
model.b = [ga1];
model.mn = [0.4];
model.md = [1 -0.6];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.0005;

M=tf(model.mn,model.md, model.TS);
G=tf(model.b,model.a, model.TS);
L=(1-M)*M;

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);
%========================================================================
% plant simul
%========================================================================
y2=lsim(G,u);
%% apply non linearity
for k=2:N
    y(k)=1*(b1*y2(k)+b3*y2(k)^3);
end

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
m_rat.n_dim   = 11;
m_rat.dim     = 11;
% to indo do
m_rat.texp    = [1 1 2 3 4 5 1 2 3 4 5];
m_rat.yu      = [1 2 2 2 2 2 2 2 2 2 2];
m_rat.regr    = [1 0 0 0 0 0 1 1 1 1 1];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r = 3 u = 2 y=1 none =0
m_rat.yplus_yur = [0 0 0 0 0 0 0 0 0 0 0];
% tels the d param
m_rat.yplus_exp = [0 0 0 0 0 0 0 0 0 0 0];
% tels the C param
m_rat.yplus_regr = [0 0 0 0 0 0 0 0 0 0 0];

m_rat.err_enable = true
%% Simulation parameters
simul=struct('N', N-1, 'nEstimates', exper, 'np', model.noise_std, 'maxError', 0.001, 'l', 100, 'diffConv', .1);

%========================================================================
% vrft
%========================================================================
for i = 1: exper
    [e y1] = f_get_vrft_e_nl(model, u, y');
    theta = f_rational_model(simul, m_rat, [u(1)], u(1:max(size(u))-1), e, 0, 0);
end

r=ones(1,N);
y=zeros(1,N);
e=zeros(1,N);
u=zeros(1,N);
u_hat=zeros(1,N);
for k=2:N
    e(k)=r(k)-y(k-1);
    u(k)=theta(1)*u(k-1)+theta(2)*e(k)+theta(3)*e(k)^2+theta(4)*e(k)^3+theta(5)*e(k)^4+theta(6)*e(k)^5+theta(7)*e(k-1)+theta(8)*e(k-1)^2+theta(9)*e(k-1)^3+theta(10)*e(k-1)^4+theta(11)*e(k-1)^5;
    u_hat(k+1)=gb1*u_hat(k)+ga1*u(k);
    y(k)=(b1*u_hat(k+1)+b3*u_hat(k)^3);
end
C=tf([0.8 -0.72],[1 -1],1)
stairs(y(2:N))
hold;
step(feedback(G*C,1))
