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
rho_size=5;
cut=1;
m=3;

% static non linearity
a1=0.5;a2=.2;b1=0.5;
b1=-0.9;b3=-0.2;

Ts=1;
exper = 1;

model.mn = [1];
model.md = [1 0];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.0005;

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);
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
m_rat.n_dim   = 4;
m_rat.dim     = 5;
% to indo do
m_rat.texp    = [1 1 1 2 1];
m_rat.yu      = [0 1 0 1 1];
m_rat.regr    = [0 0 0 1 0];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.yplus_uy = [0 0 1 1 0];
% tels the d param
m_rat.yplus_exp = [0 0 2 1 0];
% tels the C param
m_rat.yplus_regr = [0 0 1 0 0];

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
    theta = f_rational_model(simul, m_rat, yl(1:max(size(yl))-1), [yl(1)], rl)
end

nhghj
[e y] = f_get_vrft_el(model, u);


variance =var(theta);
cda=0.8; cdb=0.9; cdc=0.5; cdd=1; cde=0.36; cdf=0.7;
expect=[cda -cda*(cdb+cdc) cda*cdb*cdc (cdd+cde+cdf) -(cdd*cde+cdd*cdf+cde*cdf) cdd*cde*cdf];
ret=mean(theta);
% controller with no filter
C=tf([ret(1) ret(2) ret(3) 0], [1 -ret(4) -ret(5) -ret(6)], model.TS);
% optimal controller 
Cd=zpk([0 cdb cdc],[cdd cde cdf],cda, 1);

% costs
Jvr=f_get_vrft_Jvr(C, e, u)
Jmr=f_get_vrft_Jmr(C, model)

G=tf(model.b,model.a, model.TS);
T=feedback(C*G, 1);
step(T, M)
legend;
grid;

figure
bode(C, Cd)
legend;

