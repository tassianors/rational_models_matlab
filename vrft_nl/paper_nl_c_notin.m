% Example used in our paper: Tuning nonlinear controllers with a virtual reference approach
% Alexandre Bazanella and Tassiano Neuhaus
% 2012 - version with plant integrator
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
rho_size=7;
cut=1;
m=2;

% static non linearity
a1=.5;a2=1;b1=.25;
mn=0.4;

Ts=1;
exper = 10;

model.mn = [mn];
model.md = [1 -(1-mn)];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.05;

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
if max(y) > 10000
    error('y divergiu');
end
%========================================================================
% Controller model definition
%========================================================================
%% model parameter definition
m_rat.n_dim   = 6;
m_rat.dim     = 7;
m_rat.err_model =0;
m_rat.texp    = [1 1 1 1 1 1 1];
m_rat.yu      = [4 3 4 3 1 1 3];
m_rat.regr    = [0 0 0 1 1 1 0];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.yplus_yur = [0 0 3 3 0 3 0];
% tels the d param
m_rat.yplus_exp = [0 0 1 1 0 1 0];
% tels the C param
m_rat.yplus_regr = [0 0 1 0 0 0 0];

m_rat.err_enable = true
m_rat.err_size = 1;
%% Simulation parameters
simul=struct('N', N-1, 'nEstimates', 1, 'np', model.noise_std,'l', 100, 'verbose', true);

%========================================================================
% vrft
%========================================================================
%theta0=[mn/a2 (1-mn)/a2 mn*b1/a2 (1-mn)*b1/a2 a1/a2];
theta0=zeros(1, rho_size);

for i = 1:exper
    [e yl rl] = f_get_vrft_e_nl(model, u, y);
    ul=lsim(L,u);
    [theta(i,:) cost]=f_rational_model(simul, m_rat, [u(1)], u(1:max(size(u))-1), e, y(1:max(size(y))-1), rl);
end

mtheta=mean(theta)
vartheta=var(theta)
stdtheta=std(theta)
covtheta=cov(theta)

y=zeros(N, 1);
r=ones(N, 1);
e=zeros(N, 1);
u=zeros(N, 1);
ur=zeros(N, 1);

for k=3:N
    e(k)=r(k)-y(k-1);
    ur(k-1)=(mtheta(1)*r(k-1)^m_rat.texp(1)+mtheta(2)*y(k-1)^m_rat.texp(2)+mtheta(3)*y(k-2)^m_rat.texp(3)*r(k-1)^m_rat.yplus_yur(3)+mtheta(4)*y(k-2)^m_rat.texp(4)*y(k-1)^m_rat.yplus_yur(4)+ur(k-1)^m_rat.texp(5)+mtheta(6)*ur(k-1)^m_rat.texp(6)*y(k-1)^m_rat.yplus_yur(6))/(1+mtheta(7)*y(k-1)^m_rat.texp(7));
    u(k-1)=ur(k-1);
    y(k)=(a1*u(k-2)*y(k-1)+a2*u(k-2))/(1+b1*y(k-2)^2);
end

ym=step(M, N-2);
stairs(ym(1:N-2)-y(3:N))

grid;
title('Erro entre a saida desejada e a obtida com o controlador estimado')
xlabel('t');
ylabel('Amplitude');

figure;
stairs(ym(1:N-2))
hold;
stairs(y(3:N))


Jmr_nl=f_get_vrft_nl_Jmr(M,r(1:N-2), y(3:N)')

f_draw_elipse(theta(:,1), theta(:,2), theta0(1), theta0(2));
f_draw_elipse(theta(:,3), theta(:,4), theta0(3), theta0(4));
