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
rho_size=8;
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
    y(k)=(a1*u(k-1)*y(k-1)+a2*u(k-1))/(1+b1*y(k-2));
end
if max(y) > 100
    error('y divergiu');
end
%========================================================================
% Controller model definition
%========================================================================
%% model parameter definition
m_rat.n_dim   = rho_size-2;
m_rat.dim     = rho_size;
m_rat.err_model =0;
m_rat.texp    = [1 1 1 1 1 1 1 1];
m_rat.yu      = [3 3 3 4 4 4 3 3];
m_rat.regr    = [0 1 2 0 1 2 0 1];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.yplus_yur = [0 3 3 0 3 3 0 0];
% tels the d param
m_rat.yplus_exp = [0 1 1 0 1 1 0 0];
% tels the C param
m_rat.yplus_regr = [0 0 0 0 0 0 0 0];

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
%    function [ret cost] = f_rational_model(simul, model, ic, out_sig, in_sig, aux_sig1, aux_sig2)
    [theta(i,:) cost]=f_rational_model(simul, m_rat, [u(1)], u(1:max(size(u))-1), e, y(1:max(size(y))-1), rl);
end

mtheta=mean(theta)
vartheta=var(theta)
stdtheta=std(theta)
covtheta=cov(theta)

y2=zeros(N, 1);y=zeros(N, 1);r=ones(N, 1)*2;e=zeros(N, 1);u=zeros(N, 1);ur=zeros(N, 1);

for k=4:N
    e(k)=r(k)-y2(k-1);
    ur(k-1)=(mtheta(1)*y2(k-m_rat.regr(1))^m_rat.texp(1)*y2(k-m_rat.yplus_regr(1))^m_rat.yplus_yur(1)+    mtheta(2)*y2(k-m_rat.regr(2))^m_rat.texp(2)*y2(k-m_rat.yplus_regr(2))^m_rat.yplus_yur(2)+    mtheta(3)*y2(k-m_rat.regr(3))^m_rat.texp(3)*y2(k-m_rat.yplus_regr(3))^m_rat.yplus_yur(3)+    mtheta(4)*r(k-m_rat.regr(4))^m_rat.texp(4)*y2(k-m_rat.yplus_regr(4))^m_rat.yplus_yur(4)+    mtheta(5)*r(k-m_rat.regr(5))^m_rat.texp(5)*y2(k-m_rat.yplus_regr(5))^m_rat.yplus_yur(5)+    mtheta(6)*r(k-m_rat.regr(6))^m_rat.texp(6)*y2(k-m_rat.yplus_regr(6))^m_rat.yplus_yur(6))/(1+    mtheta(7)*y2(k-m_rat.regr(7))^m_rat.texp(7)*y2(k-m_rat.yplus_regr(7))^m_rat.yplus_yur(7)+    mtheta(8)*y2(k-m_rat.regr(8))^m_rat.texp(8)*y2(k-m_rat.yplus_regr(8))^m_rat.yplus_yur(8));    u(k-1)=ur(k-1);
    y(k-1)=(a1*u(k-1)*y(k-1)+a2*u(k-1))/(1+b1*y(k-2));
    y2(k)=y2(k-1)-y(k);
end

N_plot=20;
ym=step(M, N-2);
stairs(ym(1:N_plot-2)-y2(3:N_plot))

grid;
title('Erro entre a saida desejada e a obtida com o controlador estimado')
xlabel('t');
ylabel('Amplitude');

figure;
stairs(ym(1:N_plot-2))
hold;
stairs(y2(3:N_plot))


Jmr_nl=f_get_vrft_nl_Jmr(M,r(1:N-2), y2(3:N)')

f_draw_elipse(theta(:,1), theta(:,2), theta0(1), theta0(2));
f_draw_elipse(theta(:,3), theta(:,4), theta0(3), theta0(4));
