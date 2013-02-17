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
m=2

% static non linearity
a1=.5;a2=1;b1=.25;
mn=0.4;

Ts=1;
exper = 100

model.mn = [mn];
model.md = [1 -(1-mn)];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.005;

M=tf(model.mn,model.md, model.TS);

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);
%u=f_get_square_signal(N);
%========================================================================
% plant simul
%========================================================================
w=zeros(N, 1);
for k=3:N
    w(k)=(a1*u(k-1)*w(k-1)+a2*u(k-1))/(1+b1*w(k-2));
end
% we have a integrator after plant here
y=lsim(tf([1 0],[1 -1], model.TS), w);

%========================================================================
% Controller model definition
%========================================================================
%% model parameter definition
m_rat.n_dim            = rho_size-2;
m_rat.dim              = rho_size;
m_rat.error_model_dim  = 0;
m_rat.error_in_account = true

m_rat.a_exp            = [1 1 1 1 1 1 1 1];
m_rat.a_signal_type    = [4 4 4 3 3 3 3 3];
m_rat.a_regress        = [0 1 0 0 0 0 0 1];
m_rat.b_signal_type    = [0 0 3 0 3 3 0 0];
m_rat.b_exp            = [0 0 1 0 1 1 0 0];
m_rat.b_regress        = [0 0 2 0 1 2 0 0];

%% Simulation parameters
simul=struct('N', N-1, 'nEstimates', 1, 'np', model.noise_std,'l', 100, 'verbose', true);

%========================================================================
% vrft
%========================================================================
AA=0.4;BB=0.1;CC=0.5;
theta0=[AA BB -BB -AA -BB BB CC -CC];

for i = 1:exper
    [e yl rl] = f_get_vrft_e_nl(model, u, y);
    %    function [ret cost] = f_rational_model(simul, model, ic, out_sig, in_sig, aux_sig1, aux_sig2)
    [theta(i,:) cost]=f_rational_model(simul, m_rat, [u(1)], u(1:max(size(u))-1), e, y(1:max(size(y))-1), rl);
end

mtheta=mean(theta)
vartheta=var(theta);
stdtheta=std(theta);
covtheta=cov(theta)

y=zeros(N, 1);r=ones(N, 1);e=zeros(N, 1);u=zeros(N, 1);w=zeros(N, 1);
for k=f_model_get_max_regressor(m_rat)+2:N+f_model_get_min_regressor(m_rat)
    e(k)=r(k)-y(k-1);
    u(k)=f_y_model_k(k-1, e, y, r, mtheta,  m_rat);
    w(k)=(a1*u(k)*w(k-1)+a2*u(k))/(1+b1*w(k-2));
    y(k)=w(k)+y(k-1);
end

N_plot=20;
ym=step(M, N_plot-2);
stairs(ym(1:N_plot-2)-y(3:N_plot))

grid;
title('Erro entre a saida desejada e a obtida com o controlador estimado')
xlabel('Time');
ylabel('Amplitude');

figure;
stairs(ym(1:N_plot-2), 'r')
hold;
stairs(y(3:N_plot), 'g')

Jmr_nl=f_get_vrft_nl_Jmr(M,r(1:N-2), y(3:N)')

f_draw_elipse(theta(:,1), theta(:,2), theta0(1), theta0(2));
%f_draw_elipse(theta(:,3), theta(:,4), theta0(3), theta0(4));
