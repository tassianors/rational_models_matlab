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
P=path;
path(P,'../functions/vrft')

% simulation parameters
rho_size=4;
cut=1;
m=2;

% static non linearity
a1=.5;a2=1;b1=.25;
mn=0.4;

Ts=1;
exper = 100;

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
%========================================================================
% Controller model definition
%========================================================================

%% model parameter definition
% model example 
%y(k) = (y(k-a1)^b1)*(y(k-c1)^d1)+...+(y(k-an)^bn)*(y(k-cn)^dn)+(y(k-ua1)^ub1)*(u(k-uc1)^ud1)+...+(y(k-uan)^ubn)*(u(k-ucn)^udn)
%       1+(y(k-an1)^bn1)*(y(k-cn1)^dn1)+...+(y(k-am)^bm)*(y(k-cm)^dm)+(y(k-ua1)^ub1)*(u(k-uc1)^ud1)+...+(y(k-uan)^ubn)*(u(k-ucn)^udn)
m_rat.n_dim   = 3;
m_rat.dim     = 4;
m_rat.error_model_dim =0;
% to indo do
m_rat.a_exp    = [1 1 2 1];
m_rat.a_signal_type      = [4 4 3 3];
m_rat.a_regress    = [0 0 1 0];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.b_signal_type = [0 3 3 0];
% tels the d param
m_rat.b_exp = [0 2 1 0];
% tels the C param
m_rat.b_regress = [0 1 0 0];

m_rat.error_in_account = true
%% Simulation parameters
simul = struct('N', N-1, 'nEstimates', 1, 'np', model.noise_std, 'l', 100, 'verbose', true);

%========================================================================
% vrft
%========================================================================
%theta0=[mn/a2 (1-mn)/a2 mn*b1/a2 (1-mn)*b1/a2 a1/a2];
theta0=zeros(1, rho_size);

for i = 1:exper
    [e yl rl] = f_get_vrft_e_nl(model, u, y);
    ul=lsim(L,u);
    [theta(i,:) cost] = f_rational_model(simul, m_rat, [u(1)], u(1:max(size(u))-1), e, y(1:max(size(y))-1), rl);
end

mtheta=mean(theta)
vartheta=var(theta)
stdtheta=std(theta)
covtheta=cov(theta)

y=zeros(N, 1);
r=ones(N, 1);
e=zeros(N, 1);
u=zeros(N, 1);

for k=3:N
    e(k)=r(k)-y(k-1);
    %u(k-1)=(mtheta(1)*r(k-1)+mtheta(2)*y(k-1)+mtheta(3)*y(k-2)^2*r(k-1)+mtheta(4)*y(k-2)^2*y(k-1))/(1+mtheta(5)*y(k-1));
    u(k-1)=(mtheta(1)*r(k-1)+mtheta(2)*y(k-2)^2*r(k-1)+mtheta(3)*y(k-2)^2*y(k-1))/(1+mtheta(4)*y(k-1));
    y(k)=(a1*u(k-1)*y(k-1)+a2*u(k-1))/(1+b1*y(k-2)^2);
end

ym=step(M, N-2);
stairs(ym-y(2:N))

grid;
title('Erro entre a saida desejada e a obtida com o controlador estimado')
xlabel('t');
ylabel('Amplitude');

Jmr_nl=f_get_vrft_nl_Jmr(M,r(1:N-1), y(2:N)')

f_draw_elipse(theta(:,1), theta(:,2), theta0(1), theta0(2));
f_draw_elipse(theta(:,3), theta(:,4), theta0(3), theta0(4));
