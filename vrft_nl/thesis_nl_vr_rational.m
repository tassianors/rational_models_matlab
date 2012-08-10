% VRFT method using a non-linear model
% real system: y(k)=0.9y(k-1)+tanh(u(t-1))
%================================
close all; clear all;
clc;
P=path;
path(P,'../functions')
f_set_path();

% simulation parameters
theta_size = 5;
m          = 2;
exper      = 100;
cut        = 1;
%================================
% init variables
theta = zeros(exper, theta_size);
[u N] = f_get_prbs(m);
%================================
% Real system parameters
a1=.5; a2=1; b1=.25; mn=0.4;
theta0 = [mn/a2 (1-mn)/a2 mn*b1/a2 (1-mn)*b1/a2 a1/a2];
%================================
model.mn         = [mn];
model.md         = [1 -(1-mn)];
model.TS         = 1;
model.delay      = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std  = 0.1;

Td = tf(model.mn, model.md, model.TS);
%================================
% plant simul
%================================
y = zeros(N, 1);
for k=3:N
    y(k) = (a1*u(k-1)*y(k-1) + a2*u(k-1)) / (1+ b1*y(k-2)^2);
end
%================================
% Controller model definition
%================================
c_model.n_dim      = 4;
c_model.dim        = 5;
c_model.err_model  = 0;
c_model.texp       = [1 1 1 2 1];
c_model.yu         = [4 3 4 3 3];
c_model.regr       = [0 0 0 1 0];
c_model.yplus_yur  = [0 0 3 3 0];
c_model.yplus_exp  = [0 0 2 1 0];
c_model.yplus_regr = [0 0 1 0 0];
c_model.err_enable = true;
c_model.err_size   = 1;

%% Simulation parameters
simul = struct('N', N-1, 'nEstimates', 1, 'np', model.noise_std, 'l', 100, 'verbose', true);
%================================
% vrft
%================================
for i = 1:exper
    [e yl rl] = f_get_vrft_e_nl(model, u, y);
    [theta(i,:) cost] = f_rational_model(simul, c_model, [u(1)], u(1:max(size(u))-1), e, y(1:max(size(y))-1), rl);
end

mtheta   = mean(theta)
vartheta = var(theta)
stdtheta = std(theta)
covtheta = cov(theta)

% output signal at ideal controller
uc0 = f_y_model([0], e, y, rl, theta0, c_model);
uc  = f_y_model([0], e, y, rl, mtheta, c_model);

% VRFT cost of identification
Jvr_nl = f_get_vrft_nl_Jvr(uc0, uc)
%================================
% Plotting
y=zeros(N, 1); r=ones(N, 1); e=zeros(N, 1); u=zeros(N, 1);

for k=3:N
    e(k)   = r(k)-y(k-1);
    u(k-1) = (mtheta(1)*r(k-1) +mtheta(2)*y(k-1) +mtheta(3)*y(k-2)^2*r(k-1) +mtheta(4)*y(k-2)^2*y(k-1)) / (1 + mtheta(5)*y(k-1));
    y(k)   = (a1*u(k-1)*y(k-1) +a2*u(k-1))/ (1 + b1*y(k-2)^2);
end

ytd = step(Td, N-2);
stairs(ytd-y(2:N));
grid;
title('Erro entre a saida desejada e a obtida com o controlador estimado')
xlabel('t');
ylabel('Amplitude');

f_draw_elipse(theta(:,1),theta(:,2),theta0(1),theta0(2));
f_draw_elipse(theta(:,3),theta(:,4),theta0(3),theta0(4));
