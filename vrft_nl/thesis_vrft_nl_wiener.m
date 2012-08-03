% VRFT method using a non-linear model
% real system: y(k)=0.9y(k-1)+tanh(u(t-1))
%================================
close all; clear all;
clc;
f_set_path();

% simulation parameters
rho_size = 8;
m=1;
exper = 100;

% static non linearity - true system parameters
b1 = 1.5; b3 = 0.2; ga1 = 0.5; gb1 = 0.9;

% model structure
model.a          = [1 -gb1]; 
model.b          = [ga1];
model.mn         = [0.4];
model.md         = [1 -0.6];
model.TS         = 1;
model.delay      = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std  = 0.005;

Td     = tf(model.mn, model.md, model.TS);
G      = tf(model.b, model.a, model.TS);
Filter = tf([1],[1 -1], model.TS);
theta  = zeros(exper, rho_size);
[u N]  = f_get_prbs(m);
%================================
% plant simul
%================================
y2=lsim(G,u);
%% apply non linearity
for k=1:N
    y(k)=b1*y2(k)+b3*y2(k)^3;
end
%================================
% Controller model definition
%================================
%% model parameter definition
m_rat.n_dim      = rho_size;
m_rat.dim        = rho_size;
m_rat.texp       = [1 1 2 3 1 1 2 3];
m_rat.regr       = [0 0 0 0 1 1 1 1];
m_rat.yu         = [3 4 4 4 3 4 4 4];
m_rat.yplus_yur  = zeros(1,rho_size);
m_rat.yplus_exp  = zeros(1,rho_size);
m_rat.yplus_regr = zeros(1,rho_size);
m_rat.err_enable = true;

%% Simulation parameters
simul=struct('N', N-2, 'nEstimates', 1, 'np', model.noise_std, 'l', 100, 'verbose', true);

%================================
% vrft
%================================
for i = 1: exper
    % get virtual signals
    [e y1 rl]  = f_get_vrft_e_nl(model, u, y');
    % Apply filter
    el         = lsim(Filter, e);
    theta(i,:) = f_rational_model(simul, m_rat, [u(1)], u(1:max(size(u))-2), el(2:max(size(e))), rl(2:max(size(e))), y);
end

m_theta   = mean(theta);
std_theta = std(theta);
var_theta = var(theta);
cov_theta = cov(theta);

%================================
% Step simulation
%================================
r = ones(1,N); y = zeros(1,N); e = zeros(1,N);  u2 = zeros(1,N); w = zeros(1,N); 

for k=2:N
    % get controller input signal
    e(k) = r(k)-y(k-1);
    uu   = 0;
    for kk=1:rho_size
        if kk<=rho_size/2
            uu = uu+mtheta(kk)*e(k)^(kk);    
        else
            uu = uu+mtheta(kk)*e(k-1)^(kk-rho_size/2);                
        end
    end
    u2(k)  = u2(k-1)+uu;
    w(k+1) = gb1*w(k)+ga1*u2(k);
    y(k)   = (b1*w(k+1)+b3*w(k)^3);
end

% get VRFT costs
f_get_vrft_nl_Jmr(Td,r(1:N-1), y(2:N))
f_get_vrft_nl_Jvr(u, u2')

%================================
% Ploting
%================================
stairs(y(2:N), '-r')
hold;
step(Td);
ploting = step(Td);
v=zeros(1,N);
v2=zeros(1,N);

for k=1:N
    aux=0;
    for kk=1:rho_size/2
        aux=aux+(mtheta(kk)/0.8)*e(k)^(kk);    
    end
    v(k)=aux;
end
%================================
figure;
plot(v(2:N), w(3:N+1))
grid;
xlabel('v(t)');
ylabel('\omega(t)');
title('v(t) \times \omega(t)');
%================================
figure;
stairs(v(2:max(size(ploting))), 'r')
hold;
stairs(w(3:max(size(ploting))+1))
grid;
xlabel('t');
title('Resposta de v(t) e \omega(t) a um sinal do tipo degrau');
legend('v(t)','\omega(t)');


