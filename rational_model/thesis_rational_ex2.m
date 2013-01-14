% Rational Model
% example extracted from : rational model identification - zhu and billings 1991
%===================================
close all; clear all;
clc;
% get PATH
f_set_path();
%start time measurement
tic;

% static non linearity - true system values
a1=.2;a2=0.1;a3=1;
b1=1;b2=1

% simulation parameters
rho_size=4;
% number of prbs signal periods
m=5;
exper = 100;
theta = zeros(exper, rho_size);

% get PRBS signal
[in_sig N]=f_get_prbs(m);
%===================================
% plant simul
%===================================
y = zeros(N, 1);
for k=3:N
    y(k)=(a1*y(k-1)+a2*y(k-1)*in_sig(k-1)+a3*in_sig(k-1)) / (in_sig1+b1*y(k-1)^2+b2*y(k-2)^2);
end
%===================================
% Controller model definition
%===================================
model.n_dim      = 2;
model.dim        = 4;
model.error_model_dim  = 0;
model.a_exp       = [1 1 2 2];
model.a_signal_type         = [1 2 1 1];
model.a_regress       = [1 1 1 2];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
model.b_signal_type  = [0 0 0 0];
% tels the d param
model.b_exp  = [0 0 0 0];
% tels the C param
model.b_regress = [0 0 0 0];
model.error_in_account = true

%% Simulation parameters
simul=struct('N', N, 'nEstimates', exper, 'np', 0.005, 'l', 100, 'verbose', true);

% calc theta #exper times
theta = f_rational_model(simul, model, [y(1)], y, in_sig, 0, 0);

% get theta information
m_theta = mean(theta);
v_theta = var(theta);
c_theta = cov(theta);

%% Compare results of real system and simulated
y1=zeros(N, 1);
for k=3:N
    y1(k)=(m_theta(1)*y1(k-1)+ m_theta(2)*in_sig(k-1)) /(1+ m_theta(3)*y1(k-1)^2+ m_theta(4)*y1(k-2)^2);
end

err=0;
for k=1:N
    err = err+abs(y(k)-y1(k))^2;
end

% get step response of the system
in_sig = ones(N,1);
y1 = zeros(N, 1);
y  = zeros(N, 1);
for k=3:N
    y(k)  = (a1*y(k-1)+ a2*y(k-1)*in_sig(k-1)+ a3*in_sig(k-1)) / (1+b1*y(k-1)^2+ b2*y(k-2)^2);
    y1(k) = (m_theta(1)*y1(k-1) + m_theta(2)*in_sig(k-1)) / (1+m_theta(3)*y1(k-1)^2+ m_theta(4)*y1(k-2)^2);
end
stairs(y(1:40), 'r');
hold;
stairs(y1(1:40), 'b')
% Stop time measurement
elapsed=toc;
