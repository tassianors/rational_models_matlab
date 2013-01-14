%====================
%% DC to DC converter
% Example from Aguirre book, chapter 10.4
% 
%====================
close all; clear all;
clc;

% Add path to find functions
f_set_path();

%% model parameter definition
model.n_dim      = 3;
model.dim        = 5;
model.a_exp       = [0 3 2 1 2];
model.a_signal_type         = [1 1 1 1 1];
model.a_regress       = [1 1 1 1 1];
model.b_signal_type  = [0 0 0 0 0];
model.b_exp  = [0 0 0 0 0];
model.b_regress = [0 0 0 0 0];
model.error_model_dim  = 0;
model.error_in_account = true;

% number of Experiments
exper=100;
%% Simulation parameters
simul = struct('N', 200, 'nEstimates', 1, 'np', 0.005, 'verbose', true, 'l', 100); 

% variable init
rho_size=size(model.a_exp,2);
theta = zeros(exper, rho_size);
y = zeros(simul.N,1);
% Initial Condition
y(1)=31.5;

%% Real system - variables, from Aguirre
a=2.6204; b=99.875; c=1417.1; d=46.429;

% Simulation of real system
for k=max(abs(model.a_regress))+1:simul.N
    y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1));
end

for i = 1: exper
    [theta(i,:) cost] = f_rational_model(simul, model, [y(1)], y, y, 0, 0);
end

% get theta information: mean and covariance
m_theta = mean(theta);
c_theta = cov(theta);
cost;

% get output signal from model estimative
y_estim = f_y_model([y(1)], y, 0, 0, m_theta, model);

% plot map grafics: y x y+1
[a b] = f_aguirre_plot_map(y_estim, 0);
[c d] = f_aguirre_plot_map(y, 0);
plot(a, b, 'b.');
hold;
plot(c, d, 'r.');
% Plot covariance matriz to be used in latex
f_plot_matrix_std(c_theta);
