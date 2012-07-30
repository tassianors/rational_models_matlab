%% CC-CC BUCK
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

%% model parameter definition
model.n_dim   = 3;
model.dim     = 5;
model.texp    = [0 3 2 1 2];
model.yu      = [1 1 1 1 1];
model.regr    = [1 1 1 1 1];
model.yplus_yur = [0 0 0 0 0];
model.yplus_exp = [0 0 0 0 0];
model.yplus_regr = [0 0 0 0 0];
model.err_model   = 0;
model.err_enable = true;
model.err_size = 1;

exper=100;
rho_size=5;
theta = zeros(exper, rho_size);

%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 1, 'np', 0.005, 'verbose', true, 'l', 100); 
y=zeros(simul.N,1);
y(1)=31.5;%24.5;

%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;

% Simulation of real system
for k=max(abs(model.regr))+1:simul.N
    y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1));
end

for i = 1: exper
    [theta(i,:) cost] = f_rational_model(simul, model, [y(1)], y, y, 0, 0);
end

m_theta = mean(theta)
c_theta = cov(theta)
cost

y3=f_y_model([y(1)], y, 0,0,m_theta, model);
[a b] = f_aguirre_plot_map(y3, 0);
[c d] = f_aguirre_plot_map(y, 0);

plot(a, b, 'b.');
hold;
plot(c, d, 'r.');
f_plot_matrix_std(c_theta)



