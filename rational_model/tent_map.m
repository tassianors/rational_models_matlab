%% Tent MAP
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

%% model parameter definition
model.n_dim   = 3;
model.dim     = 5;
model.texp    = [0 2 1 1 2];
model.yu      = [1 1 1 1 1];
model.regr    = [1 1 1 1 1];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r = 0 u = 2 y=1 none =0
model.yplus_yur = [0 0 0 0 0];
% tels the d param
model.yplus_exp = [0 0 0 0 0];
% tels the C param
model.yplus_regr = [0 0 0 0 0];
model.err_model   = 0;
model.err_enable   = true;

%% Simulation parameters
simul=struct('N', 500, 'nEstimates', 1, 'np', 0, 'maxError', 0.01, 'l', 100, 'diffConv', 0.100); 

%% Real system - variables
a=1.999; b=0.5;
theta0=[0.02608 -1.325 1.325 -2.416 2.416];

%Simulation of real system
y=zeros(simul.N, 1);
for k=max(abs(model.regr))+1:simul.N
    y(k)=1-a*abs(y(k-1)-b);
end

f_aguirre_plot_map(y, 1)

theta = f_rational_model(simul, model, [y(1)], y, zeros(size(y)), 0,0)

y2=f_y_model([0], y, 0, 0, theta, model);
y3=f_y_model([0], y, 0, 0, theta0, model);

f_aguirre_plot_map(y2, 2)
f_aguirre_plot_map(y3, 3)

