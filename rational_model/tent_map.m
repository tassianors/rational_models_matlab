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
model.a_exp    = [0 2 1 1 2];
model.a_signal_type      = [1 1 1 1 1];
model.a_regress    = [1 1 1 1 1];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% r = 0 u = 2 y=1 none =0
model.b_signal_type = [0 0 0 0 0];
% tels the d param
model.b_exp = [0 0 0 0 0];
% tels the C param
model.b_regress = [0 0 0 0 0];
model.error_model_dim   = 0;
model.error_in_account   = true;

%% Simulation parameters
exper=10;
simul=struct('N', 500, 'nEstimates', 1, 'np', 0.001, 'l', 100, 'verbose', true); 

%% Real system - variables
a=1.999; b=0.5;
theta0=[0.02608 -1.325 1.325 -2.416 2.416];

%Simulation of real system
y=zeros(simul.N, 1);
for k=max(abs(model.a_regress))+1:simul.N
    y(k)=1-a*abs(y(k-1)-b);
end


for i = 1: exper
    [theta(i,:) cost] = f_rational_model(simul, model, [y(1)], y, zeros(size(y)), 0,0);
end

m_theta = mean(theta)
theta0
c_theta = cov(theta)
cost

y2=f_y_model([0], y, 0, 0, m_theta, model);
y3=f_y_model([0], y, 0, 0, theta0, model);

[a b]=f_aguirre_plot_map(y2, 0);
[c d]=f_aguirre_plot_map(y, 0);
[e f]=f_aguirre_plot_map(y3, 0);


plot(a, b, 'bo');
hold;
plot(c, d, 'g+');
% Create xlabel
xlabel({'t'});
% Create ylabel
ylabel({'t+1'});
grid;
figure(2);
plot(a, b, 'bo');
hold;
plot(e, f, 'g+');
% Create xlabel
xlabel({'t'});
% Create ylabel
ylabel({'t+1'});
grid;
