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

init=23.25;
N_it=48;
step_size=0.25;

theta = zeros(exper, rho_size);
result = zeros(N_it,rho_size+2);

%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 1, 'np', 0.005, 'l', 100, 'verbose', false); 
y=zeros(simul.N,1);

%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;

for t=1:N_it
    y(1)=init;
    % Simulation of real system
    for k=max(abs(model.regr))+1:simul.N
        y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1));
    end
    
    for i = 1: exper
        [theta(i,:) cost]= f_rational_model(simul, model, [y(1)], y, y, 0, 0);
    end
    
    m_theta = mean(theta);
    result(t,3:size(result,2))=m_theta;
    result(t,2)=cost;
    result(t,1)=init;

%c_theta = cov(theta)
str = sprintf('%s \t %s \t %s \t & %s \t & %s \t & %s \t & %s \\\\ ',num2str(result(t,1),4), num2str(result(t,2),5), num2str(result(t,3)), num2str(result(t,4), '%10.5e'), num2str(result(t,5), '%10.5e'), num2str(result(t,6),'%10.5e'), num2str(result(t,7),'%10.5e'));
disp(str);
init=init+step_size;
end
