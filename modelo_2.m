%% Aguirre 10.4
close all; clear all;
clc
%% model
num_dim=3;
den_dim=1;
model_dim=num_dim+den_dim;
mod_texp    = [2 1 1 3];
mod_yu      = [1 1 0 1];
mod_regr    = [1 2 1 2];
%%
% number of points in simulation
N=100;
% number of simulations
M=10;
%noise power
np=0.0005;
%% initialization variables
y=zeros(N, 1);
u=ones(N, 1);
% initial conditions
y(1)=rand(1)*0.1;
y(2)=rand(1)*0.1;

%% Real system - variables
a1=.3;
a2=-2;
a3=1.5;
b1=-.2;

for m=1:M
    %% Simulation of real system
    for k=max(abs(mod_regr))+1:N
        y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(1-b1*y(k-2)^3)+rand(1)*np;
    end

    psi = f_get_psi(y, u, model_dim, num_dim,  mod_texp, mod_yu, mod_regr);
    theta(1,:)=(psi'*psi)\(psi'*y);
    
    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1+np*2, model_dim);
    % we can't have a precision bigger than the noise power
    while (max(abs(err)) > np*2)
        yc=f_y_model([y(1) y(2)], u, theta(l,:), num_dim, mod_texp, mod_yu, mod_regr);
        if l ==1
            f_plot_y_y1(yc);
        end
        %% step 2 -  calc the variance
        v(l)=cov(y-yc);
        [PHY phy]=f_get_phy(y, model_dim, num_dim, mod_texp, mod_yu, mod_regr);
        
        theta(l+1,:)=(psi'*psi-v(l)*PHY)\ (psi'*y-v(l)*phy);
        delta(l,:)=theta(l+1,:)-theta(l,:);
        err(1,:)=delta(l);
        % to be used in graphic plotting
        nna(m)=theta(l+1,1);
        nnb(m)=theta(l+1,2);
        nnc(m)=theta(l+1,3);
        nda(m)=theta(l+1,4);
        l=l+1;
    end
    theta
    delta
    v'
    %f_plot_y_y1(yc);
end %J

f_draw_elipse(nna, nnb)
f_draw_elipse(nnc, nda)
