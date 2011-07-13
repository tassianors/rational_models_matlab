%% Aguirre 10.4
close all; clear all;
clc
%% model
num_dim=3;
den_dim=2;
model_dim=num_dim+den_dim;
mod_texp    = [0 3 2 1 2];
mod_yu      = [1 1 1 1 1];
mod_regr    = [1 1 1 1 1];
%%
% number of points on simulation
N=100;
% number of iteration on simulation
L=10;
% number of simulations
M=10;
%noise power
np=0.0005;
%% initialization variables
y=zeros(N, 1);
u=ones(N, 1);
v=zeros(L,1);
theta=zeros(L+1,model_dim);

y(1)=28;

%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;
f=8.658;g=0.001223;h=-0.0441;r=-0.08381;s=0.001766;
a1=.3;a2=-2;a3=1.5;
b1=-.2;b2=1.3;

for m=1:M
    
    %% Simulation of real system
    for k=2:N
        %y(k)=(a1*y(k-1)^3+a2*u(k-1)^2+a3*y(k-1))/(1+b1*y(k-1)^3+b2*y(k-1)^2);%+rand(1)*np;
        %y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(b-y(k-2)^3)+rand(1)*np;
        y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1))+rand(1)*np;
        %y(k)=(f+g*y(k-1)^3+h*y(k-1)^2)/(1+r*y(k-1)+s*y(k-1)^2)+rand(1)*np;
    end

    psi = f_get_psi(y, u, model_dim, num_dim,  mod_texp, mod_yu, mod_regr);
    theta(1,:)=(psi'*psi)\(psi'*y);
    
    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1+np*2, model_dim);
    % we can't have a precision bigger than the noise power
    while (max(abs(err)) > np*2)
        yc=f_y_model(y(1), u, theta(l,:), num_dim, mod_texp, mod_yu, mod_regr);
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
        ndb(m)=theta(l+1,5);
        l=l+1;
    end
    theta
    delta
    v(l-1)
    %f_plot_y_y1(yc);
end %J

f_draw_elipse(nda, ndb)
