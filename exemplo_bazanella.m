%% Aguirre 10.4
close all; clear all;
clc
%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;
f=8.658;g=0.001223;h=-0.0441;r=-0.08381;s=0.001766;
a1=.3;a2=-2;a3=1.5;
b1=-.2;b2=1.3;
% auxiliary variables
% number of points on simulation
N=100;
% number of iteration on simulation
L=10;
% number of simulations
M=10;
%noise power
np=0.0005;
% initialization variables
y=zeros(N, 1);
nn_dim=3;
nd_dim=2;
nn_nd_dim=nn_dim+nd_dim;
theta=zeros(L+1,nn_nd_dim);
phy=zeros(nn_nd_dim,1);
psi=zeros(N,nn_nd_dim);
%% model
t_exp=  [0 3 2 1 2];
yu=     [1 1 1 1 1];
yu_regr=[1 1 1 1 1];
%%
v=zeros(L,1);
u=ones(N, 1);
y(1)=28;

for m=1:M
    
    %% Simulation of real system
    for k=2:N
        %y(k)=(a1*y(k-1)^3+a2*u(k-1)^2+a3*y(k-1))/(1+b1*y(k-1)^3+b2*y(k-1)^2);%+rand(1)*np;
        %y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(b-y(k-2)^3)+rand(1)*np;
        y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1))+rand(1)*np;
        %y(k)=(f+g*y(k-1)^3+h*y(k-1)^2)/(1+r*y(k-1)+s*y(k-1)^2)+rand(1)*np;
    end
   % plot_y_y1(y);
    plot(y,'.');
    psi = f_get_psi(y, t_exp, nn_nd_dim, nn_dim);
    %inv(psi'*psi)*psi'*y
    theta(1,:)=(psi'*psi)\(psi'*y);
    
    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1,nn_nd_dim);
    % we can't have a precision bigger than the noise power
    while (max(abs(err)) > np*2)
        yc=f_y_model(y(1), N, u, theta(l,:), nn_dim, t_exp, yu, yu_regr);
        %% step 2 -  calc the variance
        v(l)=cov(y-yc);
        %step 4
        PHY=f_get_PHY(nn_nd_dim, nn_dim, y, t_exp);
        
        for i=nn_dim+1:nn_nd_dim
            aux=0;
            for k=1:N
                aux=aux-(y(k)^t_exp(i));
            end
            phy(i,1)=aux;
        end
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
    f_plot_y_y1(yc);
end %J
f_draw_elipse(nda, ndb)

