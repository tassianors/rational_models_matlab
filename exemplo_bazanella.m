%% Aguirre 10.4
close all; clear all;
clc;
%% model parameter definition
model.n_dim   = 3;
model.dim     = 5;
model.texp    = [0 3 2 1 2];
model.yu      = [1 1 1 1 1];
model.regr    = [1 1 1 1 1];
%% Simulation parameters
simul=struct('N', 100, 'nEstimates', 20, 'np', 0.0005); 

%% initialization variables
y=zeros(simul.N, 1);
u=ones(simul.N, 1);
% initial conditions
y(1)=28;

%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;

for m=1:simul.nEstimates
    %% Simulation of real system
    for k=max(abs(model.regr))+1:simul.N
        y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1))+rand(1)*np;
    end

    psi = f_get_psi(y, u, model);
    theta(1,:)=(psi'*psi)\(psi'*y);
    
    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1+simul.np*2, model.dim);
    % we can't have a precision bigger than the noise power
    while (max(abs(err)) > simul.np*2)
        yc=f_y_model([y(1) y(2)], u, theta(l,:), model);
        if l == 1
            f_plot_y_y1(yc);
        end
        %% step 2 -  calc the variance
        v(l)=cov(y-yc);
        [PHY phy]=f_get_phy(y, model);
        
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
    v'
    %f_plot_y_y1(yc);
end %J

f_draw_elipse(nda, ndb)
