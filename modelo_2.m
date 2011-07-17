%% Aguirre 10.4
close all; clear all;
clc;
%% model parameter definition
model.n_dim   = 3;
model.dim     = 4;
model.texp    = [2 1 1 3];
model.yu      = [1 1 0 1];
model.regr    = [1 2 1 2];
model.err_model   = 0;
%% Simulation parameters
simul=struct('N', 100, 'nEstimates', 20, 'np', 0.0005); 

%% initialization variables
y=zeros(simul.N, 1);
yc=y;
u=ones(simul.N, 1);
% initial conditions
y(1)=0;
y(2)=0;

%% Real system - variables
a1=.3;
a2=-2;
a3=1.5;
b1=-.2;

for m=1:simul.nEstimates
    %% Simulation of real system
    for k=max(abs(model.regr))+1:simul.N
        y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(1-b1*y(k-2)^3)+rand(1)*simul.np;
    end
    if m==2
        model.err_model   = 1;
    end
    psi = f_get_psi(y, yc, u, model);
    % only after the first estimative, calc using the error model

    if m ==2 && size(psi,2) ~= size(theta,2) 
        %m>1 && model.err_model == 1
        % enlarge the matrix.
        theta(m, model.dim+model.err_model)=0;
        delta(m, model.dim+model.err_model)=0;
    end
    theta(1,:)=(psi'*psi)\(psi'*y)

    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1+simul.np*2, model.dim);
    % we can't have a precision bigger than the err_model power
    while (max(abs(err)) > simul.np*2)
        yc=f_y_model([y(1) y(2)], u, theta(l,:), model);
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
        l=l+1;
    end
    theta
    delta
    v'
    %f_plot_y_y1(yc);
end %J

f_draw_elipse(nna, nnb)
f_draw_elipse(nnc, nda)
