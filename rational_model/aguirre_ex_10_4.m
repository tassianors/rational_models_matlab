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
enable=true;

use_iv = false


%% Simulation parameters
simul=struct('N', 300, 'nEstimates', 10, 'np', 0, 'maxError', 0.001, 'l', 200, 'diffConv', 100); 

%% Real system - variables
a=2.6204; b=99.875; c=1417.1; d=46.429;
expected=[8.658 0.001223 -0.0441 -0.08381 0.001766];

for m=1:simul.nEstimates
    clear theta delta v y yc u psi PHI z_iv l err v_diff;
	%% initialization variables
	y=zeros(simul.N, 1);
	yc=y;
    u=f_get_square_signal(simul.N);
	y(1)=23.2;%+rand(1)*.03;
	
    model.err_model = 0;
    % Simulation of real system
%      for k=max(abs(model.regr))+1:simul.N
%          y(k)=d*exp(22-y(k-1))+ ((a*y(k-1)^2-b*y(k-1)+c)/y(k-1));
%      end
    y = f_aguirre_get_model_output(model, simul, expected,y(1));
%     
%     f_aguirre_plot_map(y, m)
%     f_aguirre_plot_map(y2, m+1)
	% set randon noise
  	y=f_get_wnoise(y, 0.001);
	
    psi = f_get_psi(y, yc, u, 0, model);
    z_iv = zeros(size(psi));
    for t=6:simul.N
        % auxiliary instrument z
        z_iv(t, 5)=u(t-1);
        z_iv(t, 4)=u(t-2);
        z_iv(t, 3)=u(t-3);
        z_iv(t, 2)=u(t-4);
        z_iv(t, 1)=u(t-5);
    end
    if use_iv == false
        theta(1,:)=(psi'*psi)\(psi'*y)
    else
        theta(1,:)=(z_iv'*psi)\(z_iv'*y);
    end

    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1, model.dim);
    v_diff=simul.diffConv+1;
    
    % here I got the result from the first estimative
    yc=f_y_model(y(1), u, 0, theta(l,:), model);
        
    % we can't have a precision bigger than the err_model power
    while ((max(abs(err)) > simul.maxError || abs(v_diff) > simul.diffConv) && l < simul.l)
        
            
        % only after the first estimative, calc using the error model
        if l == 1 && enable == true
            model.err_model = 1;
            % enlarge the matrix
            theta(l, model.dim+model.err_model)=0;
            delta(l, model.dim+model.err_model)=0;
        end
    
        %% step 2 -  calc the variance between original signal and estimated one
        v(l)=cov(y-yc);
        
        % check if variance is converging
		if l > 1
			v_diff = abs(v(l)-v(l-1));
		else
			v_diff=v(l);
        end
        
        % here I got the phi and Phy matrix
        psi = f_get_psi(y, yc, u, 0, model);
        [PHI phi]=f_get_phi(y, model);
        
        % calculating the first aproximation (overwrite the first)
        if use_iv == false
            theta(l,:)=(psi'*psi-v(l)*PHI)\ (psi'*y-v(l)*phi);
        else
            z_iv = zeros(size(psi));
            for t=size(psi,2)+1:simul.N
                % auxiliary instrument z
                for g = 1: size(psi,2)
                    z_iv(t, size(psi,2)+1-g)=u(t-g);
                end
            end
            theta(l,:)=(z_iv'*psi+(l)*PHI)\ (z_iv'*y+v(l)*phi);
        end
        
        if l > 1
            delta(l,:)=theta(l,:)-theta(l-1,:);
           	clear err;
            err(1,:)=delta(l,:);
            if max(size(err)) > model.dim
                err(1,model.dim+model.err_model)=0;
            end
        end

        % to be used in graphic plotting
        nna(m)=theta(l,1);
        nnb(m)=theta(l,2);
        nnc(m)=theta(l,3);
        nda(m)=theta(l,4);
        ndb(m)=theta(l,5);
        yc=f_y_model(y(1), u, 0, theta(l,:), model);
        l=l+1;
        
        
    end
    theta
    delta;
    v';
end %J
result = f_y_model(y(1), u, 0, theta(size(theta, 1),:), model);
f_aguirre_plot_map(result, m+1);


f_draw_elipse(nna, nnb, expected(1), expected(2));
f_draw_elipse(nda, ndb, expected(4), expected(5));