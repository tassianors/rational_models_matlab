%% Aguirre 10.4
close all; clear all;
clc;
%% model parameter definition
model.n_dim   = 3;
model.dim     = 4;
model.a_exp    = [2 1 1 3];
model.a_signal_type      = [1 1 2 1];
model.a_regress    = [1 2 1 2];
% u = 2 y=1 none =0
model.b_signal_type = [0 0 0 0];
% tels the d param
model.b_exp = [0 0 0 0];
% tels the C param
model.b_regress = [0 0 0 0];
model.error_model_dim   = 0;
enable=true;
%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 30, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1); 

%% Real system - variables
a1=.3; a2=-2; a3=1.5; b1=.2;

for m=1:simul.nEstimates
    clear theta delta v;
	%% initialization variables
	y=zeros(simul.N, 1);
	yc=y;
	u=ones(simul.N, 1);
	% initial conditions
	y(1)=0;
	y(2)=0;
    model.error_model_dim = 0;
    %% Simulation of real system
    for k=max(abs(model.a_regress))+1:simul.N
        y(k)=(a1*y(k-1)^2+a2*y(k-2)+a3*u(k-1))/(1+b1*y(k-2)^3);
	end
	% set randon noise
    y=f_get_wnoise(y, simul.np);
    psi = f_get_psi(y, yc, u, 0, model);
    theta(1,:)=(psi'*psi)\(psi'*y);

    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1, model.dim);
    v_diff=simul.diffConv+1;
    % we can't have a precision bigger than the err_model power
   while ((max(abs(err)) > simul.maxError || abs(v_diff) > simul.diffConv) && l < simul.l)
        yc=f_y_model([y(1) y(2)], u, 0, theta(l,:), model);
    
        % only after the first estimative, calc using the error model
        if l == 2 && enable == true
            model.error_model_dim = 1;
            % enlarge the matrix
            theta(l, model.dim+model.error_model_dim)=0;
            delta(l, model.dim+model.error_model_dim)=0;
        end
    
        %% step 2 -  calc the variance
        v(l)=cov(y-yc);
		if l > 1
			v_diff = v(l)-v(l-1);
		else
			v_diff=v(l);
		end

        psi = f_get_psi(y, yc, u, 0, model);
        [PHI phi]=f_get_phi(y, model);
        
        theta(l+1,:)=(psi'*psi-v(l)*PHI)\ (psi'*y-v(l)*phi);
        delta(l,:)=theta(l+1,:)-theta(l,:);
        clear err;
        err(1,:)=delta(l,:);
        if max(size(err)) > model.dim
            err(1,model.dim+model.error_model_dim)=0;
        end
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

f_draw_elipse(nna, nnb, a1, a2);
f_draw_elipse(nnc, nda, a3, b1);
