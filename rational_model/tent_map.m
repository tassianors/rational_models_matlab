%% Tent MAP
close all; clear all;
clc;
P=path;
path(P,'../functions');
%% model parameter definition
model.n_dim   = 3;
model.dim     = 5;
model.texp    = [0 2 1 1 2];
model.yu      = [1 1 1 1 1];
model.regr    = [1 1 1 1 1];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
model.yplus_uy = [0 0 0 0 0];
% tels the d param
model.yplus_exp = [0 0 0 0 0];
% tels the C param
model.yplus_regr = [0 0 0 0 0];

model.err_model   = 0;
enable=true;
%% Simulation parameters
simul=struct('N', 200, 'nEstimates', 4, 'np', 0.5, 'maxError', 0.1, 'l', 100, 'diffConv', 100); 

%% Real system - variables
a=1.999; b=0.5;
expected=[0.02608 -1.325 1.325 -2.416 2.416];

for m=1:simul.nEstimates
    clear theta delta v;
	%% initialization variables
	y=zeros(simul.N, 1);
	yc=y;
	u=ones(simul.N, 1);
	y(1)=0;
	
    model.err_model = 0;
    % Simulation of real system
     for k=max(abs(model.regr))+1:simul.N
         y(k)=1-a*abs(y(k-1)-b);
     end
%     y2 = f_aguirre_get_model_output(model, simul, expected, y(1))
%     
%     f_aguirre_plot_map(y, m)
%     f_aguirre_plot_map(y2, m+1)
	% set randon noise
	%y=y+y.*+rand(simul.N,1)*(mean(y)/200*simul.np);
	
    psi = f_get_psi(y, yc, u, model);
    theta(1,:)=(psi'*psi)\(psi'*y);

    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1, model.dim);
    v_diff=simul.diffConv+1;
    % we can't have a precision bigger than the err_model power
    while ((max(abs(err)) > simul.maxError || abs(v_diff) > simul.diffConv) && l < simul.l)
        yc=f_y_model(y(1) , u, theta(l,:), model);
    
        % only after the first estimative, calc using the error model
        if l == 2 && enable == true
            model.err_model = 1;
            % enlarge the matrix
            theta(l, model.dim+model.err_model)=0;
            delta(l, model.dim+model.err_model)=0;
        end
    
        %% step 2 -  calc the variance
        v(l)=cov(y-yc);
		if l > 1
			v_diff = v(l)-v(l-1);
		else
			v_diff=v(l);
		end

        psi = f_get_psi(y, yc, u, model);
        [PHY phy]=f_get_phy(y, model);
        
        theta(l+1,:)=(psi'*psi-v(l)*PHY)\ (psi'*y-v(l)*phy);
        delta(l,:)=theta(l+1,:)-theta(l,:);
        clear err;
        err(1,:)=delta(l,:);
        if max(size(err)) > model.dim
            err(1,model.dim+model.err_model)=0;
        end
        % to be used in graphic plotting
        nna(m)=theta(l+1,1);
        nnb(m)=theta(l+1,2);
        nnc(m)=theta(l+1,3);
        nda(m)=theta(l+1,4);
        ndb(m)=theta(l+1,5);
        l=l+1;
    end
    theta
    delta;
    v';
end %J
result = f_aguirre_get_model_output(model, simul, theta(size(theta, 1),:), y(1));
f_aguirre_plot_map(result, m+1);

f_draw_elipse(nna, nnb, expected(1), expected(2));
f_draw_elipse(nda, ndb, expected(4), expected(5));