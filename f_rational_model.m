function ret = f_rational_model(simul, model, yuser, y_init, u)
%% DOC
%
%%

f_check_model(model);
if simul.N ~= size(yuser,1)
    error('y must be a [1 N] array]');
end
ret=zeros(simul.nEstimates, model.dim);

for m=1:simul.nEstimates
    clear theta delta v;
    %% initialization variables
    y=zeros(simul.N, 1);
    yc=y;
    model.err_model = 0;
    % set randon noise
    y=f_get_wnoise(yuser, simul.np);
    psi = f_get_psi(y, yc, u, model);
    theta(1,:)=(psi'*psi)\(psi'*y);
    
    %% here we got the first estimative, now we start the loop
    l=1;
    err=ones(1, model.dim);
    v_diff=simul.diffConv+1;
    % we can't have a precision bigger than the err_model power
    while ((max(abs(err)) > simul.maxError || abs(v_diff) > simul.diffConv) && l < simul.l)
        yc=f_y_model(y_init, u, theta(l,:), model);
        
        % only after the first estimative, calc using the error model
        if l == 2 && model.err_enable == true
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
            v_diff = v(l);
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
        for o=1:model.dim
            ret(m,o)=theta(l+1,o);
        end
        l=l+1;
    end
end %m
end %function