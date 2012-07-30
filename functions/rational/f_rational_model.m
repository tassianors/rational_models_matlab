function [ret cost]= f_rational_model(simul, model, ic, out_sig, in_sig, aux_sig1, aux_sig2)
%% DOC
%
%%
tic
f_check_model(model);

if simul.N ~= size(out_sig,1)
    error('y must be a [1 N] array]');
end

%init return variable with zeros
ret=zeros(simul.nEstimates, model.dim);

for m=1:simul.nEstimates
    clear theta v;
    %% initialization variables: 
    y = zeros(simul.N, 1);
    yc = y;
    model.err_model = 0;
    l=1;

    % set randon noise
    y = f_get_wnoise(out_sig, simul.np);
    psi = f_get_psi(y, yc, in_sig, aux_sig1, aux_sig2, model);
    theta(m,:) = (psi'*psi)\(psi'*y);
    
    %% here we got the first estimative, now we start the loop
    while (l < simul.l)
        yc = f_y_model(ic, in_sig, aux_sig1, aux_sig2, theta(l,:), model);
        
        % only after the first estimative, calc using the error model
        if l == 2 && model.err_enable == true
            model.err_model = 1;
            % enlarge the matrix
            theta(l, model.dim+model.err_model) = 0;
            delta(l, model.dim+model.err_model) = 0;
        end
        
        %% step 2 -  calc the variance
        v(l)=cov(y-yc);

        psi = f_get_psi(y, yc, in_sig, aux_sig1, aux_sig2, model);
        [PHY phy] = f_get_phy(y, model);
        
        theta(l+1,:) = (psi'*psi-v(l)*PHY)\(psi'*y-v(l)*phy);
        
        costJy(l)=sum(0.5*abs(y-yc).^2)/simul.N;
        if l > 1 && costJy(l)-costJy(l-1) > 0
            if simul.verbose == true
                str = sprintf('During iteration %d we checked that the cost rised from %s to %s', l, num2str(costJy(l-1), 5), num2str(costJy(l), 5));
                disp(str);
            end
            break;
        end
        if simul.verbose == true
            str = sprintf('Iteration %d: Cost %s covariance %s ', l, num2str(costJy(l), 5), num2str(v(l), 5));
            disp(str)
        end
        l=l+1;
    end
    % discarding the last estimative if it increased the error
    for o = 1: model.dim
        ret(m,o) = theta(l-1,o);
    end
    cost=costJy(l-1);
    % finishing the execution.
    total_time=toc;    
    if l >= simul.l
        error('procedure diverged and also took %4.2f s to complete',total_time)
        ret = zeros(simul.nEstimates, model.dim);
        cost=0;
        return
    end
    if simul.verbose == true
        str = sprintf('Execution finished in %d iterations [%4.2f seconds]: Final cost was %s\n', l, total_time, num2str(costJy(l-1), 5));
        disp(str)
    end

end %m
end %function