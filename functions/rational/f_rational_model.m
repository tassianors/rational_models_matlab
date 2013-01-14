function [ret cost] = f_rational_model(simul, model, ic, out_sig, in_sig, aux_sig1, aux_sig2)
%====================
%% Function to perform the identification of a system using a recursive algorithm.
%
% simul :: Simulation parameter. Must have:
%       N - size of data array from system
%       np - Amount of noise power to inject on each simulation
%       l - max acceptable number of algorithm iterations. This is used as a time out to avoid infinity loop.
%       verbose - used as debug to show some useful identification information 
%
% model :: Model specification
%       See #f_check_model function for further more information
% 
% ic    :: Initial system conditions
% 
% out_sig :: Output signal of system under identification
% in_sig  :: Input signal of system under identification
% aux_sig1 :: Auxiliary signal one that can be used or not in identification procedure
% aux_sig2 :: Auxiliary signal two that can be used or not in identification procedure
%====================

% start time measurement (used to know how much time has been spend in this procedure)
tic
f_check_model(model);

if simul.N ~= size(out_sig,1)
    error('out signal must be a [1 N] array]');
end

if simul.N ~= size(in_sig,1)
    error('in signal must be a [1 N] array]');
end
if simul.N ~= size(aux_sig1,1)
    error('aux signal must be a [1 N] array]');
end
if simul.N ~= size(aux_sig2,1)
    error('aux signal must be a [1 N] array]');
end

% init return variable with zeros
ret = zeros(simul.nEstimates, model.dim);

for m=1:simul.nEstimates
    clear theta v;
    %% init some variables at each estimative
    y = zeros(simul.N, 1);
    yc = y;
    model.error_model_dim = 0;
    l=1;

    % add a random noise with some defined noise power
    y = f_get_wnoise(out_sig, simul.np);
    
    psi = f_get_psi(y, yc, in_sig, aux_sig1, aux_sig2, model);
    theta(m,:) = (psi'*psi)\(psi'*y);
    
    % here we got the first estimative, now we start the loop until cost function stops to reduce.
    while (l < simul.l)
        % get output signal from model last estimative
        yc = f_y_model(ic, in_sig, aux_sig1, aux_sig2, theta(l,:), model);
        
        % only after the first estimative, calc using the error model
        if l == 2 && model.error_in_account == true
            model.error_model_dim = 1;
            % enlarge some matrix
            theta(l, model.dim+model.error_model_dim) = 0;
            delta(l, model.dim+model.error_model_dim) = 0;
        end
        
        %% step 2 -  calc covariance
        v(l)=cov(y-yc);
        
        % based on user output signal and last estimative, get new psi
        psi = f_get_psi(y, yc, in_sig, aux_sig1, aux_sig2, model);
        [PHI phi] = f_get_phi(y, model);
        
        % update theta for next estimative, if needed.
        theta(l+1,:) = (psi'*psi-v(l)*PHI)\(psi'*y-v(l)*phi);
        
        % Calc Jy cost function
        costJy(l)=sum(0.5*abs(y-yc).^2)/simul.N;
        % Check if cost is still decreasing, if not stop
        if l > 1 && costJy(l)-costJy(l-1) > 0
            if simul.verbose == true
                str = sprintf('At iteration %d we checked that the cost raised from %s to %s', l, num2str(costJy(l-1), 5), num2str(costJy(l), 5));
                disp(str);
            end
            break;
        end
        if simul.verbose == true
            str = sprintf('Iteration %d: Cost was %s covariance was %s ', l, num2str(costJy(l), 5), num2str(v(l), 5));
            disp(str)
        end
        l=l+1;
    end
    % discarding the last estimative, unless we run out by timeout
    for o = 1: model.dim
        ret(m,o) = theta(l-1,o);
    end
    % return final cost obtained.
    cost=costJy(l-1);
    
    % finishing time measurement.
    total_time=toc;    
    % check timeout occurence
    if l >= simul.l
        error('Procedure diverged and also took %4.2f s to complete',total_time)
        ret = zeros(simul.nEstimates, model.dim);
        cost=0;
        return
    end
    if simul.verbose == true
        str = sprintf('Execution finished in %d iterations [total time spend: %4.2f seconds]. Final cost was %s\n', l, total_time, num2str(costJy(l-1), 5));
        disp(str)
    end

end %m
end %function
