function psi = f_get_psi(out_sig, out_prev, in_sig, aux_sig1, aux_sig2, m)
%====================
%% Get Psi based on model structure
% out_sig:: output system data
% in_sig:: input system data 
% aux_sig1:: input aux signal 
% aux_sig2:: input aux signal 
% m:: model structure
%====================

% check parameters
f_check_model(m);

if m.err_model && max(size(out_sig)) ~= max(size(out_prev))
    error('out_sig and last out_sig (out_prev) have to have the same size');
else
    err=out_sig-out_prev;
end

N = max(size(out_sig));
psi = zeros(N, m.dim+m.err_model);

for i=f_model_get_max_regressor(m)+1:N
    for j=1:m.dim+m.err_model
        % err_model
        if j> m.dim
           % fill the err_model colunm.
           psi(i,j)=err(i);
           continue;
        end
        
        if m.yu(j) == 1
            yu=out_sig(i-abs(m.regr(j)))^m.texp(j);
        elseif m.yu(j) == 2
            yu=in_sig(i-abs(m.regr(j)))^m.texp(j);
        elseif m.yu(j) == 3
            yu=aux_sig1(i-abs(m.regr(j)))^m.texp(j);
        elseif m.yu(j) == 4
            yu=aux_sig2(i-abs(m.regr(j)))^m.texp(j);
        else
            error('not supported option');
        end
        % non linearity is yu^a*y^b
        yu2=1;
        if m.yplus_yur(j) == 1
            yu2=out_sig(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        elseif m.yplus_yur(j) == 2
            yu2=in_sig(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        elseif m.yplus_yur(j) == 3
            yu2=aux_sig1(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        elseif m.yplus_yur(j) == 4
            yu2=aux_sig2(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        end
        if j<=m.n_dim
            psi(i,j)=yu*yu2;
        else
            % here we should alway use y(i) (equation 10.44 aguirre)
            psi(i,j)=-yu*yu2*out_sig(i);
        end

    end
end
end
