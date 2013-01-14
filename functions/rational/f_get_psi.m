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

if m.error_model_dim && max(size(out_sig)) ~= max(size(out_prev))
    error('out_sig and last out_sig (out_prev) have to have the same size');
else
    err=out_sig-out_prev;
end

N = max(size(out_sig));
psi = zeros(N, m.dim+m.error_model_dim);

for i=f_model_get_max_regressor(m)+1:N
    for j=1:m.dim+m.error_model_dim
        % err_model
        if j> m.dim
           % fill the err_model colunm.
           psi(i,j)=err(i);
           continue;
        end
        
        if m.a_signal_type(j) == 1
            yu=out_sig(i-abs(m.a_regress(j)))^m.a_exp(j);
        elseif m.a_signal_type(j) == 2
            yu=in_sig(i-abs(m.a_regress(j)))^m.a_exp(j);
        elseif m.a_signal_type(j) == 3
            yu=aux_sig1(i-abs(m.a_regress(j)))^m.a_exp(j);
        elseif m.a_signal_type(j) == 4
            yu=aux_sig2(i-abs(m.a_regress(j)))^m.a_exp(j);
        else
            error('not supported option');
        end
        % non linearity is yu^a*y^b
        yu2=1;
        if m.b_signal_type(j) == 1
            yu2=out_sig(i-abs(m.b_regress(j)))^m.b_exp(j);
        elseif m.b_signal_type(j) == 2
            yu2=in_sig(i-abs(m.b_regress(j)))^m.b_exp(j);
        elseif m.b_signal_type(j) == 3
            yu2=aux_sig1(i-abs(m.b_regress(j)))^m.b_exp(j);
        elseif m.b_signal_type(j) == 4
            yu2=aux_sig2(i-abs(m.b_regress(j)))^m.b_exp(j);
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
