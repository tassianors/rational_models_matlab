function y = f_y_model_k(k, in_sig, aux_sig1, aux_sig2, theta, m)
%====================
%% Get output signal based on model structure
%
% k:: iteration
% in_sig:: input signal array
% theta:: estimative model parameters
% m:: model structure
%====================
y = 0;
num=0;
m_dim=size(theta,2);
% check parameters
f_check_model(m);
if m.dim+m.error_model_dim ~= m_dim
    error('size of theta must be the same as model dimension');
end

if f_model_get_max_regressor(m) > k
    error('k is out of bound based on Model parameters.');
end

for i=1:m.n_dim
    if m.a_signal_type(i) == 1
        yu = y(k-(m.a_regress(i)))^m.a_exp(i);
    elseif m.a_signal_type(i) == 2
        yu = in_sig(k-(m.a_regress(i)))^m.a_exp(i);
    elseif m.a_signal_type(i) == 3
        yu = aux_sig1(k-(m.a_regress(i)))^m.a_exp(i);
    elseif m.a_signal_type(i) == 4
        yu = aux_sig2(k-(m.a_regress(i)))^m.a_exp(i);
    else
        error('invalid option, just y(1) u(2) r(3) and e(4) are possible')
    end
    % non linearity is yu^a*y^b
    if m.b_signal_type(i) == 1
        yu = yu*y(k-(m.b_regress(i)))^m.b_exp(i);
    elseif m.b_signal_type(i) == 2
        yu = yu*in_sig(k-(m.b_regress(i)))^m.b_exp(i);
    elseif m.b_signal_type(i) == 3
        yu = yu*aux_sig1(k-(m.b_regress(i)))^m.b_exp(i);
    elseif m.b_signal_type(i) == 4
        yu = yu*aux_sig2(k-(m.b_regress(i)))^m.b_exp(i);
    end
    num = num+theta(i)*yu;
end
den = 1;
for i=m.n_dim+1:m.dim
    if m.a_signal_type(i) == 1
        yu = y(k-(m.a_regress(i)))^m.a_exp(i);
    elseif m.a_signal_type(i) == 2
        yu = in_sig(k-(m.a_regress(i)))^m.a_exp(i);
    elseif m.a_signal_type(i) == 3
        yu = aux_sig1(k-(m.a_regress(i)))^m.a_exp(i);
    elseif m.a_signal_type(i) == 4
        yu = aux_sig2(k-(m.a_regress(i)))^m.a_exp(i);
    end
    
    % non linearity is yu^a*y^b
    if m.b_signal_type(i) == 1
        yu = yu*y(k-abs(m.b_regress(i)))^m.b_exp(i);
    elseif m.b_signal_type(i) == 2
        yu = yu*in_sig(k-abs(m.b_regress(i)))^m.b_exp(i);
    elseif m.b_signal_type(i) == 3
        yu = yu*aux_sig1(k-abs(m.b_regress(i)))^m.b_exp(i);
    elseif m.b_signal_type(i) == 4
        yu = yu*aux_sig2(k-abs(m.b_regress(i)))^m.b_exp(i);
    end
    den=den+theta(i)*yu;
end

y= num/den;

end
