function y = f_y_model(ic, in_sig, aux_sig1, aux_sig2, theta, m)
%====================
%% Get output signal based on model structure
%
% ic:: Initial value for y, can be an array [y(1), y(2), ...]
% in_sig:: input signal array
% theta:: estimative model parameters
% m:: model structure
%====================

N = max(size(in_sig));
y = zeros(N, 1);
% init y array with initial conditions (ic)
for j=1: size(ic,2)
    y(j)=ic(j);
end
m_dim=size(theta,2);

% check parameters
f_check_model(m);
if m.dim+m.error_model_dim ~= m_dim
    error('size of theta must be the same as model dimension');
end

for k=max(max(abs(m.a_regress)), max(abs(m.b_regress)))+1:N
    y(k)=f_y_model_k(k, in_sig, aux_sig1, aux_sig2, theta, m);
end
end
