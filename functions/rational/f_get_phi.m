function [phi phim] = f_get_phi(y, m)
%====================
%% Get phi based on model structure
% y:: output system data [y1,..,yn]
% m:: model structure
%
% Return phi and PHI 
%====================

% Check if model parameter structure is valid
f_check_model(m);

% variables initialization
N    = max(size(y));
phi  = zeros(m.dim + m.error_model_dim, m.dim + m.error_model_dim);
phim = zeros(m.dim + m.error_model_dim, 1);

for i=1:m.dim
    for j=1:m.dim
        if i>m.n_dim && j>m.n_dim
            aux=0;
            for k=1:N
                aux = aux + ((y(k)^m.a_exp(j))* (y(k)^m.b_exp(j))* (y(k)^m.a_exp(i))* (y(k)^m.b_exp(i)));
            end
            phi(i,j)=aux;
        end
    end
end

for i= m.n_dim+1:m.dim
    aux=0;
    for k=1:N
        aux = aux- ((y(k)^m.a_exp(i))* (y(k)^m.b_exp(i)));
    end
    phim(i,1) = aux;
end
end %function
