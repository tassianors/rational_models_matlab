function psi = f_get_psi(y, u, m_dim, m_n_dim, m_texp, m_yu, m_regr)
%% Gets the Psi based on the model structure
% y:: output system data [y1,..,yn]
% u:: innput system data [u1,..,un]
%% Model
% m_dim:: Total dimension (num+den dimensions)
% m_n_dim:: dimension of numerator
% m_texp::exponential coef from the model
% m_yu:: one is y, 0 is u. defines witch one of the coef are dependent of y
%        and u. should be a array with the same size of #m_texp.
% m_regr:: array whith the regretion dimension (y(k-1) is 1) always
%          positive values.
%%

% check parameters
if m_n_dim > m_dim
    error('size of theta must be bigger than m_n_dim');
end

if max(size(m_yu)) < m_dim || max(size(m_regr)) < m_dim || max(size(m_texp)) < m_dim
    error('invalid parameter dimenstion');
end
%% step 1 - first estimative
N=max(size(y,1));
psi=zeros(N,m_dim);

for i=max(abs(m_regr))+1:N
    % step 3 - we can neglect the noise model (colun 6)
    for j=1:m_dim
        if m_yu(j) == 1
            yu=y(i-abs(m_regr(j)));
        else
            yu=u(i-abs(m_regr(j)));
        end
        
        if j<=m_n_dim
            psi(i,j)=yu^m_texp(j);
        else
            % here we should alway use y(i)
            psi(i,j)=-yu^m_texp(j)*y(i);
        end
    end
end
end
