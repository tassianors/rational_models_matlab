function [phy phym] = f_get_phy(y, m_dim, m_n_dim, m_texp, m_yu, m_regr)
%% Gets the phy based on the model structure
% y:: output system data [y1,..,yn]
%% Model
% m_dim:: Total dimension (num+den dimensions)
% m_n_dim:: dimension of numerator
% m_texp::exponential coef from the model
% m_yu:: one is y, 0 is u. defines witch one of the coef are dependent of y
%        and u. should be a array with the same size of #m_texp.
% m_regr:: array whith the regretion dimension (y(k-1) is 1) always
%          positive values.
%%

if m_n_dim > m_dim
    error('m_dim must be bigget than m_n_dim');
end
if max(size(m_yu)) < m_dim || max(size(m_regr)) < m_dim || max(size(m_texp)) < m_dim
    error('invalid parameter dimenstion');
end

% initialization
N=max(size(y,1));
phy=zeros(m_dim,m_n_dim);
phym=zeros(m_dim,1);

for i=1:m_dim
    for j=1:m_dim
        if i>m_n_dim && j>m_n_dim
            aux=0;
            for k=1:N
                aux=aux+((y(k)^m_texp(j))* (y(k)^m_texp(i)));
            end
            phy(i,j)=aux;
        end
    end
end
for i=m_n_dim+1:m_dim
    aux=0;
    for k=1:N
        aux=aux-(y(k)^m_texp(i));
    end
    phym(i,1)=aux;
end
end
