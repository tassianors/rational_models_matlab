function psi = f_get_psi(y, u, m)
%% Gets the Psi based on the model structure
% y:: output system data [y1,..,yn]
% u:: innput system data [u1,..,un]
%% Model
% m.dim:: Total dimension (num+den dimensions)
% m.n_dim:: dimension of numerator
% m.texp::exponential coef from the model
% m.yu:: one is y, 0 is u. defines witch one of the coef are dependent of y
%        and u. should be a array with the same size of #m.texp.
% m.regr:: array whith the regretion dimension (y(k-1) is 1) always
%          positive values.
%%

% check parameters
if m.n_dim > m.dim
    error('size of theta must be bigger than m.n_dim');
end

if max(size(m.yu)) < m.dim || max(size(m.regr)) < m.dim || max(size(m.texp)) < m.dim
    error('invalid parameter dimenstion');
end
%% step 1 - first estimative
N=max(size(y,1));
psi=zeros(N,m.dim);

for i=max(abs(m.regr))+1:N
    % step 3 - we can neglect the noise model (colun 6)
    for j=1:m.dim
        if m.yu(j) == 1
            yu=y(i-abs(m.regr(j)));
        else
            yu=u(i-abs(m.regr(j)));
        end
        
        if j<=m.n_dim
            psi(i,j)=yu^m.texp(j);
        else
            % here we should alway use y(i)
            psi(i,j)=-yu^m.texp(j)*y(i);
        end
    end
end
end
