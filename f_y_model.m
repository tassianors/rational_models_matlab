function y = f_y_model(y_init, u, theta, m_n_dim, m_texp, m_yu, m_regr)
%% Gets the y based on the model structure
% y_init:: Initial value for y, can be a array [y(1), y(2), ...]
% u:: input signal [u1,u2,u3,...]
% theta:: estimative of parameters for the model
%% Model
% m_n_dim:: Total dimension (num+den dimensions)
% m_texp::exponential coef from the model
% m_yu:: one is y, 0 is u. defines witch one of the coef are dependent of y
%        and u. should be a array with the same size of #m_texp.
% m_regr:: array whith the regretion dimension (y(k-1) is 1) always
%          positive values.
%%
N=max(size(u));
y=zeros(N, 1);
% init the y array
for j=1: size(y_init,2)
    y(j)=y_init(j);
end
m_dim=size(theta,2);

% check parameters
if m_n_dim > m_dim
    error('size of theta must be bigger than m_n_dim');
end

if max(size(m_yu)) < m_dim || max(size(m_regr)) < m_dim || max(size(m_texp)) < m_dim
    error('invalid parameter dimenstion');
end

for k=max(abs(m_regr))+1:N
    num=0;
    for i=1:m_n_dim
        if m_yu(i) == 1
            yu=y(k-abs(m_regr(i)));
        else
            yu=u(k-abs(m_regr(i)));
        end
        num=num+theta(i)*yu^m_texp(i);
    end
    den = 1;
    for i=m_n_dim+1:m_dim
        if m_yu(i) == 1
            yu=y(k-abs(m_regr(i)));
        else
            yu=u(k-abs(m_regr(i)));
        end
        den=den+theta(i)*yu^m_texp(i);
    end
    y(k)= num/den;
end
end
