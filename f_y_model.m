function y = f_y_model(y_init, N, u, theta, m_nn_dim, m_texp, m_yu, m_regr)
%% Gets the y based on the model structure
% y_init:: Initial value for y, can be a array [y(1), y(2), ...]
% N:: size of Y
% u:: input signal [u1,u2,u3,...]
% theta:: estimative of parameters for the model
%% Model
% m_nn_dim:: Total dimension (num+den dimensions)
% m_texp::exponential coef from the model
% m_yu:: one is y, 0 is u. defines witch one of the coef are dependent of y
%        and u. should be a array with the same size of #m_texp.
% m_regr:: array whith the regretion dimension (y(k-1) is 1) always
%          positive values.
%%
y=zeros(N, 1);
% init the y array
for j=1: size(y_init,2)
    y(j)=y_init(j);
end
nd_nn_dim=size(theta,2);

% check parameters
if m_nn_dim > nd_nn_dim
    error('size of theta must be bigger than m_nn_dim');
end

if max(size(m_yu)) < nd_nn_dim || max(size(m_regr)) < nd_nn_dim || max(size(m_texp)) < nd_nn_dim
    error('invalid parameter dimenstion');
end

for k=max(abs(m_regr))+1:N
    num=0;
    for i=1:m_nn_dim
        if m_yu(i) == 1
            yu=y(k-abs(m_regr(i)));
        else
            yu=u(k-abs(m_regr(i)));
        end
        num=num+theta(i)*yu^m_texp(i);
    end
    den = 1;
    for i=m_nn_dim+1:nd_nn_dim
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
