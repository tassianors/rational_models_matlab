function y = f_y_model(y_init, u, theta, m)
%% Gets the y based on the model structure
% y_init:: Initial value for y, can be a array [y(1), y(2), ...]
% u:: input signal [u1,u2,u3,...]
% theta:: estimative of parameters for the model
% m:: model
%%
N=max(size(u));
y=zeros(N, 1);
% init the y array
for j=1: size(y_init,2)
    y(j)=y_init(j);
end
m_dim=size(theta,2)
f_check_model(m);
% check parameters
if m.dim+m.err_model ~= m_dim
    error('size of theta must be the same as model dimension');
end

for k=max(abs(m.regr))+1:N
    num=0;
    for i=1:m.n_dim
        if m.yu(i) == 1
            yu=y(k-abs(m.regr(i)));
        else
            yu=u(k-abs(m.regr(i)));
        end
        num=num+theta(i)*yu^m.texp(i);
    end
    den = 1;
    for i=m.n_dim+1:m.dim
        if m.yu(i) == 1
            yu=y(k-abs(m.regr(i)));
        else
            yu=u(k-abs(m.regr(i)));
        end
        den=den+theta(i)*yu^m.texp(i);
    end
    y(k)= num/den;
end
end
