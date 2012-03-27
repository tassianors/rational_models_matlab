function y = f_y_model(y_init, u, r, theta, m)
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
m_dim=size(theta,2);
f_check_model(m);
% check parameters
if m.dim+m.err_model ~= m_dim
    error('size of theta must be the same as model dimension');
end

for k=max(abs(m.regr))+1:N
    num=0;
    for i=1:m.n_dim
        if m.yu(i) == 1
            yu=y(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 2
            yu=u(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 3
            yu=r(k-abs(m.regr(i)))^m.texp(i);
        else
            error('invalid option, just y(1) u(2) and r(3) are possible')
        end
        % non linearity is yu^a*y^b
        if m.yplus_yur(i) == 1
            yu=yu*y(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 2
            yu=yu*u(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 3
            yu=yu*r(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        end
        num=num+theta(i)*yu;
    end
    den = 1;
    for i=m.n_dim+1:m.dim
        if m.yu(i) == 1
            yu=y(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 2
            yu=u(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 3
            yu=r(k-abs(m.regr(i)))^m.texp(i);
        end
        
        % non linearity is yu^a*y^b
        if m.yplus_yur(i) == 1
            yu=yu*y(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        end
        % non linearity is yu^a*u^b
        if m.yplus_yur(i) == 2
            yu=yu*u(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        end
        if m.yplus_yur(i) == 3
            yu=yu*r(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        end
        den=den+theta(i)*yu;
    end
    y(k)= num/den;
end
end
