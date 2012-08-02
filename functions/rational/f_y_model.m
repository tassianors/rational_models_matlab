function y = f_y_model(ic, in_sig, aux_sig1, aux_sig2, theta, m)
%%%%%%%%%%%%%%%%%%%%
%% Get output signal based on model structure
%
% ic:: Initial value for y, can be an array [y(1), y(2), ...]
% in_sig:: input signal array
% theta:: estimative model parameters
% m:: model structure
%%%%%%%%%%%%%%%%%%%%

N=max(size(in_sig));
y=zeros(N, 1);
% init the y array
for j=1: size(ic,2)
    y(j)=ic(j);
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
            yu=in_sig(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 3
            yu=aux_sig1(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 4
            yu=aux_sig2(k-abs(m.regr(i)))^m.texp(i);
        else
            error('invalid option, just y(1) u(2) and r(3) are possible')
        end
        % non linearity is yu^a*y^b
        if m.yplus_yur(i) == 1
            yu=yu*y(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 2
            yu=yu*in_sig(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 3
            yu=yu*aux_sig1(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 4
            yu=yu*aux_sig2(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        end
        num=num+theta(i)*yu;
    end
    den = 1;
    for i=m.n_dim+1:m.dim
        if m.yu(i) == 1
            yu=y(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 2
            yu=in_sig(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 3
            yu=aux_sig1(k-abs(m.regr(i)))^m.texp(i);
        elseif m.yu(i) == 4
            yu=aux_sig2(k-abs(m.regr(i)))^m.texp(i);
        end
        
        % non linearity is yu^a*y^b
        if m.yplus_yur(i) == 1
            yu=yu*y(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 2
            yu=yu*in_sig(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 3
            yu=yu*aux_sig1(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        elseif m.yplus_yur(i) == 4
            yu=yu*aux_sig2(k-abs(m.yplus_regr(i)))^m.yplus_exp(i);
        end
        den=den+theta(i)*yu;
    end
%    err=0;
%    if max(size(theta)) > m.dim
%        for h=m.dim+1:max(size(theta))
%            err=err+aux_sig2(k)*theta(h);
%        end
%    end
    
    y(k)= num/den;
end
end
