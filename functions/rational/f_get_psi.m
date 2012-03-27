function psi = f_get_psi(y, yp, u, r, m)
%% Gets the Psi based on the model structure
% y:: output system data [y1,..,yn]
% u:: innput system data [u1,..,un]
% m:: model
%%

% check parameters
f_check_model(m);

%% step 1 - first estimative
N=max(size(y));
psi=zeros(N, m.dim+m.err_model);

if m.err_model && max(size(y)) ~= max(size(yp))
    error('y and last y (yp) have to have the same size');
else
    err=y-yp;
end

for i=max(abs(m.regr))+1:N
    for j=1:m.dim+m.err_model
        % err_model
        if j> m.dim
           % fill the err_model colunm.
           psi(i,j)=err(i);
           continue;
        end
        
        if m.yu(j) == 1
            yu=y(i-abs(m.regr(j)))^m.texp(j);
        elseif m.yu(j) == 2
            yu=u(i-abs(m.regr(j)))^m.texp(j);
        elseif m.yu(j) == 3
            yu=r(i-abs(m.regr(j)))^m.texp(j);
        else
            error('not supported option');
        end
        % non linearity is yu^a*y^b
        yu2=1;
        if m.yplus_yur(j) == 1
            yu2=y(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        end
        % non linearity is yu^a*u^b
        if m.yplus_yur(j) == 2
            yu2=u(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        end
        if m.yplus_yur(j) == 3
            yu2=r(i-abs(m.yplus_regr(j)))^m.yplus_exp(j);
        end
        if j<=m.n_dim
            psi(i,j)=yu*yu2;
        else
            % here we should alway use y(i) (equation 10.44 aguirre)
            psi(i,j)=-yu*yu2*y(i);
        end

    end
end
end
