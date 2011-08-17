function psi = f_get_psi(y, yp, u, m)
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
