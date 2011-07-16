function psi = f_get_psi(y, u, m)
%% Gets the Psi based on the model structure
% y:: output system data [y1,..,yn]
% u:: innput system data [u1,..,un]
% m:: model
%%

% check parameters
f_check_model(m);

%% step 1 - first estimative
N=max(size(y));
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
