function psi = f_get_psi(y, theta_exp, nn_nd_dim, nn_dim)
    %% step 1 - first estimative
    N=size(y,1);
    
    for i=2:N
        % step 3 - we can neglect the noise model (colun 6)
        for j=1:nn_nd_dim
            if j<=nn_dim
                psi(i,j)=y(i-1)^theta_exp(j);
            else 
                psi(i,j)=-y(i-1)^theta_exp(j)*y(i);
            end
        end
    end
end
