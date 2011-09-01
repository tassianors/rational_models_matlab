function [phy phym] = f_get_phy(y, m)
%% Gets the phy based on the model structure
% y:: output system data [y1,..,yn]
% m:: model
%%

f_check_model(m);

% initialization
N=max(size(y,1));
phy=zeros(m.dim+m.err_model,m.dim+m.err_model);
phym=zeros(m.dim+m.err_model,1);

for i=1:m.dim
    for j=1:m.dim
        if i>m.n_dim && j>m.n_dim
            aux=0;
            for k=1:N
                aux=aux+((y(k)^m.texp(j))*(y(k)^m.yplus_exp(j))* (y(k)^m.texp(i))*(y(k)^m.yplus_exp(i)));
            end
            phy(i,j)=aux;
        end
    end
end
for i=m.n_dim+1:m.dim
    aux=0;
    for k=1:N
        aux=aux-((y(k)^m.texp(i))*(y(k)^m.yplus_exp(i)));
    end
    phym(i,1)=aux;
end

end
