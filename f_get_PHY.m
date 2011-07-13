function PHY = f_get_PHY(nd_nn_dim, nn_dim, y, theha_exp)
N=size(y,1);
PHY=zeros(nd_nn_dim,nd_nn_dim);
for i=1:nd_nn_dim
    for j=1:nd_nn_dim
        if i>nn_dim && j>nn_dim
            aux=0;
            for k=1:N
                aux=aux+((y(k)^theha_exp(j))* (y(k)^theha_exp(i)));
            end
            PHY(i,j)=aux;
        end
    end
end
end
