function theta = f_calc_mmq_theta(m, ul, el)

phy=zeros(m.N, m.dim);

for k=abs(max(m.regr))+1:m.N
	for j=1:m.dim
		if m.eul(j) == 1
			phy(k, j)=el(k-abs(m.regr(j)));
		else
			phy(k, j)=ul(k-abs(m.regr(j)));
		end
	end
end
theta=inv(phy'*phy)*phy'*ul;
end