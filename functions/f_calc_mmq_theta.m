function theta = f_calc_mmq_theta(m, out, in)
%====================
%% Get theta parameters using LSM
% m:: Model structure
% out:: output signal system
% in:: input signal system
%====================

phi=zeros(m.N, m.dim);
out2=zeros(m.N, 1);

for k = (max(m.a_regress))+1:m.N
	for j=1:m.dim
		if m.eul(j) == 1
            phi(k, j) = in(k-(m.a_regress(j)));
		else
            phi(k, j) = out(k-(m.a_regress(j)));
		end
    end
    out2(k) = out(k);
end
theta = inv(phi'*phi)*phi'*out2;
end