function theta = f_calc_mmq_theta(m, out, in)
%====================
%% Get theta parameters using LSM
% m:: Model structure
% out:: output signal system
% in:: input signal system
%====================

phi=zeros(m.N, m.dim);
out2=zeros(m.N, 1);

for k = abs(max(m.regr))+1:m.N
	for j=1:m.dim
		if m.eul(j) == 1
			phi(k, j) = in(k-abs(m.regr(j)));
		else
			phi(k, j) = out(k-abs(m.regr(j)));
		end
    end
    out2(k) = out(k);
end
theta = inv(phi'*phi)*phi'*out2;
end