function [el yy rl] = f_get_vrft_e_nl(model, u, y)

m=tf(model.mn,model.md, model.TS);
N_orig=max(size(y));

minv=inv(m)*model.delay_func;
r=lsim(minv, y);
rl=r(model.delay+1:max(size(r)));
el=rl-y(1:N_orig-1);
yy=y(1:N_orig);
end
