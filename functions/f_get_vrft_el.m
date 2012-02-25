function [el y] = f_get_vrft_el(model, u)

g=tf(model.b,model.a, model.TS);
h=tf(model.d,model.c, model.TS);
m=tf(model.mn,model.md, model.TS);
N=max(size(u));

e=f_get_noise_signal(N, model.noise_std);
yu=lsim(g,u);
ye=lsim(h,e);
y=yu+ye;

minv=inv(m)*model.delay_func;
r=lsim(minv, y);
rl=r(model.delay+1:size(r,1));
rl(size(r,1))=0;
el=rl-y;

end
