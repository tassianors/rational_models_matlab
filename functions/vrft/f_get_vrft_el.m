function [el yy] = f_get_vrft_el(model, u)
%====================
%% Get VRFT virtual signals for linear systems
% model:: Model structure
% u:: input signal system
%
% return: Virtual signals:
%       el : controller input signal
%       yy: output signal with same size as el
%====================

% plant transfer function
g=tf(model.b,model.a, model.TS);
% noise filter transfer function
h=tf(model.d,model.c, model.TS);
% desired closed loop behavior
m=tf(model.mn,model.md, model.TS);
N_orig=max(size(u));

%Auxiliary signal MUST be 'delay_size' bigger than original u
Uaux=u;
for k=N_orig+1:N_orig+model.delay
    Uaux(k)=1;    
end

e  = f_get_noise_signal(N_orig+model.delay, model.noise_std);
yu = lsim(g,Uaux);
ye = lsim(h,e);
% get ouput signal with noise
y  = yu+ye;

minv=inv(m)*model.delay_func;
r  = lsim(minv, y);
rl = r(model.delay+1:max(size(r)));
el = rl-y(1:N_orig);
yy = y(1:N_orig);
end
