function [el yy rl] = f_get_vrft_e_nl(model, u, y)
%====================
%% Get VRFT virtual signals for Non-linear systems
% model:: Model structure
% u:: input signal system
% y:: output signal system
%
% return: Virtual signals:
%       el : controller input signal
%       yy: output signal with same size as el
%       rl: reference signal
%====================

m = tf(model.mn, model.md, model.TS);
N_orig = max(size(y));

minv = inv(m)*model.delay_func;
r = lsim(minv, y);
% shift rl. this is necessary becouse we used model.delay_func
rl = r(model.delay+1:max(size(r)));
% calculate error signal
el = rl-y(1:N_orig-1);
% truncate output signal to be sabe size as others
yy=y(1:N_orig);

end
