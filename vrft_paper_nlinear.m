%========================================================================
% Exemplo do Paper de VRFT nao linear
% real system: y(k)=0.9y(k-1)+tanh(u(t-1))
%========================================================================
close all; clear all;
clc;
P=path;
path(P,'./functions')
%========================================================================
% init
%========================================================================
model.Ts=1;
model.N=800;
%model.dim=4;
%model.regr = [1 0 1 2];
%model.eul= [0 1 1 1];

% input signal
t = 0:1:model.N-1;
u = (square(0.01*pi*t)');
%========================================================================
% real system
%========================================================================
y=zeros(model.N, 1);
r_virt=zeros(model.N, 1);
y(1)=0;
for k=2:model.N
   y(k)=0.9*y(k-1)+tanh(u(k-1)); 
end

%========================================================================
% Model reference - this behavior is exactly what we wish
% y(t)=r(t-1)
% we can use here u(t) as r(t), it is in open loop
% r_virt(k)=y(k+1)
%========================================================================
% first calculate the virtual referece, where 
for k=1:model.N-1
   r_virt(k+1)=y(k);
end

%========================================================================
% input of controler = e(t)=r(t)-y(t)
e=r_virt-y;
plot(r_virt,'b.')
hold
plot(y, 'r.')
%========================================================================
% control identification part:
% select the model
% u(t)=u(t-1)+(0.2+\theta)e(t)-\theta e(t-1)

model.dim=3;
model.regr = [1 0 1];
model.eul= [0 1 1 ];
teta=f_calc_mmq_theta(model, u, e)


%simul=struct('N', model.N, 'nEstimates', 10, 'np', 0.5, 'maxError', 0.01, 'l', 100, 'diffConv', .1);
%ret=f_rational_model(simul, model, u, [u(1)], e)
m_rat.n_dim   = 3;
m_rat.dim     = 3;
% to indo do
m_rat.texp    = [1 1 1];
m_rat.yu      = [1 0 0];
m_rat.regr    = [1 0 1];
% tels if there is some non linearity like (y(k-a)^b)*(y(k-c)^d)
% u = 2 y=1 none =0
m_rat.yplus_uy = [0 0 0];
% tels the d param
m_rat.yplus_exp = [0 0 0];
% tels the C param
m_rat.yplus_regr = [0 0 0];

m_rat.err_m_rat   = 0;
m_rat.err_enable = false
%% Simulation parameters
simul=struct('N', model.N, 'nEstimates', 10, 'np', 10, 'maxError', 0.01, 'l', 100, 'diffConv', .1);
ret = f_rational_model(simul, m_rat, u, [u(1)], e)
a1=1;
a2=0.2;
b1=0;
f_draw_elipse(ret(:,1), ret(:,2), a1, a2);
f_draw_elipse(ret(:,1), ret(:,3), a1, b1);


