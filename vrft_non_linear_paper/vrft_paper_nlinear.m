%========================================================================
% Exemplo do Paper de VRFT nao linear
% real system: y(k)=0.9y(k-1)+tanh(u(t-1))
%========================================================================
close all; clear all;
clc;
P=path;
path(P,'../functions')
%========================================================================
% init
%========================================================================
model.Ts=1;
model.N=400;
%model.dim=4;
%model.a_regress = [1 0 1 2];
%model.eul= [0 1 1 1];

% input signal
t = 0:1:model.N-1;
u =(square(0.01*pi*t)');
%u = ones(model.N, 1)*0.5;
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
% plot(r_virt,'b.')
% hold
% plot(y, 'r.')
%========================================================================
% control identification part:
% select the model
% u(t)=u(t-1)+(0.2+\theta)e(t)-\theta e(t-1)

model.dim=3;
model.a_regress = [1 0 1];
model.eul= [0 1 1 ];
teta=f_calc_mmq_theta(model, u, e)
u_out=zeros(model.N, 1);
% calc the output obtained by controler structure found
for k=2:model.N-1
    u_out(k)=teta(1)*u_out(k-1)+teta(2)*e(k)-teta(3)*e(k-1)
end
plot(u_out)
hold
plot(u)
