%% Very simple VRFT example
% PID controller in a ARX system
close all; clear all;
clc;
P=path;
path(P,'../functions')
P=path;
path(P,'../functions/signal')
P=path;
path(P,'../functions/plots')
P=path;
path(P,'../functions/vrft')

rho_size=3;

m=10;
x=0.8;
y=0.9;
Ts=1;
exper = 100;
cut=1
%% ARX system: 
% G_0(z)=z/((z-x)(z-y))
% H_0(z)=z^2/((z-0.9)(z-0.8))
% M(z)=0.4/(z-0.6)

%% Ideal Controler
% C_0(z)=(0.4(z-0.9)(z-0.8))/(z(z-1))
% y(t)=G_0(z)*u(t)+H_0(z)*e(t)

model.a = [1 -(x+y) x*y]; 
model.b = [1 0];
model.c = [1 -(x+y) x*y];
model.d = [1 0 0];
model.mn = [0.4];
model.md = [1 -0.6];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.005;

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);

C_den=[1 -1 0];
beta=[tf([1 0 0], C_den, model.TS); tf([1 0],C_den , model.TS);tf([1],C_den , model.TS)];

for i = 1: exper
    [el y]=f_get_vrft_el(model, u);
    phi2=lsim(beta, el);
    phi=phi2(cut:max(size(phi2)),:,1);
    theta(i,:)=inv(phi'*phi)*phi'*u(cut:max(size(u)));
end

expect= [0.4 -0.68 0.288];

f_draw_elipse3d(theta(:,1), theta(:,2), theta(:,3), expect(1), expect(2), expect(3));

C=tf(mean(theta),C_den, model.TS);
Jvr=f_get_vrft_Jvr(C, el, u)
Jmr=f_get_vrft_Jmr(C, model)
variance =var(theta);
variance
