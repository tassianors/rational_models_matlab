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

rho_size=3;

M=10;
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
model.noise_std = 0.01;

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(M);

C_den=[1 -1 0];
beta=[tf([1 0 0], C_den, model.TS); tf([1 0],C_den , model.TS);tf([1],C_den , model.TS)];

for i = 1: exper
    [el y]=f_get_vrft_el(model, u);
    phy2=lsim(beta, el);
    phy=phy2(cut:max(size(phy2)),:,1);
    theta(i,:)=inv(phy'*phy)*phy'*u(cut:max(size(u)));
end

expect= [0.4 -0.68 0.288];

%f_draw_elipse3d(theta(:,1), theta(:,2), theta(:,3), expect(1), expect(2), expect(3));
f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));
f_draw_elipse(theta(:,1), theta(:,3), expect(1), expect(3));

C=tf(theta(1,:),C_den, model.TS);
Cd=tf(expect,C_den, model.TS);

Jvr=f_get_vrft_Jvr(C, el, u)

%f_plot_feedback_comp(tf(model.b,model.a, model.TS), C, Cd);