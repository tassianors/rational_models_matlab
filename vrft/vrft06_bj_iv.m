% VRFT method using betas as a transfer function array
% example of luciola thesys page 90
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

% simulation parameters
rho_size=2;
m=20;
a=0.5;
b=0.85;
c=0.4;

Ts=1;
exper = 100;
cut=1;

%% BL system: 
% G_0(z)=a/(z-b)
% H_0(z)=z/(z-c)
% M(z)=0.4/(z-0.6)

%% Ideal Controler
% C_0(z)=(0.8(z-0.9))/(z-1)
% y(t)=G_0(z)*u(t)+H_0(z)*e(t)

model.a = [1 -b]; 
model.b = [a];
model.c = [1 -c];
model.d = [1 0];
model.mn = [0.4];
model.md = [1 -0.6];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.02;

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
beta=[tf([1 0],[1 -1], model.TS); tf([1],[1 -1], model.TS)];

%Instrumental Variables
IV=beta*L*(inv(M)-1);

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);
for i = 1: exper
    [el y] = f_get_vrft_el(model, u);
    [el2 y2] = f_get_vrft_el(model, u);
    ul=lsim(L,u);
    phi2=lsim(IV, y);
    instr2=lsim(IV, y2);
    phi=phi2(cut:max(size(phi2)),:,1);
    instr=instr2(cut:max(size(instr2)),:,1);
    theta(i,:)=inv(instr'*phi)*instr'*ul(cut:max(size(ul)));
end

variance =var(theta);
expect= [0.8 -0.68];
f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));

C=tf(mean(theta),[1 -1], model.TS);

Jvr=f_get_vrft_Jvr(C, el, u)
Jmr=f_get_vrft_Jmr(C, model)

variance





