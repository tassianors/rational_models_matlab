% VRFT method using betas as a transfer function array
% example of luciola thesys page 92 
% where the controller function isn't inside controller class

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
rho_size=3;
m=20;
a=0.7;
b=0.9;
c=0.3;
d=0.2;
e=0.5;

Ts=1;
exper = 100;
cut=1;

%% BL system: 
% G_0(z)=d(z-a)/((z-b)(z-e))
% H_0(z)=z/(z-c)
% M(z)=0.36z/(z-0.6)^2

%% Ideal Controler
% C_0(z)=0.8z(z-0.9)(z-0.5)
%       (z-1)(z-0.36)(z-0.7)
% y(t)=G_0(z)*u(t)+H_0(z)*e(t)

model.a = [1 -(b+e) b*e]; 
model.b = [d -a*d];
model.c = [1 -c];
model.d = [1 0];
model.mn = [0.16 0];
model.md = [1 -1.2 0.36];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.05;
C_den=[1 -1 0];

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
beta=[tf([1 0 0], C_den, model.TS); tf([1 0],C_den , model.TS);tf([1],C_den , model.TS)];

%Instrumental Variables
IV=beta*L*(inv(M)-1);

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);

for i = 1: exper
    [el y]=f_get_vrft_el(model, u);
    phy2=lsim(beta, el);
    phy=phy2(cut:max(size(phy2)),:,1);
    theta(i,:)=inv(phy'*phy)*phy'*u(cut:max(size(u)));
end

variance =var(theta);
expect= [0.8091   -0.1666   -0.3358];

CL=tf(expect,C_den, model.TS)
C=tf(mean(theta),C_den, model.TS)
Jvr=f_get_vrft_Jvr(C, el, u)
Jmr=f_get_vrft_Jmr(C, model)

JvrL=f_get_vrft_Jvr(CL, el, u)
JmrL=f_get_vrft_Jmr(CL, model)

G=tf(model.b,model.a, model.TS);
Cd=zpk([0 0.9 0.5],[1 0.36 0.7],0.8, model.TS);
TL=feedback(CL*G, 1);
T=feedback(C*G, 1);

step(M, T, TL)
%f_draw_elipse3d(theta(:,1), theta(:,2), theta(:,3), expect(1), expect(2), expect(3));
figure
bode(Cd, C, CL)
legend;

