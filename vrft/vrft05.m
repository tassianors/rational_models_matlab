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

% simulation parameters
rho_size=3;
m=10;
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
model.mn = [0.36 0];
model.md = [1 -1.2 0.36];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.005;

theta = zeros(exper, rho_size);
[u N]=f_get_prbs(m);

beta=[tf([1 0 0],[1 -1 0], model.TS); tf([1 0],[1 -1 0], model.TS);tf([1],[1 -1 0], model.TS)];

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
%L=tf([1 0],[1 0], model.TS);

for i = 1: exper
    [e y] = f_get_vrft_el(model, u);
    ul=lsim(L, u);
    el=lsim(L, e);
    phy2=lsim(beta, el);
    phy=phy2(cut:max(size(phy2)),:,1);
    theta(i,:)=inv(phy'*phy)*phy'*ul(cut:max(size(ul)));
end

var(theta)
expect= [1.1268   -0.7204   -0.3448];

G=tf(model.b,model.a, model.TS);

CL=tf(theta(1,:),[1 -1 0], model.TS);
TL=feedback(CL*G, 1);

C=tf([1.1237 -0.7222 -0.3447],[1 -1 0], model.TS);
T=feedback(C*G, 1);

Cd=zpk([0 0.9 0.5],[1 0.36 0.7],0.8, model.TS)
Td=feedback(Cd*G, 1);

step(T, TL, Td)
f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));
