%% Very simple VRFT example
close all; clear all;
clc;
P=path;
path(P,'../functions')
P=path;
path(P,'../functions/signal')
P=path;
path(P,'../functions/plots')


M=10;
x=0.8;
y=0.9;
Ts=1;
exper = 10;
cut=150;
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
model.noise_std = 0.001;

theta = zeros(exper, 3);
[u N]=f_get_prbs(M);

fil= tf([1],[0.6 0], model.TS)

for i = 1: exper
    [el y]=f_get_vrft_el(model, u);
    y2=lsim(fil,y);
    %% Controller model
    % u(t)=0.4e(t)-0.68e(t-1)-0.288e(t-2)
    mc.eul = [1 1 1];
    mc.regr = [0 1 2];
    mc.dim = 3;
    mc.N = N-cut;

    theta(i,:) = f_calc_mmq_theta(mc, u(cut:max(size(u))), y2(cut:max(size(y2))));
end

expect= [0.4 -0.68 0.288];
theta
f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));
%f_draw_elipse(theta(:,3), theta(:,4), expect(3), expect(4));
