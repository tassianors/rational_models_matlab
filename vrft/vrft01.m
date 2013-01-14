%% Very simple VRFT example
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


M=10;
x=0.8;
y=0.9;
Ts=1;
exper = 100;
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

theta = zeros(exper, 4);
[u N]=f_get_prbs(M);
for i = 1: exper

    el=f_get_vrft_el(model, u);
    %% Controller model
    % u(t)=0.4e(t)-0.68e(t-1)-0.288e(t-2)+u(t-1)
    mc.eul = [1 1 1 0];
    mc.a_regress = [0 1 2 1];
    mc.dim = 4;
    mc.N = N;

    theta(i,:) = f_calc_mmq_theta(mc, u, el);
end

expect= [0.4 -0.68 0.288 1];

f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));
f_draw_elipse(theta(:,3), theta(:,4), expect(3), expect(4));
