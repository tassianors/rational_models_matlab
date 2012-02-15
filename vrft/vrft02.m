%% tese luciola: modelo BJ

close all; clear all;
clc;
P=path;
path(P,'../functions')
P=path;
path(P,'../functions/signal')
P=path;
path(P,'../functions/plots')


M=10;
x=0.5;
y=0.9;
o=0.3
Ts=1;
exper = 100;
%% BJ system: 
% G_0(z)=x/(z-y)
% H_0(z)=z/(z-0.3)

% M(z)=0.4/(z-0.6)

%% Ideal Controler
% C_0(z)=(0.8(z-0.9))/(z-1)
% y(t)=G_0(z)*u(t)+H_0(z)*e(t)

% den of G
model.a = [1 -y];
% num of G
model.b = [x];
model.c = [1 -o];
model.d = [1 0];
model.mn = [0.4];
model.md = [1 -0.6];
model.TS = Ts;
model.delay = 1;
model.delay_func = tf([1],[1 0], model.TS);
model.noise_std = 0.002;

theta = zeros(exper, 3);
[u N]=f_get_prbs(M);
for i = 1: exper

    el=f_get_vrft_el(model, u);
    %% Controller model
    % u(t)=u(t-1)+0.8e(t)-0.72e(t-1)
    mc.eul = [0 1 1];
    mc.regr = [1 0 1];
    mc.dim = 3;
    mc.N = N;

    theta(i,:) = f_calc_mmq_theta(mc, u, el);
end

expect= [1 0.8 -0.72];

f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));
%f_draw_elipse(theta(:,3), theta(:,4), expect(3), expect(4));
