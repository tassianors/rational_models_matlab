%% Very simple VRFT example
close all; clear all;
clc;
P=path;
path(P,'../functions')

x=0.8;
y=0.9;
Ts=1;
%% ARX system: 
% G_0(z)=z/((z-x)(z-y))
% H_0(z)=z^2/((z-0.9)(z-0.8))
% M(z)=0.4/(z-0.6)

%% Ideal Controler
% C_0(z)=(0.4(z-0.9)(z-0.8))/(z(z-1))

%y(t)=G_0(z)*u(t)+H_0(z)*e(t)

model.a = [1 -(x+y) x*y];
model.b = [1 0];
model.c = [1 -(x+y) x*y];
model.d = [1 0 0];
model.mn = [0.4];
model.md = [1 -0.6];
model.TS = Ts;
model.delay = tf([1],[1 0], model.TS);

%N=250
%u=ones(N, 1);
%u=f_get_square_signal(N);
[u N]=f_get_prbs(5);
e=rand(N, 1)*0.01;
g=tf(model.b,model.a, model.TS);
h=tf(model.d,model.c, model.TS);
yu=lsim(g,u);
ye=lsim(h,e);
y=yu;%+ye;

m=tf(model.mn,model.md, model.TS);

minv=inv(m)*model.delay;
r=lsim(minv, y);
rl=r(2:size(r,1));
rl(size(r,1))=0;
el=rl-y;


%% Controller model
% u(t)=0.4e(t)-0.68e(t-1)-0.288e(t-2)+u(t-1)
mc.eul = [1 1 1 0];
mc.regr = [0 1 2 1];
mc.dim = 4;
mc.N = N;

f_calc_mmq_theta(mc, u, el)
plot(y)
expect= [0.4 -0.68 -0.288 1]