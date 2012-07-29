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
cut=1;
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

C_den=[1 -1 0];

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
beta=[tf([1 0 0], C_den, model.TS); tf([1 0],C_den , model.TS);tf([1],C_den , model.TS)];
%Instrumental Variables
IV=beta*L*(inv(M)-1);
expect= [0.4 -0.68 0.288];


%table_var=[0.10 0.08 0.06 0.05 0.04 0.01  0.008 0.005  0.003 0.001];
table_var=[0.05];
results=zeros(max(size(table_var)), 6);

str = sprintf('var \t Jvr \t Jy \t theta1 \t theta2 \t theta3');
disp(str);
for t=1:max(size(table_var))
    
    model.noise_std = table_var(t);
    theta = zeros(exper, rho_size);
    [u N]=f_get_prbs(m);
    for i = 1: exper
        [el y] = f_get_vrft_el(model, u);
        [el2 y2] = f_get_vrft_el(model, u);
        ul=lsim(L,u);
        phy2=lsim(IV, y);
        instr2=lsim(IV, y2);
        phy=phy2(cut:max(size(phy2)),:,1);
        instr=instr2(cut:max(size(instr2)),:,1);
        theta(i,:)=inv(instr'*phy)*instr'*ul(cut:max(size(ul)));
    end
    covariance =cov(theta);
    
    f_draw_elipse3d(theta(:,1), theta(:,2), theta(:,3), expect(1), expect(2), expect(3));

    C=tf(theta(1,:),C_den, model.TS);
    Jvr=f_get_vrft_Jvr(C, el, u);
    Jmr=f_get_vrft_Jmr(C, model);

    results(t,1)=table_var(t);
    results(t,2)=Jvr;
    results(t,3)=Jmr;
    results(t,4:6)=mean(theta);
    %              1          2          3          
    str = sprintf('%s \t & \t %s \t & \t %s \t & \\;[ %s \\; %s \\; %s ] \\\\ ',num2str(results(t,1)), num2str(results(t,2), '%10.5e'), num2str(results(t,3), '%10.5e'), num2str(results(t,4), 5), num2str(results(t,5),5), num2str(results(t,6),5));
    disp(str);
end
