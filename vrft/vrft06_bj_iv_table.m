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

M=tf(model.mn,model.md, model.TS);
L=(1-M)*M;
beta=[tf([1 0],[1 -1], model.TS); tf([1],[1 -1], model.TS)];

%Instrumental Variables
IV=beta*L*(inv(M)-1);
expect= [0.8 -0.68];

theta = zeros(exper, rho_size);
table_var=[0.10 0.08 0.06 0.05 0.04 0.01  0.008 0.005  0.003 0.001];
results=zeros(max(size(table_var)), 5);

iv=0;
str = sprintf('var \t Jvr \t Jy \t theta1 \t theta2 ');
disp(str);
for t=1:max(size(table_var))
    
    model.noise_std = table_var(t);
    [u N]=f_get_prbs(m);
    
    if iv == 0
    for i = 1: exper
        [el y] = f_get_vrft_el(model, u);
        phy2=lsim(beta, el);
        phy=phy2(cut:max(size(phy2)),:,1);
        theta(i,:)=inv(phy'*phy)*phy'*u(cut:max(size(u)));
    end
else
    % iv method
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
end
    %f_draw_elipse(theta(:,1), theta(:,2), expect(1), expect(2));
    C=tf(mean(theta),[1 -1], model.TS);
    Jvr=f_get_vrft_Jvr(C, el, u);
    Jmr=f_get_vrft_Jmr(C, model);

results(t,1)=table_var(t);
results(t,2)=Jvr;
results(t,3)=Jmr;
results(t,4:5)=mean(theta);
str = sprintf('\t %s \t & \t %s \t & \t %s \t & \t %s \t & \t %s \\\\ ',num2str(results(t,1)), num2str(results(t,2), '%10.5e'), num2str(results(t,3), '%10.5e'), num2str(results(t,4), 5), num2str(-results(t,5),5));
disp(str);
end



