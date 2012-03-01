clear;
P=path;
path(P,'../../functions/signal')

z = tf('z',1);

% G0 = 0.2*(z-0.6)/((z-0.9)*(z-0.5));
% H0 = 1;
% M = (1-0.6)^2*z/(z-0.6)^2;
% Cd = zpk(M/(G0-M*G0));
G0 = 0.5/(z-0.9);
H0 = z/(z-0.3);
M = (1-0.6)/(z-0.6);
Cd = zpk(M/(G0-M*G0));


%passo = 1000;

[e25 passo]=f_get_prbs(10);

%e25 = 2*prbs(12) - 1;
tempo25 = [0:1:passo-1]';
prbs_ref = [tempo25, e25];

TT=[];
TT2=[];
for b=1:10
b

seed = abs(floor(randn*1000));
r_var = 0.01;
%sim('ensaio_piD_arx');
t_amost = 1;


um = y_scope(1:passo,2);
ym = y_scope(1:passo,3);
tm = [0:1:size(ym,1)-1]';
n = size(y_scope,1);

U_fil = (z-1)*M/(1-M);
%U_fil = 0.16*z/(z-0.36);

util = lsim(U_fil,um,tm);

%y = Ctil*util + H*e
%Ctil = (k1*z^2 + k2*z + k3)/(z^3 + a1*z^2 + a2*z + a3)

% controlador pid
%beta = [1;tf([1 0],[1 -1],t_amost)];
%beta2 = [1;tf([1 -1],[1 0],t_amost)];

    DAT = IDDATA(ym,util,1);

    m=arx(DAT,[2 1 1]);
    %m2 = oe(DAT,[1 1 1])

    %m3 = inv(m2);
    [A,B,C,D,F]=polydata(m);
    teta = [B(2) A(2:3)]';
    TT=[TT teta];
end
% med = mean(TT')'
% va = var(TT')'

    Cc = 1/teta(1)*(z^2+teta(2)*z+teta(3))/(z^2-z);
    zpk(Cc)
    
    med = mean(TT')
    var = cov(TT')
%     H = z/(z+teta(3))
%     
% 
%     figure(1)  
%     step(feedback(G0*Cc,1))
%     hold
%     step(feedback(G0*Cd,1))
%     hold off
% 
% figure(2)
% TT_c = TT(1:2,:);
% hold on
% plot(TT_c(1,:),TT_c(2,:),'b*')
% plot(mean(TT_c(1,:)),mean(TT_c(2,:)),'r*')
% %axis([1.1 1.4 -1.00 -0.80])
% %plot(TT2(1,:),TT2(2,:),'m*')
% Ellipse_plot(inv(cov(TT_c'))/5.99,mean(TT_c')')
% hold off
% 
figure(4)
hold on
plot3(TT(1,:),TT(2,:),TT(3,:),'b*')
plot3(mean(TT(1,:)),mean(TT(2,:)),mean(TT(3,:)),'r*','LineWidth',2)
%axis([2.494 2.506 -1.702 -1.696 0.718 0.722])
xlabel('\rho_1','FontSize',11)
ylabel('\rho_2','FontSize',11)
zlabel('\rho_3','FontSize',11)
title('Parametros do controlador','FontSize',11)
% %plot(TT2(1,:),TT2(2,:),'m*')
Ellipse_plot(inv(cov(TT'))/7.81,mean(TT')')
hold off