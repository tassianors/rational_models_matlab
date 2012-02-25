e25 = 2*prbs(7) - 1;
tempo25 = [0:1:60]';
prbs_ref = [tempo25, e25(1:61)'];
sim('sistema_inic_maberta');

%planta G = a/(z-b)(z-c)
a = 10;
b = 0.2;
c = 0.3;

z = tf('z',1);
G = a*tf(1,[1 -b],z)*tf(1,[1 -c],z);
%rltool(G)

%pi  C = (k1z - k2)/(z-1)
k1 = 0.0219;
k2 = -k1*0.364;

C = tf([k1 -k2],[1 -1],z);

%modelo M = (1-f)(1-g)/(z-f)(z-g)
f = 0.5;
g = 0.3;

M = tf((1-f)*(1-g),[1 -(f+g) f*g],z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%metodo vrft
%passo 1
n = size(u_scope,1);


ym = y_scope(:,2);
um = u_scope(:,2);
%um = u_scope(2:62,2);

%em = e_scope(:,2);
tm = tempo25;
refm = ref_scope(:,2);

%calculo da referencia
i = 3;
for i=3:size(ym,1)
    ref_aux(i-2) = (1/0.35)*ym(i) - (0.8/0.35)*ym(i-1) + (0.15/0.35)*ym(i-2);
end
ref2 = [ref_aux, ref_aux(size(ym,1)-2), ref_aux(size(ym,1)-2)];


figure(1)
pcov(um,1)


e = ref2'-ym;


%passo 2
%filtro L(z) = (1 - M(z))M(z)Gesp(z)W
W = 1;
L1 = (1 - M)*M;
%Gesp = 1/(z-0.1);

polo = 0.7;
Gesp =1;%0.1/(z-polo);

yt = step(Gesp,n);

figure(2)
pcov(yt,1)


L = L1*1/Gesp*W;

%%Y = LSIM(SYS,U,T)
 el = lsim(L,e,tm);
 ul = lsim(L,um,tm);


%L = 1;
% el = e;
% ul = um;



%least square method
% phi(k-1) = [u(k-1) u(k-2)]

%beta = [1;tf([1 0],[1 -1],1)];

% controlador pi
beta2 = z/(z-1);
beta3 = 1/(z-1);
beta = [beta2;beta3];

% controlador otimo
% beta1 = z^2/(z^2-0.8*z-0.2);
% beta2 = z/(z^2-0.8*z-0.2);
% beta3 = 1/(z^2-0.8*z-0.2);

% beta = [beta1;beta2;beta3];


phi = lsim(beta,el,tm);
teta = inv(phi'*phi)*phi'*ul
    





