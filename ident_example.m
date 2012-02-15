%% Ident sample
% http://www.mathworks.com/products/sysid/demos.html?file=/products/demos/shipping/ident/iddemo6.html
close all; clear all;
clc;

%u = sign(randn(200,1)); % 1 inputs
%y = randn(200,1);       % 1 output
%ts = 0.1;               % The sampling interval
%z = iddata(y,u,ts);
%plot(z(:,1,1)) %Data subset with Input 1 and Output 1.

%u = z.u;   % or, equivalently u = get(z,'u');
%y = z.y;   % or, equivalently y = get(z,'y');
%set(z,'InputName',{'Voltage';'Current'},'OutputName','Speed');
%get(z)


z2 = iddata(rand(200,1),ones(200,1),0.1,'OutputName','New Output', 'InputName','New Input');
   
u = idinput([30 1 10],'sine'); % 10 periods of 30 samples
u = iddata([],u,1,'Period',30) % Making the input an IDDATA object.
m = idpoly([1 -1.5 0.7],[0 1 0.5]);  % This creates a model; see below.
y = sim(m,u,'noise') % y is the simulated response produced as an iddata object

plot(y)
m1 = armax(z2,[2 2 2 1],'maxiter',5,'search','LM'); % max 5 iterations, using the Levenberg-Marquard search direction
%m1 = armax(z2,[2 2 2 1],'Display','On');
compare(z2,m1)
[num,den]  = tfdata(m1,'v')
tfm = tf(m1)
tfm = tf(m1,'m') % 'm' for 'measured'.
   