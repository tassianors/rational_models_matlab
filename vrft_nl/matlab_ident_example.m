A = zeros(2,2,3);
B = zeros(2,1,3)
A(:,:,1) =eye(2);
A(:,:,2) = [-1.5 0.1;-0.2 1.5];
A(:,:,3) = [0.7 -0.3;0.1 0.7];
B(:,:,2) = [1;-1];
B(:,:,3) = [0.5;1.2];
m0 = idarx(A,B,1);

y=ones(1,100);
u=ones(1,100);
%[0 0 0],'linear','CustomReg',{'y1(t-1)^2','y1(t-2)*u1(t-3)'})
%m1=idnlarx(m0,Nonlinearity) 
m1=idnlarx(m0,'wavenet', 'CustmRegressors',{'y1(t-10)*u1(t-1)', 'y(t-2)*u(t-2)'})
