function [sig N] = f_get_prbs(size)
M=127
N=M*size
u=zeros(N, 1);
u(1)=0;
u(2)=1;
u(3)=1;
for j=8: N
    u(j)=rem(u(j-3)+u(j-7), 2);
end
sig= u.*2-1;
end