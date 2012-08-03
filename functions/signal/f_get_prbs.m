function [sig N] = f_get_prbs(s)
%====================
%% Get PRBS signal with order 7
% s:: Number of series repetitions
%
% Return:: 
%         sig:: PRBS Signal
%         N:: output signal length
%====================

M=127;
N=M*s;
u=zeros(N, 1);
% initialization (arbitrary value)
u(1)=0;
u(2)=1;
u(3)=1;
% get prbs signal from equation
for j=8: N
    u(j) = rem(u(j-3)+u(j-7), 2);
end
% put it with zero mean
sig= u.*2-1;
end