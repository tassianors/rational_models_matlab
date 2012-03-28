function jvr = f_get_vrft_nl_Jvr(uc0, uc)


j=0;
N=max(size(uc));
if N <= 1
    error('Parameter y must be an array!!');
end
if size(uc0) ~= size(uc)
    error('both u signals must have the same size');
end

for i=1:N
    j= j+abs(uc0(i)-uc(i))^2;    
end
jvr=j/N;
end