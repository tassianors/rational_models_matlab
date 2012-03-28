function j = f_get_vrft_jVR(C, el, u)

j=0;
N=max(size(u))
ul2=lsim(C,el);

for i=1:N
    j= j+abs(u(i)-ul2(i))^2;    
end
j
end
