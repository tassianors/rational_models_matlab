function f_plot_y_y1(vec)
% funcion util
N=size(vec,1);
yk=zeros(N-1,1);
yk1=zeros(N-1,1);
yk1(1)=vec(1);yk(1)=vec(1);
for k=2:N
    yk(k-1)=vec(k-1);
    yk1(k-1)=vec(k);
end
figure;
plot(yk, yk1, '.');
end

