function f_plot_y_y1(vec)
%% plot the array #vec, where y(k) is X and y(k+1)

N=max(size(vec));
if N <= 1
    error('Parameter must be an array!!');
end

yk=zeros(N-1,1);
yk1=zeros(N-1,1);

yk1(1)=vec(1);
yk(1)=vec(1);
for k=2:N
    yk(k-1)=vec(k-1);
    yk1(k-1)=vec(k);
end
figure;
plot(yk, yk1, '.');
end

