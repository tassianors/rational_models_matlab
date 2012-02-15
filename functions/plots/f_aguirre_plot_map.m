function a = f_aguirre_plot_map(source, index)
figure(index);
m = max(size(source));
for i = 1: m
   x(i) = source(i);
   if i < m
       y(i) = source(i+1);
   else
       y(i) = source(i);
   end
end
plot(x, y, 'b.');
end