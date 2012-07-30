function [x y] = f_aguirre_plot_map(source, index)
if index ~= 0
    figure(index);
end

m = max(size(source));
for i = 1: m
   x(i) = source(i);
   if i < m
       y(i) = source(i+1);
   else
       y(i) = source(i);
   end
end
if index ~= 0
    plot(x, y, 'b.');
end
end