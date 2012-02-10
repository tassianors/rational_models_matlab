function j = f_get_vrft_j_u(u, r, model, theta)

j=0;
N=max(size(u))

out=f_y_model(0, r, theta, model);

for k=1:N
   j=j+abs(u(k)-out(k))^2
end
j
end
