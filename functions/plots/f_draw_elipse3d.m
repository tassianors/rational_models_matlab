function [ma stda mb stdb] = f_draw_elipse3d(vetA, vetB, vetC, realA, realB, realC)
%% Plots a elipse with 95% of confiability based on data in vetA and vetB
% vetA:: data vector of variable A
% vetB:: data vector of variable B
%%
if max(size(vetA)) <= 1 || max(size(vetB)) <= 1 || max(size(vetC)) <= 1
    warning('To plot something you must pass some array, not just a value ;)');
    return
end

mark_size=8;
figure;
[N, M]=size(vetA);
if M > N
    PN=[vetA', vetB', vetC'];
else
	PN=[vetA, vetB, vetC];
end
ma=mean(vetA);
mb=mean(vetB);
mc=mean(vetC);
stda=std(vetA);
stdb=std(vetB);
stdc=std(vetC);


% from here is only to plot the estimated points
plot3(vetA, vetB, vetC, 'bo');
hold on;
plot3(ma, mb, mc, 'kp', 'MarkerSize',mark_size, 'MarkerFaceColor', 'k');

if realA ~= 0 && realB ~= 0 && realC ~= 0 
    plot3(realA, realB, realC, 'ks', 'MarkerSize',mark_size, 'MarkerFaceColor', 'k');
end

title('Estimativa dos parametros para o sistema ARX', 'FontSize',11);
xlabel('\theta_1', 'FontSize',11);
ylabel('\theta_2', 'FontSize',11);
zlabel('\theta_3', 'FontSize',11);
legend('Estimativas', 'Media', 'Real');

% chi^2 for 95% of confiability
chi = 5.991;

Ellipse_plot(inv(cov(PN))/chi,mean(PN)')
hold off
grid;
end
