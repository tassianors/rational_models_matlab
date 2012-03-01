function [ma stda mb stdb] = f_draw_elipse3d(vetA, vetB, vetC, realA, realB, realC)
%% Plots a elipse with 95% of confiability based on data in vetA and vetB
% vetA:: data vector of variable A
% vetB:: data vector of variable B
%%
figure;
[N, M]=size(vetA);
if M > N
    PN=[vetA', vetB', vetB'];
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
plot3(ma, mb, mc, 'rx');
%if realA ~= 0 && realB ~= 0
%	hold;
% 	plot(realA, realB, 'kd');
%end
title('Estimativa dos parametros para o sistema ARX')
xlabel('Estimativa Rho 1')
ylabel('Estimativa Rho 2')
zlabel('Estimativa Rho 3')
legend('Estimativas', 'Media', 'Real')
Ellipse_plot(inv(cov(PN))/7.81,mean(PN)')
hold off
grid;
end
