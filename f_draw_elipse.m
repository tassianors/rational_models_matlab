function f_draw_elipse(vetA, vetB, realA, realB)
%% Plots a elipse with 95% of confiability based on data in vetA and vetB
% vetA:: data vector of variable A
% vetB:: data vector of variable B
%%
figure;
[N, M]=size(vetA);
if M > N
	PN=[vetA', vetB'];
else
	PN=[vetA, vetB];
end
ma=mean(vetA);
mb=mean(vetB);

% from here is only to plot the estimated points
plot(vetA, vetB, 'bo');
hold;
plot(ma, mb, 'rx');
hold;
if realA ~= 0 && realB ~= 0
	hold;
 	plot(realA, realB, 'kd');
% 	hold;
% 	hold
% 	plot([realA realB], [ma mb])
end
title('Estimativa dos parametros do denominador')
xlabel('Valor da estimativa para a variavel A')
ylabel('Valor da estimativa para a variavel B')
legend('Estimativas', 'Media', 'real')

% chi^2 for 95% of confiability
chi = 5.991;
ang = linspace(0,2*pi,360)';
[avetor,SCR,avl] = princomp(PN);
Diagonal= diag(sqrt(chi*avl));
elipse=[cos(ang) sin(ang)] * Diagonal * avetor' + repmat(mean(PN), 360, 1);
line(elipse(:,1), elipse(:,2), 'linestyle', '-', 'color', 'k');
end
