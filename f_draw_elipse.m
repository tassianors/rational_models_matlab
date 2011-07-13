function f_draw_elipse(vetA, vetB)

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
title('Estimativa dos parametros do denominador')
xlabel('Valor da estimativa para a variavel A')
ylabel('Valor da estimativa para a variavel B')
legend('Estimativas', 'Media')

% chi^2 for 95% of confiability
chi = 5.991;
ang = linspace(0,2*pi,360)';
[avetor,SCR,avl] = princomp(PN);
Diagonal= diag(sqrt(chi*avl));
elipse=[cos(ang) sin(ang)] * Diagonal * avetor' + repmat(mean(PN), 360, 1);
line(elipse(:,1), elipse(:,2), 'linestyle', '-', 'color', 'k');
end
