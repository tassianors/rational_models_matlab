function jmr = f_get_vrft_Jmr(C, model)

G=tf(model.b,model.a, model.TS);
M=tf(model.mn,model.md, model.TS);
T=feedback(C*G, 1);

jmr=norm(M-T,2);
end