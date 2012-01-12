function y = f_aguirre_get_model_output(model, simul, theta, init_val)

y=zeros(simul.N, 1);
num = 0;
den = 1;

y(1)=init_val;
for k = abs(max(model.regr))+1: simul.N
    num = 0;
    den = 1;
    for d=1: model.dim
       if d <= model.n_dim
           num = num + theta(d)*y(k-abs(model.regr(d)))^(model.texp(d));
       else
           den = den + theta(d)*y(k-abs(model.regr(d)))^(model.texp(d));
       end
    end
    y(k) = num/den;
end
