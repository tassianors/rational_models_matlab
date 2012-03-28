function noise = f_get_wnoise(y, np)
%% apply white noise over y signal.

N=max(size(y));

if N <= 1
    error('Parameter y must be an array!!');
end

e=f_get_noise_signal(N, np);

noise=y+e;
end