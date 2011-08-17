function noise = f_get_wnoise(y, np)
%% apply white noise over y signal.

N=max(size(y));
if N <= 1
    error('Parameter y must be an array!!');
end
noise_= rand(N,1);
n_mean = mean(noise_);
noise = y+y.*+(noise_-n_mean/2)*(mean(y)/200*np);
end