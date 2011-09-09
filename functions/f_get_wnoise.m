function noise = f_get_wnoise(y, np)
%% apply white noise over y signal.

if max(size(y)) <= 1
    error('Parameter y must be an array!!');
end
noise_= rand(size(y, 1),size(y, 2));
n_mean = mean(noise_);
noise = y+y.*+(noise_-n_mean/2)*(mean(y)/200*np);
end