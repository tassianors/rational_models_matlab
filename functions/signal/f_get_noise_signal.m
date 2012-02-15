function noise = f_get_noise_signal(size, stdVal)
tmp_noise = randn(size,1);
tmp_std = std(tmp_noise);
noise = tmp_noise*(stdVal/tmp_std);
end