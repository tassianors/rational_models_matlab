function jmr = f_get_vrft_nl_Jmr(M, in_sig, out_sig)

ym=lsim(M, in_sig);
jmr=norm(ym-out_sig',2);
end