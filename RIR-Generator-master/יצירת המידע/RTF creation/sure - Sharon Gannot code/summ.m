cr_in = [-5.6 -2.7 -0.7   1.2  4.3  7.2 10.2 14.1];
cr_bf = [ 6.1  9.1  11.3 13.4 18.2 21.2 23.8 26.8];
cr_ot = [16.5 20.8 23.1 25.6 30.1 33.0 35.6 39.0];

figure(1)
plot(cr_in,cr_in,'-+',cr_in,cr_bf,'-*',cr_in,cr_ot,'-O')
xlabel('input SNR [dB]')
ylabel('output SNR [dB]')
grid('on')
legend('input','w/o p. proc.','w. p. proc.')
print -f1 -depsc2 '/home/nesher/elect/sharong/thesis/phd/beam_stat/ver3/snr_cr.eps'


wt_in = [ -5.9  -3.0  -1.0  3.9  7.0  9.0  11.0 13.9];
wt_bf = [  7.3  10.3  12.3 17.4 19.1 16.1  17.8 26.7];
wt_ot = [ 16.2  23.8  26.0 32.0 33.6 30.6  32.3 40.4];

figure(2)
plot(wt_in,wt_in,'-+',wt_in,wt_bf,'-*',wt_in,wt_ot,'-O')
xlabel('input SNR [dB]')
ylabel('output SNR [dB]')
grid('on')
legend('input','w/o p. proc.','w. p. proc.')
print -f2 -depsc2 '/home/nesher/elect/sharong/thesis/phd/beam_stat/ver3/snr_wt.eps'

df_in = [-5.0 -2.3 -0.5 1.4 4.3 7.3 10.2 14.1];
df_bf = [-2.0  0.8  2.5 2.4 6.9 9.9 13.2 17.5];
df_ot = [-1.2  6.3 12.6 14.2 18.5 22.4 26.4 31.6];

uc_in = [-5.7   -2.7   -0.8    1.2    4.2    7.2   10.2 14.1];
uc_bf = [-2.3    0.6   2.7    4.8    7.9   10.9   13.8   17.8];
uc_ot = [-1.6    4.7   11.8   15.8   21.0   24.8   28.1 32.2];
