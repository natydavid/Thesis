disp('RESULTS_fr_cr_-10.mat')

SIGNAL =  'Timit' % 'Timit' 'Timit_br'
NOISE =   'Dir_Crown' % 'Dir_White_Gauss''Dir_Crown''Diffuse''Ucorr'
M = 5;  % No. of Microphones
SNR =    5*ones(M,1); 
snr = 100;


[z,T_sg,vst] = construct(SIGNAL,NOISE,SNR,snr,M);
[yO,h,g] = algorithm(z,T_sg,vst);
zz = z(1,:);
%yMX = ncanSH1(yO');
%zMX = ncanSH1(zz');
%save RESULTS_fr_cr_-10.mat zz h g yO yMX

% --------------------------------------------------
