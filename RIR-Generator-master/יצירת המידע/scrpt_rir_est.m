
% chirp signals

[sig_source,fs]=audioread('wav_files_source\chirp_log_iter3_old.wav');
%pt [sig_source,fs]=audioread('wav_files_source\chirp3win_src.wav');
%[sig_source,fs]=audioread('wav_files_source\chirp_log_14kHz_iter3.wav');
% sig_record=audioread('recordings\340ms\340ms_loudspkr4_chirp.wav');
sig_record=audioread('2017-01-15__10-34-08.wav');
%sig_record=audioread('..\160317-six\2016-03-17__14-03-30__rir__chirp.wav');

% parameters
J=size(sig_record,2);
iter_num=3;
N=length(sig_source);
sig_record=sig_record(1:N,:);

% more parameters
Niter=N/iter_num;
Nfft=2*Niter;

% estimate the RIR for each iteration
H=zeros(Nfft,J);
for p=1:iter_num
    idxs=((p-1)*Niter+1:p*Niter).';
    sig_source_p=sig_source(idxs);
    sig_record_p=sig_record(idxs,:);
    H=H+fft(sig_record_p,Nfft)./repmat(fft(sig_source_p,Nfft),[1,J]);
end
load('misc\bp_filt_blackman_4000.mat'); % BP - from 80 to 7900 Hz, linear phase
H=real(ifft(H.*repmat(fft(Num,Nfft),[1,J])));

% evaluate the T60 using EDC
t60=zeros(J,1); for j=1:J, figure; t60(j)=edc(H(1:2*fs,j),fs,.144,.19); end

fprintf('%.3f +- %.3fs\n', mean(t60), std(t60));
