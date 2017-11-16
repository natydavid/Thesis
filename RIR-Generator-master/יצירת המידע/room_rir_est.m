function [ rir ] = room_rir_est( sig_source,sig_record,fs,n )
% room_rir_est estimates the rir using the source and a record of the source after going through a system,
% Expected Input:
% sig_source -  a matrix of dim X num_of_samples X channels (optional)
% sig_record -  same as sig_source
if length(size(sig_source))==3
    permute_order = [ 2 1 3];
else
    permute_order = [ 2 1];
end
flag = 0;
if size(sig_source,2)>size(sig_source,1), sig_source = permute(sig_source,permute_order); flag=1;   end;
if size(sig_record,2)>size(sig_record,1), sig_record = permute(sig_record,permute_order); end;

if ~exist('n','var')
    n = 3000;
end
order = 4000;
f0 = 200;
f1 = 7900;
w0 = f0/fs;
w1 = f1/fs;
% w0 = 0.0050;
% w1 = 0.4938;

Num = fir1(order,[w0 w1],[],blackman(order+1))'; %% BP using blackman window
% load('bp_filt_blackman_4000.mat'); % BP - from 80 to 7900 Hz, linear phase
% parameters
if length(size(sig_source))==3
    J=size(sig_record,3);
else
    J=1;  
end

if size(sig_record,2)>size(sig_source,2)
    sig_source = repmat(sig_source,1,size(sig_record,2),J);
end

iter_num=3;
N=size(sig_source,1);
sig_record = sig_record(1:N,:,:);
num_of_smp = size(sig_source,2);


% more parameters
Niter=floor(N/iter_num);
Nfft=2*Niter;

% estimate the RIR for each iteration
H=zeros(Nfft,num_of_smp,J);
for p=1:iter_num
    idxs=((p-1)*Niter+1:p*Niter).';
    sig_source_p=sig_source(idxs,:,:);
    sig_record_p=sig_record(idxs,:,:);
    H=H+fft(sig_record_p,Nfft)./repmat(fft(sig_source_p,Nfft),[1,1,J]);
end
rir=real(ifft(H.*repmat(fft(Num,Nfft),[1,num_of_smp,J])));

if flag == 1
    rir = permute(rir,permute_order);
end

% undoing filter's delay
rir = circshift(rir,-floor(order/2),2);

% throwing past RT60 samples 
rir = rir(:,1:n,:);

% max_val = max(rir,2);
% 
% fs
% t60 = edc(rir(3,:),fs)
% evaluate the T60 using EDC
% t60=zeros(J,1); for j=1:J, figure; t60(j)=edc(rir(5,1:2*fs,j),fs,.144,.19); end

% fprintf('%.3f +- %.3fs\n', mean(t60), std(t60));
end

