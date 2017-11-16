function [ RTF ] = room_RTF_est_gevd( signal,noise,Nfft,Jump,ref_mic )
% room_RTF_est_gevd recievs a noisy signal and a segment with only noise
% Inputs:
% sig_source -  a matrix of dim X num_of_samples X channels (optional)
% sig_record -  same as sig_source

M = size(signal,1);

Rvv = calcR_mul(noise,Nfft,Jump);
Rss = calcR_mul(signal,Nfft,Jump);
vv = gevd_freq_mul(Rss,Rvv);
RTF_freq_domain = squeeze(vv(:,1,:,:)./repmat(vv(ref_mic,1,:,:),[M,1,1,1]));

RTF_freq_domain(ref_mic,:,:) = [];
RTF_time_domain = ifft(RTF_freq_domain,[],2);
RTF = fftshift(RTF_time_domain,2);
RTF = permute(RTF,[3 2 1]);
end

