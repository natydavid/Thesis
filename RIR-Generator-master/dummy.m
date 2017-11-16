clearvars
close all

%% Mode
room_vol =  '20[m^3]'; %'20[m^3]' '40[m^3]'
only_white_noise_flag = 1;

if 0
    load('1_to_4_ratio_rectangle_cubical_noise_RIR.mat','room_rir_list','room_list');
    noise_rir  = room_rir_list{2,2};
    clear room_rir_list
    load('1_to_4_ratio_rectangle_cubical.mat','room_rir_list','room_list');
    signal_rir = room_rir_list{2,2};
    clear room_rir_list
end

switch room_vol
    case '20[m^3]'
        load signal_and_noise_RIRs_20[m^3].mat
    case '40[m^3]'
        load signal_and_noise_RIRs_40[m^3].mat
end

if only_white_noise_flag
    spatial_noise_dB_under_signal = 0;
else
    spatial_noise_dB_under_signal = 10;
end
white_noise_dB_under_spatial = 50;

speech_file = 'SA1.wav';
[speech_anechoic,fs] = audioread(speech_file);
speech_anechoic = speech_anechoic';
speech_anechoic_len = length(speech_anechoic);

%% create recieved speech

speech1 = conv2(1,speech_anechoic,signal_rir1);
% speech1 = speech_anechoic;
speech2 = conv2(1,speech_anechoic,signal_rir2);
% speech2 = speech_anechoic;
speech_len = length(speech1);

speech_power = sum(speech1.^2)/speech_len;
speech_anechoic_power = sum(speech_anechoic.^2)/speech_anechoic_len;

%% creating spatial noise

% noise compared to anechoic speech
spatial_noise_anechoic_power = 10^(-spatial_noise_dB_under_signal/10)*speech_anechoic_power;
spatial_noise_anechoic_std = sqrt(spatial_noise_anechoic_power);
spatial_noise_anechoic_before_RIR = normrnd(0,spatial_noise_anechoic_std,[1 speech_anechoic_len]);

% noise compared to speech
spatial_noise_power = 10^(-spatial_noise_dB_under_signal/10)*speech_power;
spatial_noise_std = sqrt(spatial_noise_power);
spatial_noise_before_RIR = normrnd(0,spatial_noise_std,[1 speech_anechoic_len]);

%recieved noise
spatial_noise1 = conv2(1,spatial_noise_before_RIR,noise_rir1);
spatial_noise2 = conv2(1,spatial_noise_before_RIR,noise_rir2);

spatial_noise_anechoic1 = conv2(1,spatial_noise_anechoic_before_RIR,noise_rir1);
spatial_noise_anechoic2 = conv2(1,spatial_noise_anechoic_before_RIR,noise_rir2);

% white noise
white_noise_power = 10^(-white_noise_dB_under_spatial/10)*spatial_noise_anechoic_power;
white_noise_std = sqrt(white_noise_power);
white_noise = normrnd(0,white_noise_std,[2 length(speech1)]);


% speech and noise
if only_white_noise_flag
    noise1 = white_noise(1,:);
    noise2 = white_noise(2,:);
else
    noise1 = spatial_noise_anechoic1 + white_noise(1,:);
    noise2 = spatial_noise_anechoic2 + white_noise(2,:);
end

speech_and_noise1 = speech1 + noise1;
speech_and_noise2 = speech2 + noise2;

% ax1 = subplot(1,3,1);
% plot(speech_and_noise1);
% ax2 = subplot(1,3,2);
% plot(spatial_noise_anechoic1);
% ax3 = subplot(1,3,3);
% plot(speech1);
%
% linkaxes([ax3,ax2,ax1],'y');
%
% figure
% plot(speech_and_noise1);
% hold on
% plot(speech_and_noise2);

%% RTF estimate

Nfft = 2048; % length(noise_rir1)/2;
Jump = Nfft/2;

% noise PSD Est

% noise_end_idx = 2600;
% noise_for_psd_est = [speech_and_noise1(1:noise_end_idx) ; speech_and_noise2(1:noise_end_idx)];

% Estimation of RTF of speech with noise
noise_for_psd_est = [noise1 ; noise2];

if 1
    temp_noise = cat(3,noise_for_psd_est,noise_for_psd_est);
    speech_for_RTF_est = [speech_and_noise1 ; speech_and_noise2];
    temp_speech = cat(3,speech_for_RTF_est,speech_for_RTF_est);
    Rss_1 = calcR_mul(temp_speech,Nfft,Jump);

    Rvv_1 = calcR_mul(temp_noise,Nfft,Jump);
    Rss = calcR(speech_for_RTF_est,Nfft,Jump);
    Rvv = calcR(noise_for_psd_est,Nfft,Jump);
    vv = gevd_freq(Rss,Rvv);
    vv_1 = gevd_freq_mul(Rss_1,Rvv_1);
%     room_RTF_est_gevd(temp_speech,temp_noise,Nfft,Jump,1);
end

Rvv = calcR(noise_for_psd_est,Nfft,Jump);
Rss = calcR([speech_and_noise1 ; speech_and_noise2],Nfft,Jump);


vv = gevd_freq(Rss,Rvv);
RTF = squeeze(vv(:,1,:)./repmat(vv(1,1,:),[2,1,1]));


% Estimation of RTF of clean speech
RT60_len = Nfft/2;
clean_rtf= room_RTF_est( speech1,speech1,fs,RT60_len );

%% plotting

if 1
    figure('WindowStyle','docked');
    plot(signal_rir1);
    title('RIR')
    
    figure('WindowStyle','docked');
    ax1 = subplot(1,3,1);
    plot(speech1);
    title('Clean signal -  only speech')
    ax2 = subplot(1,3,2);
    plot(noise1);
    title('Noise - Directional and Spatial white noise')
    ax3 = subplot(1,3,3);
    plot(speech_and_noise1);
    title('Signal - speech and noise')
    linkaxes([ax3,ax2,ax1],'y');
    
    figure('WindowStyle','docked');
    plot(clean_rtf)
    title('RTF extracted from clean signal using basic RTF extraction algorithm')
    
    figure('WindowStyle','docked');
    plot(fftshift(ifft(RTF(2,:))))
    title('RTF extracted from noisy signal using gevd_freq')
    
end

