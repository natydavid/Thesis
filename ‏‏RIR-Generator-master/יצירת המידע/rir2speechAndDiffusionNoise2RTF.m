function [ room_speech_RTF_list ] = rir2speechAndDiffusionNoise2RTF( room_rir_list,room_noise_rir_list,room_list,speech_file,save_name,nworkers )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Initalizations
rir_sz = size(room_rir_list);
room_speech_RTF_list = cell(rir_sz);
spatial_noise_dB_under_signal = 10;
white_noise_dB_under_spatial = 10;
ref_mic = 1;

vol_cls = rir_sz(1);
shp_cls = rir_sz(2);
if length(rir_sz)==3
    beta_num = rir_sz(3);
else
    beta_num = 1;
end

Fs = 16e3;
M = 2;
sample_size = roomSampleSize(room_rir_list);

[speech_anechoic,fs] = audioread(speech_file);
speech_anechoic = speech_anechoic';

%% noise power calculation

% spatial noise
speech_anechoic_len = length(speech_anechoic);
speech_anechoic_power = sum(speech_anechoic.^2)/speech_anechoic_len;
spatial_noise_anechoic_power = 10^(-spatial_noise_dB_under_signal/10)*speech_anechoic_power;
spatial_noise_anechoic_std = sqrt(spatial_noise_anechoic_power);
spatial_noise_anechoic_before_RIR = normrnd(0,spatial_noise_anechoic_std,[1 speech_anechoic_len]);

% white noise
white_noise_power = 10^(-white_noise_dB_under_spatial/10)*spatial_noise_anechoic_power;
white_noise_std = sqrt(white_noise_power);

%%
for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_speech_RTF_list{vlvl,spsp,bb} = cell(size(room_rir_list{vlvl,spsp,bb}));
        end
    end
end

max_min_flag = 1 ; % min - 0 , max - 1
if max_min_flag
    num_samples = 0;
else
    num_samples = inf;
end

for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            if max_min_flag
                num_samples = max(size( room_rir_list{vlvl,spsp,bb}{1,1},2),num_samples);
            else
                num_samples = min(size( room_rir_list{vlvl,spsp,bb}{1,1},2),num_samples);
            end
        end
    end
end
%tracking parameters
totalIter = vol_cls*shp_cls*beta_num;
n=0;
Nfft = num_samples/2;
for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            n=n+1;
            fprintf('Volume#: %d/%d\t Shape#: %d/%d\t Beta#: %d/%d\n',...
                vlvl,vol_cls,...
                spsp,shp_cls,...
                bb,beta_num);
            
            [num_spkrPos,num_of_dis] = size( room_speech_RTF_list{vlvl,spsp,bb});
            rir_len = size( room_rir_list{vlvl,spsp,bb}{1,1},2);
            
            for ss = 1:num_spkrPos
                fprintf('speaker#: %d/%d \t', ss,num_spkrPos);
                tim = tic;
                
                %             Diffusion noise creation
                
                convoled_signal_len = rir_len + speech_anechoic_len - 1;
                nSamples = [sample_size{vlvl,spsp,bb,ss,:}];
                nSamples = nSamples(1:2:end)/2;
                diffusedNoise = roomDiffusionNoise( M, convoled_signal_len, nSamples, num_of_dis, white_noise_power, Fs );
                aa=1;
                if exist('nworkers','var')
                    thisRirs  =room_rir_list{vlvl,spsp,bb}(ss,:);
                    thisNoiseRirs  =room_noise_rir_list{vlvl,spsp,bb}(ss,:);
                    
                    temp_speech_RTF = cell(1,num_of_dis);
                    parfor (dd = 1:num_of_dis,nworkers)
                        if ~isempty(thisRirs{dd})
                            mic1 = thisRirs{dd}(1:2:end,:);
                            mic2 = thisRirs{dd}(2:2:end,:);
                            mic1_spatial_noise = thisNoiseRirs{dd}(1:2:end,:);
                            mic2_spatial_noise = thisNoiseRirs{dd}(2:2:end,:);
                            
                            temp_speech1 = conv2(1,speech_anechoic,mic1);
                            temp_speech2 = conv2(1,speech_anechoic,mic2);
                            temp_spatial_noise1 = conv2(1,spatial_noise_anechoic_before_RIR,mic1_spatial_noise);
                            temp_spatial_noise2 = conv2(1,spatial_noise_anechoic_before_RIR,mic2_spatial_noise);
                            white_noise = normrnd(0,white_noise_std,[size(temp_speech1) 2]);
                            noise1 = temp_spatial_noise1 + white_noise(:,:,1) + diffusedNoise{dd}(:,:,1);
                            noise2 = temp_spatial_noise2 + white_noise(:,:,2) + diffusedNoise{dd}(:,:,2);
                            
                            noise = cat(3,noise1,noise2);
                            noise = permute(noise,[3 2 1]);
                            
                            temp_speech1 = temp_speech1 + noise1;
                            temp_speech2 = temp_speech2 + noise2;
                            
                            signal = cat(3,temp_speech1,temp_speech2);
                            signal = permute(signal,[3 2 1]);
                            
                            temp_speech_RTF{dd} = room_RTF_est_gevd( signal,noise,Nfft,Nfft/2,ref_mic );
                        end
                    end
                    room_speech_RTF_list{vlvl,spsp,bb}(ss,:) = temp_speech_RTF;
                else
                    for dd = 1:num_of_dis
                        if ~isempty(room_rir_list{vlvl,spsp,bb}{ss,dd})
                                                       
                            mic1 = room_rir_list{vlvl,spsp,bb}{ss,dd}(1:2:end,:);
                            mic2 = room_rir_list{vlvl,spsp,bb}{ss,dd}(2:2:end,:);
                            mic1_spatial_noise = room_noise_rir_list{vlvl,spsp,bb}{ss,dd}(1:2:end,:);
                            mic2_spatial_noise = room_noise_rir_list{vlvl,spsp,bb}{ss,dd}(2:2:end,:);
                            
                            temp_speech1 = conv2(1,speech_anechoic,mic1);
                            temp_speech2 = conv2(1,speech_anechoic,mic2);
                            temp_spatial_noise1 = conv2(1,spatial_noise_anechoic_before_RIR,mic1_spatial_noise);
                            temp_spatial_noise2 = conv2(1,spatial_noise_anechoic_before_RIR,mic2_spatial_noise);
                            white_noise = normrnd(0,white_noise_std,[size(temp_speech1) 2]);
                            noise1 = temp_spatial_noise1 + white_noise(:,:,1) + diffusedNoise{dd}(:,:,1);
                            noise2 = temp_spatial_noise2 + white_noise(:,:,2) + diffusedNoise{dd}(:,:,2);
                            
                            noise = cat(3,noise1,noise2);
                            noise = permute(noise,[3 2 1]);
                            
                            temp_speech1 = temp_speech1 + noise1;
                            temp_speech2 = temp_speech2 + noise2;
                            
                            signal = cat(3,temp_speech1,temp_speech2);
                            signal = permute(signal,[3 2 1]);
                            
                            room_speech_RTF_list{vlvl,spsp,bb}{ss,dd} = room_RTF_est_gevd( signal,noise,Nfft,Nfft/2,ref_mic );
                        end
                    end
                end
                
                tim = toc(tim);
                fprintf('[elaps = %.2f s]\n',tim);
            end
            
            prec = (n/totalIter)*100;
            fprintf('Progress: %d%s\n',prec,'%');
        end
    end
end



if exist('save_name','var')
    save(save_name,'room_speech_RTF_list');
    room_file = matfile(save_name,'Writable',true);
    room_file.room_list = room_list;
    room_file.fs = fs;
end

end

