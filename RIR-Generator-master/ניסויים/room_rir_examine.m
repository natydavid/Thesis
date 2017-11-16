
% clc
clearvars
% close all
dbstop if error
save_flag = 0;
saveFig = 0 ;

rir_est_flag = 0;
rir_RTF_flag = 0;
rir_feat_flag = 0;
rir_speech_RTF_flag = 0;
rir_speech_and_noise_RTF_flag = 1;

speech_file = 'SA1.wav';
[speech_anechoic,fs] = audioread(speech_file);
speech_anechoic = speech_anechoic';

if rir_speech_and_noise_RTF_flag
    
    spatial_noise_dB_under_signal = 10;
    white_noise_dB_under_spatial = 10;
    ref_mic = 1;

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
end

%% Going to path
rooms_data = 0; % saving data in the new computer and not on drive
if rooms_data
    rooms_data_path = which('rooms_data_path.m');
    [rooms_data_path,~,~] = fileparts(rooms_data_path);
    cd(rooms_data_path);
else
    rir_generator_path = which('rir_generator_path.m');
    [rir_generator_path,~,~] = fileparts(rir_generator_path);
    cd(rir_generator_path);
end
%% parameters
dim = 3;
fs = 16e3;
CDfactor = [ 0.1 1 2 3];
exp_num = [  4 ];
scnd_mic_dis = 0.05;

if length(exp_num)>1
    double_exp_flag = 1;
else
    double_exp_flag = 0;
end

% experiment 1
beta_pick = 1;
vol_pick = [1 4];
shp_pick = 1;

vargin.mode = 'create_room'
vargin.exp_num = exp_num(1);
vargout = experimentParameters(vargin);

if double_exp_flag
    beta = vargout.beta(beta_pick);
    room_shp = vargout.room_shp(shp_pick,:);
    room_vol = vargout.room_vol(vol_pick);
else
    beta = vargout.beta;
    room_shp = vargout.room_shp;
    room_vol = vargout.room_vol;
end

if double_exp_flag
    
    % experiment 2
    beta_pick = 1;
    vol_pick = 1;
    shp_pick = 1;
    
    
    vargin.mode = 'create_room'
    vargin.exp_num = exp_num(2);
    vargout = experimentParameters(vargin);
    
    beta_tmp = vargout.beta(beta_pick);
    room_shp_tmp = vargout.room_shp(shp_pick,:);
    room_vol_tmp = vargout.room_vol(vol_pick);
    
    % unite experiment
    beta = unique([beta beta_tmp]);
    room_shp = unique([room_shp; room_shp_tmp],'rows');
    room_vol = unique([room_vol room_vol_tmp]);
    
end

%% rooms creation
room_list = rooms_creation( room_vol,room_shp,beta );

if rooms_data
    rooms_data_path = which('rooms_data_path.m');
    [rooms_data_path,~,~] = fileparts(rooms_data_path);
    cd(rooms_data_path);
else
    rir_generator_saved_data = which('rir_generator_saved_data.m');
    [rir_generator_saved_data,~,~] = fileparts(rir_generator_saved_data);
    cd(rir_generator_saved_data);
end




room_sz = size(room_list);
room_points_list = cell(room_sz);
room_rir_list = cell(room_sz);
if rir_speech_and_noise_RTF_flag
    room_noise_rir_list = room_rir_list;
end

vol_cls = length(room_vol);
shp_cls = size(room_shp,1);
beta_num = length(beta);

vol_picked = 1;
lens = [];



for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            if rir_RTF_flag || rir_speech_RTF_flag || rir_speech_and_noise_RTF_flag
                room_points_list{vlvl,spsp,bb} = roomSimulationPoints_examine( room_list(vlvl,spsp,bb),CDfactor,scnd_mic_dis );
            else
                room_points_list{vlvl,spsp,bb} = roomSimulationPoints_examine( room_list(vlvl,spsp,bb),CDfactor );
            end
        end
    end
end

%tracking parameters
totalIter = vol_cls*shp_cls*beta_num;
n=0;
max_samples = 0;

for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            n=n+1;
            fprintf('rir generator:\n');
            fprintf('Volume#: %d/%d\t Shape#: %d/%d\t Beta#: %d/%d\n',...
                vlvl,vol_cls,...
                spsp,shp_cls,...
                bb,beta_num);
            room_rir_list{vlvl,spsp,bb} = generic_rir_generator( room_points_list{vlvl,spsp,bb} );
            
            if rir_speech_and_noise_RTF_flag
                room_noise_rir_list{vlvl,spsp,bb} = generic_noise_rir_generator( room_points_list{vlvl,spsp,bb} );
            end
            
            num_of_examine = size(room_rir_list{vlvl,spsp,bb}{1},1);
            if rir_RTF_flag
                mic1 = room_rir_list{vlvl,spsp,bb}{1}(1:num_of_examine/2,:);
                mic2 = room_rir_list{vlvl,spsp,bb}{1}(num_of_examine/2+1:end,:);
                room_rir_list{vlvl,spsp,bb}{1} = room_RTF_est( mic1,mic2,fs );
            elseif rir_speech_RTF_flag || rir_speech_and_noise_RTF_flag
                max_samples = max(size( room_rir_list{vlvl,spsp,bb}{1,1},2),max_samples);
            end
            prec = (n/totalIter)*100;
            fprintf('Progress: %d%s\n',prec,'%');
        end
    end
end


%%

if rir_speech_RTF_flag || rir_speech_and_noise_RTF_flag
    n=0;
    Nfft = max_samples/2;
    RT60_len = max_samples;
    for vlvl = 1:vol_cls
        for spsp = 1:shp_cls
            for bb = 1:beta_num
                n=n+1;
                fprintf('RTF generator:\n');
                fprintf('Volume#: %d/%d\t Shape#: %d/%d\t Beta#: %d/%d\n',...
                    vlvl,vol_cls,...
                    spsp,shp_cls,...
                    bb,beta_num);
                num_of_examine = size(room_rir_list{vlvl,spsp,bb}{1},1);
                mic1 = room_rir_list{vlvl,spsp,bb}{1}(1:num_of_examine/2,:);
                mic2 = room_rir_list{vlvl,spsp,bb}{1}(num_of_examine/2+1:end,:);
                temp_speech1 = conv2(1,speech_anechoic,mic1);
                temp_speech2 = conv2(1,speech_anechoic,mic2);
                
                if rir_speech_and_noise_RTF_flag
                    mic1_spatial_noise = room_noise_rir_list{vlvl,spsp,bb}{1}(1:2:end,:);
                    mic2_spatial_noise = room_noise_rir_list{vlvl,spsp,bb}{1}(2:2:end,:);
                    temp_spatial_noise1 = conv2(1,spatial_noise_anechoic_before_RIR,mic1_spatial_noise);
                    temp_spatial_noise2 = conv2(1,spatial_noise_anechoic_before_RIR,mic2_spatial_noise);
                    white_noise = normrnd(0,white_noise_std,[size(temp_speech1) 2]);
                    noise1 = temp_spatial_noise1 + white_noise(:,:,1);
                    noise2 = temp_spatial_noise2 + white_noise(:,:,2);
                    
                    noise = cat(3,noise1,noise2);
                    noise = permute(noise,[3 2 1]);
                    
                    temp_speech1 = temp_speech1 + noise1;
                    temp_speech2 = temp_speech2 + noise2;
                    
                    signal = cat(3,temp_speech1,temp_speech2);
                    signal = permute(signal,[3 2 1]);
                    room_rir_list{vlvl,spsp,bb}{1} = room_RTF_est_gevd( signal,noise,Nfft,Nfft/2,ref_mic );
                else
                    room_rir_list{vlvl,spsp,bb}{1} = room_RTF_est( temp_speech1,temp_speech2,fs,RT60_len );
                end
                
                prec = (n/totalIter)*100;
                fprintf('Progress: %d%s\n',prec,'%');
            end
        end
    end
end

if rir_feat_flag
    room_rir_list = rir2feat( room_rir_list );
end

if rir_est_flag
    room_rir_list = rir2speech2estRir( room_rir_list,room_list,speech_file );
end

[dis_to_plot,~] = size(room_rir_list{1}{1});
global Ts
Ts = 1/16000;

for ii = 1:dis_to_plot
    
    if ii<=length(CDfactor)
        fig_name = [num2str(CDfactor(ii)) '-critical distance ' ];
    else
        fig_name = 'Distance of 1 meter ';
    end
    n=0;
    %     figure
    figure('Name',fig_name,'WindowStyle','docked');
    for vlvl = 1:vol_cls
        for spsp = 1:shp_cls
            for bb = 1:beta_num
                roomDisc = room_list(vlvl,spsp,bb).roomDisc;
                n=n+1;
                subplot(totalIter,1,n);
                this_rir = room_rir_list{vlvl,spsp,bb}{1}(ii,:);
                t = (1:length(this_rir))*Ts;
                if rir_feat_flag
                    plot(this_rir);
                else
                    plot(t,this_rir);
                end
                %                 xlim([0 0.3]);
                %                 ylim([-0.02 0.12]);
                title([fig_name 'at room ' roomDisc ])
                xlabel('time[s]');
                ylabel('Amp');
                
            end
        end
    end
    linkaxes
    
    if saveFig
        savefig([fig_name '.fig']);
    end
end

if rir_feat_flag
    for ii = 1:dis_to_plot
        
        if ii<=length(CDfactor)
            fig_name = [num2str(CDfactor(ii)) '-critical distance ' ];
        else
            fig_name = 'Distance of 1 meter ';
        end
        this_rir = [];
        figure('Name',fig_name,'WindowStyle','docked');
        for vlvl = 1:vol_cls
            for spsp = 1:shp_cls
                for bb = 1:beta_num
                    this_rir = [this_rir; room_rir_list{vlvl,spsp,bb}{1}(ii,:)];
                end
            end
        end
        plot(this_rir(:,40:end)');
        title(fig_name)
        xlabel('feature');
        ylabel('Amp');
        if saveFig
            savefig([fig_name '.fig']);
        end
    end
end

if save_flag
    if double_exp_flag
        save_name =  ['exp_' num2str(exp_num(1)) '_' num2str(exp_num(2)) '_rir_generator_examine.mat'];
    else
        vargout = experimentParameters();
        exp_name = vargout.exp_list{exp_num};
        save_name =  [exp_name '_rir_generator_examine.mat'];
    end
    room_file = matfile(save_name,'Writable',true);
    room_file.room_list = room_list;
    room_file.room_rir_list = room_rir_list;
    room_file.room_points_list = room_points_list;
end


