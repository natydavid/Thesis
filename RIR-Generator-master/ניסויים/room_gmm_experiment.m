% clc
clearvars
close all

mySeed = 10;
rng(mySeed);

% flags 

train_phase = 1;
test_phase = 1;
normalize_data = 1;
rir_est_flag = 0;
rir_speech_RTF_flag = 0;
rir_RTF_flag = 0;
rir_feat_flag = 1;
log_flag = 0;
global feat_flag speech_RTF_flag
feat_flag = rir_feat_flag;
speech_RTF_flag = rir_speech_RTF_flag;
squeeze_pad_flag = 1; % 0 - none , 1 - squeeze ,2 - zero pad
simulator_name = 'rir_generator'; % 'rir_generator' 'MCRoomSim'


%% Going to path
rooms_data = 1; % saving data in the new computer and not on drive
if rooms_data
    rooms_data_path = which('rooms_data_path.m');
    [rooms_data_path,~,~] = fileparts(rooms_data_path);
    cd(rooms_data_path);
else
    rir_generator_path = which('rir_generator_path.m');
    [rir_generator_path,~,~] = fileparts(rir_generator_path);
    cd(rir_generator_path);
end

%% Experiment cases


exp_vargin.mode = 'gmm_experiment';
exp_vargin.exp_num = 12;
exp_vargin.feat_flag = rir_feat_flag;
exp_vargin.RTF_flag = rir_RTF_flag;
exp_vargin.speech_RTF_flag = rir_speech_RTF_flag;
vargout = experimentParameters(exp_vargin);

exp_disc = vargout.exp_disc;
load_exp = vargout.load_exp;
save_name = vargout.save_name;


if rir_est_flag
    exp_disc = [exp_disc ' - using rir estimation'];
    
    feat_idx = strfind(load_exp,'_feat');
    load_exp = [load_exp(1:feat_idx-1) '_rir_est' load_exp(feat_idx:end)];
    feat_idx = strfind(save_name,'_feat');
    save_name = [save_name(1:feat_idx-1) '_rir_est' save_name(feat_idx:end)];
end


if strcmp(simulator_name,'MCRoomSim')
    exp_disc = [exp_disc ' - using MCRoomSim simulstor'];
    
    feat_idx = strfind(load_exp,'_feat');
    load_exp = [load_exp(1:feat_idx-1) '_MCRSim' load_exp(feat_idx:end)];
    feat_idx = strfind(save_name,'_feat');
    save_name = [save_name(1:feat_idx-1) '_MCRSim' save_name(feat_idx:end)];
end


%% load room expriment

load(load_exp);
vol = [1 2];
shp = [1 ];
beta = 1;
if rir_feat_flag
    [room_feat_list,room_list] = pickRooms( room_list,room_feat_list,vol,shp,beta );
elseif rir_RTF_flag
    [room_feat_list,room_list] = pickRooms( room_list,room_RTF_list,vol,shp,beta );
elseif rir_speech_RTF_flag
    [room_feat_list,room_list] = pickRooms( room_list,room_speech_RTF_list,vol,shp,beta );
end

full_speaker = 1:20;
vargin.room_rir_list = room_feat_list;
vargin.nmix = 4;
vargin.normalize_data = normalize_data;
vargin.speaker = 1:15;
vargin.CD = 1:4;
vargin.save_name  = save_name;
vargin.squeeze_pad_flag  = squeeze_pad_flag;

%% train phase
if train_phase
    gmm_models = room_gmm( vargin );
end

%% test phase
if test_phase

    load(save_name);
    
    vargin_test.room_rir_list = room_feat_list(:);
    vargin_test.room_list = room_list(:);
    vargin_test.normalize_data = normalize_data;
    vargin_test.gmm_models = gmm_models(:);
    test_speaker_idx = ~ismember(full_speaker,vargin.speaker);
    vargin_test.speaker = full_speaker(test_speaker_idx);
    vargin_test.CD = vargin.CD;
    vargin_test.exp_disc = exp_disc;
    vargin_test.log_flag = log_flag;
    vargin_test.squeeze_pad_flag = squeeze_pad_flag;
    room_gmm_test(vargin_test);
end