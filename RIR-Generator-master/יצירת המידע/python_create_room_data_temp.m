clc
clearvars
close all


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

%% load room expriment
global feat_flag

exp_mode = 'speech_and_noise_RTF_flag'; % 'feat_flag' 'speech_RTF_flag' 'speech_and_noise_RTF_flag' 'examine_flag' 

feat_flag = 0; % use feature or rirs
speech_RTF_flag = 0;
speech_and_noise_RTF_flag = 0;
examine_flag = 0;

eval([exp_mode '=1;']);

no_test_train_split = 1;
simulator_name = 'rir_generator'; % 'rir_generator' 'MCRoomSim'
RTF_peek_samples = 200;

%% Experiment picks
exp_num = [  7 ];
vol_list = {[ 1 2 3 4  ];[ 1 2 3]};
shp_list = {[1 ];[1]};
beta_list = {[1],[1]};

room_rir_list_total = [];
room_list_total = [];

number_of_experiments = length(exp_num);


for ee = 1:number_of_experiments
    
    %% Experiment cases
    
    exp_vargin.mode = 'python_create';
    exp_vargin.exp_num = exp_num(ee);
    exp_vargin.feat_flag = feat_flag;
    exp_vargin.speech_RTF_flag = speech_RTF_flag;
    exp_vargin.speech_and_noise_RTF_flag = speech_and_noise_RTF_flag;
    exp_vargin.examine_flag = examine_flag;
    exp_vargin.simulator_name = simulator_name;
    vargout = experimentParameters(exp_vargin);
    
    load_exp = vargout.load_exp;
    exp_name = vargout.exp_name;
    
    load(load_exp);
    vargin.binary = 0;
    if feat_flag
        vargin.room_rir_list = room_feat_list(:);
        vargin.binary = 0;
    elseif speech_RTF_flag
        vargin.room_rir_list = room_speech_RTF_list;
        clear room_speech_RTF_list;
    else
        vargin.room_rir_list = room_rir_list;
    end
    
    vol = vol_list{ee};
    shp = shp_list{ee};
    beta = beta_list{ee};
    
    [vargin.room_rir_list,room_list] = pickRooms( room_list,vargin.room_rir_list,vol,shp,beta );
    
    room_rir_list_total = [room_rir_list_total ; vargin.room_rir_list(:)];
    room_list_total = [room_list_total ; room_list(:)];
end

vargin.room_rir_list = room_rir_list_total;
vargin.room_list = room_list_total;
vargin.normalize_data = 0;
vargin.num_of_chunks = 1;
vargin.relative_part = 1; % taking a relative part of the samples due to memory limitations
vargin.same_amount_flag = 1; %making sure all rooms will have the same amount of samples
vargin.train_speaker_pos = 1:15;
vargin.CD  = [ 1 2 3 4];
if speech_RTF_flag
    vargin.RTF_peek_samples = RTF_peek_samples;
end

if vargin.binary
    vargin.save_name = [exp_name '_pythonData_binary'];
elseif  examine_flag
    vargin.save_name = [exp_name '_' simulator_name '_pythonData_examine'];
elseif no_test_train_split
    vargin.save_name = [exp_name '_pythonData_no_split'];
else
    vargin.save_name = [exp_name '_pythonData'];
end

if rooms_data
    python_saved_data_path = which('python_room_saved_data.m');
    [python_saved_data_path,~,~] = fileparts(python_saved_data_path);
else
    python_saved_data_path = which('python_saved_data_path.m');
    [python_saved_data_path,~,~] = fileparts(python_saved_data_path);
end
vargin.save_path = python_saved_data_path;

if examine_flag
    data = python_orgnize_data_examine( vargin );
elseif no_test_train_split
    data = python_orgnize_data_no_split( vargin );
else
    [train_data, test_data] = python_orgnize_data( vargin );
end