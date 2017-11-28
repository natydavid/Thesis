clc 
clearvars
close all
dbstop if error
speech_file = 'SA1.wav';
par_calc = 1;

if par_calc
    nworkers = 12;
    nworkers = min(nworkers, feature('NumCores'));
    
    p = gcp('nocreate')
    if isempty(p)
        parpool(nworkers)
    end
end

%% Going to path
rooms_data = 1; % saving data in the new computer and not on drive
if rooms_data
    rooms_data_path = which('rooms_data_path.m');
    [rooms_data_path,~,~] = fileparts(rooms_data_path);
    cd(rooms_data_path);
else
    rir_generator_path = which('rir_generator_saved_data.m');
    [rir_generator_path,~,~] = fileparts(rir_generator_path);
    cd(rir_generator_path);
end

%% Experiment cases
vargout = experimentParameters();
exp_list = vargout.exp_list;
exp_num = 3;

exp_name = exp_list{exp_num};

if strcmp(exp_name,'vol_basic_cubic')
    load_exp = 'vol_basic_cubic_data.mat';
    save_name = 'vol_basic_cubic_data_speech.mat';
else
    load_exp = [exp_name '.mat'];
    save_name = [exp_name '_speech.mat'];
end

%% load room expriment
load(load_exp);

vol = -1;
shp = 1 ;
beta = 1;
[room_rir_list,room_list] = pickRooms( room_list,room_rir_list,vol,shp,beta );

if par_calc
    room_speech_list = rir2speech( room_rir_list,room_list,speech_file,save_name,nworkers );
else
    room_speech_list = rir2speech( room_rir_list,room_list,speech_file,save_name );
end

