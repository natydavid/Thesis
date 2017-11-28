% clc
% clearvars
close all
dbstop if error

if 0
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

if exist('vargin','var')
    exp_num = vargin.exp_num;
else
    exp_num = 18;
end

rir_est_flag = 0;
simulator_name = 'rir_generator'; % 'rir_generator' 'MCRoomSim'


exp_name = exp_list{exp_num};

if rir_est_flag
    exp_name = [exp_name '_rir_est'];
end

if strcmp(simulator_name,'MCRoomSim')
    exp_name = [exp_name '_MCRSim'];
end


if strcmp(exp_name,'vol_basic_cubic')
    load_exp = 'vol_basic_cubic_data.mat';
    save_name = 'vol_basic_cubic_data_feat.mat';
else
    load_exp = [exp_name '.mat'];
    save_name = [exp_name '_feat.mat'];
end

% switch exp_name
%     case 'vol_basic_cubic' % 2 Vol , cubic shape , same beta
%         load_exp = 'vol_basic_cubic_data.mat';
%         save_name = 'vol_basic_cubic_data_feat.mat';
%         
%     case 'small_Vs_concertHall_cubic' % small office Vs concert hall, cubic shape, same beta
%         load_exp = 'small_Vs_concertHall_cubic.mat';
%         save_name = 'small_Vs_concertHall_cubic_feat.mat';
%         
%     case '1_to_4_ratio_rectangle' % small office Vs concert hall, cubic shape, same beta
%         load_exp = '1_to_4_ratio_rectangle.mat';
%         save_name = '1_to_4_ratio_rectangle_feat.mat';
%         
%     otherwise
%         load_exp = 'rooms_data.mat';
%         save_name = 'rooms_data_feat.mat';
%         
% end


%% load room expriment
load(load_exp);
% room_feat_list = roomQuantiz( room_rir_list );

if rir_est_flag
    room_feat_list = rir2feat( room_rir_est_list,save_name,room_list );
else
    room_feat_list = rir2feat( room_rir_list,save_name,room_list );
end

