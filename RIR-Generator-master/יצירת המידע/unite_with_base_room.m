clc
clearvars
close all
dbstop if error


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
vargout = experimentParameters();
exp_list = vargout.exp_list;
base_num = 12;
add_num = 13;
rir_est_flag = 0;

base_name = exp_list{base_num};
add_name = exp_list{add_num};

if rir_est_flag
    base_name = [ base_name '_rir_est'];
    add_name = [ add_name '_rir_est'];
end


%% uniting rirs
if strcmp(base_name,'vol_basic_cubic')
    load_exp = 'vol_basic_cubic_data.mat';
else
    load_exp = [base_name '.mat'];
end
load(load_exp);


vol_flag = 1;
shp_flag = 0;

%% volume dealing

if vol_flag
    
    base_taken = [1];
    shp_taken = 2;
    
    % distilling first rirs
    room_list_a = room_list(base_taken,shp_taken);
    if ~rir_est_flag
        room_points_list_a = room_points_list(base_taken,shp_taken);
    end
    
    if rir_est_flag
        room_rir_est_list_a = room_rir_est_list(base_taken,shp_taken);
        vars = {'room_list','room_points_list','room_rir_est_list'};
    else
        room_rir_list_a = room_rir_list(base_taken,shp_taken);
        vars = {'room_list','room_points_list','room_rir_list'};
    end
    
    clear(vars{:})
    
    if strcmp(add_name,'vol_basic_cubic')
        load_exp = 'vol_basic_cubic_data.mat';
    else
        load_exp = [add_name '.mat'];
    end
    load(load_exp);
    
    room_file = matfile(load_exp,'Writable',true);
    room_file.room_list = [room_list_a; room_list];
    
    if rir_est_flag
        room_file.room_rir_est_list = [room_rir_est_list_a; room_rir_est_list];
    else
        room_file.room_points_list = [room_points_list_a; room_points_list];
        room_file.room_rir_list = [room_rir_list_a; room_rir_list];
    end
    
    %% uniting featurs
    if strcmp(base_name,'vol_basic_cubic')
        load_exp = 'vol_basic_cubic_data_feat.mat';
    else
        load_exp = [base_name '_feat.mat'];
    end
    load(load_exp);
    
    % distilling first rirs
    room_feat_list_a = room_feat_list(base_taken,shp_taken);
    room_list_a = room_list(base_taken,shp_taken);
    
    
    vars = {'room_feat_list','room_list'};
    clear(vars{:})
    
    if strcmp(add_name,'vol_basic_cubic')
        load_exp = 'vol_basic_cubic_data_feat.mat';
    else
        load_exp = [add_name '_feat.mat'];
    end
    load(load_exp);
    
    room_file = matfile(load_exp,'Writable',true);
    room_file.room_list = [room_list_a; room_list];
    room_file.room_feat_list = [room_feat_list_a; room_feat_list];
    
end

%% shape dealing

if shp_flag
      base_taken = [1 ];
        vol_taken = 1;
    % distilling first rirs
    room_list_a = room_list(vol_taken,base_taken);
    room_points_list_a = room_points_list(vol_taken,base_taken);
    room_rir_list_a = room_rir_list(vol_taken,base_taken);
    
    vars = {'room_list','room_points_list','room_rir_list'};
    clear(vars{:})
    
    if strcmp(add_name,'vol_basic_cubic')
        load_exp = 'vol_basic_cubic_data.mat';
    else
        load_exp = [add_name '.mat'];
    end
    load(load_exp);
    
    room_file = matfile(load_exp,'Writable',true);
    room_file.room_list = [room_list_a room_list];
    room_file.room_points_list = [room_points_list_a room_points_list];
    room_file.room_rir_list = [room_rir_list_a room_rir_list];
    
    
    %% uniting featurs
    if strcmp(base_name,'vol_basic_cubic')
        load_exp = 'vol_basic_cubic_data_feat.mat';
    else
        load_exp = [base_name '_feat.mat'];
    end
    load(load_exp);
    
    % distilling first rirs
    room_feat_list_a = room_feat_list(vol_taken,base_taken);
    room_list_a = room_list(vol_taken,base_taken);
    
    
    vars = {'room_feat_list','room_list'};
    clear(vars{:})
    
    if strcmp(add_name,'vol_basic_cubic')
        load_exp = 'vol_basic_cubic_data_feat.mat';
    else
        load_exp = [add_name '_feat.mat'];
    end
    load(load_exp);
    
    room_file = matfile(load_exp,'Writable',true);
    room_file.room_list = [room_list_a room_list];
    room_file.room_feat_list = [room_feat_list_a room_feat_list];
end