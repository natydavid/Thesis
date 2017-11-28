clc
clearvars
close all
dbstop if error
par_calc = 1;

seed = 10;
rng(seed);

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
exp_num = 4;
simulator_name = 'rir_generator'; % 'rir_generator' 'MCRoomSim'


exp_name = exp_list{exp_num};

if strcmp(exp_name,'vol_basic_cubic')
    load_exp = 'vol_basic_cubic_data.mat';
    save_name = 'vol_basic_cubic_data_noise_RIR.mat';
else
    load_exp = [exp_name '.mat'];
    save_name = [exp_name '_noise_RIR2.mat'];
end

%% load room expriment
load(load_exp,'room_points_list','room_list');

vol = [1 2];
shp = [ 2] ;
beta = 1;
room_points_list = room_points_list(vol,shp,beta);
room_list = room_list(vol,shp,beta);
room_sz = size(room_list);
room_noise_rir_list = cell(room_sz);

vol_cls = length(vol);
shp_cls = length(shp);
beta_num = length(beta);

%tracking parameters
totalIter = vol_cls*shp_cls*beta_num;
n=0;

for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            n=n+1;
            fprintf('Volume#: %d/%d\t Shape#: %d/%d\t Beta#: %d/%d\n',...
                vlvl,vol_cls,...
                spsp,shp_cls,...
                bb,beta_num);
            if par_calc
                switch simulator_name
                    case 'rir_generator'
                        room_noise_rir_list{vlvl,spsp,bb} = generic_noise_rir_generator( room_points_list{vlvl,spsp,bb},nworkers );
                    case 'MCRoomSim'
                        room_noise_rir_list{vlvl,spsp,bb} = generic_MCRoomSim_generator( room_points_list{vlvl,spsp,bb} ,nworkers);
                end
            else
                switch simulator_name
                    case 'rir_generator'
                        room_noise_rir_list{vlvl,spsp,bb} = generic_noise_rir_generator( room_points_list{vlvl,spsp,bb});
                    case 'MCRoomSim'
                        room_noise_rir_list{vlvl,spsp,bb} = generic_MCRoomSim_generator( room_points_list{vlvl,spsp,bb});
                end
            end
            prec = (n/totalIter)*100;
            fprintf('Progress: %d%s\n',prec,'%');
        end
    end
end

save( save_name,'room_list');

room_file = matfile(save_name,'Writable',true);

room_file.room_points_list = room_points_list;
room_file.room_rir_list = room_noise_rir_list;

