
clc
clearvars
close all

par_calc = 0;

if par_calc
    nworkers = 12;
    nworkers = min(nworkers, feature('NumCores'));
    
    p = gcp('nocreate')
    if isempty(p)
        parpool(nworkers)
    end
end
%% Going to path
rooms_data = 1; % saving data in the new computer and not on google drive
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
simulator_name = 'rir_generator'; % 'rir_generator' 'MCRoomSim'

vargin.mode = 'create_room'
vargin.exp_num = 4;
vargout = experimentParameters(vargin);

beta = vargout.beta;
room_shp = vargout.room_shp;
room_vol = vargout.room_vol;
save_file = vargout.save_file;


if strcmp(simulator_name,'MCRoomSim')
    save_file = [save_file(1:end-4) '_MCRSim' save_file(end-3:end)]
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


vol_cls = length(room_vol);
shp_cls = size(room_shp,1);
beta_num = length(beta);
lens = [];



for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_points_list{vlvl,spsp,bb} = roomSimulationPoints( room_list(vlvl,spsp,bb) );
        end
    end
end

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
                        room_rir_list{vlvl,spsp,bb} = generic_rir_generator( room_points_list{vlvl,spsp,bb},nworkers );
                    case 'MCRoomSim'
                        room_rir_list{vlvl,spsp,bb} = generic_MCRoomSim_generator( room_points_list{vlvl,spsp,bb} ,nworkers);
                end
            else
                switch simulator_name
                    case 'rir_generator'
                        room_rir_list{vlvl,spsp,bb} = generic_rir_generator( room_points_list{vlvl,spsp,bb} );
                    case 'MCRoomSim'
                        room_rir_list{vlvl,spsp,bb} = generic_MCRoomSim_generator( room_points_list{vlvl,spsp,bb} );
                end
            end
            prec = (n/totalIter)*100;
            fprintf('Progress: %d%s\n',prec,'%');
        end
    end
end

save( save_file,'room_list');

room_file = matfile(save_file,'Writable',true);

room_file.room_points_list = room_points_list;
room_file.room_rir_list = room_rir_list;



if 1 %% extracting features
    feature_create_room_data
end
clearvars