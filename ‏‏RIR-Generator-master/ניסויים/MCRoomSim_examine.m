
clc
clearvars
close all
dbstop if error

%%  flags
rir_est_flag = 0;
parllel_flag = 0;
save_flag = 1;

speech_file = 'SA1.wav';

if parllel_flag
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
    rir_generator_path = which('rir_generator_path.m');
    [rir_generator_path,~,~] = fileparts(rir_generator_path);
    cd(rir_generator_path);
end
%% parameters
saveFig = 0 ;
dim = 3;
CDfactor = [ 0.1 1 2 3];
exp_num = [ 4  ]
if length(exp_num)>1
    double_exp_flag = 1;
else
    double_exp_flag = 0;
end

% experiment 1
beta_pick = 1;
vol_pick = 1;
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


vol_cls = length(room_vol);
shp_cls = size(room_shp,1);
beta_num = length(beta);

vol_picked = 1;
lens = [];



for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_points_list{vlvl,spsp,bb} = roomSimulationPoints_examine( room_list(vlvl,spsp,bb),CDfactor );
        end
    end
end

%% RIR creation

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
            if parllel_flag
                room_rir_list{vlvl,spsp,bb} = generic_MCRoomSim_generator(room_points_list{vlvl,spsp,bb}, nworkers );
            else
                room_rir_list{vlvl,spsp,bb} = generic_MCRoomSim_generator(room_points_list{vlvl,spsp,bb} );
            end
      
            prec = (n/totalIter)*100;
            fprintf('Progress: %d%s\n',prec,'%');
        end
    end
end

if rir_est_flag
    room_rir_list = rir2speech2estRir( room_rir_list,room_list,speech_file );
end

%% Plotting

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
%                 this_rir = room_rir_list{vlvl,spsp,bb}{ii}';
                this_rir = room_rir_list{vlvl,spsp,bb}{1}(ii,:);
                t = (1:length(this_rir))*Ts;
                plot(t,this_rir);
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

if save_flag    
    if double_exp_flag
        save_name =  ['exp_' num2str(exp_num(1)) '_' num2str(exp_num(2)) '_MCRoomSim_examine.mat'];
    else
        vargout = experimentParameters();
        exp_name = vargout.exp_list{exp_num};
        save_name =  [exp_name '_MCRoomSim_examine.mat'];
    end
    room_file = matfile(save_name,'Writable',true);
    room_file.room_list = room_list;
    room_file.room_rir_list = room_rir_list;
    room_file.room_points_list = room_points_list;   
end



