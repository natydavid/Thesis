function [ room_rir ] = generic_MCRoomSim_generator( room_point_list,nworkers )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Falgs
no_scattering_flag = 0;
SimDiff = true;

%% parameters
c = 340;                    % Sound velocity (m/s)
fs = 44100;                 % Sample frequency (samples/s)
Ts = 1/fs;
L = room_point_list.L;
beta = room_point_list.beta*ones(1,6);
part_of_RT60 = 0.8;


sz = size(room_point_list.points);
num_of_dis = sz(2);
num_spkrPos = sz(1);
room_rir = cell(num_spkrPos,num_of_dis);
for ss = 1:num_spkrPos
    fprintf('speaker#: %d/%d \t', ss,num_spkrPos);
    tim = tic;
    points = room_point_list.points(ss,:);
    valid = room_point_list.validMsk(ss,:);
    s = points(1).s;
    
    if exist('nworkers','var')
        
        parfor dd = 1:num_of_dis
            validMsk = valid{dd}(:);
            r = points(dd).r(:);
            r = cell2mat(r(validMsk));
            if ~isempty(r)
                
                % room setup
                Room = SetupRoom('Dim',L);
                if no_scattering_flag
                    Room.Scattering = ones(6)*0;
                end
                
                % source setup
                Sources = AddSource([],'Location',s);
                
                % receiver setup
                rec_num = size(r,1);
                Receivers = [];
                for rr = 1:rec_num
                    Receivers  = AddReceiver(Receivers,'Location',r(rr,:));
                end
                
                % option setup
                beta_tmp = max(1-Room.Absorption);
                RT60 = -1;
                Options = MCRoomSimOptions('Duration',RT60, ...
                    'Verbose', false,...
                    'SimDirect',     true,       ...
                    'SimSpec',     true,       ...
                    'Order',       [-1,-1,-1], ...
                    'SimDiff',     SimDiff, ...
                    'Fs', fs ...
                    );
                
                room_rir{ss,dd} = RunMCRoomSim(Sources,Receivers,Room,Options);
                room_rir{ss,dd}=cell2mat(room_rir{ss,dd}')';
                RT60 = calcRT60(L,beta_tmp);
                room_rir{ss,dd} = room_rir{ss,dd}(:,1:floor(part_of_RT60*RT60*fs));
            end
        end
        
    else
        
        for dd = 1:num_of_dis
            validMsk = valid{dd}(:);
            r = points(dd).r(:);
            r = cell2mat(r(validMsk));
            if ~isempty(r)
                
                % room setup
                Room = SetupRoom('Dim',L);
                if no_scattering_flag
                    Room.Scattering = ones(6)*0;
                end
                
                % source setup
                Sources = AddSource([],'Location',s);
                
                % receiver setup
                rec_num = size(r,1);
                Receivers = [];
                for rr = 1:rec_num
                    Receivers  = AddReceiver(Receivers,'Location',r(rr,:));
                end
                
                % option setup
                beta_tmp = max(1-Room.Absorption);
                RT60 = -1;
                Options = MCRoomSimOptions('Duration',RT60, ...
                    'Verbose', false,...
                    'SimDirect',     true,       ...
                    'SimSpec',     true,       ...
                    'Order',       [-1,-1,-1], ...
                    'SimDiff',     SimDiff, ...
                    'Fs', fs ...
                    );
                room_rir{ss,dd} = RunMCRoomSim(Sources,Receivers,Room,Options);
                room_rir{ss,dd}=cell2mat(room_rir{ss,dd}')';
                RT60 = calcRT60(L,beta_tmp);
                room_rir{ss,dd} = room_rir{ss,dd}(:,1:floor(part_of_RT60*RT60*fs));
            end
        end
        
    end
    
    tim = toc(tim);
    fprintf('[elaps = %.2f s]\n',tim);
end

end

