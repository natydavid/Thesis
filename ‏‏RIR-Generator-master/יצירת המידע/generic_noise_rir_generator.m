function [ room_rir ] = generic_noise_rir_generator( room_point_list,nworkers )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
c = 340;                    % Sound velocity (m/s)
fs = 16000;                 % Sample frequency (samples/s)
Ts = 1/fs;
L = room_point_list.L;
beta = room_point_list.beta*ones(1,6);
RT60 = calcRT60(L,beta(1));
n = ceil((0.8*RT60)/Ts);

%% noise position creation


roomDim = L;
dim = length(roomDim);
num_noisePos = 1;
min_gap = 0.5; % rir_generator needs a minimum gap n order to work well
roomLims = [ min_gap*ones(1,dim); roomDim-min_gap];

noisePos = (roomLims(2,:)-roomLims(1,:)).*rand(num_noisePos,dim) + roomLims(1,:);
%%
sz = size(room_point_list.points);
num_of_dis = sz(2);
num_spkrPos = sz(1);
room_rir = cell(num_spkrPos,num_of_dis);
for ss = 1:num_spkrPos
    fprintf('speaker#: %d/%d \t', ss,num_spkrPos);
    tim = tic;
    points = room_point_list.points(ss,:);
    valid = room_point_list.validMsk(ss,:);
    s = noisePos;
    
    if exist('nworkers','var')
        parfor (dd = 1:num_of_dis,nworkers)
            validMsk = valid{dd}(:);
            r = points(dd).r(:);
            r = cell2mat(r(validMsk));
            if ~isempty(r)
                room_rir{ss,dd} = rir_generator(c, fs, r, s, L, beta,n);
            end
        end
    else
        for dd = 1:num_of_dis
            validMsk = valid{dd}(:);
            r = points(dd).r(:);
            r = cell2mat(r(validMsk));
            if ~isempty(r)
                room_rir{ss,dd} = rir_generator(c, fs, r, s, L, beta,n);
            end
        end
    end
    
    tim = toc(tim);
    fprintf('[elaps = %.2f s]\n',tim);
end

end

% for dd = 1:num_of_dis
%     validMsk = valid{dd}(:);
%     r = points(dd).r(:);
%     r = cell2mat(r(validMsk));
%     size(r)
% end