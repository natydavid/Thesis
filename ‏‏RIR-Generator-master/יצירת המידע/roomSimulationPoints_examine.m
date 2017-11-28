function [ room_sim ] = roomSimulationPoints_examine( room,CDfactor, scnd_mic_dis)
% roomSimulationPoints recieves room properties and return a set of points
% ready to be handled to the rir_generator


roomDim = room.roomDim;
dim = length(roomDim);
beta = room.beta;

sim_point = struct('r',[],...
    's',[]...
    ); % all fields correspond to rir_generator inputs

% critical distance calculation
min_gap = 0;
[RT60,V] = calcRT60(roomDim,beta);
gama = 1;
critDis =  0.057*sqrt((gama*V)/(RT60));
[~,longest_dim] = max(roomDim(1:2));
roomLims = [ min_gap*ones(1,dim); roomDim-min_gap];
meanDim = mean(diff(roomLims,1)); % range to optional distances
maxDim = max(diff(roomLims,1,1)); % max optional distance
disToRec = CDfactor*critDis; %different distance between speaker and reciever
disToRec(disToRec>0.5*maxDim) = 0.5*maxDim;
disToRec = unique(disToRec);
validMsk = true(length(disToRec)+1,1);
spkrPos = roomDim*0.5;

spkrPos_rep = repmat(spkrPos,length(disToRec)+1,1);

num_of_dis = 1;
num_spkrPos = 1;

room_sim.points = repmat(sim_point,num_spkrPos,num_of_dis);
room_sim.validMsk = cell(size(room_sim.points));
room_sim.L = roomDim;
room_sim.beta = beta;

zeros_dis = zeros(length(disToRec)+1,dim);
zeros_dis(:,longest_dim) = zeros_dis(:,longest_dim)+[disToRec' ;1];
x_y_z = spkrPos_rep + zeros_dis;

if exist('scnd_mic_dis','var')
    num_of_examine = size(x_y_z,1);
    xy = scnd_mic_dis*exp(rand(num_of_examine,1)*(2*pi)*1i);
    x = real(xy);
    y = imag(xy);
    tmp = x_y_z;
    tmp(:,1) = tmp(:,1)+x;
    tmp(:,2) = tmp(:,2)+y;
    x_y_z = [x_y_z ; tmp];
%     x_y_z = cat(3,x_y_z,tmp);
%     room_sim.points = repmat(sim_point,num_of_examine,1);
%     room_sim.validMsk = cell(size(room_sim.points));
%     for ee = 1:num_of_examine
%         tmp = x_y_z(ee,:,:);
%         tmp = permute(tmp,[3 2 1]);
%         room_sim.points(ee,1).s = spkrPos;
%         room_sim.validMsk{ee,1} =  true;
%         room_sim.points(ee,1).r = {tmp};
%     end
     room_sim.points(1,1).s = spkrPos;
    room_sim.validMsk{1,1} =  true;
    room_sim.points(1,1).r = {x_y_z};

else
    room_sim.points(1,1).s = spkrPos;
    room_sim.validMsk{1,1} =  true;
    room_sim.points(1,1).r = {x_y_z};
end



end


