function [ room_sim ] = roomSimulationPoints_examine_MCRoomSim( room,CDfactor )
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
maxDim = max(diff(roomLims,longest_dim)); % max optional distance 
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

room_sim.points(1,1).s = spkrPos;
room_sim.validMsk{1,1} =  true;
room_sim.points(1,1).r = {x_y_z};


    

end


