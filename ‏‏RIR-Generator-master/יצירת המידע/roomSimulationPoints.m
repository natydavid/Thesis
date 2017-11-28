function [ room_sim ] = roomSimulationPoints( room )
% roomSimulationPoints recieves room properties and return a set of points
% ready to be handled to the rir_generator


roomDim = room.roomDim;
dim = length(roomDim);
beta = room.beta;

sim_point = struct('r',[],...
    's',[]...
    ); % all fields correspond to rir_generator inputs
          
% critical distance calculation
[RT60,V] = calcRT60(roomDim,beta);
gama = 1;
critDis =  0.057*sqrt((gama*V)/(RT60));
CDfactor = [ 0.5 1 2 3];

min_gap = 0.5; % rir_generator needs a minimum gap n order to work well
num_of_dis = 4;
spehereRes = 25;
num_spkrPos = 20;
scnd_mic_dis = 0.05; % 5[cm]

roomLims = [ min_gap*ones(1,dim); roomDim-min_gap];
meanDim = mean(diff(roomLims,1)); % range to optional distances 
maxDim = max(diff(roomLims,1)); % max optional distance 
disToRec = CDfactor*critDis; %different distance between speaker and reciever
disToRec(disToRec>maxDim) = maxDim;
% disToRec = disToRec(disToRec>=min_gap);
disToRec = unique(disToRec);

spkrPos = (roomLims(2,:)-roomLims(1,:)).*rand(num_spkrPos,dim) + roomLims(1,:);

room_sim.points = repmat(sim_point,num_spkrPos,num_of_dis);
room_sim.validMsk = cell(size(room_sim.points));
room_sim.L = roomDim;
room_sim.beta = beta;



for dd = 1:num_of_dis
    [x,y,z] = sphere(spehereRes);
    
    %second mic position randomization
    vecIdx = unidrnd((spehereRes+1)^2,1,1);
    vec = [x(vecIdx) y(vecIdx) z(vecIdx)]*scnd_mic_dis;
    
    [x,y,z] = scaleXYZ(x,y,z,disToRec,dd);
    for ss = 1:num_spkrPos
        [vec(1),vec(2)] = rotateXY(vec(1),vec(2));
        thisPos = spkrPos(ss,:);
        room_sim.points(ss,dd).s=thisPos;
        
        [r_x,r_y,r_z]=shiftXYZ(x,y,z,thisPos);
        
        validMsk = validChk(r_x,r_y,r_z,roomLims);
        room_sim.validMsk{ss,dd} = validMsk;
        
        x_y_z = cat(3,r_x,r_y,r_z);
        scnd_mic = x_y_z ;
        scnd_mic(:,:,1) = scnd_mic(:,:,1)+vec(1);
        scnd_mic(:,:,2) = scnd_mic(:,:,2)+vec(2);
        scnd_mic(:,:,3) = scnd_mic(:,:,3)+vec(3);
        x_y_z = cat(4,x_y_z,scnd_mic);
        x_y_z = permute(x_y_z,[4 3 1 2]);
        x_y_z = num2cell(x_y_z,[1 2]);
        x_y_z = permute(x_y_z,[3 4 1 2]);

        
        
%         x_y_z = permute(x_y_z,[1 3 2]);
%         x_y_z = num2cell(x_y_z,2);
%         x_y_z = permute(x_y_z,[1 3 2]);
        room_sim.points(ss,dd).r = x_y_z;
    end
end

    

end



function [x,y,z] = scaleXYZ(x,y,z,disToRec,dd)
x = x*disToRec(dd);
y = y*disToRec(dd);
z = z*disToRec(dd);
end

function [r_x,r_y,r_z]=shiftXYZ(x,y,z,thisPos)
r_x = x+thisPos(1);
r_y = y+thisPos(2);
r_z = z+thisPos(3);
end

function [validMsk] = validChk(r_x,r_y,r_z,roomLims)
Val_x = r_x >= roomLims(1,1) & r_x <= roomLims(2,1);
Val_y = r_y >= roomLims(1,2) & r_y <= roomLims(2,2);
Val_z = r_z >= roomLims(1,3) & r_z <= roomLims(2,3);
validMsk = Val_x & Val_y & Val_z;
end

function [x,y] = rotateXY(x,y)
xy = x+y*1i;
ang = rand*(2*pi);
xy = xy*exp(ang*1i);
x = real(xy);
y = imag(xy);
end


