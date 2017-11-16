function [ room ] = createRoom( vol,shp,beta )
% createRoom takes in volume of a room and its shape and returns a struct
% of room 

room = struct('roomDim',nan,... % room dimensions
              'beta',nan,... % Reverberation time (s) or reflection coeff
              'roomDisc',[]... % room discription: [x,y,z] - surface area - Vol[m^3] 
              );
dim = length(shp);
          
set_idxs = shp<0;  %set size indexs
set_val  = shp(set_idxs);
set_num = sum(set_idxs);
x_y_z = zeros(1,dim);
srfc_area = 0;
roomDisc = [];

%% calculating dimensions
base_val = nthroot(vol/prod(abs(shp)),3-set_num);
x_y_z(set_idxs) = -set_val;
temp = shp*base_val;
x_y_z(~set_idxs) = temp(~set_idxs);
areaCalc
roomDisc = sprintf('[%.2f,%.2f,%.2f] - A:%.2f - V:%.2f',x_y_z(1),x_y_z(2),x_y_z(3),srfc_area,vol);



room.roomDim = x_y_z;
room.beta = beta;
room.roomDisc = roomDisc;
% room.roomSmpls = 


%% nested functions
    function areaCalc
        for d1 = 1:dim
            for d2 = dim:-1:d1+1
                srfc_area = srfc_area+2*(x_y_z(d1)*x_y_z(2));
            end
        end
    end

end


% base_val = nthroot(vol/prod(shp),3-set_num);
% x_y_z = shp*base_val;
% areaCalc
% roomDisc = sprintf('[%d,%d,%d] - A:%d - V:%d',x_y_z(1),x_y_z(2),x_y_z(3),srfc_area,vol);

