function [ room_list ] = rooms_creation( room_vol,room_shp,beta )

%% rooms creation

room = struct('roomDim',nan,... % room dimensions
              'beta',nan,... % Reverberation time (s) or reflection coeff
              'roomDisc',[]... % room discription: [x,y,z] - surface area - Vol[m^3] 
              );

%% Going to path
rir_generator_path = which('rir_generator_path.m');
[rir_generator_path,~,~] = fileparts(rir_generator_path);
cd(rir_generator_path);

%% parameters

dim = 3;
vol_cls = length(room_vol);
shp_cls = size(room_shp,1);
beta_num = length(beta);
room_cls = shp_cls*vol_cls;
room_list = repmat(room,vol_cls,shp_cls,beta_num);


for vlvl = 1:vol_cls
    for spsp = 1:shp_cls
        for bb = 1:beta_num
            room_list(vlvl,spsp,bb) = createRoom(room_vol(vlvl),room_shp(spsp,:),beta(bb));
            
        end
    end
end
end


