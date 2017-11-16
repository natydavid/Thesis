function [ room_rir_list,room_list,room_points_list ] = pickRooms( room_list,room_rir_list,vol,shp,beta,room_points_list )
room_sz = size(room_list);
if vol == -1
    vol = 1:room_sz(1);
end
if shp == -1
    shp = 1:room_sz(2);
end
if beta == -1
    beta = 1:room_sz(3);
end

room_list = room_list(vol,shp,beta);
room_rir_list = room_rir_list(vol,shp,beta);

if exist('room_points_list','var')
    room_points_list = room_points_list(vol,shp,beta);
end

end

