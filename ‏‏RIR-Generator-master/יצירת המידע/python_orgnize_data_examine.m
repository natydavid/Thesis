function  [ data ]= python_orgnize_data_examine( vargin )
% this function receives rooms' rirs and creates a gmm model for each room
% room_rirs is a rooms X 1 cell array
% room is speakerPointsInRoom X CD cell array
% input:
% vargin.room_rir_list -  room rirs
% vargin.train_speaker_pos - a vector indicating which speaker points should be used for training
% vargin.CD - a vector indicating which CD to take
% vargin.save_name - name to save the model
% vargin.save_path - path to save the model


room_rir_list = vargin.room_rir_list;
room_list = vargin.room_list ;
normalize_data = vargin.normalize_data;
rooms_num = size(room_rir_list,1);


%% zero padding the data

room_rir_list = room_zero_pad( room_rir_list );

%% transfering all data to one matrix to deliver to algorithm

room_matrixs = cell(rooms_num,1);

for rr = 1:rooms_num
    room_matrixs{rr} = room_rir_list{rr}{1};
end

label = [];
feature = [];
target_names = {};
for rr = 1:rooms_num
    
    points_num = size(room_matrixs{rr},1);
    label = [label ; (rr-1)*ones(points_num,1)];
    feature = [feature ; room_matrixs{rr}];
    target_names = [target_names room_list(rr).roomDisc];
    
end

if isfield(vargin,'normalize_data')
        if normalize_data
            vargout = normalizeFeature( feature );
            feature = vargout.feature;
            min_val = vargout.min_val;
            max_val = vargout.max_val;
        end
end

points = length(label);

if isfield(vargin,'binary')
    if vargin.binary
        precntile = 20;
        precThresh = prctile(abs(feature),precntile,2);
%         precThresh = prctile(abs(feature(:)),precntile);
        bin_feature = feature;
        lesOrEq = bsxfun(@le,abs(bin_feature),precThresh);
        bin_feature(lesOrEq)=0; 
        grtr = bsxfun(@gt,abs(bin_feature),precThresh);
        bin_feature(grtr)=1;
        data.feature = bin_feature;
    else
        
        data.feature = feature;
    end
else
    data.feature = feature;
end
data.label = label;
data.points = points;
data.target_names = target_names;



%% saving models 

if isfield(vargin,'save_path')
    cd(vargin.save_path);
end

if isfield(vargin,'save_name')
    save_name = [vargin.save_name '.mat'];
    save(save_name,'data');
end

end

