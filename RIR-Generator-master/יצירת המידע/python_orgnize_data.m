function  [train_data, test_data]= python_orgnize_data( vargin )
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
train_speaker_pos = vargin.train_speaker_pos;
CD = vargin.CD ;
room_list = vargin.room_list ;
normalize_data = vargin.normalize_data;
rooms_num = size(room_rir_list,1);
room_sz = size(room_rir_list{1});
speaker_num = room_sz(1);
CD_num = room_sz(2);

all_speaker = 1:speaker_num;
test_speaker_pos = find(~ismember(all_speaker,train_speaker_pos));

%% zero padding the data

room_rir_list = room_squeeze( room_rir_list );

%% transfering all data to one matrix to deliver to algorithm

% train data
room_matrixs = cell(rooms_num,1);

for rr = 1:rooms_num
    for ss = train_speaker_pos
        for cc = CD
            if isempty(room_matrixs{rr})
                room_matrixs{rr} = room_rir_list{rr}{ss,cc};
            else
                room_matrixs{rr} = [room_matrixs{rr}; room_rir_list{rr}{ss,cc}];
            end
        end
    end
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
        if vargin.normalize_data
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
        train_data.feature = bin_feature;
    else
        
        train_data.feature = feature;
    end
else
    train_data.feature = feature;
end
train_data.label = label;
train_data.points = points;
train_data.target_names = target_names;


% test data

room_matrixs = cell(rooms_num,1);

for rr = 1:rooms_num
    for ss = test_speaker_pos
        for cc = CD
            if isempty(room_matrixs{rr})
                room_matrixs{rr} = room_rir_list{rr}{ss,cc};
            else
                room_matrixs{rr} = [room_matrixs{rr}; room_rir_list{rr}{ss,cc}];
            end
        end
    end
end

label = [];
feature = [];
for rr = 1:rooms_num
    
    points_num = size(room_matrixs{rr},1);
    label = [label ; (rr-1)*ones(points_num,1)];
    feature = [feature ; room_matrixs{rr}];
    
end

if isfield(vargin,'normalize_data')
        if vargin.normalize_data
            feature = normalizeFeature( feature,min_val,max_val );
        end
end

points = length(label);
if isfield(vargin,'binary')
    if vargin.binary
        bin_feature = feature;
        precThresh = prctile(abs(feature),precntile,2);
        lesOrEq = bsxfun(@le,abs(bin_feature),precThresh);
        bin_feature(lesOrEq)=0; 
        grtr = bsxfun(@gt,abs(bin_feature),precThresh);
        bin_feature(grtr)=1;
        test_data.feature = bin_feature;
        
    else
        test_data.feature = feature;
    end
else
    test_data.feature = feature;
end
test_data.label = label;
test_data.points = points;

%% saving models 

if isfield(vargin,'save_path')
    cd(vargin.save_path);
end

if isfield(vargin,'save_name')
    save_name = [vargin.save_name '.mat'];
    save(save_name,'train_data','test_data');
end

end

