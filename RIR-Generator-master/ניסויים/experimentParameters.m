function [ vargout ] = experimentParameters( vargin )

exp_list = {'vol_basic_cubic'... % #1 2 Vol , cubic shape , same beta
    'small_Vs_concertHall_cubic'... % #2 small office Vs concert hall, cubic shape, same beta
    '1_to_4_ratio_rectangle'... % #3
    '1_to_4_ratio_rectangle_cubical'...% #4
    'vol_limit_check'...% #5
    'vol_limit_check_plus'...% #6
    'vol_limit_check_combine'...% #7
    'vol_limit_check__1_9__2'...% #8
    'vol_limit_check__2_5__3'...% #9
    'shp_limit_check__2_5_to_3_5_2nd_dim'...% #10
    'shp_limit_check__8_2nd_dim'...% #11
    '1_to_2_ratio_rectangle_cubical_base_40[m^3]'...% #12
    'vol_limit_check__2_5__3_base_40[m^3]'...% #13
    'vol_limit_check__1_1__1_4_base_40[m^3]'...% #14
    'vol_limit_check__1_5__1_9_base_40[m^3]'...% #15
    '1_to_2_ratio_rectangle_cubical_base_80[m^3]'...% #16
    'vol_limit_check__2_5__3_base_80[m^3]'...% #17
    'bug_test_base_40[m^3]_1_1__1_3'...% #18
    };


if ~exist('vargin','var')
    vargout.exp_list = exp_list;
    return;
end
exp_name = exp_list{vargin.exp_num(1)};

mode = vargin.mode;
vargout.exp_list = exp_list;

switch mode
    case 'gmm_experiment'
        feat_flag = vargin.feat_flag; % use feature or rirs
        RTF_flag = vargin.RTF_flag; % use RTF
        speech_RTF_flag = vargin.speech_RTF_flag; % use RTF of speech
        exp_disc = {'2 Vol with 1 to 2 ratio , cubic shape , same beta'... % #1 2 Vol , cubic shape , same beta
            'small_Vs_concertHall_cubic'... % #2 small office Vs concert hall, cubic shape, same beta
            '1_to_4_ratio_rectangle'... % #3
            '1_to_4_ratio_rectangle_cubical'...% #4
            'Changing Vol with respect to first one, rectangle shape , same beta'...% #5
            'Changing Vol with respect to first one, rectangle shape , same beta'...% #6
            'Changing Vol with respect to first one, rectangle shape , same beta'...% #7
            'Changing Vol with respect to first one 1.9,2, rectangle shape , same beta'...% #8
            'Changing Vol with respect to first one 2.5,3, rectangle shape , same beta'...% #9
            '1 Vol ,checking changes in the 2nd dim with respect to the 1st and 3rd dim 2, 2.5, 3, 3.5 , same beta'...% #10
            '1 Vol ,checking changes in the 2nd dim with respect to the 1st and 3rd dim 8 times , same beta'...% #11
            '1_to_2_ratio_rectangle_cubical_base_40[m^3]-vol_check'...% #12
            'Changing Vol with respect to first one 2.5,3, rectangle shape , same beta_base_40[m^3]'...% #13
            'Changing Vol with respect to first one 1.1-1.4, rectangle shape , same beta_base_40[m^3]'...% #14
            'Changing Vol with respect to first one 1.5-1.9, rectangle shape , same beta_base_40[m^3]'...% #15
            '1_to_2_ratio_rectangle_cubical_base_80[m^3]-vol_check'...% #16
            'Changing Vol with respect to first one 2.5,3, rectangle shape , same beta_base_80[m^3]'...% #17
            'bug_test_base_40[m^3]_1_1__1_3'...% #18
            };
        
        if feat_flag
            vargout.exp_disc = exp_disc{vargin.exp_num};
            vargout.load_exp = [exp_name '_feat.mat'];
            vargout.save_name = [exp_name '_feat_gmm_model.mat'];
        elseif RTF_flag
            vargout.exp_disc = [exp_disc{vargin.exp_num} ' - using RTF'];
            vargout.load_exp = [exp_name '_RTF.mat'];
            vargout.save_name = [exp_name '_RTF_gmm_model.mat'];
        elseif speech_RTF_flag
            vargout.exp_disc = [exp_disc{vargin.exp_num} ' - using speech RTF'];
            vargout.load_exp = [exp_name '_speech_RTF.mat'];
            vargout.save_name = [exp_name '_speech_RTF_gmm_model.mat'];
        else
            vargout.exp_disc = exp_disc{vargin.exp_num};
            vargout.load_exp = [exp_name '.mat'];
            vargout.save_name = [exp_name '_gmm_model.mat'];
        end
    case 'create_room'
        %% defining volume classes
        
        room_vol(1) = 2.6*3.07*2.5; %standard kids living room
        room_vol(2) = room_vol(1)*2; % big room
        room_vol(3) = room_vol(1)*4; % hall
        room_vol(4) = 42; % small office
        room_vol(5) = 18200; % concert hall
        
        
        
        %% defining shapes classes
        % We defined room ration between x,y,z. if one of the numbers is negative
        % it means it is not ration but a set size
        room_shp(1,:) = [1,1,1]; % cubical
        room_shp(2,:)  = [1,5,-2.5]; % long
        room_shp(3,:)  = [1,5,2]; % long
        room_shp(4,:) = [1,2.5,-2.5]; % rectangle,fixed hight
        room_shp(5,:) = [1,2,1]; % rectangle
        
        %% beta values
        
        beta = [0.5 0.8 0.9];
        
        %% Experiment cases
        switch exp_name
            case 'vol_basic_cubic' % 2 Vol , cubic shape , same beta
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(1,:);
                vargout.room_vol = room_vol(1:2);
                vargout.save_file = 'vol_basic_cubic_data.mat';
            case 'small_Vs_concertHall_cubic' % small office Vs concert hall, cubic shape, same beta
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(1,:);
                vargout.room_vol = room_vol(4:5);
                vargout.save_file = 'small_Vs_concertHall_cubic.mat';
            case '1_to_4_ratio_rectangle' % small office Vs concert hall, cubic shape, same beta
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                vargout.room_vol = room_vol([1 3]);
                vargout.save_file = '1_to_4_ratio_rectangle.mat';
            case '1_to_4_ratio_rectangle_cubical' % small office Vs concert hall, cubic shape, same beta
                vargout.beta = 0.85;
                vargout.room_shp = room_shp([1 5],:);
                vargout.room_vol = room_vol([1 3]);
                vargout.save_file = '1_to_4_ratio_rectangle_cubical.mat';
            case 'vol_limit_check' % #5
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(1);
                res_start = 1.1;
                res = 0.1;
                res_chk_num = 4;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = [room_vol room_vol*res];
                vargout.save_file = 'vol_limit_check.mat';
            case 'vol_limit_check_plus' % #6
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(1);
                res_start = 1.5;
                res = 0.2;
                res_chk_num = 2;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = [room_vol room_vol*res];
                vargout.save_file = 'vol_limit_check_plus.mat';
            case 'vol_limit_check__1_9__2' % #8
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(1);
                res_start = 1.9;
                res = 0.1;
                res_chk_num = 2;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = room_vol*res;
                vargout.save_file = 'vol_limit_check__1_9__2.mat';
            case 'vol_limit_check__2_5__3' % #9
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(1);
                res_start = 2.5;
                res = 0.5;
                res_chk_num = 2;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = room_vol*res;
                vargout.save_file = 'vol_limit_check__2_5__3.mat';
            case 'shp_limit_check__2_5_to_3_5_2nd_dim' % #10
                vargout.beta = 0.85;
                base_shp = room_shp(1,:);
                room_vol = room_vol(1);
                res_start = 2.5;
                res = 0.5;
                res_chk_num = 3;
                res = res_start+((1:res_chk_num)-1)*res;
                room_shp_rep = repmat(base_shp,res_chk_num,1);
                room_shp_rep(:,2) = res';
                vargout.room_vol = room_vol;
                vargout.room_shp = room_shp_rep;
                vargout.save_file = [exp_name '.mat'];
            case 'shp_limit_check__8_2nd_dim' % #11
                vargout.beta = 0.85;
                base_shp = room_shp(1,:);
                room_vol = room_vol(1);
                res_start = 8;
                res = 0.5;
                res_chk_num = 1;
                res = res_start+((1:res_chk_num)-1)*res;
                room_shp_rep = repmat(base_shp,res_chk_num,1);
                room_shp_rep(:,2) = res';
                vargout.room_vol = room_vol;
                vargout.room_shp = room_shp_rep;
                vargout.save_file = [exp_name '.mat'];
            case '1_to_2_ratio_rectangle_cubical_base_40[m^3]' % #12
                vargout.beta = 0.85;
                vargout.room_shp = room_shp([1 5],:);
                vargout.room_vol = room_vol([2 3]);
                vargout.save_file = [exp_name '.mat'];
            case 'vol_limit_check__2_5__3_base_40[m^3]' % #13
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(2);
                res_start = 2.5;
                res = 0.5;
                res_chk_num = 2;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = room_vol*res;
                vargout.save_file = [exp_name '.mat'];
            case 'vol_limit_check__1_1__1_4_base_40[m^3]' % #14
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(2);
                res_start = 1.1;
                res = 0.1;
                res_chk_num = 4;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = room_vol*res;
                vargout.save_file = [exp_name '.mat'];
            case 'vol_limit_check__1_5__1_9_base_40[m^3]' % #15
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(2);
                res_start = 1.5;
                res = 0.1;
                res_chk_num = 4;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = room_vol*res;
                vargout.save_file = [exp_name '.mat'];
            case '1_to_2_ratio_rectangle_cubical_base_80[m^3]' % #16
                vargout.beta = 0.85;
                vargout.room_shp = room_shp([1 5],:);
                vargout.room_vol = [room_vol(3) room_vol(3)*2] ;
                vargout.save_file = [exp_name '.mat'];
            case 'vol_limit_check__2_5__3_base_80[m^3]' % #17
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(3);
                res_start = 2.5;
                res = 0.5;
                res_chk_num = 2;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = room_vol*res;
                vargout.save_file = [exp_name '.mat'];
            case 'bug_test_base_40[m^3]_1_1__1_3' % 2 Vol , cubic shape , same beta
                vargout.beta = 0.85;
                vargout.room_shp = room_shp(5,:);
                room_vol = room_vol(2);
                res_start = 1.1;
                res = 0.1;
                res_chk_num = 3;
                res = res_start+((1:res_chk_num)-1)*res;
                vargout.room_vol = [room_vol room_vol*res];
                vargout.save_file = [exp_name '.mat'];
            otherwise
                vargout.save_file = 'rooms_data.mat';
        end
    case 'feature_create'
        if strcmp(exp_name,'vol_basic_cubic')
            vargout.load_exp = 'vol_basic_cubic_data.mat';
            vargout.save_name = 'vol_basic_cubic_data_feat.mat';
        else
            vargout.load_exp = [exp_name '.mat'];
            vargout.save_name = [exp_name '_feat.mat'];
        end
    case 'python_create'
        examine_flag = vargin.examine_flag;
        feat_flag = vargin.feat_flag; % use feature or rirs
        speech_RTF_flag = vargin.speech_RTF_flag; % use RTF of speech
        speech_and_noise_RTF_flag = vargin.speech_and_noise_RTF_flag; % use RTF of speech
        if examine_flag
            simulator_name = vargin.simulator_name;
            exp_num = vargin.exp_num;
            if length(exp_num)>1
                vargout.load_exp =  ['exp_' num2str(exp_num(1)) '_' num2str(exp_num(2)) '_' simulator_name '_examine.mat'];
                vargout.exp_name = exp_name;
            else
                vargout.load_exp =  [exp_name '_' simulator_name '_examine.mat'];
                vargout.exp_name = exp_name;
            end
        else
            if feat_flag
                exp_name = [exp_name '_feat'];
                vargout.exp_name = exp_name;
                vargout.load_exp = [exp_name '.mat'];
            elseif speech_RTF_flag
                exp_name = [exp_name '_speech_RTF'];
                vargout.exp_name = exp_name;
                vargout.load_exp = [exp_name '.mat'];
            elseif speech_and_noise_RTF_flag
                exp_name = [exp_name '_speech_and_noise_RTF'];
                vargout.exp_name = exp_name;
                vargout.load_exp = [exp_name '.mat'];   
            else
                vargout.exp_name = exp_name;
                vargout.load_exp = [exp_name '.mat'];
            end
        end
        
        
end

end
