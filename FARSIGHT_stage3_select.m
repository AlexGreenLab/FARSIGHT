%{ 
    FARSIGHT Design Code v. 1.0

    Design Stage 3: Compilation, scoring, and selection of top FARSIGHT designs

    Copyright (c) 2025 Alexander A. Green/Department of Biomedical Engineering, Boston University
    This project is licensed under an Academic Open Source License - see LICENSE.txt file for details
    Contact: aagreen@bu.edu
%} 

addpath('lib');
base_file_dir = [pwd,'/NUPACK_base_designs'];
base_mfile_dir = pwd;
input_aptamer_set;

combined_str_data = readtable('NUPACKdesigns/AAB_combined_seq_design_data.csv');
combined_data = readmatrix('NUPACKdesigns/AAB_combined_seq_design_nums.csv');
% combined_str_labels = combined_str_data(1,:);
% combined_str_data(1,:) = [];
num_designs = 8;
target_set = readcell('mutant_target_input.csv');
target_set(1,:) = [];
output_dir = 'FARSIGHT_design_output';
[~,~,~] = mkdir(output_dir);
combined_design_output = [];
for c0 = 1:size(aptamer_info_set,1)
    for c1 = 1:size(target_set,1)
        index_set = zeros(size(combined_str_data,1),1);
        base_name_string = sprintf('%s_%s',aptamer_info_set{c0,1},target_set{c1,1});
        for c2 = 1:size(combined_str_data,1)
            temp = strfind(combined_str_data.Name{c2},base_name_string);
            if ~isempty(temp) && temp == 1
                index_set(c2,1) = 1;
            end
        end
        indices = find(index_set == 1);
        sub_combined_str_data = combined_str_data(indices,:);
        sub_combined_data = combined_data(indices,:);
        
        %Scoring designs

        %Target ddgINT_OFF values of approximately -0.5, ideally [-1.5,0.5]
        ddgINT_OFF0 = -0.5; %optimal value for ddg(INT-OFF)
        ddgOFF_INT_score = abs(sub_combined_str_data.ddg_INT_OFF_ - 0.5);
        
        %Target ddgON_INT values <= -1.5. Empirically designs between -15
        %and -4 performed well with preference for designs around -10.
        ddgON_INT0 = -10;
        ddgON_INT_score = (sub_combined_str_data.ddg_ON_INT_ <= -1.5) .* abs(sub_combined_str_data.ddg_ON_INT_ - ddgON_INT0)/6 + (sub_combined_str_data.ddg_ON_INT_ > -1.5)*100;

        %FARSIGHT defect scores. Eliminate designs where FARSIGHT or
        %aptamerDefect with correct target is greater than 0.5.
        total_defect_score = (sub_combined_str_data.FARSIGHTDefect <= 0.5 | sub_combined_str_data.aptamerDefectCT <= 0.5) .* (sub_combined_str_data.FARSIGHTDefect + sub_combined_str_data.aptamerDefectCT) + ...
            (sub_combined_str_data.FARSIGHTDefect > 0.5 | sub_combined_str_data.aptamerDefectCT > 0.5)*100;
        
        % Add together components to generate an overall design score.
        % Lower scores are better.
        total_design_score = ddgOFF_INT_score + ddgON_INT_score + total_defect_score;
        sub_combined_str_data.total_design_score = total_design_score;
        [~,indices] = sort(total_design_score);
        sub_combined_str_data = sub_combined_str_data(indices,:);
        index = length(base_name_string) + 10;
        for c2 = 1:size(sub_combined_str_data,1)
            curr_name = sub_combined_str_data.Name{c2};
            left_name = curr_name(1:index);
            right_name = curr_name(index+1:end);
            new_name = sprintf('%srank%03d_%s',left_name,c2,right_name);
            sub_combined_str_data.Name{c2} = new_name;
        end
        writetable(sub_combined_str_data,sprintf('%s/%s_full_design_info.csv',output_dir,base_name_string))
        combined_design_output = [combined_design_output;sub_combined_str_data(1:num_designs,:)];
    end
end
writetable(combined_design_output,sprintf('%s/FARSIGHT_top_designs.csv',output_dir));