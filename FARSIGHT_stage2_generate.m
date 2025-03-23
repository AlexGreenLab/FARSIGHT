%{ 
    FARSIGHT Design Code v. 1.0

    Design Stage 2: Generation and analysis of FARSIGHT designs

    Copyright (c) 2025 Alexander A. Green/Department of Biomedical Engineering, Boston University
    This project is licensed under an Academic Open Source License - see LICENSE.txt file for details
    Contact: aagreen@bu.edu
%} 

addpath('lib');
Ns = char(zeros(1,10000)+'N');
base_prevented = ['T06';'A06';'C04';'G04';'K10';'M10';'R10';'S10';'W10';'Y10'];
prevent_strand_indices = [];

strand_indices = [1,5];
duplicate_indices = [1];
key_defect_indices = [1,1;5,1]; %col1: index; col2: weighting

T7_primer = 'GCGCTAATACGACTCACTATAGGG';
primer_Tm = 57;
primer_index = 5;

L3_len = 6;
complex_linker_len = 8;
complex_linker_seq = char(zeros(1,complex_linker_len)+'A');

%stop_min is the stop condition applied to the design. stop_min can be
%raised or lowered depending on the time the design process takes and
%what the user deems is an acceptable defect level for the designs returned
stop_min = 26;%24;%30;%32;
stop_max = 36;
repeat_counter_max = 2;

%load information on the designs that will be generated
design_info_set = readcell('design_info/FARSIGHT_design_v1_design_info_set.csv');
design_info_nums = readmatrix('design_info/FARSIGHT_design_v1_design_info_nums.csv');

base_file_dir = [pwd,'/NUPACK_base_designs'];
home_dir = pwd;

[~,~,~] = mkdir('NUPACKdesigns');
cd('NUPACKdesigns');

% number of sequences to generate for each design blueprint
num_designs_overall = 2;
design_info_set = [design_info_set];

rng('shuffle');

design_info_set0 = design_info_set;
design_info_set_original = design_info_set;
design_info_nums0 = design_info_nums;
max_aptamer_num = max(design_info_nums(:,end-2));
max_mutation_num = max(design_info_nums(:,end-1));

[~,indices] = sort(rand(size(design_info_set,1),1));
design_info_set = design_info_set0(indices,:);
design_info_nums = design_info_nums0(indices,:);
NUPACKdesigns_dir = pwd;

% This loop cycles through designs of FARSIGHTs with different secondary
% structure blueprints. It also increases the stop condition and prevented
% sequences in case NUPACK design does not return a result.
total_designs_needed = size(design_info_set,1)*num_designs_overall;
while 1
    num_complete = 0;
    cd(NUPACKdesigns_dir);
    remove_list = [];
    for c1 = 1:size(design_info_set,1)
        file_name = design_info_set{c1,1};
        file_list = dir(sprintf('%s/*.npo',file_name));
        num_complete = num_complete + min(length(file_list),num_designs_overall);
        if min(length(file_list)) >= num_designs_overall
            remove_list(end+1,1) = c1;
        end
    end
    fprintf('***** %d out of %d designs are complete (%0.1f%%) *****\n',num_complete,total_designs_needed,num_complete/total_designs_needed*100);
    if num_complete >= total_designs_needed
        break;
    end
    design_info_set_temp = design_info_set;
    design_info_set_temp(remove_list,:) = [];

    for c1 = 1:min(size(design_info_set_temp,1),20)
        result = [];
        cd(NUPACKdesigns_dir);
        file_name = design_info_set_temp{c1,1};
        [prevent_string0,prevent_array0] = generatePreventDNASeq(base_prevented,'');
        [~,~,~] = mkdir(file_name);
        cd(file_name);
        temp_file_list = dir('*.npo');
        design_number = 1;
        design_number_overall = length(temp_file_list) + 1;
        design_number_offset = 1;   
        while ~isempty(dir(sprintf('%s_%d.np',file_name,design_number_offset)))
            design_number_offset = design_number_offset + 1;
        end
        design_number_offset = design_number_offset - 1;
        for prevent_extra = 0:1:4
            for prevent_extra2 = prevent_extra
                for stop_num = stop_min:2:stop_max
                    while design_number_overall <= num_designs_overall
                        design_file_name = sprintf('%s_%d.np',file_name,design_number+design_number_offset);
                        fid = fopen(design_file_name,'w');
                        fid2 = fopen(sprintf('%s/%s_base.txt',base_file_dir,file_name),'r');
                        temp_string = fgetl(fid2);
                        while ~isnumeric(temp_string)
                            fprintf(fid,'%s\n',temp_string);
                            temp_string = fgetl(fid2);
                        end
                        fclose(fid2);
                        fprintf(fid,'\n');
                        fprintf(fid,'prevent L1, L2, L3, L4 = ');
                        temp_prevent_extra_set = [prevent_extra+zeros(1,4),prevent_extra2+zeros(1,20)];
                        for c2 = 1:size(prevent_array0,1)
                            if c2 ~= size(prevent_array0,1)
                                fprintf(fid,'%s%d, ',prevent_array0{c2}(1),length(prevent_array0{c2})+temp_prevent_extra_set(c2));
                            else
                                fprintf(fid,'%s%d\n',prevent_array0{c2}(1),length(prevent_array0{c2})+temp_prevent_extra_set(c2));
                            end
                        end
                        fprintf(fid,'stop[%%] = %d\n',stop_num);
                        fclose(fid);
                        temp_prevent_array = [];
                        for cx2 = 1:size(prevent_array0,1)
                            temp_prevent_array(cx2,1:3) = sprintf('%s%02d',prevent_array0{cx2}(1),length(prevent_array0{cx2})+prevent_extra);
                        end
                        temp_prevent_array = char(temp_prevent_array);
                        fprintf('Running %s out %d of %d (%d,%d,%d)...\n',design_file_name,design_number,num_designs_overall,stop_num,prevent_extra,prevent_extra2);
                        repeat_counter = 1;
                        while repeat_counter <= repeat_counter_max
                            result = call_multitubedesign(design_file_name,temp_prevent_array(1:end,:),prevent_strand_indices,1);
                            if isempty(result)
                                break;
                            end
                            fprintf('\t%s %s (%d,%d,%d).\n',file_name,result,stop_num,prevent_extra,prevent_extra2);
                            repeat_counter = repeat_counter + 1;
                        end
                        if isempty(result)
                            design_number = design_number + 1;
                            design_number_overall = design_number_overall + 1;
                            fprintf('\t%s design complete.\n',design_file_name);
                        else
                            break;
                            fprintf('\tERROR: %s did not result in any sequences. Adjusting design file.\n',design_file_name);
                        end
                    end
                    if design_number_overall > num_designs_overall
                        break;
                    end
                end
            end
            if design_number_overall > num_designs_overall
                break;
            end
        end
    end
end

design_info_set = design_info_set0;
design_info_nums = design_info_nums0;
while 1
    [~,indices] = sort(rand(size(design_info_set,1),1));
    design_info_set = design_info_set(indices,:);
    design_info_set0 = design_info_set;
    keep_indices = ones(size(design_info_set0,1),1);
    file_list = dir(sprintf('zz_*_final_designs.csv'));
    file_list0 = file_list;

    missing_count = 0;
    complete_count = 0;
    for c1 = 1:size(design_info_set,1)
        file_name = design_info_set{c1,1};
        fprintf('Processing %s designs (%d of %d):\n',file_name,c1,size(design_info_set,1));
        full_file_name = sprintf('zz_%s_final_designs.csv',file_name);
        flag = 0;
        for c2 = 1:length(file_list)
            if issame2(full_file_name,file_list(c2).name)
                flag = c2;
                keep_indices(c1) = 0;
                complete_count = complete_count + 1;
            end
        end
        file_list(c2) = [];
        if flag == 0
            missing_count = missing_count + 1;
        end
        if missing_count > 400
            break;
        end
    end

    design_info_set = design_info_set0(find(keep_indices == 1),:);
    output_set = {};
    combined_data = [];
    combined_str_data = {};
    if isempty(design_info_set)
        break;
    end
    for c1 = 1:min(200,size(design_info_set,1))
        file_name = design_info_set{c1,1};

        full_hpin_ON_struc = design_info_set{c1,2};
        full_hpin_INT_struc = design_info_set{c1,3};
        full_hpin_OFF_struc = design_info_set{c1,4};
        full_aptamer_struc = design_info_set{c1,5};
        full_FARSIGHT_struc = design_info_set{c1,6};
        len_aptamer_struc = length(find(full_aptamer_struc ~= 'N'));
        fprintf('Processing %s designs (%d of %d):\n',file_name,c1,size(design_info_set,1))
        if ~isempty(dir(sprintf('zz_%s_final_designs.csv',file_name)))
            fprintf('\tDesign already processed.\n');
        elseif ~isempty(dir(file_name))
            cd(file_name);
            sub_seq_set = {};
            sub_total_defect = [];
            file_list = dir(sprintf('%s_*.npo',file_name));
            for c2 = 1:length(file_list)
                [seq_output,temp_defect] = parseNUPACKfile2SeqInfo(file_list(c2).name);
                if c2 == 1
                    sub_seq_set(1,1:length(strand_indices)) = seq_output(strand_indices);
                    sub_seq_set{1,end+1} = 'Total Defect';
                    sub_seq_set{1,end+1} = 'Selected Defect';
                end
                temp_row = seq_output(strand_indices,2)';
                temp_row{1,end+1} = num2str(temp_defect);
                temp_select_defect = 0;
                for c3 = 1:size(key_defect_indices,1)
                    temp_select_defect = temp_select_defect + str2num(seq_output{key_defect_indices(c3,1),end-1})*key_defect_indices(c3,2);
                end
                temp_select_defect = temp_select_defect/sum(key_defect_indices(:,2));
                temp_row{1,end+1} = num2str(temp_select_defect);
                sub_total_defect(end+1,:) = [temp_defect,temp_select_defect];
                sub_seq_set = [sub_seq_set;temp_row];
            end
            cd ../

            %*********** *********** *********** *********** 
            %check for duplicates in probe sequences
            remove_list = [];
            for c2 = 1:length(duplicate_indices)
                curr_index = duplicate_indices(c2);
                for c3 = 2:size(sub_seq_set,1)
                    for c4 = c3+1:size(sub_seq_set,1)
                        if issame2(sub_seq_set{c3,curr_index},sub_seq_set{c4,curr_index})
                            remove_list(end+1,1) = c4;
                        end
                    end
                end
            end
            sub_seq_set(remove_list,:) = [];
            sub_total_defect(remove_list-1,:) = [];
            fprintf('\tDuplicate designs removed = %d.\n',length(remove_list));
    
            sub_seq_set{1,1} = 'FARSIGHT RNA sequence';
            sub_seq_set{1,2} = 'Trigger RNA sequence';
            insert_set = {'Short FARSIGHT RNA sequence','Joined RNA sequence (energy balance)','Mutant Trigger A','Mutant Trigger B','Mutant Trigger C','Mutant Trigger D'};
            primer_Tm = 57;
            primer_index = 5;
            a_len = design_info_nums(c1,1);
            b_len = design_info_nums(c1,2);
            c_len = design_info_nums(c1,3);
            d_len = design_info_nums(c1,4);
            e_len = design_info_nums(c1,5);
            f_len = design_info_nums(c1,6);
            shorten_len = c_len + d_len + L3_len;
            mutant_pos_set = [...
                3+a_len-1;...
                3+a_len+b_len-1;...
                3+a_len+b_len+c_len-1;...
                3+a_len+b_len+c_len+d_len-1;...
                ];
            for c2 = 2:size(sub_seq_set,1)
                curr_seq = rna2dna2(sub_seq_set{c2,1});
                curr_seq0 = curr_seq;
                curr_seq = [T7_primer(1:end-3),concatDNA(T7_primer(end-2:end),curr_seq)];
                FARSIGHT_seq = curr_seq;
                anti_FARSIGHT_seq = d_revcomp(curr_seq);
    
                curr_seq = curr_seq0(1:end-shorten_len);
                short_FARSIGHT_RNA_seq = dna2rna2(curr_seq);
                short_FARSIGHT_seq = [T7_primer(1:end-3),concatDNA(T7_primer(end-2:end),curr_seq)];
    
                curr_seq = rna2dna2(sub_seq_set{c2,2});
                curr_seq = [T7_primer(1:end-3),concatDNA(T7_primer(end-2:end),curr_seq)];
                anti_trigger_seq = d_revcomp(curr_seq);
    
                temp_FARSIGHT_seq = sub_seq_set{c2,1};
                temp_trigger_seq = sub_seq_set{c2,2};
                joined_seq = [temp_trigger_seq,complex_linker_seq,temp_FARSIGHT_seq];
    
                mutant_trigger_set = {};
                base_trigger_seq = sub_seq_set{c2,2};
                for c3 = 1:length(mutant_pos_set)
                    mutant_pos = mutant_pos_set(c3);
                    new_trigger_seq = base_trigger_seq;
                    new_trigger_seq(mutant_pos) = r_revcomp(new_trigger_seq(mutant_pos));
                    mutant_trigger_set{1,c3} = new_trigger_seq;
                end
    
                insert_set(end+1,:) = [short_FARSIGHT_RNA_seq,joined_seq,mutant_trigger_set];
            end
            extra_labels = {'defect ON','defect INT','defect OFF','min complex defect','aptamer defect CT','aptamer defect mutA','aptamer defect mutB','aptamer defect mutC','aptamer defect mutD','FARSIGHT defect','Predicted deltaG','dg ON','dg INT','dg OFF','ddg(ON - INT)','ddg(ON - OFF)','ddg(INT-OFF)','Mean dg','Std. dg'};
            base_label = {'Name','Total Defect','Selected Defect'};
            full_label = [base_label,extra_labels];
            temp_calc_table = [];
            for c2 = 2:size(sub_seq_set,1)
                fprintf('\tAssessing design %d of %d...\n',c2-1,size(sub_seq_set,1)-1);
                FARSIGHT_seq = sub_seq_set{c2,1};
                trigger_seq = sub_seq_set{c2,2};
                curr_seq = [trigger_seq,complex_linker_seq,FARSIGHT_seq];
    
                defect_ON_full = checkDefectArb1999(curr_seq,full_hpin_ON_struc);
                defect_INT_full = checkDefectArb1999(curr_seq,full_hpin_INT_struc);
                defect_OFF_full = checkDefectArb1999(curr_seq,full_hpin_OFF_struc);
                defect_FARSIGHT_full = checkDefectArb1999(FARSIGHT_seq,full_FARSIGHT_struc);
    
                dg_ON_full = checkStructureFreeEnergyRNA1999(curr_seq,full_hpin_ON_struc);
                dg_INT_full = checkStructureFreeEnergyRNA1999(curr_seq,full_hpin_INT_struc);
                dg_OFF_full = checkStructureFreeEnergyRNA1999(curr_seq,full_hpin_OFF_struc);
    
                [deltaG_full,~] = computeRNAdeltaG1999(curr_seq);
    
                mean_dg = mean([dg_ON_full,dg_INT_full,dg_OFF_full]);
                std_dg = std([dg_ON_full,dg_INT_full,dg_OFF_full]);
    
                matrix_full_aptamer_struc = convertStrucString2PairProbTable(full_aptamer_struc);
                pp_curr_seq = computePairProbTableRNA1999(curr_seq);
                aptamer_struc_defect_min = 1-sum(sum(matrix_full_aptamer_struc.*pp_curr_seq))./len_aptamer_struc;
    
                base_trigger_seq = sub_seq_set{c2,2};
                mutant_aptamer_struc_defect_min_set = [];
                %Generate a set of mutant triggers compared to the correct
                %target to assess the sensitivity of the FARSIGHT to target
                %sequence changes.
                for c3 = 1:length(mutant_pos_set)
                    mutant_pos = mutant_pos_set(c3);
                    new_trigger_seq = base_trigger_seq;
                    new_trigger_seq(mutant_pos) = r_revcomp(new_trigger_seq(mutant_pos));
                    mut_curr_seq = [new_trigger_seq,complex_linker_seq,FARSIGHT_seq];
                    pp_curr_seq = computePairProbTableRNA1999(mut_curr_seq);
                    mutant_aptamer_struc_defect_min_set(1,c3) = 1-sum(sum(matrix_full_aptamer_struc.*pp_curr_seq))./len_aptamer_struc;
                end
                temp_calc_table(c2-1,:) = [defect_ON_full,defect_INT_full,defect_OFF_full,min([defect_ON_full,defect_INT_full,defect_OFF_full]),aptamer_struc_defect_min,mutant_aptamer_struc_defect_min_set,defect_FARSIGHT_full,...
                    deltaG_full,dg_ON_full,dg_INT_full,dg_OFF_full,dg_ON_full-dg_INT_full,dg_ON_full-dg_OFF_full,dg_INT_full-dg_OFF_full,mean_dg,std_dg];
            end
            temp_calc_str_table = {};
            for c2 = 1:size(temp_calc_table,1)
                for c3 = 1:size(temp_calc_table,2)
                    temp_calc_str_table{c2,c3} = num2str(temp_calc_table(c2,c3));
                end
            end
            sub_total_defect = [sub_total_defect,temp_calc_table];
            sub_seq_set = [sub_seq_set,insert_set,[extra_labels;temp_calc_str_table]];
    
            %just sort designs at this stage based on st. dev. of dg values of the ON,
            %OFF, and INT states
            [~,indices] = sort(sub_total_defect(:,end));
            sub_total_defect0 = sub_total_defect;
            sub_seq_set0 = sub_seq_set;
    
            sub_total_defect = sub_total_defect(indices,:);
            sub_seq_set = sub_seq_set([1;indices+1],:);
            temp_set = sub_seq_set(:,3:4);
            sub_seq_set(:,3:4) = [];
            calc_len = size(temp_calc_table,2);
            sub_seq_set = [sub_seq_set(:,1:end-calc_len),temp_set,sub_seq_set(:,end-calc_len+1:end)];
            writecell(sub_seq_set,sprintf('zz_%s_final_designs.csv',file_name));
        end
    end
end
design_info_set = design_info_set_original;

output_set = {};
extra_labels = {'defect ON','defect INT','defect OFF','min complex defect','aptamer defect CT','aptamer defect mutA','aptamer defect mutB','aptamer defect mutC','aptamer defect mutD','FARSIGHT defect','Predicted deltaG','dg ON','dg INT','dg OFF','ddg(ON - INT)','ddg(ON - OFF)','ddg(INT-OFF)','Mean dg','Std. dg'};
base_label = {'Name','Total Defect','Selected Defect'};
full_label = [base_label,extra_labels];

combined_data = [];
combined_str_data = {};
for c1 = 1:size(design_info_set,1)
       
    file_name = design_info_set{c1,1};
    
    full_hpin_ON_struc = design_info_set{c1,2};
    full_hpin_INT_struc = design_info_set{c1,3};
    full_hpin_OFF_struc = design_info_set{c1,4};
    full_aptamer_struc = design_info_set{c1,5};
    full_FARSIGHT_struc = design_info_set{c1,6};
    len_aptamer_struc = length(find(full_aptamer_struc ~= 'N'));
    fprintf('Processing %s designs:\n',file_name)
    if ~isempty(dir(sprintf('zz_%s_final_designs.csv',file_name)))
        sub_seq_set = readcell(sprintf('zz_%s_final_designs.csv',file_name));
        sub_total_defect = [];
        for c2 = 2:size(sub_seq_set,1)
            for c3 = 9:size(sub_seq_set,2)
                sub_total_defect(c2-1,c3-8) = sub_seq_set{c2,c3};
            end
        end
        left_temp_data = {'Name'};
        for c2 = 2:size(sub_seq_set,1)
            left_temp_data{c2,1} = sprintf('%s_N%02d',file_name,c2-1);
        end
        sub_seq_set = [left_temp_data,sub_seq_set];
        if c1 == 1
            combined_str_data = sub_seq_set;
        else
            combined_str_data = [combined_str_data;sub_seq_set(2:end,:)];
        end
        combined_data = [combined_data;[zeros(size(sub_total_defect,1),1)+c1,sub_total_defect]];
    end
end
writecell(combined_str_data,'AAB_combined_seq_design_data.csv');
table_labels = full_label;
table_labels{1} = 'Design Index';
combined_table = array2table(combined_data,"VariableNames",table_labels);
writetable(combined_table,'AAB_combined_seq_design_nums.csv');

cd ../;