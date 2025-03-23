%{ 
    FARSIGHT Design Code v. 1.0

    Design Stage 1: Definition of parameters for FARSIGHT designs to be generated

    Copyright (c) 2025 Alexander A. Green/Department of Biomedical Engineering, Boston University
    This project is licensed under an Academic Open Source License - see LICENSE.txt file for details
    Contact: aagreen@bu.edu
%} 

addpath('lib');
T = 37;
Ns = char(zeros(1,1000)+'N');
dots = char(zeros(1,1000)+'.');

design_info_set = {};
design_info_nums = [];

dock_len = 21;
L1_len = 10;
L4_len = L1_len;
L2_len = L1_len*2;
L3_len = 6;

e_len = 3;
complex_linker_len = 8;
complex_linker_seq = char(zeros(1,complex_linker_len)+'A');

% mut_pos_range = [0,1,2,3];
mut_pos_range = [0,1]; %limit mutation position to just two locations to reduce design time

% define domain length ranges
% forward and reverse toehold 1: domains a and f
% forward and reverse toehold 2: domains b and d
domain_range_set_labels = {'a','b','c','d','e','f','M','L1','L2','L3','L4'};
domain_range_set = [];
for a_len = 3:4
    for f_len = a_len-1:a_len+1
        for M_len = mut_pos_range + a_len
            for b_len = 3:5
                for d_len = b_len-2:b_len-1
                    for c_len = 12:12 %10-b_len:2:14-b_len
                        domain_range_set(end+1,:) = [a_len,b_len,c_len,d_len,e_len,f_len,M_len,L1_len,L2_len,L3_len,L4_len];
                    end
                end
            end
        end
    end
end

input_aptamer_set;

for c1 = 1:size(aptamer_info_set,1)
    aptamer_info_set{c1,3} = DotParens2DUnotation(aptamer_info_set{c1,3});
end

input_RNA_table = readcell('mutant_target_input.csv');
input_RNA_table(1,:) = [];

%define the secondary structures and sequences of the designs to be
%generated using different combinations of domain lengths. 
design_info_set = cell(size(aptamer_info_set,1)*size(input_RNA_table,1)*size(domain_range_set,1),18);
design_info_nums = zeros(size(aptamer_info_set,1)*size(input_RNA_table,1)*size(domain_range_set,1),9);
counter = 1;
for cn2 = 1:size(aptamer_info_set,1)
    aptamer_core_name = aptamer_info_set{cn2,1};
    aptamer_core_seq = dna2rna2(aptamer_info_set{cn2,2});
    aptamer_core_strucDU = aptamer_info_set{cn2,3};
    for cn1 = 1:size(input_RNA_table,1)
        targ_name = input_RNA_table{cn1,1};
        WT_seq0 = input_RNA_table{cn1,2};
        SN_seq0 = input_RNA_table{cn1,3};
        fprintf('%s: Generating structures for %s (%d of %d)...\n',aptamer_core_name,targ_name,cn1,size(input_RNA_table,1));
        dim1 = size(input_RNA_table,1);
        dim2 = size(domain_range_set,1);
        for c0 = 1:size(domain_range_set,1)
            name_block = [];
            for c1 = 1:length(domain_range_set_labels)
                eval(sprintf('%s_len = domain_range_set(c0,c1);',domain_range_set_labels{c1}));
                if c1 == 3
                    name_block = sprintf('%s_%s%02d',name_block,domain_range_set_labels{c1},domain_range_set(c0,c1));
                elseif c1 <= 7
                    name_block = sprintf('%s_%s%d',name_block,domain_range_set_labels{c1},domain_range_set(c0,c1));
                end
            end
            name_block(1) = '';
            targ_len = e_len + d_len + c_len + b_len + a_len + L4_len + dock_len;
            mut_pos_fwd = targ_len - M_len - dock_len - L4_len + 1;
            index = find(WT_seq0(1:min(end,length(SN_seq0))) ~= SN_seq0(1:min(end,length(WT_seq0))));
            index_start = index(1) - mut_pos_fwd + 1;
            index_end = index_start + targ_len - 1;

            WT_seq = WT_seq0(index_start:index_end);
            SN_seq = SN_seq0(index_start:index_end);

            SN_seq_star = r_revcomp(SN_seq);
            dock_dom = SN_seq_star(1:dock_len);
            a_dom = SN_seq_star(dock_len+L4_len+1:dock_len+L4_len+a_len);
            b_dom = SN_seq_star(dock_len+L4_len+a_len+1:dock_len+L4_len+a_len+b_len);
            c_dom = SN_seq_star(dock_len+L4_len+a_len+b_len+1:dock_len+L4_len+a_len+b_len+c_len);
            d_dom = SN_seq_star(dock_len+L4_len+a_len+b_len+c_len+1:dock_len+L4_len+a_len+b_len+c_len+d_len);
            e_dom = SN_seq_star(dock_len+L4_len+a_len+b_len+c_len+d_len+1:dock_len+L4_len+a_len+b_len+c_len+d_len+e_len);
            L4_dom = SN_seq_star(dock_len+1:dock_len+L4_len);
            assert(~isempty(strfind(SN_seq,r_revcomp([a_dom,b_dom,c_dom,d_dom,e_dom]))))

            %assume complex will be joined by A8 sequence
            ON_struc = DUnotation2DotParens(sprintf('D%d (U%d D%d (U%d U3) U%d) U%d D%d (%s) U%d',e_len+d_len+c_len+b_len+a_len,...
                L4_len,dock_len,complex_linker_len,L1_len,f_len*2+L2_len+e_len+d_len,c_len+b_len,aptamer_core_strucDU,d_len+L3_len+d_len+c_len));
            OFF_struc = DUnotation2DotParens(sprintf('U%d D%d (U%d U3) U%d D%d U%d %s U%d D%d U%d',e_len+d_len+c_len+b_len+a_len+L4_len,...
                dock_len,complex_linker_len,L1_len+a_len,b_len+c_len+d_len+e_len+f_len,L2_len,aptamer_core_strucDU,b_len,c_len+d_len,L3_len));

            %includes unwinding of first stem, but not the aptamer formation
            INT_struc = DUnotation2DotParens(sprintf('D%d (U%d D%d (U%d U3) U%d) U%d %s U%d D%d U%d',e_len+d_len+c_len+b_len+a_len,...
                L4_len,dock_len,complex_linker_len,L1_len,f_len*2+L2_len+e_len+d_len+c_len+b_len,aptamer_core_strucDU,b_len,c_len+d_len,L3_len));

            %includes only the aptamer secondary structure
            aptamer_struc = DUnotation2DotParens(sprintf('N%d N%d N%d N%d N3 N%d N%d D%d (%s) N%d',(e_len+d_len+c_len+b_len+a_len)*2,...
                L4_len,dock_len*2,complex_linker_len,L1_len,f_len*2+L2_len+e_len+d_len,c_len+b_len,aptamer_core_strucDU,d_len+L3_len+d_len+c_len));

            FARSIGHT_struc = DUnotation2DotParens(sprintf('U3 U%d U%d U%d D%d U%d %s U%d D%d U%d',dock_len,L1_len,a_len,b_len+c_len+d_len+e_len+f_len,L2_len,aptamer_core_strucDU,b_len,c_len+d_len,L3_len));

            design_name = sprintf('%s_%s_FARSIGHT_%s',aptamer_core_name,targ_name,name_block);
            design_info_set((cn2-1)*dim1*dim2+(cn1-1)*dim2+c0,:) = {design_name,ON_struc,INT_struc,OFF_struc,aptamer_struc,FARSIGHT_struc,SN_seq,WT_seq,SN_seq0,WT_seq0,dock_dom,a_dom,b_dom,c_dom,d_dom,e_dom,L4_dom,aptamer_core_seq};
            design_info_nums((cn2-1)*dim1*dim2+(cn1-1)*dim2+c0,:) = [a_len,b_len,c_len,d_len,e_len,f_len,cn2,cn1,c0];%,findstr(full_seq,'AUG')];
        end
    end
end
indices = find(design_info_nums(:,7) == 0);
design_info_nums(indices,:) = [];
design_info_set(indices,:) = [];

[~,~,~] = mkdir('design_info');
writecell(design_info_set,'design_info/FARSIGHT_design_v1_design_info_set.csv');
design_info_table = table();
design_info_table.a_len = design_info_nums(:,1);
design_info_table.b_len = design_info_nums(:,2);
design_info_table.c_len = design_info_nums(:,3);
design_info_table.d_len = design_info_nums(:,4);
design_info_table.e_len = design_info_nums(:,5);
design_info_table.f_len = design_info_nums(:,6);
writetable(design_info_table,'design_info/FARSIGHT_design_v1_design_info_nums.csv');

base_dir = 'NUPACK_base_designs';
[~,~,~] = mkdir(base_dir);
cd(base_dir);

%Output the design information to files that will later be used to run
%NUPACK design function on to generate FARSIGHT sequences
domain_input = {'a','b','c','d','e','f','L1','L2','L3','L4','dock'}';
for c1 = 1:size(design_info_set,1)
    base_name = design_info_set{c1,1};
    if mod(c1,100) == 0
        fprintf('Outputting %s design (%d of %d)...\n',base_name,c1,size(design_info_set,1));
    end
    a_len = design_info_nums(c1,1);
    b_len = design_info_nums(c1,2);
    c_len = design_info_nums(c1,3);
    d_len = design_info_nums(c1,4);
    e_len = design_info_nums(c1,5);
    f_len = design_info_nums(c1,6);
    aptamer_index = design_info_nums(c1,7);
    
    aptamer_core_name = aptamer_info_set{aptamer_index,1};
    aptamer_core_seq = dna2rna2(aptamer_info_set{aptamer_index,2});
    aptamer_core_strucDU = aptamer_info_set{aptamer_index,3};

    domain_set = {};
    for c2 = 1:length(domain_input)
        domain_len = eval(sprintf('%s_len',domain_input{c2,1}));
        if issame2(domain_input{c2,1},'dock')
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,11}};
        elseif domain_input{c2,1} == 'a'
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,12}};
        elseif domain_input{c2,1} == 'b'
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,13}};
        elseif domain_input{c2,1} == 'c'
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,14}};
        elseif domain_input{c2,1} == 'd'
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,15}};
        elseif domain_input{c2,1} == 'e'
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,16}};
        elseif sum(domain_input{c2,1} == 'L4') == 2
            domain_set(end+1,:) = {domain_input{c2,1},design_info_set{c1,17}};
        else
            domain_set(end+1,:) = {domain_input{c2,1},Ns(1:domain_len)};
        end
    end
    domain_set(end+1,:) = {'preGGG','GGG'};...
    domain_set(end+1,:) = {'x',aptamer_core_seq};...

    strand_set = {...
        'FARSIGHT_s','preGGG dock L1 a b c d e f L2 f* e* d* c* b* x b c d L3 d* c*';...
        'trigger_s','e* d* c* b* a* L4* dock*';...
        'left_s','preGGG dock L1 a b c d e f L2';...
        'middle_s','L2 f* e* d* c* b*';...
        'right_s','x b c d L3 d* c*';...
        };
    complex_set = {...
        'FARSIGHT','FARSIGHT_s';...
        'trigger','trigger_s';...
        'left','left_s';...
        'middle','middle_s';...
        'right','right_s';...
        'trigger_FARSIGHT','trigger_s FARSIGHT_s';...
        };
    structure_set = {...
        'FARSIGHT',sprintf('U3 U%d U%d U%d D%d U%d U%d U%d D%d U%d',dock_len,L1_len,a_len,b_len+c_len+d_len+e_len+f_len,L2_len,length(aptamer_core_seq),b_len,c_len+d_len,L3_len);...
        'trigger',sprintf('U%d',e_len+d_len+c_len+b_len+a_len+L4_len+dock_len);...
        'left',sprintf('U3 U%d',dock_len+L1_len+a_len+b_len+c_len+d_len+e_len+f_len+L2_len);...
        'middle',sprintf('U%d',L2_len+f_len+e_len+d_len+c_len+b_len);...
        'right',sprintf('U%d U%d D%d U%d',length(aptamer_core_seq),b_len,c_len+d_len,L3_len);...
        'trigger_FARSIGHT',sprintf('D%d (U%d D%d (+ U3) U%d) U%d D%d (%s) U%d',e_len+d_len+c_len+b_len+a_len,...
            L4_len,dock_len,L1_len,f_len*2+L2_len+e_len+d_len,c_len+b_len,aptamer_core_strucDU,d_len+L3_len+d_len+c_len);...
            };

    design_file_name = sprintf('%s_base.txt',base_name);
    fid = fopen(design_file_name,'w');
    fprintf(fid,'material = rna1999\n');
    fprintf(fid,'temperature = %d\n',T);
    fprintf(fid,'trials = 10\n');
    fprintf(fid,'allowwobble = true\n\n');

    for c2 = 1:size(domain_set,1)
        fprintf(fid,'domain %s = %s\n',domain_set{c2,1},domain_set{c2,2});
    end
    fprintf(fid,'\n');
    for c2 = 1:size(strand_set,1)
        fprintf(fid,'strand %s = %s\n',strand_set{c2,1},strand_set{c2,2});
    end
    fprintf(fid,'\n');
    for c2 = 1:size(complex_set,1)
        fprintf(fid,'complex %s = %s\n',complex_set{c2,1},complex_set{c2,2});
    end
    fprintf(fid,'\n');
    for c2 = 1:size(structure_set,1)
        fprintf(fid,'%s.structure = %s\n',structure_set{c2,1},structure_set{c2,2});
    end
    fprintf(fid,'\n');
    fclose(fid);
end
cd ../