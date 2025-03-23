function [deltaG,dot_struc] = loadNupackMFEsingleStrand(file_name)
%function [deltaG,dot_struc] = loadNupackMFEsingleStrand(file_name)

fid = fopen(file_name,'r');
%move to complex data
search_string = '% %%%%%%%%%%%%%%%%%%%%';
temp_string = fgetl(fid);
while ischar(temp_string)
    if length(findstr(temp_string,search_string)) ~= 0
        break;
    end
    temp_string = fgetl(fid);
end
temp_string = fgetl(fid);
deltaG = str2num(fgetl(fid));
dot_struc = fgetl(fid);
fclose(fid);