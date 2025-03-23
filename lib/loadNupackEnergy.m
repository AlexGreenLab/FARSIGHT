function energy = loadNupackEnergy(FileName)
%function energy = loadNupackEnergy(FileName)

out_seq = 1;
fid = fopen(FileName,'r');
if fid == -1
    return;
end
temp_string = fgetl(fid);
energy = [];
while ~isnumeric(temp_string)
    if isempty(temp_string)
        temp_string = fgetl(fid);
        continue;
        energy = [];
        fclose(fid);
        return;
    end
    if temp_string(1) ~= '%' && temp_string(1) ~= '*'
        energy = str2num(temp_string);
        fclose(fid);
        return;
    end
    if temp_string(1) == '*' && ~isempty(findstr(temp_string,'Error'))
        energy = NaN;
        fclose(fid);
        return;
    end
    temp_string = fgetl(fid);
end
fclose(fid);
