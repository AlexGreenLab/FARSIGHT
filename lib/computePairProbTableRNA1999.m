function pair_prob_table = computePairProbTableRNA1999(rseq)
%function pair_prob_table = computePairProbTableRNA1999(rseq)

if length(rseq) == 0
    struc_out = char(zeros(1,10)+'(');
    return 
end
filename = randseq(8);
fid = fopen([filename,'.in'],'w');
fprintf(fid,'%s',rseq);
fclose(fid);
a = 1;
counter = 1;
while a ~= 0 && counter < 5
    a = system(sprintf('pairs -material rna1999 -T 37 %s',filename));
    if a ~= 0
        continue;
    end
    pair_prob_table = loadNupackPairProbTable([filename,'.ppairs']);
    if isempty(pair_prob_table)
        fprintf('computePairProbTableRNA1999: pair_prob_table is empty.\n');
        a = 1;
        continue;
    end
    counter = counter + 1;
end

%pair_prob_table = loadNupackPairProbTableMulti([filename,'.ppairs'],length(findstr(rseq,'+')));
% fid = fopen([filename,'.ppairs'],'r');
% if fid == -1
%     pair_prob_table = 0;
%     return;
% end
% temp_string = fgetl(fid);
% while length(temp_string) ~= 0
%     temp_string = fgetl(fid);
% end
% temp_string = fgetl(fid); % % line
% temp_string = fgetl(fid); % length
% temp_string = fgetl(fid); % deltaG
% struc_out = fgetl(fid);
% fclose(fid);
system(sprintf('rm %s.in %s.ppairs',filename,filename));