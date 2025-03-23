function dot_bracket = DUnotation2DotBracket(du_string)
%function dot_bracket = DUnotation2DotBracket(du_string)

du_string(findstr(du_string,' ')) = [];
left_string = [];
right_string = [];
c1 = 1;
while c1 <= length(du_string)
    if du_string(c1) == '('
        left_pos = c1+1;
        bracket_num = 1;
        for c2 = c1+1:length(du_string)
            if du_string(c2) == '('
                bracket_num = bracket_num + 1;
            elseif du_string(c2) == ')'
                if bracket_num == 1
                    right_pos = c2 - 1;
                    break;
                else
                    bracket_num = bracket_num - 1;
                end
            end
        end
        left_string = [left_string,DUnotation2DotBracket(du_string(left_pos:right_pos)),right_string];
        right_string = [];
        c1 = right_pos + 2;
    elseif du_string(c1) == '+'
        left_string = [left_string,'+',right_string];
        right_string = [];
        c1 = c1 + 1;
    else
        curr_letter = du_string(c1);
        c1 = c1 + 1;
        num_string = du_string(c1);
        c1 = c1 + 1;
        while c1 <= length(du_string) && ~isempty(str2num(du_string(c1)))
            num_string = [num_string,du_string(c1)];
            c1 = c1 + 1;
        end
        curr_number = str2num(num_string);
        if curr_letter == 'D'
            left_string = [left_string,char(zeros(1,curr_number)+'(')];
            right_string = [char(zeros(1,curr_number)+')'),right_string];
        elseif curr_letter == 'U'
            left_string = [left_string,char(zeros(1,curr_number)+'.'),right_string];
            right_string = [];
        elseif curr_letter == 'N'
            left_string = [left_string,char(zeros(1,curr_number)+'N'),right_string];
            right_string = [];
        end
    end
end
dot_bracket = left_string;
        

% 
% counter = 1;
% while du_string(counter) ~= 'D' || du_string(counter) ~= 'U' || du_string(counter) ~= '+'
%     counter = counter + 1;
% end
% curr_letter = du_string(counter);
% if curr_letter = '+';
%     dot_bracket = '+';
%     return;
% end
% counter = counter + 1;
% number_string = [];
% while ~isempty(str2num(du_string(counter)))
%     num_string = [num_string,du_string(counter)];
%     counter = counter + 1;
% end
% curr_number = str2num(num_string);
% if curr_letter == 'D'
%     dot_bracket = [char(zeros(1,curr_number)+'('),DUnotation2DotBracket(du_string(counter+end))