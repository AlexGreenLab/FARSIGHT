function du_string = DotBracket2DUnotation(dot_bracket)
%function du_string = DotBracket2DUnotation(dot_bracket)

dot_bracket(findstr(dot_bracket,' ')) = [];
left_string = [];
c1 = 1;
possible_letters = {'D','U','+','N'};
du_string = [];
while c1 <= length(dot_bracket)
    if dot_bracket(c1) == '.'
        num_dots = 1;
        for c2 = c1+1:length(dot_bracket)
            if dot_bracket(c2) == '.'
                num_dots = num_dots + 1;
            else
                break;
            end
        end
        left_string = [left_string,'U',num2str(num_dots),' '];
        c1 = c1 + num_dots;
    elseif dot_bracket(c1) == 'N'
        num_dots = 1;
        for c2 = c1+1:length(dot_bracket)
            if dot_bracket(c2) == 'N'
                num_dots = num_dots + 1;
            else
                break;
            end
        end
        left_string = [left_string,'N',num2str(num_dots),' '];
        c1 = c1 + num_dots;
    elseif dot_bracket(c1) == '+'
        left_string = [left_string,'+ '];
        c1 = c1 + 1;
    elseif dot_bracket(c1) == '('
        num_brackets = 1;
        for c2 = c1+1:length(dot_bracket)
            if dot_bracket(c2) == '('
                num_brackets = num_brackets + 1;
            else
                break;
            end
        end
        max_num_brackets = num_brackets;
        consecutive_brackets = 0;
        for c2 = c1+max_num_brackets:length(dot_bracket)
            if dot_bracket(c2) == ')'
                consecutive_brackets = consecutive_brackets + 1;
                num_brackets = num_brackets - 1;
                if num_brackets == 0
                    consecutive_brackets = min(max_num_brackets,consecutive_brackets);
                    left_string = [left_string,'D',num2str(consecutive_brackets),' '];
                    new_du_string = DotBracket2DUnotation(dot_bracket(c1+consecutive_brackets:c2-consecutive_brackets));
                    string_pos = [];
                    for c3 = 1:length(possible_letters)
                        string_pos = [string_pos,findstr(new_du_string,possible_letters{c3})];
                    end
                    if length(string_pos) == 1
                        left_string = [left_string,new_du_string];
                    else
                        left_string = [left_string,'(',new_du_string(1:end-1),') '];
                    end
                    c1 = c2 + 1;
                    break;
                end
            elseif dot_bracket(c2) == '('
                consecutive_brackets = 0;
                num_brackets = num_brackets + 1;
            else
                consecutive_brackets = 0;
            end
        end
    end
end
du_string = left_string;
