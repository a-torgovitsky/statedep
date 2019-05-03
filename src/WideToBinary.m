%###############################################################################
% WideToBinary.m
%
% Y is a "wide format" array of 0's and 1's.
% YStr is this same sequence of 0's and 1's in a string.
% YInt is the integer representation of this binary string.
%###############################################################################
function [YInt YStr] = WideToBinary(Y)
    assert(all(ismember(Y(:), [0, 1])));
    for i = 1:1:size(Y,1)
        YStr{i} = sprintf('%d', Y(i,:));

        % 01/28/17
        % For some idiotic long-forgotten reason I decided to start my binary
        % numbering at 1 instead of 0. Too much trouble to change it now.
        YInt(i) = bin2dec(YStr{i}) + 1;
    end
end
