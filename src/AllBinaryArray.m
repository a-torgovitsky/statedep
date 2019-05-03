%###############################################################################
% AllBinaryArray
%
% Return in wide format all binary sequences of length Len.
%###############################################################################
function [YWide] = AllBinaryArray(Len)
    YWide = IntToBinaryArray((0:1:(2^Len - 1)), Len);
end

%###############################################################################
% IntToBinaryArray
%
% YInt is the integer representation for binary strings of length Len
% YWideStr converts these to strings.
% YWide converts these to a "wide format" array.
%###############################################################################
function [YWide] = IntToBinaryArray(YInt, Len)
    YWideStr = dec2bin(YInt, Len);
    YWide = YWideStr - '0'; % Converts strings to arrays
end
