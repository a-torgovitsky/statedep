%###############################################################################
% MatchPattern
%
% A "pattern" is a wide array sequence with some 0's and 1's but other
% elements that are not 0 and 1
% The elements that are 0 or 1 are "fixed"
% The elements that are not 0 or 1 are "free"
% The routine generates all possible binary sequences over the
% "free" elements, keeping the fixed elements as what they were in Pat.
% Then return this in both Int and Wide format.
%###############################################################################
function [Int Wide] = MatchPattern(Pat)
    L1 = length(Pat);
    IdxFix = find(ismember(Pat, [0 1]));
    IdxFree = setdiff((1:1:L1), IdxFix);
    L2 = length(IdxFree);

    Wide = -1*ones(2^L2, L1);
    Wide(:, IdxFix) = repmat(Pat(IdxFix), size(Wide, 1), 1);
    Wide(:, IdxFree) = AllBinaryArray(L2);
    assert(all(ismember(Wide(:), [0 1])));
    Int = WideToBinary(Wide);
    Int = Int(:);
end

