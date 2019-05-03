%###############################################################################
% CFHN
%###############################################################################
function [CFHN] = ComputeCFHNBounds(Settings, Data)
    for d = 0:1:1
        tot = 0;
        for t = 1:1:Settings.T
            % Find all sequences of y such that
            % Y(s) = (1-d) for all s < t
            % Y(t) = d
            % Y(t+1) = 1
            % sum up their probabilities
            % sum up over all t
            % this then becomes the lower bound
            % see pg. 553 of CFHN

            YPat = -1*ones(1, 1 + Settings.T);
            YPat(1:(t-1)) = (1-d);
            YPat(t) = d;
            YPat(t+1) = 1;

            [~, YWide] = MatchPattern(YPat);
            Match = ismember(Data.Y, YWide, 'rows');

            tot = tot + sum(Match);
        end
        CFHN.ASFLB(d+1) = tot/Settings.N;

        % Find all sequences of y such that
        % Y(s) = (1-d) for all s <= (T-1)
        % the probability over these sequences is the width of the
        % bounds -- see pg. 553 of CFHN
        YPat = -1*ones(1, 1 + Settings.T);
        YPat(1:Settings.T) = (1-d);
        [~, YWide] = MatchPattern(YPat);
        Match = ismember(Data.Y, YWide, 'rows');
        CFHN.ASFUB(d+1) = CFHN.ASFLB(d+1) + sum(Match)/Settings.N;
    end
    % Implied bounds on ATE
    CFHN.ATELB = CFHN.ASFLB(2) - CFHN.ASFUB(1);
    CFHN.ATEUB = CFHN.ASFUB(2) - CFHN.ASFLB(1);

    % Adding monotonicity does not affect lower bound in untreated state
    CFHN.ASFLBMon(1) = CFHN.ASFLB(1);

    % Monotonicity increases lower bound in treated state
    % by P[Y = (0,0,...,0,1)]
    YPat = -1*ones(1, 1 + Settings.T);
    YPat(1:Settings.T) = 0;
    YPat(Settings.T+1) = 1;
    [~, YWide] = MatchPattern(YPat);
    Match = ismember(Data.Y, YWide, 'rows');
    CFHN.ASFLBMon(2) = CFHN.ASFLB(2) + sum(Match)/Settings.N;

    % Adding monotonicity does not affect upper bound in treated state
    CFHN.ASFUBMon(2) = CFHN.ASFUB(2);

    % New upper bound is non-monotonicity lower bound
    % plus P[Y = (1,1,...,1)]
    YPat = -1*ones(1, 1 + Settings.T);
    YPat(1:Settings.T) = 1;
    YPat(1:(Settings.T+1)) = 1;
    [~, YWide] = MatchPattern(YPat);
    Match = ismember(Data.Y, YWide, 'rows');
    CFHN.ASFUBMon(1) = CFHN.ASFLB(1) + sum(Match)/Settings.N;

    % Implied ATE bounds
    CFHN.ATELBMon = CFHN.ASFLBMon(2) - CFHN.ASFUBMon(1);
    CFHN.ATEUBMon = CFHN.ATEUB;
end
