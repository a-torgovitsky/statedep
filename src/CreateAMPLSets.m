%###############################################################################
% CreateAMPLSets
%
% Create set definitions in AMPL.
%###############################################################################
function [] = CreateAMPLSets(ampl, Settings, Data)

TicTotal = tic;

if (Settings.Noise >= 1)
    disp(repmat('=', 1, Settings.DisplaySepLen));
    disp('Started creating set definitions...')
end

T = Settings.T;
ampl.getParameter('T').setValues(T);

%###############################################################################
% All sequences of length 1..T
% Could be done in AMPL, but it is helpful later on here as well
%###############################################################################
for t = 1:1:(T+1)
    AllYSeqsWide{t} = AllBinaryArray(t);
    AllYSeqsInt{t} = WideToBinary(AllYSeqsWide{t});
    FillAMPLSet(ampl, 'YSEQS', [t], AllYSeqsInt{t}(:)');
end

%###########################################################################
% Fill in YHAT
% I found that not doing this first could lead to errors with the
% Matlab AMPL API since other sets are indexed based on it.
%###########################################################################
TicYHat = tic;
YHat = unique(Data.Y, 'rows');
YHatInt = WideToBinary(YHat);

FillAMPLSet(ampl, 'YHAT', [], YHatInt(:)');
clear YHatInt;
TimeYHat = toc(TicYHat);

%##########################################################################
% U_OBSEQ[y] and U_HAT
%
% Fix an observed sequence y.
% The pattern for U sequences that need to be summed over
% to match observational equivalence for y is:
%   u_{t_{0} = y_{0}
%   if y_{0} = 0 then u_{1}(0) = y_{1}
%                else u_{1}(1) = y_{1}
%   if y_{1} = 0 then u_{2}(0) = y_{2}
%                else u_{2}(1) = y_{2}
%   and so on
%
% NOTE: The convention for u sequences in the code differs
% from that in the paper.
% In the code it is:
%   (u(0), u_{1}(0),...,u_{T}(0),
%          u_{1}(1),...,u_{T}(1))
%
% UHAT is then the union of U_OBSEQ[y] for all y
%##########################################################################
TicUHat = tic;
for y = 1:1:size(YHat, 1)
    YSeq = YHat(y,:);
    UPatObsEq = BuildUPatternToMatchY(YSeq, T);

    M = MatchPattern(UPatObsEq);
    UObsEqInt(y,:) = M(:)';
    YSeqInt = WideToBinary(YSeq);
    UObsEqIdx(y,:) = [YSeqInt];
end
UHatInt = sort(UObsEqInt(:)); % Keep this around--needed later
FillAMPLSet(ampl, 'UHAT', [], UHatInt(:)');
FillAMPLSet(ampl, 'U_OEQ', UObsEqIdx, UObsEqInt);
clear UObsEqInt UObsEqIdx;
TimeUHat = toc(TicUHat);

%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
% PARAMETERS
%###############################################################################
%###############################################################################
%###############################################################################
%###############################################################################
TicParams = tic;
if isempty(Settings.OptPeriod)
    Settings.OptPeriod = T + 1;
else
    if ~isnumeric(Settings.OptPeriod) ...
        | (Settings.OptPeriod < 1) | (Settings.OptPeriod > T + 1)
        error('OptPeriod must be a number between 1 and (T+1).')
    end
end

%###################################################################
% U_PSD
% Pattern for PSD is:
%   u_{t}(0) = 0
%   u_{t}(1) = 1
%###################################################################
for t = 1:1:T
    UPat = -1*ones(1, 1 + 2*T);
    UPat(1 + t) = 0;
    UPat(1 + T + t) = 1;
    UPSDInt(t,:) = MatchPattern(UPat);
    UList = sort(UPSDInt(t,:));
    FillAMPLSet(ampl, 'U_PSD', [t], UList(:)');
end
clear UPSDInt;

%*******************************************************************************
% U_AE1
% Pattern is:
%   u_{t}(1) = 1
% so a special case of above
%*******************************************************************************
for t = 1:1:T
    UPat = -1*ones(1, 1 + 2*T);
    UPat(1 + T + t) = 1;
    UAE1Int(t,:) = MatchPattern(UPat);
    UList = sort(UAE1Int(t,:));
    FillAMPLSet(ampl, 'U_AE1', [t], UList(:)');
end
clear UAE1Int;

%###################################################################
% U_NSD
% Pattern for NSD is:
%   u_{t}(0) = 1
%   u_{t}(1) = 0
%###################################################################
for t = 1:1:T
    UPat = -1*ones(1, 1 + 2*T);
    UPat(1 + t) = 1;
    UPat(1 + T + t) = 0;
    UNSDInt(t,:) = MatchPattern(UPat);
    UList = sort(UNSDInt(t,:));
    FillAMPLSet(ampl, 'U_NSD', [t], UList(:)');
end
clear UNSDInt;

%*******************************************************************************
% U_AE0
% Pattern is:
%   u_{t}(0) = 1
% so a special case of above
%*******************************************************************************
for t = 1:1:T
    UPat = -1*ones(1, 1 + 2*T);
    UPat(1 + t) = 1;
    UAE0Int(t,:) = MatchPattern(UPat);
    UList = sort(UAE0Int(t,:));
    FillAMPLSet(ampl, 'U_AE0', [t], UList(:)');
end
clear UAE0Int;

%###############################################################################
% PSD_G0, PSD_G00
%
% PSD_G0
% Loop over t
%   For each t, find all Y, such that Y_{t} = 0.
%   Pass this set to AMPL as Y_G0[t].
%   Loop over y in Y_G0[t]
%       For each y in Y_G0[t], find all U = u that could generate y
%       Pass this set to AMPL as U_PSD_G0_DEN
%       Then further restrict these U to those with U_{t}(0) = 0, U_{t}(1) = 1
%       Pass this set to AMPL as U_PSD_G0_NUM
%
% PSD_G00
%   Start with all Y such that Y_{t} = 0 from above.
%   Restrict this set to all Y such that also Y_{t-1} = 0.
%   Then repeat the same steps as above with this set.
%###############################################################################
for t = 1:1:T
    % Given Y_{t} = 0
    YPat = -1*ones(1, 1 + T);
    YPat(t+1) = 0;
    [YG0Int YG0Wide] = MatchPattern(YPat);
    FillAMPLSet(ampl, 'Y_G0', [t], YG0Int(:)');
    UG0Num = FindConditionalPSDSequences(YG0Wide, T, t);
    FillAMPLSet(ampl, 'U_PSD_G0_NUM', [t], UG0Num(:)');

    % Given Y_{t} = 0, Y_{t-1} = 0
    I = find(YG0Wide(:,1+(t-1)) == 0);
    YG00Int = YG0Int(I);
    YG00Wide = YG0Wide(I,:);
    FillAMPLSet(ampl, 'Y_G00', [t], YG00Int(:)');
    UG00Num = FindConditionalPSDSequences(YG00Wide, T, t);
    FillAMPLSet(ampl, 'U_PSD_G00_NUM', [t], UG00Num(:)');
end

%###############################################################################
% PSD_G1, PSD_G11
%
% Almost identical to the above
%###############################################################################
for t = 1:1:T
    % Given Y_{t} = 1
    YPat = -1*ones(1, 1 + T);
    YPat(t+1) = 1;
    [YG1Int YG1Wide] = MatchPattern(YPat);
    FillAMPLSet(ampl, 'Y_G1', [t], YG1Int(:)');
    UG1Num = FindConditionalPSDSequences(YG1Wide, T, t);
    FillAMPLSet(ampl, 'U_PSD_G1_NUM', [t], UG1Num(:)');

    % Given Y_{t} = 1, Y_{t-1} = 1
    I = find(YG1Wide(:,1+(t-1)) == 1);
    YG11Int = YG1Int(I);
    YG11Wide = YG1Wide(I,:);
    FillAMPLSet(ampl, 'Y_G11', [t], YG11Int(:)');
    UG11Num = FindConditionalPSDSequences(YG11Wide, T, t);
    FillAMPLSet(ampl, 'U_PSD_G11_NUM', [t], UG11Num(:)');
end

TimeParams = toc(TicParams);

%###############################################################################
% ST
% Sequences required to compute ST
%
% U_ST_EQUATE is all binary sequences of length 2(DimST + 1)
% This is the set of all sequences that
%  (U_{t}^{m}(0), U_{t}^{m}(1)) can assume
% where m = DimST
%
% For each u in U_ST_EQUATE,
%   U_ST[t,u] is the set of all sequences such that U_{t}^{m} = u
% This set gets intersected with UHAT then passed to AMPL.
%###############################################################################
TicST = tic;
if Settings.Assumption_ST
    DimST = Settings.Assumption_DimST;
    MaxDimST = max(0, T-2);
    if ~((DimST >= 0) & (DimST <= MaxDimST))
        error('DimST is %d, but should be between 0 and %d.',...
            DimST, MaxDimST);
    end

    ampl.getParameter('DIMST').setValues(DimST);
    USTEquateWide = AllBinaryArray(2*(DimST + 1));
    USTEquateInt = WideToBinary(USTEquateWide);
    USTEquateInt = sort(USTEquateInt);
    FillAMPLSet(ampl, 'U_ST_EQUATE', [], USTEquateInt(:)');

    for u = 1:1:size(USTEquateWide, 1)
        USeq = USTEquateWide(u,:);
        USeqInt = WideToBinary(USeq);

        for t = 1:1:(T - DimST)
            UPat = -1*ones(1, 1 + 2*T);

            % First half of USeq gets assigned to (U_{t}(0),...,U_{t+m}(0))
            % Second half gets assigned to (U_{t}(1),...,U_{t+m}(1))
            UPat((1 + t):(1 + t + DimST)) ...
                = USeq(1:(DimST + 1));
            UPat((1 + T + t):(1 + T + t + DimST))...
                = USeq((DimST + 2):end);

            USTInt = MatchPattern(UPat);
            USTInt = sort(USTInt);

            FillAMPLSet(ampl, 'U_ST', [t USeqInt], USTInt(:)');
        end
    end
    clear USTEquateWide USTEquateInt;

    Idx = [];
    Val = [];
    for t = 1:1:(T - DimST)
        for tt = 1:1:(T-DimST)
            Idx = [Idx; [t tt]];
            if  (tt == t)
                Val = [Val; 0];
            else
                Val = [Val; Settings.Assumption_SigmaST];
            end
        end
    end
    vSIGMAST = ampl.getParameter('SIGMAST');
    vSIGMAST.setValues(Idx, Val);
end

TimeST = toc(TicST);

%###############################################################################
% DSC
% Sequences required to compute DSC
%
% U_DSC[t,s,d] is the set of all binary sequences such that
%   P[U_{t}(d) = 1, U_{s}(d) = 1]
% for t ~= s, and d in 0..1
%###############################################################################
TicDSC = tic;
if Settings.Assumption_DSC
    USeqIdx = [];
    USeqInt = [];
    for t = 1:1:T
        for s = 1:1:T
            if (s == t)
                continue;
            end

            for d = 0:1:1
                UPat = -1*ones(1, 1 + 2*T);
                UPat(1 + d*T + t) = 1;
                UPat(1 + d*T + s) = 1;

                UDSCInt = MatchPattern(UPat);
                UDSCInt = sort(UDSCInt);

                USeqIdx = [USeqIdx; [t s d]];
                USeqInt = [USeqInt; UDSCInt(:)'];
            end
        end
    end
    FillAMPLSet(ampl, 'U_DSC', USeqIdx, USeqInt);
    clear UPat UDSCInt UIdx USeqInt count;
end
TimeDSC = toc(TicDSC);

%###############################################################################
% MTS
% Sequences required to compute MTS
%
% Write the conditional probability as Numerator/Denominator. Then:
%   - Find all possible denominator sequences in the condition
%   - Find all sequences of Y that need to be summed over to get the
%     probability of any given denominator sequence.
%   - Find all sequences of U that need to be summed over to get the
%     probability of any given numerator sequence.
%###############################################################################
TicMTS = tic;
DimMTS = Settings.Assumption_DimMTS;
if (rem(DimMTS, 1) ~= 0)
    error('Assumption_DimMTS should be an integer, but is %7.5f.', DimMTS);
end
if ~ismember(DimMTS, [2:1:T])
    error([ 'Assumption_DimMTS is %d, but should be between 2'...
            ' and T = %d.'], DimMTS, T);
end

%###############################################################################
% A ``conditioning sequence'' refers to
%    (Y_{t-2},...,Y_{t-q})
% but does not include Y_{t-1}
% for q = 2,...,DimMTS
%
% The maximum length of a conditioning sequence that starts at time t
% is (DimMTS - 1) if t - DimMTS >= 0
% and is (t-1) if t - DimMTS < 0
%###############################################################################
MaxLen = nan(T,1);
for t = 2:1:T
    if (t - DimMTS) >= 0
        MaxLen(t,1) = DimMTS - 1;
    else
        MaxLen(t,1) = t - 1;
    end
end
ampl.getParameter('Y_MTS_LENMAX').setValues((2:1:t)', MaxLen(2:1:t));

for t = 2:1:T
for q = 1:1:MaxLen(t)
    YMTSSumInt = [];
    YMTSSumWide = [];
    YMTSSumIdx = [];

    for ytm1 = 0:1:1
        % For each conditioning sequence y' = (ytm1, y)
        % find all full sequences yy of length T
        % such that Pr[y'] = sum {yy} Pr[yy]
        % This set gets summed over for the denominator term
        for y = 1:1:size(AllYSeqsWide{q}, 1)
            YPat = -1*ones(1, 1 + T);

            % Note that position t of YPat corresponds to Y_{t-1}
            % since MATLAB only does base-1 numbering
            YPat(t) = ytm1;
            YPat((t-q):(t-1)) = AllYSeqsWide{q}(y,:);
            [YMTSSumInt(y,:) YMTSSumWide{y}] = MatchPattern(YPat);
            YMTSSumIdx(y,:) = [t ytm1 q AllYSeqsInt{q}(y)];
        end
        FillAMPLSet(ampl, 'Y_MTS_DENOM_SUM', YMTSSumIdx, YMTSSumInt);

        % Now for each d, and each conditioning sequence y' = (ytm1,y),
        % find the set of u to sum over in the numerator.
        % The way I do this is to take the sequences yy of length T
        % that were determined above to be the ones that would get summed
        % over to get conditioning sequence y'.
        % Truncate these sequences yy at time (t-1).
        % Keep only the unique sequences.
        % Then use each one to create a pattern for U sequences.
        % Also, add to this pattern that U_{t}(d) = 1
        % Get all such U's
        % Repeat for every unique yy
        for y = 1:1:size(AllYSeqsWide{q}, 1)
            YMTSSumTrunc = YMTSSumWide{y}(:,1:t); % time 0 to time t-1
            YMTSSumTrunc = unique(YMTSSumTrunc, 'rows');

            MTSNumerInt = cell(2,1);
            for yy = 1:1:size(YMTSSumTrunc, 1)
                YYSeq = YMTSSumTrunc(yy,:);

                UPatBase = BuildUPatternToMatchY(YYSeq, T);

                for d = 0:1:1
                    UPat = UPatBase;
                    UPat(1 + d*T + t) = 1;

                    UList = MatchPattern(UPat);

                    MTSNumerInt{d+1} = [MTSNumerInt{d+1}; UList(:)];
                end
            end
            for d = 0:1:1
                FillAMPLSet(ampl, 'U_MTS_NUMER',...
                    [t d ytm1 q AllYSeqsInt{q}(y)], MTSNumerInt{d+1}');
            end
        end
    end
end
end
TimeMTS = toc(TicMTS);

%###############################################################################
% TIV
%###############################################################################
TicTIV = tic;
if Settings.Assumption_TIV
    for t = 1:1:T
        for r = 0:1:(t-1)
            for y = 1:1:size(AllYSeqsWide{r+1}, 1)
                YSeq = AllYSeqsWide{r+1}(y,:);
                UTIVPatBase = BuildUPatternToMatchY(YSeq, T);

                for u0 = 0:1:1
                    for u1 = 0:1:1
                        UTIVPat = UTIVPatBase;
                        UTIVPat(1 + t) = u0;
                        UTIVPat(1 + T + t) = u1;
                        UTIVList = MatchPattern(UTIVPat);
                        UTIVIdx = [t r u0 u1 AllYSeqsInt{r+1}(y)];
                        FillAMPLSet(ampl, 'U_TIV', UTIVIdx, UTIVList(:)');
                    end
                end
            end
        end
    end
end
TimeTIV = toc(TicTIV);

if (Settings.Noise >= 1)
    disp(sprintf('Finished creating set definitions in %5.3f seconds:',...
        toc(TicTotal)));
    fmt = '\t%10s in %5.3f seconds.';
    disp(sprintf(fmt, 'YHat', TimeYHat));
    disp(sprintf(fmt, 'UHat', TimeUHat));
    disp(sprintf(fmt, 'Parameters', TimeParams));
    disp(sprintf(fmt, 'ST', TimeST));
    disp(sprintf(fmt, 'DSC', TimeDSC));
    disp(sprintf(fmt, 'MTS', TimeMTS));
    disp(sprintf(fmt, 'TIV', TimeTIV));
    disp(repmat('=', 1, Settings.DisplaySepLen));
end

end
%###############################################################################
%###############################################################################
%###############################################################################

%###############################################################################
% BuildUPatternToMatchY
%
% Input is a Y sequences of some length and model length T which determines
% the length of UPat.
% The Y is assumed to start from Y_{0}, so (Y_{0}, Y_{1},...,Y_{s})
% for some s.
%
% Create the pattern of all sequences U that could generate this Y.
%   The first component of U is Y_{0}, so that must match directly.
%   Then proceed iteratively:
%       if Y_{0} = 0, then U_{1}(0) = Y_{1}
%       else U_{1}(1) = Y_{1}
%       if Y_{1} = 0, then U_{2}(0) = Y_{2}
%       else U_{2}(1) = Y_{2}
%       and so on
%###############################################################################
function [UPat] = BuildUPatternToMatchY(YSeq, T)
    UPat = -1*ones(1, 1 + 2*T);
    assert(length(YSeq) <= (1 + T));

    UPat(1) = YSeq(1);
    for t = 2:1:length(YSeq);
        if (YSeq(t-1) == 0)
            UPat(1 + (t-1)) = YSeq(t);
        else
            UPat(1 + T + (t-1)) = YSeq(t);
        end
    end
end

%###############################################################################
% FindConditionalPSDSequences
%
% YCond is a matrix, each row is a Y sequence.
% Loop over each Y sequence.
%   Construct all U sequences that could have generated that Y sequence.
%   Then limit these U to those that also have U_{t}(0) = 0, U_{t}(1) = 1.
%   This becomes UNum (numerator).
%###############################################################################
function UNum = FindConditionalPSDSequences(YCond, T, t)
    UNum = [];
    for y = 1:1:size(YCond, 1)
        YSeq = YCond(y,:);
        UPat = BuildUPatternToMatchY(YSeq, T);
        [UInt UWide] = MatchPattern(UPat);

        I = find((UWide(:,1+t) == 0).*(UWide(:,1+T+t) == 1));
        UInt = UInt(I);
        UNum = [UNum; UInt(:)];
    end
end

%###############################################################################
% FillAMPLSet
%
% SetName is the name (a string) in AMPL used to define the set.
% Each row of Idx is an index of the set. (Pass Idx = [] if no indexing.)
% Each row of Val are the values to be set to this index of the set.
%
% Note that Val(i,:)' needs to be a column or AMPL will say "too many members"
%###############################################################################
function FillAMPLSet(ampl, SetName, Idx, Val)
    a = ampl.getSet(SetName);
    if isempty(Idx)
        assert(min(size(Val,1)) == 1);
        a.setValues(Val(:));
    else
        for i = 1:1:size(Idx, 1)
            Idx(i,:);
            s = a.get(Idx(i,:));
            s.setValues(Val(i,:)');
        end
    end
end
