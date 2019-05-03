%###############################################################################
% BuildConfidenceRegions
%###############################################################################
function CR = BuildConfidenceRegions(ampl, Settings, Data, Bounds)

% Initial brackets
WorstUBLeft = Bounds(:,1);
WorstLBRight = Bounds(:,2);

if (length(Settings.ParametersToTest) == 0)
    warning('Called BuildConfidenceRegions with ParametersToTest empty.');
end

for p = 1:1:length(Settings.ParametersToTest)
Settings.ActiveParam = Settings.ParametersToTest(p);
RejectList = zeros(length(Settings.LevelsCR), 0, length(Settings.Tests));
PointList = [];

for j = 1:1:length(Settings.Tests)
    for a = 1:1:length(Settings.LevelsCR)
        Settings.ActiveLevel = Settings.LevelsCR(a);
        Settings.ActiveTest = Settings.Tests(j);

        if (Settings.Noise >= 1)
            disp(repmat('-', 1, Settings.DisplaySepLen));
            disp(sprintf([  'Constructing %-.0f%% confidence region with '...
                            'test %s.'],...
                    (100 - 100*Settings.ActiveLevel),...
                    Settings.ActiveTest{:}));
        end

        if ((j == 1) & (a == 1))
            LeftBracket(1) = 0;
            LeftBracket(2) = WorstUBLeft(p,1);
            RightBracket(1) = WorstLBRight(p,1);
            RightBracket(2) = 1;
        else
            % Determine bracket based on saved rejection info thus far
            % Find simulations where the current test/size was rejected (I1)
            % and where it was not rejected (I0)
            I0 = find(RejectList(...
                        Index(Settings.ActiveLevel, Settings.LevelsCR),...
                        :, j) == 0);
            I1 = find(RejectList(...
                        Index(Settings.ActiveLevel, Settings.LevelsCR),...
                        :, j) == 1);

            % Set the upper bound on the left bracket to the smallest point that
            % has not been rejected by test j at level a unless this point is
            % larger than the initial bracketing point.
            LeftBracket(2) = min([PointList(I0); WorstUBLeft(p,1)]);
            % Set the lower bound on the left bracket to the largest point that
            % has been rejected by test j at level a, but which is smaller than
            % LeftBracket(2). If no such point, then set it to 0.;
            PLI1 = PointList(I1);
            LT = find(PLI1 <= LeftBracket(2));
            if ~isempty(LT)
                LeftBracket(1) = max(PLI1(LT));
            else
                LeftBracket(1) = 0;
            end

            % Do the symmetric thing for the right bracket
            RightBracket(1) = ...
                max([PointList(I0); WorstLBRight(p,1)]);
            GT = find(PLI1 >= RightBracket(1));
            if ~isempty(GT)
                RightBracket(2) = min(PLI1(GT));
            else
                RightBracket(2) = 1;
            end
        end

        if (Settings.Noise >= 1)
            disp(sprintf(['Parameter is ' Settings.ActiveParam{:} '.']));
            disp(sprintf('\tLeft endpoint', LeftBracket));
        end

        % Passing back settings so keep track of points that have been tested
        % already.
        [LeftEndpoint PointList RejectList] = BracketCREndpoint(ampl,...
            LeftBracket(2), LeftBracket(1),...
            Settings, Data, PointList, RejectList);

        if (Settings.Noise >= 1)
            disp(sprintf('\tRight endpoint', RightBracket));
        end
        [RightEndpoint PointList RejectList] = BracketCREndpoint(ampl,...
            RightBracket(1), RightBracket(2),...
            Settings, Data, PointList, RejectList);

        if (Settings.Noise >= 1)
            disp(sprintf(['Final %-.0f%% confidence region for %s '...
                           'is [%7.6f, %7.6f]'],...
                           (100 - 100*Settings.ActiveLevel),...
                           Settings.ActiveParam{:},...
                           LeftEndpoint, RightEndpoint));
        end

        CR(p,j,a,:) = [LeftEndpoint RightEndpoint];
    end
end
end
end

%###############################################################################
% BracketCREndpoints
%
% Determine either a lower or upper endpoint of a confidence region through
% a bracketing procedure.
% "In" is a point that is known to be within the confidence region.
% "Out" is a point that is known to be outside of the confidence region.
% Depending on whether In is smaller than Out or not the script determines where
% this is a left endpoint or right endpoint, then proceeds with bracketing.
%
% PointList and RejectList keep track of all points that have been tested
% so far and whether they were rejected or not at the levels and tests
% for which confidence interval are to be built.
%###############################################################################
function [Endpoint PointList RejectList] = ...
    BracketCREndpoint(ampl, In, Out, Settings, Data, PointList, RejectList)
%###############################################################################

    if (In < Out) % Right-hand side bracket
        RightBracket = 1;
        dirstr = 'right';
        LB = In;
        UB = Out;
    else % Left-hand bracket
        RightBracket = 0;
        dirstr = 'left';
        LB = Out;
        UB = In;
    end

    numfmt = '%8.6f';
    bracketfmt = ['[' numfmt ', ' numfmt ']'];
    collenwide = length(sprintf(bracketfmt, pi, pi)) + 4;
    colfmtwide = sprintf('%%-%ds', collenwide);
    collen = round(collenwide/2);
    colfmt = sprintf('%%-%ds', collen);
    rowfmt1 = ['\t' colfmtwide colfmt]; % First half of a row line
    rowfmt2 = [colfmt colfmt colfmt '\n']; % Second half of a row line
    if (Settings.Noise >= 1)
        fprintf([rowfmt1 rowfmt2], 'Bracket', 'Midpoint', 'TS', 'CV', 'Reject');
    end

    while (UB - LB) > Settings.BracketTol
        t = (UB + LB)/2; % Test the midpoint
        PointList = [PointList; t];

        if (Settings.Noise >= 1)
            fprintf(rowfmt1, sprintf(bracketfmt, LB, UB), sprintf(numfmt, t));
        end

        [TS PValue CV Reject] ...
            = TestListOfPoints(ampl, Settings, Data, t, Settings.LevelsCR);
        RejectList(:, size(RejectList,2) + 1, :) = Reject;

        CurReject = ...
            Reject( Index(Settings.ActiveLevel, Settings.LevelsCR),...
                    Index(Settings.ActiveTest, Settings.Tests));
        CurCV = CV( Index(Settings.ActiveLevel, Settings.LevelsCR),...
                     Index(Settings.ActiveTest, Settings.Tests));

        if (Settings.Noise >= 1)
            fprintf(rowfmt2, sprintf(numfmt, TS), sprintf(numfmt, CurCV),...
                     sprintf('%d', CurReject));
        end
        diary off; diary on; % Flush

        % Adjust bracket depending on whether this was a right or left bracket
        % and depending on whether there was a rejection at the midpoint
        if RightBracket
            if CurReject
                UB = t;
            else
                LB = t;
            end
        else
            if CurReject
                LB = t;
            else
                UB = t;
            end
        end
    end
    Endpoint = (LB + UB)/2;
end
