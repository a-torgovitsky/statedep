%###############################################################################
% TestListOfPoints
%
% Output:
%   TS: test statistics at each element in Points
%   PValue: smallest level at which one would reject for each element in Points
%
%   if Levels was passed then also return:
%   CV:     critical value for
%           each test in Settings.Tests,
%           each element in Points,
%           and at each level in Levels
%   Reject: binary rejection indicator to correspond with CV
%###############################################################################
function [TS PValue CV Reject] ...
    = TestListOfPoints(ampl, Settings, Data, Points, Levels)
%###############################################################################
    TestUniverse = {'CNS', 'SS'};
    if ~all(ismember(Settings.Tests, TestUniverse))
        error('Invalid list of tests.')
    end

    if isempty(Points)
        warning(['Points is empty. Skipping this routine.'],...
                 Settings.TestAListOfPoints);
    end

    Settings.InitialSeed = round(Settings.InitialSeed);
    if (Settings.InitialSeed <= 0)
        error('InitialSeed is not a positive integer.')
    end
    Settings.B = round(Settings.B);
    if (Settings.B <= 0)
        error('Number of bootstrap replications is not a positive integer.')
    end

    if ~isnumeric(Settings.SSExp) | (Settings.SSExp <= 0) | (Settings.SSExp > 1)
        error('SSExp must be a number between 0 and 1.')
    end

    TS = zeros(1, length(Points));
    TS = ComputeTestStatistics(ampl, Points, Settings);
    TS = TS(:)';

    % If TS was basically 0, then no point in doing the test
    % since you know you're not going to reject.
    % Indicate this by setting BSStat = +Inf, so critical values
    % will also be +Inf. Then skip testing these points.
    BSStat = zeros(Settings.B, length(Points), length(Settings.Tests));
    IdxContinue = find(TS > Settings.SkipTestingTol);
    IdxPass = find(TS <= Settings.SkipTestingTol);
    BSStat(:,IdxPass,:) = +Inf;
    PointsContinue = Points(IdxContinue);
    Settings.SavedTS = TS(IdxContinue); % This gets used in CNS

    if length(PointsContinue) > 0
        if ismember('SS', Settings.Tests)
            t = Index('SS', Settings.Tests);
            BSStat(:,IdxContinue,t) = SolveBootstrapProblems(ampl,...
                PointsContinue, 'SS', Settings, Data);
        end

        if ismember('CNS', Settings.Tests)
            t = Index('CNS', Settings.Tests);
            BSStat(:,IdxContinue,t) = SolveBootstrapProblems(ampl,...
                PointsContinue, 'CNS', Settings, Data);
        end
    end

    % Compute p-values, i.e. one minus the quantile of the largest CV for which
    % one would still get a rejection
    PValue = -1*ones(length(Points), length(Settings.Tests));
    for t = 1:1:length(Settings.Tests)
        for j = 1:1:length(Points)
            PValue(j,t) = 1 - mean(TS(j) > ...
                BSStat(:,j,t) + Settings.RejectTol);
        end
    end

    % If Levels was passed, then return also matrices of critical values
    % and rejection at the specified levels
    CV = [];
    Reject = [];
    if exist('Levels')
    if ~isempty(Levels)
        % Compute critical values (quantiles)
        %   note that Matlab's built-in does interpolation
        %   which is why I am writing my own code for this
        CV = zeros(length(Levels), length(Points), length(Settings.Tests));
        for t = 1:1:length(Settings.Tests)
            for p = 1:1:length(Points)
                CV(:,p,t) = ...
                    ComputeQuantile(BSStat(:,p,t), 1 - Levels);
            end
        end

        % Determine rejection
        Reject = zeros(size(CV));
        for t = 1:1:length(Settings.Tests)
            for j = 1:1:length(Points)
                for a = 1:1:length(Levels)
                    Reject(a,j,t) = (TS(j) > ...
                        CV(a,j,t) + Settings.RejectTol);
                    RejectCheck = ...
                        (PValue(j,t) <= Levels(a) + Settings.RejectTol);

                    if RejectCheck ~= Reject(a,j,t)
                        error('Something is wrong with CV or PValues.');
                    end
                end
                assert(issorted(Reject(:,j,t)));
            end
        end
    end
    end
end

%###############################################################################
% SolveBootstrapProblems
%
% Solve bootstrap problem for test "Type" at all points in the vector Points.
% Return:
%   A vector of B bootstrap statistics for each point in Points
%###############################################################################
function [BSStat] = SolveBootstrapProblems(ampl, Points, Type, Settings, Data)
    AcceptedTypes = {'SS', 'CNS'};
    assert(ismember(Type, AcceptedTypes));

    if strcmp(Type, 'SS')
        ChangeOptimizationProblem(ampl, Settings, 'Criterion');
        ResampleSize = round(Settings.N^Settings.SSExp);
        WithReplacement = 0;
        CriterionName = 'minCriterion';
        IDStrStub = 'SolveBootstrapProblems (SS)';
        FlagCNS = 0;
    else % CNS
        ChangeOptimizationProblem(ampl, Settings, 'CNS');
        ResampleSize = Settings.N;
        WithReplacement = 1;
        CriterionName = 'minCriterion_CNS';
        IDStrStub = 'SolveBootstrapProblems (CNS)';
        FlagCNS = 1;
    end

    B = Settings.B;
    BSStat = zeros(B, length(Points));

    % Save sample quantities that are used in the resampling procedures.
    % Note that the AMPL Q variable itself gets overwritten with bootstrap
    % draws which is why the need for Q_Sample.
    ampl.eval('let {y in YHAT} Q_Sample[y] := Q[y];');
    CriterionHat = ampl.getParameter('CriterionHat');

    for b = 1:1:B
        % Draw a bootstrap sample with replacement and apply to AMPL
        DataBS = ResampleData(Data, ResampleSize,...
            WithReplacement, b + Settings.InitialSeed);
        UpdateAMPLData(ampl, Settings, DataBS);

        for t = 1:1:length(Points)
            ampl.getParameter('Fix').setValues(Points(t));

            if FlagCNS
                CriterionHat.setValues(Settings.SavedTS(t));
            end

            if Settings.NoisyOptimization
                eval('ampl.solve');
            else
                evalc('ampl.solve');
            end
            SolveResult = ampl.getValue('solve_result');

            % Identifier for printing output
            IDStr = [IDStrStub ' '];
            if isfield(Settings, 'MCPoints')
                IDStr = [IDStr 'm = ' int2str(Settings.CurrentSim) ', '];
            end
            IDStr = [IDStr 'b = ' int2str(b), ' t = ' num2str(Points(t))];

            % If return code is not solved (or ``solved?'')
            % then keep increasing the tolerance
            % until we get one or we exceed some maximum tolerance.
            if (isempty(strfind(SolveResult, 'solved')))
                [BSStat(b,t) SolveResult] = ...
                    OptimizeWithHigherTolerance(...
                        ampl, CriterionName, IDStr, Settings);
                % Restore original tolerance
                SetTolerance(ampl, Settings.FeasTolDefault);
            else
                BSStat(b,t) = SafelyGetObjective(ampl, CriterionName);
            end

            ErrorCheckOptimization(ampl, IDStr, 1);
        end
    end

    % Restore the original data to AMPL
    UpdateAMPLData(ampl, Settings, Data);
end

%*******************************************************************************
% OptimizeWithHigherTolerance
%
% This is a routine called when the previous optimization failed.
% It adjusts the solve tolerance from its default level up to
% a maximum of Settings.FeasTolMax in multiplicative steps of
% factor Settings.FeasTolStepFactor.
%
% If Settings.FeasTolMax has been hit and there's still no good solve
% then throw an error and stop the program.
%*******************************************************************************
function [Solution SolveResult] =...
    OptimizeWithHigherTolerance(ampl, Criterion, IDStr, Settings)

    ToleranceCurrent = Settings.FeasTolDefault;
    SolveResult = 'infeasible';
    while (isempty(strfind(SolveResult, 'solved')))

        if (ToleranceCurrent > Settings.FeasTolMax)
            InfoStr = [IDStr '\n'...
                'Quitting in OptimizeWithHigherTolerance:\n' ...
                '\t solve_result = %s\n' ...
                '\t solve_result_num = %d.'
            ];
            disp(sprintf(InfoStr, SolveResult, SolveResultNum));
            display(...
                [Criterion ': '...
                 'Problem still not solved and tolerance has' ...
                 ' gone above the maximum allowed.']);
            display('Continuing for now...')
            break;
        end

        ToleranceCurrent = ...
            ToleranceCurrent*Settings.FeasTolStepFactor;
        SetTolerance(ampl, ToleranceCurrent);
        if Settings.NoisyOptimization
            eval('ampl.solve');
        else
            evalc('ampl.solve');
        end
        SolveResult = ampl.getValue('solve_result');
        SolveResultNum = ampl.getValue('solve_result_num');
    end
    Solution = ampl.getValue(Criterion);
    % Restore original tolerance
    SetTolerance(ampl, Settings.FeasTolDefault);
end
