%*******************************************************************************
% MonteCarlo
%*******************************************************************************
function MonteCarlo(DPOSettings, MCSettingsIn)
    % Default MC Settings
    MCSettings.M = 500;
    MCSettings.NMultiplier = 1;
    MCSettings.InitialSeed = 3131;
    MCSettings.ProgressFrequency = 10;
    MCSettings.PrintCols = 8;

    % Replace with user input
    if exist('MCSettingsIn')
        MCSettings = UpdateStruct(MCSettings, MCSettingsIn, 1);
    end

    % Fill in any defaults for DPOSettings since some may be used below
    DPOSettings.GetDefaultSettings = 1;
    [~, DPOSettings] = DPO(DPOSettings);

    % Some error checking for running a MC
    if DPOSettings.TestAListOfPoints & (length(DPOSettings.Parameters) > 1)
        warning(['DPOSettings.Parameters has more than 1 element.'...
                 ' Only the first parameter will be recorded.']);
    end
    DPOSettings.ParametersToTest = {DPOSettings.Parameters{1}};
    DPOSettings.BuildConfidenceRegions = 0;
    DPOSettings.RunMisspecificationTest = 0;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run DPO once with original (DGP) data, and record true bounds
    %   Add the endpoints to the points to test
    %   Then also record these points
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    TestOptionSave = DPOSettings.TestAListOfPoints;
    DPOSettings.TestAListOfPoints = 0;
    DPOSettings.Noise = 1;
    disp('Running with DGP data to find bounds...');
    [DGPResults DGPSettings DGPData] = DPO(DPOSettings);

    if (DGPResults.MinCriterion > 0)
        warning('DGP is infeasible with MinCriterion = %10.7f.',...
                DGPResults.MinCriterion);
    end

    RecordBounds(DGPSettings.Parameters, DGPResults.Bounds, 'TrueBounds.out');

    DPOSettings.TestAListOfPoints = TestOptionSave;
    if DPOSettings.TestAListOfPoints
        DPOSettings.PointsToTest{1} ...
            = sort([DPOSettings.PointsToTest{1} DGPResults.Bounds(1,:,1)]);
        fid = fopen('TestPoints.out', 'w');
        fprintf(fid, '%6.5f ', DPOSettings.PointsToTest{1}(:)');
        fprintf(fid, '\n');
        fclose(fid);
    else
        DPOSettings.PointsToTest{1} = sort([DGPResults.Bounds(1,:,1)]);
        DPOSettings.Tests = {};
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Prepare files for MC output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FileBounds{1} = fopen('EstimatedLB.out','w');
    FileBounds{2} = fopen('EstimatedUB.out','w');

    header = '';
    SFmt = sprintf('%%%ds', MCSettings.PrintCols);
    for j = 1:1:length(DPOSettings.Parameters)
        header = [header sprintf(SFmt, DPOSettings.Parameters{j})];
        if j < length(DPOSettings.Parameters)
            header = [header ','];
        else
            header = [header '\n'];
        end
    end
    fprintf(FileBounds{1}, header);
    fprintf(FileBounds{2}, header);

    if DPOSettings.TestAListOfPoints
        FilePValue = fopen('PValue.out', 'w');
        for j = 1:1:length(DPOSettings.Tests)
            for a = 1:1:length(DPOSettings.LevelsTestList)
                FileCV(a,j) = ...
                    fopen(['CriticalValue_A'...
                            int2str(DPOSettings.LevelsTestList(a)*100)...
                            '_' DPOSettings.Tests{j} '.out'], 'w');
                FileRej(a,j) = ...
                    fopen(['Reject_A' ...
                            int2str(DPOSettings.LevelsTestList(a)*100) ...
                            '_' DPOSettings.Tests{j} '.out'], 'w');
                FileRejProb(a,j) = ...
                    fopen(['RejProb_A' ...
                            int2str(DPOSettings.LevelsTestList(a)*100) ...
                            '_' DPOSettings.Tests{j} '.out'], 'w');
            end
        end
    end
    FileMinCriterion = fopen('MinCriterion.out', 'w');
    FileTS = fopen('TestStatistic.out', 'w');
    FileTimes = fopen('Times.out', 'w');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Define some variables for formatting output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    PFmt = sprintf('%%%d.%df', MCSettings.PrintCols, MCSettings.PrintCols-3);
    PFmtPtsVec = '';
    PFmtPtsVecInt = '';
    for p = 1:1:length(DPOSettings.PointsToTest{1})
        PFmtPtsVec = [PFmtPtsVec ' ' PFmt];
        PFmtPtsVecInt = [PFmtPtsVecInt ' %d'];
    end
    PFmtPtsVec = [PFmtPtsVec '\n'];
    PFmtPtsVecInt = [PFmtPtsVecInt '\n'];

    PFmtBoundsVec = '';
    for p = 1:1:length(DPOSettings.Parameters)
        PFmtBoundsVec = [PFmtBoundsVec PFmt];
        if p < length(DPOSettings.Parameters)
            PFmtBoundsVec = [PFmtBoundsVec ','];
        else
            PFmtBoundsVec = [PFmtBoundsVec '\n'];
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run MC
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DPOSettings.Noise = 0;
    N = round(DGPSettings.N*MCSettings.NMultiplier); % Sample size for MC
    disp('Beginning Monte Carlo simulation.')
    for m = 1:1:MCSettings.M
        if (mod(m,MCSettings.ProgressFrequency) == 0)
            disp(sprintf('Starting replication %d.', m));
        end
        diary off; diary on; % Flush output
        tic; % Start timing

        % Redraw data
        Seed = MCSettings.InitialSeed + m;
        Data = ResampleData(DGPData, N, 1, Seed);

        % Run DPO
        [Results(m), ~, ~] = DPO(DPOSettings, Data);

        % Record
        fprintf(FileBounds{1}, [PFmtBoundsVec], Results(m).Bounds(:,1));
        fprintf(FileBounds{2}, [PFmtBoundsVec], Results(m).Bounds(:,2));
        if DPOSettings.TestAListOfPoints
            fprintf(FilePValue, PFmtPtsVec, Results(m).PValue{1}(:));

            for j = 1:1:length(DPOSettings.Tests)
                for a = 1:1:length(DPOSettings.LevelsTestList)
                    fprintf(FileCV(a,j),...
                            PFmtPtsVec,...
                            Results(m).CV{1}(a,:,j));
                    fprintf(FileRej(a,j),...
                            PFmtPtsVecInt,...
                            Results(m).Reject{1}(a,:,j));
                    RejCount = zeros(size(Results(1).Reject{1}(a,:,j)));
                    for mm = 1:1:m
                        RejCount =    RejCount ...
                                    + Results(mm).Reject{1}(a,:,j);
                    end
                    RejProb = RejCount/m;
                    fprintf(FileRejProb(a,j), PFmtPtsVec, RejProb);
                end
            end
        end
        fprintf(FileMinCriterion, [PFmt '\n'], Results(m).MinCriterion);
        fprintf(FileTS, PFmtPtsVec, Results(m).TS{1}(:));
        fprintf(FileTimes, [PFmt '\n'], toc/60);
    end
    fclose('all');
    disp(repmat('=', 1, DPOSettings.DisplaySepLen));
    disp(repmat('=', 1, DPOSettings.DisplaySepLen));
    disp(repmat('=', 1, DPOSettings.DisplaySepLen));
    disp('All done with Monte Carlo.')
    disp(repmat('=', 1, DPOSettings.DisplaySepLen));
    disp(repmat('=', 1, DPOSettings.DisplaySepLen));
    disp(repmat('=', 1, DPOSettings.DisplaySepLen));
end
