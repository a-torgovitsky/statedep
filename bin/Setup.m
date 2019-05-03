%*******************************************************************************
% Setup actions that are common to all run files
%*******************************************************************************
function [CleanUpPath CleanUpLog SaveDir] = Setup(ScriptTag, SaveDir)

    OriginalPath = pwd;
    CleanUpPath = onCleanup(@()cd(OriginalPath));

    evalc(['run ../cfg/Config']); % Load user-defined paths, etc.

    if ~exist(ResultsPath, 'dir')
        error('ResultsPath %s does not exist.', ResultsPath);
    end

    % Create a directory path name
    FlagNoSaveDir = 1;
    if exist('SaveDir')
        if ~isempty(SaveDir)
            SaveDir = [ResultsPath SaveDir];
            FlagNoSaveDir = 0;
        end
    end
    if FlagNoSaveDir
        SaveDir = MakeTargetDirName(ResultsPath, ScriptTag, 1);
    end

    % Create a new directory, copy everything there and cd to it
    SaveDir = CreateSaveDirectory(SaveDir);

    % Create a results subdirectory and move there
    if ~exist('./results', 'dir')
        mkdir('results');
    end
    cd('results');

    % Add paths to make all scripts available
    addpath(genpath(SaveDir));

    % Setup AMPL
    if ~exist(AMPLAPISetupPath, 'file')
        error('Path %s for AMPL API SetUp.m file does not exist.',...
                AMPLAPISetupPath);
    else
        evalc(['run ' AMPLAPISetupPath]);
    end

    % Start log and create destructor
    diary off;
    diary('Log.txt');
    CleanUpLog = onCleanup(@()diary('off'));
end

%*******************************************************************************
% MakeTargetDirName
%*******************************************************************************
function [TargetDir, tag] = MakeTargetDirName(BASEDIR, STUB, Interactive, tag)
    if ~exist('Interactive')
        Interactive = 1;
    end

    if Interactive
        tag = input('Do you want to add a tag to the directory?\n', 's');
    else
        if ~exist('tag')
            tag = [];
        end
    end
    if ~isempty(tag)
        if ~strcmp(tag(1), '-')
            tag = ['-' tag];
        end
    end
    if ~strcmp(BASEDIR(length(BASEDIR)), '/')
        BASEDIR = [BASEDIR '/'];
    end
    TargetDir = ...
        [BASEDIR STUB tag datestr(now, '-mmddyy-HHMMSS')];
end

%*******************************************************************************
% CreateSaveDirectory
%
% Create a directory with a current copy of the code (so that further changes
% don't affect this simulation) and where we can store results.
% Create a subdirectory in that directory to store results.
% And change directories to that subdirectory.
%*******************************************************************************
function [SaveDir] = CreateSaveDirectory(SaveDir)
    if isempty(SaveDir)
        error('SaveDir is an empty string!');
    end
    if ~exist(SaveDir, 'dir')
        disp(sprintf('Saving in new directory %s.', SaveDir));
        success = mkdir(SaveDir);
        if ~success
            error('Something is wrong; error creating directory.')
        end
        copyfile('../*', SaveDir);
    elseif (exist(SaveDir, 'dir') == 7)
        disp(sprintf('Saving in existing directory %s.', SaveDir));
        disp(sprintf('\t Directory already exists; not copying code.'));
    end
    cd(SaveDir);
end
