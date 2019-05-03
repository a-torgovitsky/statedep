function OriginalPath = CreateResultsDir(SaveDir, ResultsSubdir)
    if ~exist('SaveDir', 'var')
        error('Must pass in SaveDir.')
    else
        if isempty(SaveDir)
            error('SaveDir should not be empty.')
        end
    end

    if ~exist('ResultsSubdir', 'var')
        ResultsSubdir = '';
    end
    if isempty(ResultsSubdir)
        ResultsSubdir = CreateNewNumberedDirName('.', '');
    end

    FullResultsDirPath = fullfile(SaveDir, 'results', ResultsSubdir);
    if ~exist(FullResultsDirPath, 'dir')
        mkdir(FullResultsDirPath)
    else
        error('Danger -- final results directory already exists. Quitting.')
    end
    OriginalPath = cd(FullResultsDirPath);
end

function DirName = CreateNewNumberedDirName(DirBase, DirStub)
    DirName = '';
    DirNumber = [];
    for i = 0:1:999
        DirName = fullfile(DirBase, sprintf([DirStub '%03d'], i + 1));
        if ~exist(fullfile(DirName), 'dir')
            DirNumber = i;
            break;
        end
    end
end
