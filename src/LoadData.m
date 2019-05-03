%###############################################################################
% LoadData
%
% Load data into Matlab from long-form tab separated file.
% The columns in this file should be:
%   ID t Y X1 X2 X3 ...
% times should be numbered starting from 0
%
% Data should be sorted by ID then t
%###############################################################################
function [Settings Data] = LoadData(Settings)
    if exist(Settings.DataPath, 'file') ~= 2
        error('Could not find data file %s.', Settings.DataPath);
    end
    if (Settings.Noise >= 1)
        disp(sprintf('Loading long-form tab-delimited data from %s.',...
            Settings.DataPath));
    end
    DataIn = csvimport(Settings.DataPath, 'delimiter', '\t');
    ColHeaders = DataIn(1,:);
    DataIn = DataIn(2:end,:); % Toss out the header
    Settings.N = length(unique(cell2mat(DataIn(:,1))));
    if ~(Settings.N > 0)
        error('Something is wrong; sample size is %d.', Settings.N);
    end

    MaxT = max(cell2mat(DataIn(:,2)));
    if isempty(Settings.T)
        Settings.T = MaxT;
    end
    if (Settings.T > MaxT)
        error('Settings.T is larger than the largest T in data of %d', MaxT);
    end
    if (Settings.T < 1)
        error('Settings.T is smaller than 1.');
    end
    if ~(mod(Settings.T, 1) == 0)
        error('Settings.T is not an integer.');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % LOAD Y DATA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Data.Y = nan(Settings.N, Settings.T+1);
    IdxT = find(cell2mat(DataIn(:,2)) <= Settings.T);
    Data.Y = reshape(cell2mat(DataIn(IdxT, 3)), Settings.T + 1, Settings.N)';
    if ~(all(ismember(Data.Y(:), [0 1])))
        error('Something is wrong: Data.Y is not all in {0, 1}.');
    end

    %###########################################################################
    % HARDCODED
    %###########################################################################
    Data.Age = reshape(cell2mat(DataIn(IdxT, 4)), Settings.T + 1, Settings.N)';
end
