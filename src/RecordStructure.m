%*******************************************************************************
% RecordStructure
%
% Save the Settings structure in a file.
%*******************************************************************************
function RecordStructure(Settings, Outfilename)
    fid = fopen(Outfilename, 'wt');
    PrintStructure(Settings, fid);
    fclose(fid);
end

%*******************************************************************************
% PrintStructure
%
% Recurse through a structure and print its contents to fid
%*******************************************************************************
function PrintStructure(s, fid)
    fields = repmat(fieldnames(s), numel(s), 1);
    values = struct2cell(s);

    for i = 1:1:length(fields)
        if isstruct(s.(fields{i}))
            fprintf(fid, '-----------------\n');
            fprintf(fid, '%s\n', fields{i});
            fprintf(fid, '-----------------\n');
            PrintStructure(s.(fields{i}), fid); % Recursion
        else
            if isnumeric(values{i})
                values{i} = num2str(values{i});
            end
            stringout = [fields{i} ': '];
            if iscell(values{i})
                for j = 1:1:length(values{i})
                    strvalues{j} = num2str(values{i}{j});
                end
                if length(values{i}) > 0
                    stringout = [stringout strjoin(strvalues)];
                end
            else
                stringout = [stringout values{i}];
            end

            fprintf(fid, '%s\n', stringout);
        end
    end
end

