%*******************************************************************************
% RecordAssumptions
%*******************************************************************************
function AssumptionString = RecordAssumptions(AssumptionStr, Outfilename)
    fid = fopen(Outfilename, 'w');

    for j = 1:1:length(AssumptionStr)
        fprintf(fid, '%s\n', AssumptionStr{j});
    end
    fclose(fid);
end
