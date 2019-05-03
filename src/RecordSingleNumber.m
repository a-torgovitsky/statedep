%*******************************************************************************
% RecordSingleNumber
%*******************************************************************************
function RecordSingleNumber(Number, Label, Outfilename)
    fid = fopen(Outfilename, 'wt');
    fprintf(fid, [Label ' %10.8f'], Number);
    fclose(fid);
end
