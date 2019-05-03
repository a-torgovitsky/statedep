%*******************************************************************************
% RecordBounds
%*******************************************************************************
function RecordBounds(Parameters, Bounds, Outfilename)
    fid = fopen(Outfilename, 'w');
    for p = 1:1:length(Parameters)
        fprintf(fid, '%-10s %8.5f %8.5f\n',...
                Parameters{p}, Bounds(p,1), Bounds(p,2));
    end
    fclose(fid);
end

