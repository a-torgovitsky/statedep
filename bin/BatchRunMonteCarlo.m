function [] = BatchRunMonteCarlo(SaveDir, SimNumber)

    errstr = 'Need to pass nonempty SaveDir for this routine.';
    if ~exist('SaveDir', 'var')
        error(errstr);
    else
        if isempty(SaveDir)
            error(errstr);
        end
    end
    errstr = 'Need to pass nonempty SimNumber (!= 2) for this routine.';
    if ~exist('SimNumber', 'var')
        error(errstr);
    else
        if isempty(SimNumber)
            error(errstr);
        end
        if SimNumber == 2
            error(errstr);
        end
    end


    sbase = ['!matlab -nodesktop -nosplash -singleCompThread'...
         ' -r "RunMonteCarlo(''%s'', %d, %d)" &'];

    %###########################################################################
    % HARDCODING
    %###########################################################################
    NMULTIPLIERLIST = [.5 1 2];
    for n = 1:1:length(NMULTIPLIERLIST)
        s = sprintf(sbase, SaveDir, SimNumber, NMULTIPLIERLIST(n));
        disp(s);
        eval(s);
        pause(10);
    end
end
