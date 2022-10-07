load("PVTmanipulations.mat")


ROWS_PER_SUB = 3; %stores number of rows with screen visits excluded

% group into three tables, ignore pre-ride visits
placRuns = pvt(strcmp(string(pvt.DRUG), 'PPREDRUG') | strcmp(string(pvt.DRUG), 'PPOSTDRUG') | strcmp(string(pvt.DRUG), 'PPOSTRIDE'),:);
chlorRuns = pvt(strcmp(string(pvt.DRUG), 'CPREDRUG') | strcmp(string(pvt.DRUG), 'CPOSTDRUG') | strcmp(string(pvt.DRUG), 'CPOSTRIDE'),:);
clepRuns = pvt(strcmp(string(pvt.DRUG), 'CEPREDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTRIDE'),:);

% % sort each group by phase of drug trial
% sortedPlacRuns = sortrows(placeboRuns, 'DRUG');
% sortedChlorRuns = sortrows(chlorRuns, 'DRUG');
% sortedClepRuns = sortrows(clepRuns, 'DRUG');

% make table of differences of rxn time ALL_MEAN from:
%   1. pre-drug to pre-ride
%   2. pre-ride to post-ride
%   3. pre-drug to post-ride


meanDiffsPlac = make_diff_table(placRuns, 'ALL_MEAN', ROWS_PER_SUB)

% subjectRow = table("SUBJECT", "DRUG_PRERIDE", "PRERIDE_POSTRIDE", "DRUG_POSTRIDE");
% startRow = placRuns(1*3-2, :);
% subjectRow.SUBJECT = startRow.SUBJECT;



% make_diff_table constructs a table that stores the three differences
% between pre-drug,pre-ride,post-ride results of an input stat

function t = make_diff_table(tableIn, tableInVar, rows_per_sub)
    
    %preallocate table
    subjects = height(tableIn)/rows_per_sub;
    sz = [subjects 4];
    varTypes = ["string","double","double","double"];
    varNames = ["SUBJECT","DRUG_PRERIDE","PRERIDE_POSTRIDE", "DRUG_POSTRIDE"];
    t = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);

    for subject = 1:subjects
        
        %proc one patient
        startNum = subject*3-2;

        startRow = tableIn(startNum, :);
        subjectRow.SUBJECT = startRow.SUBJECT;

        %pre/post reference the chair ride
        drugMean = tableIn(startNum, :).(tableInVar);
        preMean = tableIn(startNum+1, :).(tableInVar);
        postMean = tableIn(startNum+2, :).(tableInVar);
    
    end

end

% meanDiffsChlor
% meanDiffsClep

% % concat sorted means across three drug phases
% sPlacMeans = renamevars(sortedPlacRuns(:, 7), 'ALL_MEAN', 'PLAC_MEAN');
% sChlorMeans = renamevars(sortedChlorRuns(:, 7), 'ALL_MEAN', 'CHLOR_MEAN');
% sClepMeans = renamevars(sortedClepRuns(:, 7), 'ALL_MEAN', 'CLEP_MEAN');
% 
% % join into one table and convert to matrix
% sortedMeans = [sPlacMeans sChlorMeans sClepMeans];
% sortedMeansArray = table2array(sortedMeans)
% 
% p = friedman(sortedMeansArray, 18)