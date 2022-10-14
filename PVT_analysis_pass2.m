load("PVTmanipulations.mat")


ROWS_PER_SUB = 3; %stores number of rows with screen visits excluded
SUBCOUNT = 18;

% group into three tables, ignore pre-ride visits
placRuns = pvt(strcmp(string(pvt.DRUG), 'PPREDRUG') | strcmp(string(pvt.DRUG), 'PPOSTDRUG') | strcmp(string(pvt.DRUG), 'PPOSTRIDE'),:);
chlorRuns = pvt(strcmp(string(pvt.DRUG), 'CPREDRUG') | strcmp(string(pvt.DRUG), 'CPOSTDRUG') | strcmp(string(pvt.DRUG), 'CPOSTRIDE'),:);
clepRuns = pvt(strcmp(string(pvt.DRUG), 'CEPREDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTRIDE'),:);

% make table of differences of rxn time ALL_MEAN from:
%   1. pre-drug to pre-ride
%   2. pre-ride to post-ride
%   3. pre-drug to post-ride
meanDiffsPlac = make_diff_table(placRuns, 'ALL_MEAN', ROWS_PER_SUB);
meanDiffsChlor = make_diff_table(chlorRuns, 'ALL_MEAN', ROWS_PER_SUB);
meanDiffsClep = make_diff_table(clepRuns, 'ALL_MEAN', ROWS_PER_SUB);


%concat of each difference of stats
placDrPr = renamevars(meanDiffsPlac(:, 2), 'DRUG_PRERIDE', 'PLAC_DRPR');
chlorDrPr = renamevars(meanDiffsChlor(:, 2), 'DRUG_PRERIDE', 'CHLOR_DRPR');
clepDrPr = renamevars(meanDiffsClep(:, 2), 'DRUG_PRERIDE', 'CLEP_DRPR');

diffComp_DRPR_TABLE = [placDrPr chlorDrPr clepDrPr];


placPrPo = renamevars(meanDiffsPlac(:, 3), 'PRERIDE_POSTRIDE', 'PLAC_PRPO');
chlorPrPo = renamevars(meanDiffsChlor(:, 3), 'PRERIDE_POSTRIDE', 'CHLOR_PRPO');
clepPrPo = renamevars(meanDiffsClep(:, 3), 'PRERIDE_POSTRIDE', 'CLEP_PRPO');

diffComp_PRPO_TABLE = [placPrPo chlorPrPo clepPrPo];


placDrPo = renamevars(meanDiffsPlac(:, 4), 'DRUG_POSTRIDE', 'PLAC_DRPO');
chlorDrPo = renamevars(meanDiffsChlor(:, 4), 'DRUG_POSTRIDE', 'CHLOR_DRPO');
clepDrPo = renamevars(meanDiffsClep(:, 4), 'DRUG_POSTRIDE', 'CLEP_DRPO');

diffComp_DRPO_TABLE = [placDrPo chlorDrPo clepDrPo];

%array-ify

diffComp_DRPR = table2array(diffComp_DRPR_TABLE);
diffComp_PRPO = table2array(diffComp_PRPO_TABLE);
diffComp_DRPO = table2array(diffComp_DRPO_TABLE);

%See if any diffs are sig
%PROG 10/7: "Must have at least 2 rows and columns error" got more than
%that for sure
p1 = friedman(diffComp_DRPR, SUBCOUNT);
p2 = friedman(diffComp_PRPO, SUBCOUNT);
p3 = friedman(diffComp_DRPO, SUBCOUNT);




% make_diff_table constructs a table that stores the three differences
% post-preride, preride-predrug, post-predrug
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

        %ID, pre/post reference the chair ride
        subID = tableIn(startNum, :).SUBJECT;
        drugMean = tableIn(startNum, :).(tableInVar);
        preMean = tableIn(startNum+1, :).(tableInVar);
        postMean = tableIn(startNum+2, :).(tableInVar);
    
        t(subject,:).SUBJECT = subID;
        t(subject,:).DRUG_PRERIDE = preMean - drugMean;
        t(subject,:).PRERIDE_POSTRIDE = postMean - preMean;
        t(subject,:).DRUG_POSTRIDE = postMean - drugMean;
        
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