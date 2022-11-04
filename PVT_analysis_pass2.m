load("PVTmanipulations.mat")

ROWS_PER_SUB = 3; %stores number of rows with screen visits excluded
SUBCOUNT = 18; %number of study subjects/number of columns
NUM3DIFFS = 3; %number of columns for Friedman's test
NUM2DIFFS = 2;
CRIT_CHI2_2DF = 5.991; % critical chi2 val for df = 2 (comparing 3 columns)
CRIT_CHI2_1DF = 3.841; % critical chi2 val for df = 1 (comparing 2 columns)


testres1 = studyColumn3diffs('ALL_MEAN', pvt, ROWS_PER_SUB, NUM3DIFFS, CRIT_CHI2_2DF);
testres2 = studyColumn2diffs('ALL_MEAN', pvt, ROWS_PER_SUB, NUM2DIFFS, CRIT_CHI2_1DF);



% takes a statistic from column names and calculates if there
% is a significant difference across each patient for 3 differences:
    %   1. pre-drug to pre-ride
    %   2. pre-ride to post-ride
    %   3. pre-drug to post-ride
function results = studyColumn3diffs(statistic, pvt, ROWS_PER_SUB, NUM3DIFFS, CRIT_CHI2_2DF)
    % group into three tables, ignore pre-ride visits
    placRuns = pvt(strcmp(string(pvt.DRUG), 'PPREDRUG') | strcmp(string(pvt.DRUG), 'PPOSTDRUG') | strcmp(string(pvt.DRUG), 'PPOSTRIDE'),:);
    chlorRuns = pvt(strcmp(string(pvt.DRUG), 'CPREDRUG') | strcmp(string(pvt.DRUG), 'CPOSTDRUG') | strcmp(string(pvt.DRUG), 'CPOSTRIDE'),:);
    clepRuns = pvt(strcmp(string(pvt.DRUG), 'CEPREDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTRIDE'),:);
    
    % make table of differences of the statistic from:
    %   1. pre-drug to pre-ride
    %   2. pre-ride to post-ride
    %   3. pre-drug to post-ride
    meanDiffsPlac = make_3diff_table(placRuns, statistic, ROWS_PER_SUB);
    meanDiffsChlor = make_3diff_table(chlorRuns, statistic, ROWS_PER_SUB);
    meanDiffsClep = make_3diff_table(clepRuns, statistic, ROWS_PER_SUB);
    
    
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
    %PROG 11/4: "Must have at least 2 rows and columns error" from 10/7
    %           fixed, was using number of rows rather than cols
    
    p1 = friedman(diffComp_DRPR, NUM3DIFFS);
    p2 = friedman(diffComp_PRPO, NUM3DIFFS);
    p3 = friedman(diffComp_DRPO, NUM3DIFFS);
    
    % initialize results cell array to 0s == False
    results = {0, 0, 0};
    
    if p1 > CRIT_CHI2_2DF
        results{1} = 1;
    end
    if p2 > CRIT_CHI2_2DF
        results{2} = 1;
    end
    if p3 > CRIT_CHI2_2DF
        results{3} = 1;
    end
end

% takes a statistic from column names and calculates if there
% is a significant difference across each patient for 2 differences:
    %   1. pre-drug to pre-ride
    %   2. pre-drug to post-ride
function results = studyColumn2diffs(statistic,pvt, ROWS_PER_SUB, NUM2DIFFS, CRIT_CHI2_1DF)
    % group into three tables, ignore pre-ride visits
    placRuns = pvt(strcmp(string(pvt.DRUG), 'PPREDRUG') | strcmp(string(pvt.DRUG), 'PPOSTDRUG') | strcmp(string(pvt.DRUG), 'PPOSTRIDE'),:);
    chlorRuns = pvt(strcmp(string(pvt.DRUG), 'CPREDRUG') | strcmp(string(pvt.DRUG), 'CPOSTDRUG') | strcmp(string(pvt.DRUG), 'CPOSTRIDE'),:);
    clepRuns = pvt(strcmp(string(pvt.DRUG), 'CEPREDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTRIDE'),:);
    
    % make table of differences of the statistic from:
    %   1. pre-drug to pre-ride
    %   2. pre-drug to post-ride
    meanDiffsPlac = make_2diff_table(placRuns, statistic, ROWS_PER_SUB);
    meanDiffsChlor = make_2diff_table(chlorRuns, statistic, ROWS_PER_SUB);
    meanDiffsClep = make_2diff_table(clepRuns, statistic, ROWS_PER_SUB);
    
    
    %concat of each difference of stats
    placDrPr = renamevars(meanDiffsPlac(:, 2), 'DRUG_PRERIDE', 'PLAC_DRPR');
    chlorDrPr = renamevars(meanDiffsChlor(:, 2), 'DRUG_PRERIDE', 'CHLOR_DRPR');
    clepDrPr = renamevars(meanDiffsClep(:, 2), 'DRUG_PRERIDE', 'CLEP_DRPR');
    
    diffComp_DRPR_TABLE = [placDrPr chlorDrPr clepDrPr];
    
    placDrPo = renamevars(meanDiffsPlac(:, 3), 'DRUG_POSTRIDE', 'PLAC_DRPO');
    chlorDrPo = renamevars(meanDiffsChlor(:, 3), 'DRUG_POSTRIDE', 'CHLOR_DRPO');
    clepDrPo = renamevars(meanDiffsClep(:, 3), 'DRUG_POSTRIDE', 'CLEP_DRPO');
    
    diffComp_DRPO_TABLE = [placDrPo chlorDrPo clepDrPo];
    
    %array-ify
    
    diffComp_DRPR = table2array(diffComp_DRPR_TABLE);
    diffComp_DRPO = table2array(diffComp_DRPO_TABLE);
    
    %See if any diffs are sig
    %PROG 11/4: "Must have at least 2 rows and columns error" from 10/7
    %           fixed, was using number of rows rather than cols
    
    p1 = friedman(diffComp_DRPR, NUM2DIFFS);
    p2 = friedman(diffComp_DRPO, NUM2DIFFS);
    
    % initialize results cell array to 0s == False
    results = {0, "NOT EXAMINED", 0};
    
    if p1 > CRIT_CHI2_1DF
        results{1} = 1;
    end
    if p2 > CRIT_CHI2_1DF
        results{3} = 1;
    end
end
% maybe need to remove "PRERIDE_POSTRIDE". Just the predrug to preride
% indicates drug effects, predrug to postride with drug runs compared
% to control show how well performance loss is mitigated.

% make_diff_table constructs a table that stores the three differences
% post-preride, preride-predrug, post-predrug
% between pre-drug,pre-ride,post-ride results of an input stat
function t = make_3diff_table(tableIn, tableInVar, rows_per_sub)
    
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

%same as make_3diff_table, but without preride-postride diff
function t = make_2diff_table(tableIn, tableInVar, rows_per_sub)
    
    %preallocate table
    subjects = height(tableIn)/rows_per_sub;
    sz = [subjects 3];
    varTypes = ["string","double","double"];
    varNames = ["SUBJECT","DRUG_PRERIDE","DRUG_POSTRIDE"];
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
        t(subject,:).DRUG_POSTRIDE = postMean - drugMean;
        
    end

end