% load("Updates 11-4.mat")

% v2.0 PVT significant differences
% Patrick Howard
% Space Medicine Innovations Lab, Dartmouth Hitchcock Medical Center

% work in bonferroni or tukey posthoc testing 

% This script isolates significant relationships between reaction time
% aggregate statistics. Applicable to a study comparing reaction times
% under the effects of motion sickness with management from
% chlorpheniramine, chlorpheniramine+ephedrine (Chlorphedra), and placebo
% treatment. P-value cutoff can be set here:
P_CUT = .05;

% Each statistic is individually evaluated. If two phases of the study have
% a significant difference for a given statistic, it is reported as:
%   [{STATISTIC} {PHASE1} {PHASE2} {P-VALUE} {SIGCOUNT} {TOTCOUNT}]
% where SIGCOUNT is the number of significant relationships observed up to
% that point, and TOTCOUNT is the total number of relationships tested
% (including forward/backward pairs, i.e. a:b and b:a are distinct). 
% 
% The reported p-value between two significantly different phases Phase1 
% and Phase2 is the minimum of Phase1 vs. Phase2 and Phase2 vs. Phase1


%%%%%% Set relevant constants %%%%%%%
ROWS_PER_SUB = 3; %stores number of rows with screen visits excluded
SUBCOUNT = 18; %number of study subjects/number of columns
NUM3DIFFS = 3; %number of columns for Friedman's test
NUM2DIFFS = 2;
NUMPHASES = 9; %number of timepoints measured
REPS = 1; %reps for Friedman's parameter
CRIT_CHI2_2DF = 5.991; % critical chi2 val for df = 2 (comparing 3 columns)
CRIT_CHI2_1DF = 3.841; % critical chi2 val for df = 1 (comparing 2 columns)
set(0,'DefaultFigureVisible','off') %turns off excessive graphics

%Specify which groups are being assessed
stats = {'ALL_MEAN', 'ALL_MED', 'SLOW_MEAN', 'FAST_MEAN', 'IALL_MEAN', 'IALL_MED'};
%stats = {'ALL_MEAN'};

%contains stats mapped to their results
results_alldiff = containers.Map();
results_3diff = containers.Map();
results_2diff = containers.Map();

%build pvt table
pvt = readtable("pvtDataNov11.csv");

%init analysis
res = studyColumnAll(stats, pvt, ROWS_PER_SUB, NUMPHASES, P_CUT);


% Convert cell to a table and use first row as variable names
resTable = cell2table(res);
resTable.Properties.VariableNames = ["Statistic" "Phase1" "Phase2" "P-Value" "Num Result" "Total Rels"];
 
% Write the table to a CSV file
writetable(resTable,'allStatsSignificant.csv')

function results = studyColumnAll(stats, pvt, ROWS_PER_SUB, NUMPHASES, P_CUT)
    
    stats_tables = containers.Map();

    %build table for each statistic
    for idx = 1:numel(stats)
        stat = stats{idx};     
        stats_tables(stat) = flattenData(stat, pvt, ROWS_PER_SUB);
    end
    
    statistics = keys(stats_tables);
    tables = values(stats_tables);
    sigMap = containers.Map();
    totalRels = 0; %counts relationships tested

    %2d array indicating if study phases had a sig difference
    for statIdx = 1:numel(statistics)
        
        statArray = {-1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1;
                     -1, -1, -1, -1, -1, -1, -1, -1, -1};

        stat = statistics{statIdx};
        table = tables{statIdx};

        %map phase numbers to their col names, more easily iterable
        phaseMap = containers.Map({1,2,3,4,5,6,7,8,9}, ...
            {'PPREDRUG','PPOSTDRUG','PPOSTRIDE', ...
            'CPREDRUG','CPOSTDRUG','CPOSTRIDE',...
            'CEPREDRUG','CEPOSTDRUG','CEPOSTRIDE'});
        

        % repeat comps in loops below (a:b and b:a) guarantee lowest pval

        for phase1 = 1:NUMPHASES
            for phase2 = 1:NUMPHASES
                
                if phase1 ~= phase2
                    
                    %count another tested relationship
                    totalRels = totalRels + 1;

                    %pull data columns of each phase in consideration
                    studyArray = table2array(table(:,{phaseMap(phase1),phaseMap(phase2)}));
                    p = friedman(studyArray, 1, 'off');
                    
                    %if significant, record relationship
                    if p <= P_CUT
                        
                        %find curr vals
                        curr_p_f = cell2mat(statArray(phase1, phase2));
                        curr_p_b = cell2mat(statArray(phase2, phase1));

                        %update both to lowest if currently -1, or not -1
                        %but greater than new p-value
                        if or(curr_p_f == -1, and(curr_p_f ~= -1, p < curr_p_f))
                            statArray(phase1, phase2) = num2cell(p);
                        end

                        if or(curr_p_b == -1, and(curr_p_b ~= -1, p < curr_p_b))
                            statArray(phase2, phase1) = num2cell(p);
                        end

                    end                   
    
                end
            end
        end
        
        %add results to the global results
        sigMap(stat) = statArray;

    end

    % init sig relationship counter, results 
    numSig = 0;
    results = cell([1 6]);

    %for each stat, search for significant relationships, and print
    for statIdx = 1:numel(statistics)
        stat = statistics{statIdx};
        sigArray = sigMap(stat);
        
        for row = 1:NUMPHASES
            for col = 1:NUMPHASES
            
                if col ~= row
                    
                    %If significant (!= -1) report result
                    if cell2mat(sigArray(row, col)) ~= -1 
                        numSig = numSig + 1;
                        report = [stat phaseMap(row) phaseMap(col) sigArray(row, col) numSig totalRels];
                        results(numSig,:) = report;
                    end

                end

            end
        end

    end

end






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


function t = flattenData(stat, pvt, rows_per_sub)
    
    statIsolated = pvt(:, {'SUBJECT', 'DRUG', stat});
    
    statUnstacked = unstack(statIsolated, stat, 'DRUG');

    statSorted = sortrows(statUnstacked,{'SUBJECT'});

    statFinal=statSorted(:,{'PPREDRUG','PPOSTDRUG','PPOSTRIDE','CPREDRUG','CPOSTDRUG','CPOSTRIDE',...
    'CEPREDRUG','CEPOSTDRUG','CEPOSTRIDE'});
    
    t = statFinal; 

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



%old code for running 3diff/2diff analysis

% % for each stat, complete the 3diff and 2diff analysis
% for idx = 1:numel(stats)
%     
%     stat = stats{idx};
%     
%     results_alldiff(stat) = studyColumnAll()
%     results_3diff(stat) = studyColumn3diffs(stat, pvt, ROWS_PER_SUB, NUM3DIFFS, CRIT_CHI2_2DF);
%     results_2diff(stat) = studyColumn2diffs(stat, pvt, ROWS_PER_SUB, NUM2DIFFS, CRIT_CHI2_1DF);
% 
% end
% 
% %return results of 3diff analysis
% k = keys(results_3diff) ;
% val = values(results_3diff) ;
% for i = 1:length(results_3diff)
%  [k{i} val{i}];
% end
% 
% %return results of 2diff analysis
% k = keys(results_2diff) ;
% val = values(results_2diff) ;
% for i = 1:length(results_2diff)
%  [k{i} val{i}];
% end