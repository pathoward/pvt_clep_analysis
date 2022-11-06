% Do analysis for PVT data from Chlorphedra project
 
% Author: Jay Buckey
% Created: 31 March 2020
% Last modified: 31 March 2020

% Revised by Patrick Howard
% Revision: 4 November 2022
% Notes: changed Repeated Measures anova to Friedman's nonparametric

%clear %Clears workspace of previous variables
close all hidden; %Closes figures that may be open from previous analysis
% Add paths to key functions
%addpath('/Users/jaybuckey/Documents/MATLAB/Key functions');%Adds the functions folder to the path

PVT=readtable('PVTALL combined data.xls','ReadRowNames',false,'Sheet','PVT Combined');

% Sort by subject and drug
PVT=sortrows(PVT,{'SUBJECT','DRUG'});

% Variable for repeats (not currently used)
Repeats= table([1 2 3 4 5 6 7 8 9]','VariableNames',{'Measurements'});


PVTALL_MEAN=PVT(:,{'SUBJECT','DRUG','ALL_MEAN'});

PVTALL_MEAN=unstack(PVTALL_MEAN,'ALL_MEAN','DRUG');

PVTALL_MEAN=PVTALL_MEAN(:,{'PPREDRUG','PPOSTDRUG','PPOSTRIDE','CPREDRUG','CPOSTDRUG','CPOSTRIDE',...
    'CEPREDRUG','CEPOSTDRUG','CEPOSTRIDE'});

Time=[1 2 3 4 5 6 7 8 9];
rmPVTALL_MEAN = fitrm(PVTALL_MEAN,'PPREDRUG-CEPOSTRIDE ~ 1','WithinDesign',Time);
ranovatblPVTALL_MEAN = ranova(rmPVTALL_MEAN);
TukeyFactPVTALL_MEAN=multcompare(rmPVTALL_MEAN,'Time');
BonferroniPVTALL_MEAN=multcompare(rmPVTALL_MEAN,'Time','ComparisonType','bonferroni');
MargmeanFactall=margmean(rmPVTALL_MEAN,'Time');
disp(ranovatblPVTALL_MEAN);
figure('Name','ALL_MEAN')
plot(rmPVTALL_MEAN);

PVTALL_MED=PVT(:,{'SUBJECT','DRUG','ALL_MED'});

PVTALL_MED=unstack(PVTALL_MED,'ALL_MED','DRUG');

PVTALL_MED=PVTALL_MED(:,{'PPREDRUG','PPOSTDRUG','PPOSTRIDE','CPREDRUG','CPOSTDRUG','CPOSTRIDE',...
    'CEPREDRUG','CEPOSTDRUG','CEPOSTRIDE'});

Time=[1 2 3 4 5 6 7 8 9];
rmPVTALL_MED = fitrm(PVTALL_MED,'PPREDRUG-CEPOSTRIDE ~ 1','WithinDesign',Time);

ranovatblPVTALL_MED = ranova(rmPVTALL_MED);
TukeyFactPVTALL_MED=multcompare(rmPVTALL_MED,'Time');
BonferroniPVTALL_MED=multcompare(rmPVTALL_MED,'Time','ComparisonType','bonferroni');
MargmeanFactall=margmean(rmPVTALL_MED,'Time');
disp(ranovatblPVTALL_MED);
figure('Name','ALL_MED')
plot(rmPVTALL_MED);


PVTSLOW_MEAN=PVT(:,{'SUBJECT','DRUG','SLOW_MEAN'});

PVTSLOW_MEAN=unstack(PVTSLOW_MEAN,'SLOW_MEAN','DRUG');

PVTSLOW_MEAN=PVTSLOW_MEAN(:,{'PPREDRUG','PPOSTDRUG','PPOSTRIDE','CPREDRUG','CPOSTDRUG','CPOSTRIDE',...
    'CEPREDRUG','CEPOSTDRUG','CEPOSTRIDE'});

Time=[1 2 3 4 5 6 7 8 9];
rmPVTSLOW_MEAN = fitrm(PVTSLOW_MEAN,'PPREDRUG-CEPOSTRIDE ~ 1','WithinDesign',Time);
ranovatblPVTSLOW_MEAN = ranova(rmPVTSLOW_MEAN);
TukeyFactPVTSLOW_MEAN=multcompare(rmPVTSLOW_MEAN,'Time');
BonferroniPVTSLOW_MEAN=multcompare(rmPVTSLOW_MEAN,'Time','ComparisonType','bonferroni');
MargmeanFactall=margmean(rmPVTSLOW_MEAN,'Time');
disp(ranovatblPVTSLOW_MEAN);
figure('Name','SLOW_MEAN')
plot(rmPVTSLOW_MEAN);

PVTFAST_MEAN=PVT(:,{'SUBJECT','DRUG','FAST_MEAN'});

PVTFAST_MEAN=unstack(PVTFAST_MEAN,'FAST_MEAN','DRUG');

PVTFAST_MEAN=PVTFAST_MEAN(:,{'PPREDRUG','PPOSTDRUG','PPOSTRIDE','CPREDRUG','CPOSTDRUG','CPOSTRIDE',...
    'CEPREDRUG','CEPOSTDRUG','CEPOSTRIDE'});

Time=[1 2 3 4 5 6 7 8 9];
rmPVTFAST_MEAN = fitrm(PVTFAST_MEAN,'PPREDRUG-CEPOSTRIDE ~ 1','WithinDesign',Time);
ranovatblPVTFAST_MEAN = ranova(rmPVTFAST_MEAN);
TukeyFactPVTFAST_MEAN=multcompare(rmPVTFAST_MEAN,'Time');
BonferroniPVTFAST_MEAN=multcompare(rmPVTFAST_MEAN,'Time','ComparisonType','bonferroni');
MargmeanFactall=margmean(rmPVTFAST_MEAN,'Time');
disp(ranovatblPVTFAST_MEAN);
figure('Name','FAST_MEAN')
plot(rmPVTFAST_MEAN);