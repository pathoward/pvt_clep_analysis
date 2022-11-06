load("PVTmanipulations.mat");

%group into three tables, ignore pre-ride visits
placRuns = pvt(strcmp(string(pvt.DRUG), 'PPREDRUG') | strcmp(string(pvt.DRUG), 'PPOSTDRUG') | strcmp(string(pvt.DRUG), 'PPOSTRIDE'),:);
chlorRuns = pvt(strcmp(string(pvt.DRUG), 'CPREDRUG') | strcmp(string(pvt.DRUG), 'CPOSTDRUG') | strcmp(string(pvt.DRUG), 'CPOSTRIDE'),:);
clepRuns = pvt(strcmp(string(pvt.DRUG), 'CEPREDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTDRUG') | strcmp(string(pvt.DRUG), 'CEPOSTRIDE'),:);

% sort each group by phase of drug trial
sortedPlacRuns = sortrows(placeboRuns, 'DRUG');
sortedChlorRuns = sortrows(chlorRuns, 'DRUG');
sortedClepRuns = sortrows(clepRuns, 'DRUG');

%concat sorted means across three drug phases
sPlacMeans = renamevars(sortedPlacRuns(:, 7), 'ALL_MEAN', 'PLAC_MEAN');
sChlorMeans = renamevars(sortedChlorRuns(:, 7), 'ALL_MEAN', 'CHLOR_MEAN');
sClepMeans = renamevars(sortedClepRuns(:, 7), 'ALL_MEAN', 'CLEP_MEAN');

%join into one table and convert to matrix
sortedMeans = [sPlacMeans sChlorMeans sClepMeans];
sortedMeansArray = table2array(sortedMeans)

%test add

p = friedman(sortedMeansArray, 18)