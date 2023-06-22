% 12 Apr 2023 NPG ripped this code from the specvis.m script to create a
% loadable, preprocessed absorbance datafile. 

clear
clc
close all

dataPath = '../datasets/UV-Vis/';
specInfo = readtable('../datasets/UV-Vis/Sample Table.csv');
baseline = readtable('../datasets/UV-Vis/100% or 0 Absorbance Baseline.Correction.Raw.csv');
l = baseline.nm;
specData = table('Size', [length(l),length(specInfo.Description)+1],...
    'VariableNames', ['l';specInfo.SampleID],...
    'VariableTypes', repelem("double",length(specInfo.Spike)+1,1));
specData.l(:) = l;

for ii=1:size(specData,2)-1
    temp = readtable(horzcat(dataPath,specInfo.FileName{ii}), 'Range', "B:B");
    specData(:,ii+1) = temp(:,:);
end

repGroups = findgroups(specInfo.Matrix, specInfo.Spike, specInfo.timePoint, specInfo.isQC);
[groupInd, ia, ib] = unique(unique(repGroups, 'stable'));
meanNames = ["blank0", "MQ0","MQs0","ASW0","ASWs0","VSW0","VSWs0",...
    "blank12","ASW12","ASWs12","VSW12","VSWs12"];
horzmean = @(x) mean(x,2);
preMean = splitapply(horzmean,specData{:,2:end},repGroups');
preMean = preMean(:,ib); 

specmeans = array2table(preMean, 'VariableNames', meanNames);
specmeans.l = l;
means_BLsub = array2table(specmeans{:,1:12} - baseline.A,'VariableNames',meanNames);
means_BLsub.l = l;

diffmeans = table(l);
diffmeans.VSW = specmeans.VSW12 - specmeans.VSW0;
diffmeans.VSWs = specmeans.VSWs12 - specmeans.VSWs0;
diffmeans.ASW = specmeans.ASW12 - specmeans.ASW0;
diffmeans.ASWs = specmeans.ASWs12 - specmeans.ASWs0;
diffmeans.blanks = specmeans.blank12 - specmeans.blank0;
diffmeans.t0s = specmeans.MQs0 - specmeans.MQ0;
diffmeans.t0A = specmeans.ASW0 - specmeans.MQ0;
diffmeans.t0V = specmeans.VSW0 - specmeans.MQ0;
diffmeans.t0As = specmeans.ASWs0 - specmeans.MQ0;
diffmeans.t0Vs = specmeans.VSWs0 - specmeans.MQ0;

clear ii premean ia ib preMean temp

Napierian = @(x) 2.303.*x./0.1; % Converts to Napierian absorbance using a 0.1 m cuvette (units m^-1)

a_specmeans = specmeans;     a_specmeans{:,2:end-1} = Napierian(a_specmeans{:,2:end-1});
a_means_BLsub = means_BLsub; a_means_BLsub{:,2:end-1} = Napierian(a_means_BLsub{:,2:end-1});
a_diffmeans = diffmeans;     a_diffmeans{:,2:end-1} = Napierian(a_diffmeans{:,2:end-1});

save("../datasets/absorbance.mat")