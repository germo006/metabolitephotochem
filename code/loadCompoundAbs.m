% 30 Jun 2023 NPG 
% This is for reading in the measured absorbance spectra of the compounds I
% measured. 

clear
clc
close all

dataPath = '../datasets/temp/MeasuredSpectra/';
specInfo = readtable('../datasets/temp/MeasuredSpectra/Sample Table.csv');
specInfo.FileName = [string(specInfo{:,"SampleID"}) + '.Sample.Raw.csv'];
baseline = readtable('../datasets/temp/MeasuredSpectra/100% or 0 Absorbance Baseline.Correction.Raw.csv');
l = baseline.nm;
specData = table('Size', [length(l),length(specInfo.Description)+1],...
    'VariableNames', ['l';specInfo.SampleID],...
    'VariableTypes', repelem("double",length(specInfo.SampleID)+1,1));
specData.l(:) = l;

for ii=1:size(specData,2)-1
    temp = readtable(horzcat(dataPath,specInfo.FileName{ii}), 'Range', "B:B");
    specData(:,ii+1) = temp(:,:);
end

% I'm going to take this dataset into Excel, actually. 

save("../datasets/MeasuredSpectraCompounds.mat")