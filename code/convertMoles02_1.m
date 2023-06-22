function [mtabData_pM,MaxStd_pM,ratingFlags,mtabElem,LOQ_pM] =...
    convertMoles02(negTransitions, posTransitions, mtabNames, mtabData,...
    MaxStd, LOQ)

%% Converting to molar basis
% This process involves reloading the transition list files and matching
% names to molecular weights. For DT5, I had to manually correct many, many
% names. That's why the data file loaded at the beginning contains two
% variables for mtabNames.
% NPG 30 Dec 2022: Editing the code to V02 to use the LOQ calculated from
% replicate standards. 

posInfo = readtable(posTransitions);
posInfo(posInfo.isParent == 0,:) = [];
MWp = table(posInfo.CompoundName, posInfo.StdMW, posInfo.rating123, ...
    posInfo.C, posInfo.O, posInfo.N, 'VariableNames', ...
    {'CompoundName','StdMW','Rating', 'C', 'O', 'N'});
negInfo = readtable(negTransitions);
negInfo(negInfo.isParent == 0,:) = [];
MWn = table(negInfo.CompoundName, negInfo.StdMW, negInfo.rating123, ...
    negInfo.C, negInfo.O, negInfo.N, 'VariableNames', ...
    {'CompoundName','StdMW','Rating', 'C', 'O', 'N'});
MW = [MWp;MWn];
clear MWp MWn
MW = unique(MW, 'rows');
clear posInfo negInfo

% Making both compound name columns into strings and removing the neg/pos
% identifier.
MW.CompoundName = string(MW.CompoundName);
mtabNamesAgnostic = strrep(strrep(string(mtabNames), ' pos', ''),' neg','');

% Time to index where each unique molecule is found in mtabNames.
[~, iNames] = ismember(mtabNamesAgnostic, MW.CompoundName);
iNames(iNames==0) = [];
% Use those indices to sort MW values.
if sum(mtabNamesAgnostic == MW.CompoundName(iNames)) ~= length(mtabNames)
    disp("name mismatch")
    return
end
MWtoConvert = MW.StdMW(iNames);
ratingFlags = MW.Rating(iNames);
% convert from pg/mL to pM
mtabData_pM = mtabData.*1000./MWtoConvert;
MaxStd_pM = MaxStd.*1000./MWtoConvert;

LOQ_pM = LOQ.*1000./MWtoConvert;
mtabElem = table(MW.C(iNames), MW.N(iNames), MW.O(iNames),'VariableNames',...
    {'C','N','O'});

end