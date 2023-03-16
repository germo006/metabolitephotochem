%% Noah Germolus 9 Jan 2023
% Here I will attempt to preprocess and visualize spectrophotometry data
% from BCP2. 

clear
clc
close all

dataPath = '../../UV-vis/NPG_BCP2/';
specInfo = readtable('../../UV-vis/NPG_BCP2/Sample Table.csv');
baseline = readtable('../../UV-vis/NPG_BCP2/100% or 0 Absorbance Baseline.Correction.Raw.csv');
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

means = array2table(preMean, 'VariableNames', meanNames);
means.l = l;
means_BLsub = array2table(means{:,1:12} - baseline.A,'VariableNames',meanNames);
means_BLsub.l = l;

diffmeans = table(l);
diffmeans.VSW = means.VSW12 - means.VSW0;
diffmeans.VSWs = means.VSWs12 - means.VSWs0;
diffmeans.ASW = means.ASW12 - means.ASW0;
diffmeans.ASWs = means.ASWs12 - means.ASWs0;
diffmeans.blanks = means.blank12 - means.blank0;
diffmeans.t0s = means.MQs0 - means.MQ0;
diffmeans.t0A = means.ASW0 - means.MQ0;
diffmeans.t0V = means.VSW0 - means.MQ0;
diffmeans.t0As = means.ASWs0 - means.MQ0;
diffmeans.t0Vs = means.VSWs0 - means.MQ0;

%% Advanced difference testing. 
% First, I will excise the spectra for which I can make t0/t12 comparisons
% and subtract the blank MQ baseline. 

t0 = specData{:,12:20} - specData.blank1;
t0grp = findgroups(repGroups(11:19));
t12 = specData{:,25:end} - specData.blank2;
t12grp = findgroups(repGroups(24:end));

% Now calculate the mean difference between t0 and t12.
diff_comp_means = splitapply(horzmean, t12, t12grp') - ...
    splitapply(horzmean, t0, t0grp');
diff_comp_means = [diffmeans(:,1), array2table(diff_comp_means, "VariableNames",["sASW", "VSW", "sVSW"])];

% Now take those sample groups and do a t-test. 
h_t = splitapply(@ttest_spec, t0, t12, t0grp');
h_t = logical(h_t);

if 0
    figure
    h1a = scatter(diff_comp_means.l(~h_t(:,1)), diff_comp_means.sASW(~h_t(:,1)),...
        3, [0.2, 0.2, 0.2], "filled", "o", "HandleVisibility","on");

    hold on
    h2a = scatter(diff_comp_means.l(~h_t(:,2)), diff_comp_means.VSW(~h_t(:,2)),...
        3, [0.2, 0.2, 0.2], "filled", "o", "HandleVisibility","off");
    h3a = scatter(diff_comp_means.l(~h_t(:,3)), diff_comp_means.sVSW(~h_t(:,3)),...
        3, [0.2, 0.2, 0.2], "filled", "o", "HandleVisibility","off");
    h1 = scatter(diff_comp_means.l(h_t(:,1)), diff_comp_means.sASW(h_t(:,1)),...
        5, chainsaw{1}, "filled", "o");
    h2 = scatter(diff_comp_means.l(h_t(:,2)), diff_comp_means.VSW(h_t(:,2)),...
        5, chainsaw{3}, "filled", "o");
    h3 = scatter(diff_comp_means.l(h_t(:,3)), diff_comp_means.sVSW(h_t(:,3)),...
        5, chainsaw{4}, "filled", "o");
    legend("Change Not Significant", "sASW", "VSW","sVSW")
    xlim([249,801])
    xlabel("\lambda, nm")
    ylabel("\Delta_{abs}")
end

%% Loading irradiance, Calculating Light Absorption
% See my giant word doc outlining what I'm doing here, but basically I need
% to evaluate the following equation:

% Wabs = @(I,A) 1.308e-7.*I.*A;

% ...for each sample. 
% So let's load that irradiance data I interpolated in loadradiometry.m
load("irradiance.mat")

% Somewhat unfortunately, the radiometer had a limited range. So we'll
% truncate the calculation to that range. 
% Oh, also the order of the wavelengths is flipped in the radiometer
% measurements vs. the UV-Vis. The indexing that follows will take care of
% that.
ia = find(ismember(l, l_SunTest));

Wabs_sASW = zeros(length(l_SunTest), 5);
Wabs_sVSW = Wabs_sASW;

for ii = 1:5
Wabs_sASW(:,ii) = Wcyl(flip(E_ASW(ii,:)'), means_BLsub.ASWs0(ia(ii)), 7);
Wabs_sVSW(:,ii) = Wcyl(flip(E_VSW(ii,:)'), means_BLsub.VSWs0(ia(ii)), 7);
end

[X,Y] = meshgrid(l_SunTest,1:5);
figure
subplot(2,1,1)
contourf(X,Y,Wabs_sASW'./1000)
hold on
text(675, 4.7, "sASW", "Color","w", "FontWeight","bold")
ylabel("Position")
c1 = colorbar(gca, "eastoutside");
c1.Label.String = "W_{abs}, \mumol photons s^{-1} L";
subplot(2,1,2)
contourf(X,Y,Wabs_sVSW'./1000)
hold on
text(675, 4.7, "sVSW", "Color","w", "FontWeight","bold")
xlabel("\lambda, nm")
ylabel("Position")
c2 = colorbar(gca, "eastoutside");
c2.Label.String = "W_{abs}, \mumol photons s^{-1} L";

%%
setDefaultFigs
load("AlbumMaps.mat", "chainsaw")

figure
% subplot(3,1,1)
% plot(means,"l", meanNames, 'LineWidth', 2)
% xlabel("\lambda, nm")
% ylabel("Absorbance")
% xlim([250,800])
% legend(meanNames, "NumColumns",2)
% title("Raw Spectra")

subplot(2,1,1)
h = plot(diffmeans, "l", ["t0s","t0A", "t0V", "t0As", "t0Vs"], 'LineWidth', 2);
h(1).Color = chainsaw{5};
h(2).Color = chainsaw{2};
h(3).Color = 0.7.*chainsaw{1};
h(4).Color = chainsaw{2};
h(5).Color = 0.7.*chainsaw{1};
h(4).LineStyle = "--";
h(5).LineStyle = "--";
xlabel("\lambda, nm")
ylabel("Absorbance")
xlim([250,800])
legend(["MQs", "ASW", "VSW", "ASWs", "VSWs"])
title("Initial Differences (- MQ blank)")

subplot(2,1,2)
h = plot(diffmeans, "l", ["blanks","ASWs", "VSW", "VSWs"], 'LineWidth', 2);
h(1).Color = chainsaw{5};
h(2).Color = chainsaw{2};
h(3).Color = 0.7.*chainsaw{1};
h(4).Color = 0.7.*chainsaw{1};
h(2).LineStyle = "--";
h(4).LineStyle = "--";
xlabel("\lambda, nm")
ylabel("Absorbance")
xlim([250,800])
legend(["blanks", "ASWs", "VSW","VSWs"])
title("Final Differences t_{12} - t_0")

