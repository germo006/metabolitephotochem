% Noah Germolus 23 July 2022
% This is an analysis script for looking at the first photochem BC
% experiment.
% 1 Jan 2023 Reformatting for BCP2

%% Loading Data and Adding Information
clear
close all
clc
load("..\datasets\BCP2.2023.03.16.mat")
load("AlbumMaps.mat")
mtabNames = strrep(mtabNames,'′','''');
if ~exist("../graphs", "dir")
    mkdir("../graphs");
end

% We will use the sInfo variable a lot, to categorize and sort. So, we'll
% need to add columns for things like time point, type of sample, and notes
% for errors. 
mtabData_exp = mtabData_pM;
sInfo.cName{string(sInfo.cName) == 'BCP2_54_ASW_blank_1'} = 'BCP2_54_ASW_t5_ctrl_1';

sInfo.matrix = strings(length(sInfo.cName),1);
sInfo.timePoint = strings(length(sInfo.cName),1);
sInfo.sType = strings(length(sInfo.cName),1);
sInfo.rep = strings(length(sInfo.cName),1);

for i = 1:length(sInfo.cName)
    splitInfo = split(sInfo.cName(i),"_");
    if length(splitInfo) == 6
        sInfo.matrix(i) = splitInfo{3};
        sInfo.timePoint(i) = splitInfo{4};
        sInfo.sType(i) = splitInfo{5};
        sInfo.rep(i) = splitInfo{6};
    else
        sInfo.matrix(i) = splitInfo{3};
        sInfo.timePoint(i) = splitInfo{4};
        sInfo.sType(i) = "sample";
        sInfo.rep(i) = splitInfo{5};
    end
end

sInfo.timePoint = string(sInfo.timePoint);
sInfo.matrix = string(sInfo.matrix);
times = datetime(["02-12-2022 11:50:00";...
    "02-12-2022 13:11:00";...
    "02-12-2022 14:11:00";...
    "02-12-2022 16:14:00";...
    "02-12-2022 20:13:00";...
    "02-12-2022 23:50:00"], "InputFormat","dd-MM-uuuu HH:mm:ss");
timepts = "t"+ (0:5)';
sInfo.times = datetime(zeros(length(sInfo.cName),3));
for ii=1:6
    sInfo.times(sInfo.timePoint == timepts(ii)) = ...
        times(ii);
end

sInfo.plotGroups = findgroups(sInfo.matrix, sInfo.timePoint, sInfo.sType);

% groupmeans = array2table(splitapply(@means,mtabData_exp',sInfo_exp.plotGroups),'VariableNames',mtabNames);
% groupstd = array2table(splitapply(@stds,mtabData_exp',sInfo_exp.plotGroups)./sqrt(3),'VariableNames',mtabNames);
grpmeans = splitapply(@means,mtabData_exp',sInfo.plotGroups);
grpstd = splitapply(@stds,mtabData_exp',sInfo.plotGroups)./sqrt(3);

% singleRowInfo = sInfo_exp(cell2mat(sInfo_exp.rep)=='1', :);
% iTimesVSW = (singleRowInfo.matrix=="VSW" & singleRowInfo.sType=="sample");
% iTimesASW = (singleRowInfo.matrix=="ASW" & singleRowInfo.sType=="sample");
% iCV = (singleRowInfo.matrix=="VSW" & singleRowInfo.sType=="ctrl");
% iCA = (singleRowInfo.matrix=="ASW" & singleRowInfo.sType=="ctrl");
% iBV = (singleRowInfo.matrix=="VSW" & singleRowInfo.sType=="blank");
% iBA = (singleRowInfo.matrix=="ASW" & singleRowInfo.sType=="blank");

iTimesVSW = (sInfo.matrix=="VSW" & sInfo.sType=="sample");
iTimesASW = (sInfo.matrix=="ASW" & sInfo.sType=="sample");
iCV = (sInfo.matrix=="VSW" & sInfo.sType=="ctrl");
iCA = (sInfo.matrix=="ASW" & sInfo.sType=="ctrl");
iBV = (sInfo.matrix=="VSW" & sInfo.sType=="blank");
iBA = (sInfo.matrix=="ASW" & sInfo.sType=="blank");

%% This is a short diversion to look at flow cytometry data. 

if 0
    f = FCM_Vis;
    saveas(f, "../graphs/flowcyto.png", "png")
    close(f)
    clear f
end

%%
tBuffer = range(times)*.05;
tRange = [min(times)-tBuffer, max(times)+tBuffer];
tRange_g = rem(datenum(tRange), 1)*24;
tRange_g(2) = tRange_g(2) + 24;
tRange_g = tRange_g - tRange_g(1);
Durations = hours(times-min(times));
if ~exist("../graphs/png/", "dir")
    mkdir("../graphs/png/");
end
if ~exist("../graphs/logx/", "dir")
    mkdir("../graphs/logx/");
end
if 0
    for ii=1:length(mtabNames)
        f = figure("Visible","off");
        G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
        G2 = findgroups(sInfo.matrix(iTimesVSW), sInfo.timePoint(iTimesVSW), sInfo.sType(iTimesVSW));
        m1 = splitapply(@means,mtabData_exp(ii,iTimesASW)',G1);
        m1(isnan(m1)) = 0;
        m2 = splitapply(@means,mtabData_exp(ii,iTimesVSW)',G2);
        m2(isnan(m2)) = 0;
        ba = mtabData_exp(ii,iBA); ba(isnan(ba)) = 0;
        bv = mtabData_exp(ii,iBV); bv(isnan(bv)) = 0;
        ca = mtabData_exp(ii,iCA); ca(isnan(ca)) = 0;
        cv = mtabData_exp(ii,iCV); cv(isnan(cv)) = 0;
        
        EB1 = errorbar(Durations, m1,...
            splitapply(@stds,mtabData_exp(ii,iTimesASW)',G1),...
            'Color', shallowseas{1}, "LineWidth",1.5); hold on;
        EB2 = errorbar(Durations, m2,...
            splitapply(@stds,mtabData_exp(ii,iTimesVSW)',G2),...
            'Color', shallowseas{3}, "LineWidth",1.5);
        EB3 = scatter(repelem(Durations(1),1,3), ba, 30, shallowseas{2},...
            "filled", "Marker","^");
        EB4 = scatter(repelem(Durations(1),1,6), bv, 30,'k',...
            "filled", "Marker","v");
        EB5 = scatter(repelem(Durations(6),1,3), ca,30,shallowseas{4},...
            'Filled');
        EB6 = scatter(repelem(Durations(6),1,3), cv,30,shallowseas{5},...
            'Filled');
        ax = gca;
        set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
        title(mtabNames(ii))
        ylabel("Concentration, pM")
        xlabel("Time (h)")
        xlim(tRange_g)
        ylim([0,max(mtabData_exp(ii,:))])
        plot([0,13],[LOQ_pM(ii),LOQ_pM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
        plot([0,13],[MaxStd_pM(ii),MaxStd_pM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
        legend(["ASW Exp",...
            "VSW Exp",...
            "ASW Blank",...
            "VSW Blank",...
            "ASW Ctrl",...
            "VSW Ctrl",...
            "LOQ",...
            "MaxStd"], "Location","eastoutside")
        
        saveas(f, "../graphs/png/"+mtabNames(ii)+"_photo1.png", "png")
        close(f)
    end
    
    for ii=1:length(mtabNames)
        f = figure("Visible","off");
        G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
        G2 = findgroups(sInfo.matrix(iTimesVSW), sInfo.timePoint(iTimesVSW), sInfo.sType(iTimesVSW));
        m1 = splitapply(@means,mtabData_exp(ii,iTimesASW)',G1);
        m1(isnan(m1)) = 0;
        m2 = splitapply(@means,mtabData_exp(ii,iTimesVSW)',G2);
        m2(isnan(m2)) = 0;
        ba = mtabData_exp(ii,iBA); ba(isnan(ba)) = 0;
        bv = mtabData_exp(ii,iBV); bv(isnan(bv)) = 0;
        ca = mtabData_exp(ii,iCA); ca(isnan(ca)) = 0;
        cv = mtabData_exp(ii,iCV); cv(isnan(cv)) = 0;
        
        EB1 = errorbar(Durations, m1,...
            splitapply(@stds,mtabData_exp(ii,iTimesASW)',G1),...
            'Color', shallowseas{1}, "LineWidth",1.5); hold on;
        EB2 = errorbar(Durations, m2,...
            splitapply(@stds,mtabData_exp(ii,iTimesVSW)',G2),...
            'Color', shallowseas{3}, "LineWidth",1.5);
        EB3 = scatter(repelem(Durations(1),1,3), ba, 30, shallowseas{2},...
            "filled", "Marker","^");
        EB4 = scatter(repelem(Durations(1),1,6), bv, 30,'k',...
            "filled", "Marker","v");
        EB5 = scatter(repelem(Durations(6),1,3), ca, 30, shallowseas{4},...
            'Filled');
        EB6 = scatter(repelem(Durations(6),1,3), cv, 30, shallowseas{5},...
            'Filled');
        ax = gca;
        set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
        title(mtabNames(ii))
        ylabel("Concentration, pM")
        xlabel("Time (h)")
        xlim(tRange_g)
        ylim([0,max(mtabData_exp(ii,:))])
        set(ax, "xscale", "log")
        plot([0,13],[LOQ_pM(ii),LOQ_pM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
        plot([0,13],[MaxStd_pM(ii),MaxStd_pM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
        legend(["ASW Exp",...
            "VSW Exp",...
            "ASW Blank",...
            "VSW Blank",...
            "ASW Ctrl",...
            "VSW Ctrl",...
            "LOQ",...
            "MaxStd"], "Location","eastoutside")
        saveas(f, "../graphs/logx/"+mtabNames(ii)+"_photo1.png", "png")
        close(f)
    end
end

%% Glutamine, Glutamate, and GABA
% Figured I'd make a figure. 

GluBAnames = {"glutamic acid pos", "glutamine neg", "GABA pos"};
f = figure("Visible","off");
for ii=1:length(GluBAnames)
    ind = find(ismember(mtabNames, GluBAnames{ii}));
    G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
    m1 = splitapply(@means,mtabData_exp(ind,iTimesASW)',G1)./1000;
    mn = m1-m1(1);
    mn(isnan(m1)) = 0;
    errorbar(Durations, mn,...
        splitapply(@stds,mtabData_exp(ind,iTimesASW)',G1)./1000,...
        'Color', chainsaw{ii+1}, "LineWidth",1.5); hold on;
    hold on
end

ax = gca;
set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
ylabel("C-C_0 (nM)")
xlabel("Time (h)")
xlim(tRange_g)
legend(["glutamate",...
    "glutamine",...
    "GABA"], "Location","northeast")
saveas(f, "../graphs/GLUBA_photo1.png", "png")

%% Histidine, Aspartate, Asparagine
% Figured I'd make a figure. 

HAAnames = {"histidine pos", "asparagine neg", "aspartate pos"};
f = figure("Visible","on");
for ii=1:length(HAAnames)
    ind = find(ismember(mtabNames, HAAnames{ii}));
    G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
    m1 = splitapply(@means,mtabData_exp(ind,iTimesASW)',G1)./1000;
    mn = m1-m1(1);
    mn(isnan(m1)) = 0;
    errorbar(Durations, mn,...
        splitapply(@stds,mtabData_exp(ind,iTimesASW)',G1)./1000,...
        'Color', chainsaw{ii+1}, "LineWidth",1.5); hold on;
    hold on
end

ax = gca;
set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
ylabel("C-C_0 (nM)")
xlabel("Time (h)")
xlim(tRange_g)
legend(["histidine",...
    "asparagine",...
    "aspartate"], "Location","northeast")
saveas(f, "../graphs/HAA_photo1.png", "png")


%% Which metabolites changed? 
% Here's the part where we compare the initial and final time points, as
% well as the controls. 

i0V = (sInfo.matrix=="VSW" & sInfo.sType=="sample" & sInfo.timePoint=="t0");
i0A = (sInfo.matrix=="ASW" & sInfo.sType=="sample" & sInfo.timePoint=="t0");
i5V = (sInfo.matrix=="VSW" & sInfo.sType=="sample" & sInfo.timePoint=="t5");
i5A = (sInfo.matrix=="ASW" & sInfo.sType=="sample" & sInfo.timePoint=="t5");

diffVSW_time = ttest_spec(mtabData_pM(:,i0V),mtabData_pM(:,i5V));
diffASW_time = ttest_spec(mtabData_pM(:,i0A),mtabData_pM(:,i5A));
diffVSW_ctrl = ttest_spec(mtabData_pM(:,iCV),mtabData_pM(:,i5V));
diffASW_ctrl = ttest_spec(mtabData_pM(:,iCA),mtabData_pM(:,i5A));
diffVSW_good = diffVSW_ctrl & diffVSW_time;
diffASW_good = diffASW_ctrl & diffASW_time;

VSWChangedNames = mtabNames(diffVSW_good);
ASWChangedNames = mtabNames(diffASW_good);
ALLChangedNames = mtabNames(diffASW_good|diffVSW_good);

diffVSW_time_05 = ttest_spec(mtabData_pM(:,i0V),mtabData_pM(:,i5V), 0.05);
diffASW_time_05 = ttest_spec(mtabData_pM(:,i0A),mtabData_pM(:,i5A), 0.05);
diffVSW_ctrl_05 = ttest_spec(mtabData_pM(:,iCV),mtabData_pM(:,i5V), 0.05);
diffASW_ctrl_05 = ttest_spec(mtabData_pM(:,iCA),mtabData_pM(:,i5A), 0.05);
diffVSW_good_05 = diffVSW_ctrl_05 & diffVSW_time_05;
diffASW_good_05 = diffASW_ctrl_05 & diffASW_time_05;

VSWChangedNames_05 = mtabNames(diffVSW_good_05);
ASWChangedNames_05 = mtabNames(diffASW_good_05);
ALLChangedNames_05 = mtabNames(diffASW_good_05|diffVSW_good_05);
%% Calculating effective reaction rates.
% Here we will use the "linfit.m" function provided by Dave Glover to
% incorporate the uncertainties and data into linear regressions to get k,
% an effective reaction coefficient (first-order).
% LMAO nope I had to code my own function for a nonlinear regression, sort
% of.

ALLChangedNames = ALLChangedNames_05; % Fit all
if 1
    % Coefficients and uncertainties for ASW
    kASW = zeros(length(ALLChangedNames),1); % Rate coefficient (slope)
    bASW = zeros(length(ALLChangedNames),1); % Intercept (C0)
    rASW = zeros(length(ALLChangedNames),1); % Correlation coefficient
    kuASW = zeros(length(ALLChangedNames),1); % upper 90% CI in k
    buASW = zeros(length(ALLChangedNames),1); % upper 90% CI in b
    kbASW = zeros(length(ALLChangedNames),1); % upper 90% CI in k
    bbASW = zeros(length(ALLChangedNames),1); % upper 90% CI in b

    % Coefficients and uncertainties for VSW
    kVSW = zeros(length(ALLChangedNames),1);
    bVSW = zeros(length(ALLChangedNames),1);
    rVSW = zeros(length(ALLChangedNames),1);
    kuVSW = zeros(length(ALLChangedNames),1);
    buVSW = zeros(length(ALLChangedNames),1);
    kbVSW = zeros(length(ALLChangedNames),1);
    bbVSW = zeros(length(ALLChangedNames),1);


    k = table(ALLChangedNames, kASW, bASW, rASW, kuASW, kbASW, buASW, bbASW,...
        kVSW, bVSW, rVSW, kuVSW, kbVSW, buVSW, bbVSW);

    G1 = findgroups(sInfo.timePoint(iTimesASW));
    G2 = findgroups(sInfo.timePoint(iTimesVSW));

    for im=1:length(ALLChangedNames)
        % Retooling so I only do rates for things that I know change in a
        % relevant way
        
        ii = find(mtabNames==ALLChangedNames(im));
        

        %Trying a nM conversion to bring the regression vals into a similar
        %order of magnitude.
        m1 = splitapply(@means,mtabData_exp(ii,iTimesASW)',G1)./1000;
        m2 = splitapply(@means,mtabData_exp(ii,iTimesVSW)',G2)./1000;
        s1 = splitapply(@stds,mtabData_exp(ii,iTimesASW)',G1)./1000;
        s2 = splitapply(@stds,mtabData_exp(ii,iTimesVSW)',G2)./1000;
        if sum(isnan(m1)|m1==0)>0
            iBDL = find(isnan(m1)|m1==0);
            m1(iBDL(1)) = 1e-3;
            s1(iBDL(1)) = 1e-3;
        end
        if sum(isnan(m2)|m2==0)>0
            iBDL = find(isnan(m2)|m2==0);
            m2(iBDL(1)) = 1e-3;
            s2(iBDL(1)) = 1e-3;
        end
        missingVals1 = isnan(m1)|m1==0;
        missingVals2 = isnan(m2)|m2==0;
        m1(missingVals1) = [];
        m2(missingVals2) = [];
        s1(missingVals1) = [];
        s2(missingVals2) = [];
        tx1 = Durations; tx1(missingVals1) = [];
        tx2 = Durations; tx2(missingVals2) = [];
        % If the values dip below LOQ, set them to LOQ with stdev also LOQ.
        % I don't set them as zero because that would make the logarithm
        % busted.
        m1(m1==1e-3) = LOQ_pM(ii)./1000; s1(s1==1e-3) = LOQ_pM(ii)./1000;
        m2(m2==1e-3) = LOQ_pM(ii)./1000; s2(s2==1e-3) = LOQ_pM(ii)./1000;
        if sum(s1==0)>0
            disp("Check on the error values for "+mtabNames(ii)+" in ASW")
            s1(s1==0) = 0.05*m1(s1==0); % Set error to 5% of measurement.
        end
        if isempty(tx1) || length(m1)<3
            disp("For " + mtabNames(ii) + " there were too many missing values in ASW")
            k.bASW(im) = NaN; k.buASW(im) = NaN; k.bbASW(im) = NaN;
            k.kASW(im) = NaN; k.kuASW(im) = NaN; k.kbASW(im) = NaN;
            k.rASW(im) = NaN;
        elseif length(tx1) ~= length(m1)
            disp(tx1)
            disp(m1)
            disp("Unequal lengths for "+mtabNames(ii)+" in ASW")
            break
        else
            [coef, ~, r2, ~, ~, ccoef, bcmean] = expregress1(tx1,m1,s1);
            k.bASW(im) = bcmean(1);
            k.kASW(im) = -bcmean(2); % regression result is -k
            k.bbASW(im) = ccoef(1,1);
            k.buASW(im) = ccoef(1,2);
            k.kbASW(im) = -ccoef(2,1);
            k.kuASW(im) = -ccoef(2,2);
            k.rASW(im) = r2;
        end


        if sum(s2==0)>0
            disp("Check on the error values for "+mtabNames(ii)+" in VSW")
            s2(s2==0) = 0.05*m2(s2==0); % Set error to 5% of measurement.
        end
        if isempty(tx2) || length(m2)<3
            disp("For " + mtabNames(ii) + " there were too many missing values in VSW")
            k.bVSW(im) = NaN; k.buVSW(im) = NaN; k.bbVSW(im) = NaN;
            k.kVSW(im) = NaN; k.kuVSW(im) = NaN; k.kbVSW(im) = NaN;
            k.rVSW(im) = NaN;
        elseif length(tx2) ~= length(m2)
            disp(tx2)
            disp(m2)
            disp("Unequal lengths for "+mtabNames(ii)+" in VSW")
            break
        else
            [coef, ~, r2, ~, ~, ccoef, bcmean] = expregress1(tx2,m2,s2);
            k.bVSW(im) = bcmean(1);
            k.kVSW(im) = -bcmean(2); % regression result is -k
            k.bbVSW(im) = ccoef(1,1);
            k.buVSW(im) = ccoef(1,2);
            k.kbVSW(im) = -ccoef(2,1);
            k.kuVSW(im) = -ccoef(2,2);
            k.rVSW(im) = r2;
        end

        if strcmp("tryptophan neg", mtabNames(ii))
            %Set debug%
            % disp("stopped")
        end
        %     if strcmp(mtabNames(ii),"tryptamine pos")
        %         disp("Here it is!")
        %         disp(m1); disp(s1); disp(tx1);
        %     end
        clear coef r2 ccoef bcmean
    end


    % Degradation only! This will not take into account the possibility of
    % photoproduction (see kynurenine)
    kR = k;

    variableNames = ["Name","rate_h","intercept_nM","corrCoef","URate", "LRate","UInt", "LInt"];
    ASWrates = table(kR.ALLChangedNames, kR.kASW, kR.bASW, kR.rASW, kR.kuASW, kR.kbASW, kR.buASW, kR.bbASW,...
        'VariableNames', variableNames);
    VSWrates = table(kR.ALLChangedNames, kR.kVSW, kR.bVSW, kR.rVSW, kR.kuVSW, kR.kbVSW, kR.buVSW, kR.bbVSW,...
        'VariableNames', variableNames);

    removeASW = isnan(ASWrates.rate_h) | (imag(ASWrates.rate_h)~=0);
    ASWrates(removeASW,:) = [];
    ASWrates.rate_h = real(ASWrates.rate_h);
    ASWrates.URate = real(ASWrates.URate);
    ASWrates.LRate = real(ASWrates.LRate);
    removeVSW = isnan(VSWrates.rate_h) | (imag(VSWrates.rate_h)~=0);
    VSWrates(removeVSW,:) = [];
    VSWrates.rate_h = real(VSWrates.rate_h);
    VSWrates.URate = real(VSWrates.URate);
    VSWrates.LRate = real(VSWrates.LRate);
    fileBase = 'DegRates'; % Set this, don't mess with the automatic date system.
    today = datestr(datetime('now'),'yyyy-mm-dd_');
    NameOfFile = string(['../datasets/',today,fileBase,'.mat']);
    save(NameOfFile, "ASWrates", "VSWrates")
end

%% Now for some more graphs
load('../datasets/2023-03-16_DegRates.mat')
C = @(A,k,t) A.*exp(-k.*t);
if 1
    t = 0:0.1:12;
    durASW = repelem(Durations,3); durASW(5,:) = [];
    durVSW = repelem(Durations,3);
    tRange = [-0.5,12.5];
    mtabData_exp_nM = mtabData_exp;
    mtabData_exp_nM = mtabData_exp_nM./1000;
    if ~exist("../graphs\rateplots", "dir")
        mkdir("../graphs\rateplots");
    end
    LOQ_nM = LOQ_pM./1000;
    MaxStd_nM = MaxStd_pM./1000;

    for ii=1:length(mtabNames)
        rateIndexVSW = ismember(mtabNames(ii),VSWrates.Name);
        rateIndexASW = ismember(mtabNames(ii),ASWrates.Name);
        if rateIndexASW == 0 && rateIndexVSW == 0 
            continue
        elseif rateIndexASW == 1 && rateIndexVSW == 0 
            iASW = ismember(ASWrates.Name, mtabNames(ii));

            Aupper = C(ASWrates.UInt(iASW), ASWrates.LRate(iASW),t);
            Alower = C(ASWrates.LInt(iASW), ASWrates.URate(iASW),t);
            Ahat = C(ASWrates.intercept_nM(iASW), ASWrates.rate_h(iASW),t);

            f = figure("Visible","off");
            ba =mtabData_exp_nM(ii,iBA); ba(isnan(ba)) = 0;
            ca =mtabData_exp_nM(ii,iCA); ca(isnan(ca)) = 0;

            ub = plot(t, Aupper, "--", "LineWidth", 1.5, "Color", chainsaw{3});
            hold on
            mb = plot(t, Ahat, "--", "LineWidth", 1.5, "Color", chainsaw{5});
            lb = plot(t, Alower, "--", "LineWidth", 1.5, "Color", chainsaw{3},...
                "HandleVisibility","off");
            scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 30, "filled",...
                'Color', chainsaw{1});
            scBA = scatter(repelem(Durations(1),1,3), ba, 30, chainsaw{2},...
                "filled", "Marker","v");
            scCA = scatter(repelem(Durations(6),1,3), ca,30,chainsaw{4},...
                'filled', "Marker","^");

            ax = gca;
            set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
            title(mtabNames(ii)+ ", ASW")
            ylabel("Concentration, nM")
            xlabel("Time (h)")
            xlim(tRange)
            ylim([0,max(mtabData_exp_nM(ii,:))])
            plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
            plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
            legend(["90% CI",...
                "Predicted Decay",...
                "Observations",...
                "Blank",...
                "Dark Incubations",...
                "LOQ",...
                "High Standard"], "Location","eastoutside")

        elseif rateIndexASW == 0 && rateIndexVSW == 1
            iVSW = ismember(VSWrates.Name, mtabNames(ii));

            Vupper = C(VSWrates.UInt(iVSW), VSWrates.LRate(iVSW),t);
            Vlower = C(VSWrates.LInt(iVSW), VSWrates.URate(iVSW),t);
            Vhat = C(VSWrates.intercept_nM(iVSW), VSWrates.rate_h(iVSW),t);
            
            f = figure("Visible","off");
            bv =mtabData_exp_nM(ii,iBV); bv(isnan(bv)) = 0;
            cv =mtabData_exp_nM(ii,iCV); cv(isnan(cv)) = 0;

            ub = plot(t, Vupper, "--", "LineWidth", 1.5, "Color", chainsaw{3});
            hold on
            mb = plot(t, Vhat, "--", "LineWidth", 1.5, "Color", chainsaw{5});
            lb = plot(t, Vlower, "--", "LineWidth", 1.5, "Color", chainsaw{3},...
                "HandleVisibility","off");
            scV = scatter(durVSW,mtabData_exp_nM(ii,iTimesVSW), 30, "filled",...
                'Color', chainsaw{1});
            scBV = scatter(repelem(Durations(1),1,6), bv, 30,chainsaw{2},...
                "filled", "Marker","v");
            scCV = scatter(repelem(Durations(6),1,3), cv,30,chainsaw{4},...
                'filled', "Marker","^");

            ax = gca;
            set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
            title(mtabNames(ii)+ ", VSW")
            ylabel("Concentration, nM")
            xlabel("Time (h)")
            xlim(tRange)
            ylim([0,max(mtabData_exp_nM(ii,:))])
            plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
            plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
            legend(["90% CI",...
                "Predicted Decay",...
                "Observations",...
                "Blank",...
                "Dark Incubations",...
                "LOQ",...
                "High Standard"], "Location","eastoutside")
        else
            iVSW = ismember(VSWrates.Name, mtabNames(ii));
            iASW = ismember(ASWrates.Name, mtabNames(ii));

            Aupper = C(ASWrates.UInt(iASW), ASWrates.LRate(iASW),t);
            Alower = C(ASWrates.LInt(iASW), ASWrates.URate(iASW),t);
            Ahat = C(ASWrates.intercept_nM(iASW), ASWrates.rate_h(iASW),t);

            Vupper = C(VSWrates.UInt(iVSW), VSWrates.LRate(iVSW),t);
            Vlower = C(VSWrates.LInt(iVSW), VSWrates.URate(iVSW),t);
            Vhat = C(VSWrates.intercept_nM(iVSW), VSWrates.rate_h(iVSW),t);

            f = figure("Visible","off", "Position",[5, 5, 1200, 600], "Units","inches");
            ba =mtabData_exp_nM(ii,iBA); ba(isnan(ba)) = 0;
            bv =mtabData_exp_nM(ii,iBV); bv(isnan(bv)) = 0;
            ca =mtabData_exp_nM(ii,iCA); ca(isnan(ca)) = 0;
            cv =mtabData_exp_nM(ii,iCV); cv(isnan(cv)) = 0;

            subplot(1,2,1)
            ub = plot(t, Aupper, "--", "LineWidth", 2, "Color", chainsaw{3});
            hold on
            mb = plot(t, Ahat, "--", "LineWidth", 2, "Color", chainsaw{5});
            lb = plot(t, Alower, "--", "LineWidth", 2, "Color", chainsaw{3},...
                "HandleVisibility","off");
            scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 50, "filled",...
                'Color', chainsaw{1});
            scBA = scatter(repelem(Durations(1),1,3), ba, 50, chainsaw{2},...
                "filled", "Marker","v");
            scCA = scatter(repelem(Durations(6),1,3), ca,50,chainsaw{4},...
                'filled', "Marker","^");

            ax = gca;
            set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
            title("ASW")
            ylabel("Concentration, nM")
            xlabel("Time (h)")
            xlim(tRange)
            ylim([0,max(mtabData_exp_nM(ii,:))])
            plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
            plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
            legend(["90% CI",...
                "Predicted Decay",...
                "Observations",...
                "Blank",...
                "Dark Incubations",...
                "LOQ",...
                "High Standard"], "Location","northeast")

            subplot(1,2,2)
            ub = plot(t, Vupper, "--", "LineWidth", 2, "Color", chainsaw{3});
            hold on
            mb = plot(t, Vhat, "--", "LineWidth", 2, "Color", chainsaw{5});
            lb = plot(t, Vlower, "--", "LineWidth", 2, "Color", chainsaw{3},...
                "HandleVisibility","off");
            scV = scatter(durVSW,mtabData_exp_nM(ii,iTimesVSW), 50, "filled",...
                'Color', chainsaw{1});
            scBV = scatter(repelem(Durations(1),1,6), bv, 50,chainsaw{2},...
                "filled", "Marker","v");
            scCV = scatter(repelem(Durations(6),1,3), cv,50,chainsaw{4},...
                'Filled', "Marker","^");

            ax = gca;
            set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
            title("VSW")
            xlabel("Time (h)")
            xlim(tRange)
            ylim([0,max(mtabData_exp_nM(ii,:))])
            plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
            plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
            legend(["90% CI",...
                "Predicted Decay",...
                "Observations",...
                "Blank",...
                "Dark Incubations",...
                "LOQ",...
                "High Standard"], "Location","northeast")
        end

        saveas(f, "../graphs/rateplots/"+mtabNames(ii)+"_rate.png", "png")
        close(f)
    end
end

%% Now for the magnum opus: tryptophan quantum yields.
eps = readtable("../datasets/trpext.xlsx", "Range","C1:E502");
eps.l = eps.ReducedLambda;
eps.ep = eps.Eps_zero;
eps(:,1:3) = [];
lrange = 280:699;
eps_quant = eps.ep(ismember(eps.l, lrange));
eps_quant(eps_quant==0,:) = 1;
eps_quant = movmean(eps_quant,5,"Endpoints","shrink");


%only worth calculating where there's acutally trp absorption

load('../datasets/irradiance.mat')

EA = E_ASW(1:5,1:128);
EV = E_VSW(1:5,1:128);
lrange = lrange(1:128);
eps_quant = eps_quant(1:128);

ii = find(mtabNames=="tryptophan pos");

LOQ_nM = LOQ_pM./1000;

G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
G2 = findgroups(sInfo.matrix(iTimesVSW), sInfo.timePoint(iTimesVSW), sInfo.sType(iTimesVSW));
m1 = splitapply(@means,mtabData_exp(ii,iTimesASW)',G1)./1000;
s1 = splitapply(@stds,mtabData_exp(ii,iTimesASW)',G1)./1000;
m1(isnan(m1)|m1==0,:) = LOQ_nM(ii); s1(isnan(s1)|s1==0,:) = 0.05*LOQ_nM(ii);
m2 = splitapply(@means,mtabData_exp(ii,iTimesVSW)',G2)./1000;
s2 = splitapply(@stds,mtabData_exp(ii,iTimesVSW)',G2)./1000;
m2(isnan(m2)|m2==0,:) = LOQ_nM(ii); s2(isnan(s2)|s2==0,:) = 0.05*LOQ_nM(ii);

%% Separate section to run
Q = eps_quant./(2*max(eps_quant)); 
Qmax = ones(size(eps_quant));
Qmin = 1e-7*Q;
chiA = @(Q) chisq_quantum(m1,s1,EA,eps_quant,Q,Durations,lrange);
opts = optimset("PlotFcns", "optimplotfval", "MaxFunEvals", 3e7, "TolX", 1e-6);
QAest = fmincon(chiA,Q, [],[],[],[],Qmin,Qmax,[],opts);

% chiV = @(Q) chisq_quantum(m2,s2,EV,eps_quant,Q,Durations,lrange);
% opts = optimset("PlotFcns", "optimplotfval", "MaxFunEvals", 3e6);
% QVest = fmincon(chiV,Q, [],[],[],[],Qmin,Qmax,[],opts);
figure
% subplot(1,2,1)
plot(lrange, QAest, "LineWidth",2)
xlabel("\lambda, nm"); ylabel("\Phi_{ASW}");
xlim([280 407])
set(gca, "YScale", "log")
yyaxis right
plot(lrange, eps_quant, "LineWidth",2)
ylabel("\epsilon, M^{-1} cm^{-1}")
set(gca, "YColor", "k", "YScale", "log")
% subplot(1,2,2)
% plot(lrange, QVest, "LineWidth",2)
% xlabel("\lambda, nm"); ylabel("\Phi_{VSW}");
% xlim([280 407])

%% Estimating Decay with QY

durASW = repelem(Durations,3); durASW(5,:) = [];
% durVSW = repelem(Durations,3);
tRange = [-0.5,12.5];
mtabData_exp_nM = mtabData_exp./1000;
if ~exist("../graphs\trpPlot", "dir")
    mkdir("../graphs\trpPlot");
end
LOQ_nM = LOQ_pM./1000;
MaxStd_nM = MaxStd_pM./1000;
t = 0:0.1:12;

f = figure("Visible","on", "Position",[5, 5, 1200, 600], "Units","inches");
ba =mtabData_exp_nM(ii,iBA); ba(isnan(ba)) = 0;
bv =mtabData_exp_nM(ii,iBV); bv(isnan(bv)) = 0;
ca =mtabData_exp_nM(ii,iCA); ca(isnan(ca)) = 0;
cv =mtabData_exp_nM(ii,iCV); cv(isnan(cv)) = 0;

% subplot(1,2,1)
for ii = 1:5
    dCdt = @(t,Ct) -Ct.*int_solver(EA(ii,:)',eps_quant,QAest,lrange);
    tc = 0:1/3600:Durations(ii+1);
    [ty, yhat] = ode45(dCdt,tc,m1(1));
    ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", chainsaw{3});
    if ii~=1
        ub(ii).HandleVisibility = "off";
    end
    hold on
end
ii = find(mtabNames=="tryptophan pos");
scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 50, "filled",...
    'Color', chainsaw{1});
scBA = scatter(repelem(Durations(1),1,3), ba, 50, chainsaw{2},...
    "filled", "Marker","v");
scCA = scatter(repelem(Durations(6),1,3), ca,50,chainsaw{4},...
    'filled', "Marker","^");

ax = gca;
set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
title("ASW")
ylabel("Concentration, nM")
xlabel("Time (h)")
xlim(tRange)
ylim([0,max(mtabData_exp_nM(ii,:))])
plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
legend(["Predicted Decay",...
    "Observations",...
    "Blank",...
    "Dark Incubations",...
    "LOQ",...
    "High Standard"], "Location","northeast")

% subplot(1,2,2)
% for ii = 1:5
%     dCdt = @(t,Ct) -Ct.*int_solver(EV(ii,:)',eps_quant,QVest,lrange);
%     tc = 0:1/3600:Durations(ii+1);
%     [ty, yhat] = ode45(dCdt,tc,m2(1));
%     ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", chainsaw{3});
%     if ii~=1
%         ub(ii).HandleVisibility = "off";
%     end
%     hold on
% end
% ii = find(mtabNames=="tryptophan pos");
% scV = scatter(durVSW,mtabData_exp_nM(ii,iTimesVSW), 50, "filled",...
%     'Color', chainsaw{1});
% scBV = scatter(repelem(Durations(1),1,6), bv, 50,chainsaw{2},...
%     "filled", "Marker","v");
% scCV = scatter(repelem(Durations(6),1,3), cv,50,chainsaw{4},...
%     'Filled', "Marker","^");
% 
% ax = gca;
% set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
% title("VSW")
% xlabel("Time (h)")
% xlim(tRange)
% ylim([0,max(mtabData_exp_nM(ii,:))])
% plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
% plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
% legend(["Predicted Decay",...
%     "Observations",...
%     "Blank",...
%     "Dark Incubations",...
%     "LOQ",...
%     "High Standard"], "Location","northeast")

saveas(f, "../graphs/trpPlot/trpQY.png", "png")
close(f)

%% BINNED QY ESTIMATION
addpath("./quantum_bin")
Q = eps_quant./(2*max(eps_quant)); 
Qmax = ones(size(eps_quant));
Qlength = size(Q,1);
Qmin = 1e-7*Q;
Qbin = ones(5,1);
Qmin_bin = ones(5,1);
Qmax_bin = ones(5,1);
Qbounds = [1;20;40;65;85;128];
for ii = 1:size(Qbounds,1)-1
    Qbin(ii) = mean(Q(ii:ii+1));
    Qmax_bin(ii) = mean(Qmax(ii:ii+1));
    Qmin_bin(ii) = mean(Qmin(ii:ii+1));
end
chiA = @(Qbin) chisq_quantum_bins(m1,s1,EA,eps_quant,Qlength, Qbin,Qbounds,Durations,lrange);
opts = optimset("PlotFcns", "optimplotfval", "MaxFunEvals", 3e7, "TolX", 1e-6);
QAest = fmincon(chiA,Qbin, [],[],[],[],Qmin_bin,Qmax_bin,[],opts);

QA_breakout = ones(size(eps_quant,1),1);
for ii = 1:size(Qbounds,1)-1
    QA_breakout(Qbounds(ii):Qbounds(ii+1)) = QAest(ii);
end

figure
plot(lrange, QA_breakout, "LineWidth",2)
xlabel("\lambda, nm"); ylabel("\Phi_{ASW}");
xlim([280 407])
set(gca, "YScale", "log")
yyaxis right
plot(lrange, eps_quant, "LineWidth",2)
ylabel("\epsilon, M^{-1} cm^{-1}")
set(gca, "YColor", "k", "YScale", "log")


%% Estimating Decay 

durASW = repelem(Durations,3); durASW(5,:) = [];
tRange = [-0.5,12.5];
mtabData_exp_nM = mtabData_exp./1000;
if ~exist("../graphs\trpPlot", "dir")
    mkdir("../graphs\trpPlot");
end
LOQ_nM = LOQ_pM./1000;
MaxStd_nM = MaxStd_pM./1000;
t = 0:0.1:12;

f = figure("Visible","on", "Position",[5, 5, 1200, 600], "Units","inches");
ba =mtabData_exp_nM(ii,iBA); ba(isnan(ba)) = 0;
bv =mtabData_exp_nM(ii,iBV); bv(isnan(bv)) = 0;
ca =mtabData_exp_nM(ii,iCA); ca(isnan(ca)) = 0;
cv =mtabData_exp_nM(ii,iCV); cv(isnan(cv)) = 0;

for ii = 1:5
    dCdt = @(t,Ct) -Ct.*int_solver(EA(ii,:)',eps_quant,QA_breakout,lrange);
    tc = 0:1/60:Durations(ii+1);
    [ty, yhat] = ode45(dCdt,tc,m1(1));
    ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", chainsaw{3});
    if ii~=1
        ub(ii).HandleVisibility = "off";
    end
    hold on
end
ii = find(mtabNames=="tryptophan pos");
scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 50, "filled",...
    'Color', chainsaw{1});
scBA = scatter(repelem(Durations(1),1,3), ba, 50, chainsaw{2},...
    "filled", "Marker","v");
scCA = scatter(repelem(Durations(6),1,3), ca,50,chainsaw{4},...
    'filled', "Marker","^");

ax = gca;
set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
title("ASW")
ylabel("Concentration, nM")
xlabel("Time (h)")
xlim(tRange)
ylim([0,max(mtabData_exp_nM(ii,:))])
plot([0,13],[LOQ_nM(ii),LOQ_nM(ii)],"Color",[0.5 0.5 0.5], 'LineStyle', '--', 'LineWidth',2)
plot([0,13],[MaxStd_nM(ii),MaxStd_nM(ii)],"Color",[0.3 0.3 0.3], 'LineStyle', ':', 'LineWidth',2)
legend(["Predicted Decay",...
    "Observations",...
    "Blank",...
    "Dark Incubations",...
    "LOQ",...
    "High Standard"], "Location","northeast")

saveas(f, "../graphs/trpPlot/trpQY.png", "png")
%close(f)



