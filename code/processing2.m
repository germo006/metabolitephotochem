% Noah Germolus 23 July 2022
% This is an analysis script for looking at the first photochem BC
% experiment.
% 1 Jan 2023 Reformatting for BCP2

%% Loading Data and Adding Information
clear
close all
clc
load("..\datasets\BCP2.2023.06.06.mat")
load("../datasets/onlyGoodMtabs_10Apr2023.mat", "mtabNamesGood")
load("AlbumMaps.mat", "FPC")
col = FPC;

% Delete the metabolites we're not working with.
iDelete = ~ismember(mtabNames, mtabNamesGood);
vars = whos();
siz = 0;
for ii=1:size(vars,1)
    if strcmp(vars(ii).name, 'iDelete')
        continue
    end
    if vars(ii).size(1)==96
        siz = siz +1;
    end
end
toModify = cell(siz,2);
jj=1;
for ii=1:size(vars,1)
    if strcmp(vars(ii).name, 'iDelete')
        continue
    end
    if vars(ii).size(1)==96
        toModify{jj,1} = vars(ii).name;
        toModify{jj,2} = eval(vars(ii).name);
        toModify{jj,2}(iDelete,:) = [];
        assignin('base',toModify{jj,1},toModify{jj,2})
        jj= jj+1;
    end
end


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

grpmeans = splitapply(@means,mtabData_exp',sInfo.plotGroups);
grpstd = splitapply(@stds,mtabData_exp',sInfo.plotGroups);
grpnstd = grpstd./grpmeans;
grpAvgnstd = mean(mean(grpnstd,1, 'omitnan'),2,'omitnan');

iTimesVSW = (sInfo.matrix=="VSW" & sInfo.sType=="sample");
iTimesASW = (sInfo.matrix=="ASW" & sInfo.sType=="sample");
iCV = (sInfo.matrix=="VSW" & sInfo.sType=="ctrl");
iCA = (sInfo.matrix=="ASW" & sInfo.sType=="ctrl");
iBV = (sInfo.matrix=="VSW" & sInfo.sType=="blank");
iBA = (sInfo.matrix=="ASW" & sInfo.sType=="blank");

%% How much light did each sample absorb, anyway?
% We know irradiance varied spatially, so maybe that accounts for some
% variance in concentrations.
load("../datasets/AbsorbedPhotons.mat")
Dur_sec = seconds(times-min(times));
Dur_sec(1,:) = [];

% The variables in AbsorbedPhotons are in terms of geometry in the SunTest,
% meaning that the first column corresponds to the sample closest to the
% front, which is also the first one removed. I will multiply the values in
% each column by the number of seconds (approximately) that they remained
% in the SunTest.
AbsTot_VSW = Dur_sec.*Wabs_sVSW'.*0.012./1e6; %converts from umol h^-1 L^-1 to mol photons
AbsTot_ASW = Dur_sec.*Wabs_sASW'.*0.012./1e6;

% Okay, but how much does the absorbance vary relative to the time they've
% spent in the solar sim?
% Let's calculate a mean and then see how the different samples compare.

AbsRel_VSW = Wabs_sVSW./mean(Wabs_sVSW,2);
AbsRel_ASW = Wabs_sASW./mean(Wabs_sASW,2);

%% Ingredients for time-course plots.
tBuffer = range(times)*.05;
tRange = [min(times)-tBuffer, max(times)+tBuffer];
tRange_g = rem(datenum(tRange), 1)*24;
tRange_g(2) = tRange_g(2) + 24;
tRange_g = tRange_g - tRange_g(1);
Durations = hours(times-min(times));


%% Glutamine, Glutamate, and GABA
% Figured I'd make a figure. 

% GluBAnames = {"glutamic acid pos", "glutamine neg", "GABA pos"};
% f = figure("Visible","off");
% for ii=1:length(GluBAnames)
%     ind = find(ismember(mtabNames, GluBAnames{ii}));
%     G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
%     m1 = splitapply(@means,mtabData_exp(ind,iTimesASW)',G1)./1000;
%     mn = m1-m1(1);
%     mn(isnan(m1)) = 0;
%     errorbar(Durations, mn,...
%         splitapply(@stds,mtabData_exp(ind,iTimesASW)',G1)./1000,...
%         'Color', col{ii+1}, "LineWidth",1.5); hold on;
%     hold on
% end
% 
% ax = gca;
% set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
% ylabel("C-C_0 (nM)")
% xlabel("Time (h)")
% xlim(tRange_g)
% legend(["glutamate",...
%     "glutamine",...
%     "GABA"], "Location","northeast")
% saveas(f, "../graphs/GLUBA_photo1.png", "png")

%% Histidine, Aspartate, Asparagine
% Figured I'd make a figure. 

% HAAnames = {"histidine pos", "asparagine neg", "aspartate pos"};
% f = figure("Visible","on");
% for ii=1:length(HAAnames)
%     ind = find(ismember(mtabNames, HAAnames{ii}));
%     G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
%     m1 = splitapply(@means,mtabData_exp(ind,iTimesASW)',G1)./1000;
%     mn = m1-m1(1);
%     mn(isnan(m1)) = 0;
%     errorbar(Durations, mn,...
%         splitapply(@stds,mtabData_exp(ind,iTimesASW)',G1)./1000,...
%         'Color', col{ii+1}, "LineWidth",1.5); hold on;
%     hold on
% end
% 
% ax = gca;
% set(ax, "Box", "on", "LineWidth", 2, "FontSize", 14)
% ylabel("C-C_0 (nM)")
% xlabel("Time (h)")
% xlim(tRange_g)
% legend(["histidine",...
%     "asparagine",...
%     "aspartate"], "Location","northeast")
% saveas(f, "../graphs/HAA_photo1.png", "png")


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
            k.kASW(im) = NaN; k.kuASW(im) = NaN; k.kbASW(imr2) = NaN;
            k.rASW(im) = NaN;
        elseif length(tx1) ~= length(m1)
            disp(tx1)
            disp(m1)
            disp("Unequal lengths for "+mtabNames(ii)+" in ASW")
            break
        else
            a0 = [m1(1), 0.1];
            [coef, ~, SSR, ~, ~, ccoef, bcmean] = expregress2(tx1(2:end),m1(2:end),s1(2:end),a0,1);
            k.bASW(im) = bcmean(1);
            k.kASW(im) = -bcmean(2); % regression result is -k
            k.bbASW(im) = ccoef(1,1);
            k.buASW(im) = ccoef(1,2);
            k.kbASW(im) = -ccoef(2,1);
            k.kuASW(im) = -ccoef(2,2);
            k.rASW(im)= SSR;
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
            a0 = [m2(1), 0.1];
            [coef, ~, SSR, ~, ~, ccoef, bcmean] = expregress2(tx2(2:end),m2(2:end),s2(2:end),a0,1);
            k.bVSW(im) = bcmean(1);
            k.kVSW(im) = -bcmean(2); % regression result is -k
            k.bbVSW(im) = ccoef(1,1);
            k.buVSW(im) = ccoef(1,2);
            k.kbVSW(im) = -ccoef(2,1);
            k.kuVSW(im) = -ccoef(2,2);
            k.rVSW(im)= SSR;
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

    variableNames = ["Name","rate_h","intercept_nM","SSR","URate", "LRate","UInt", "LInt"];
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
    VSWrates.halflife_h = 0.693./VSWrates.rate_h;
    ASWrates.halflife_h = 0.693./ASWrates.rate_h;
    fileBase = 'DegRates'; % Set this, don't mess with the automatic date system.
    today = datestr(datetime('now'),'yyyy-mm-dd_');
    NameOfFile = string(['../datasets/',today,fileBase,'.mat']);
    save(NameOfFile, "ASWrates", "VSWrates")
end


%% Addition 28 Jul 2023: Wideband AQY for a few compounds. 
% This section only calculates AQY using the wavelength band for which
% molar extinction is greater than zero. The next section does not truncate
% in this way, and calculates for the full range of what was measured by
% the radiometer.

AQYType = "Truncated"; 

eps = readtable("../datasets/ExtinctionCalcs.xlsx","Sheet","FinalMolarExt","Range", "A2:F553");
lrange = 280:699';
eps_quant = eps(ismember(eps.l, lrange),:);
%eps_quant = flip(eps_quant,1);
load('../datasets/irradiance.mat')
load("../datasets/Napierian.mat")
addpath('quantum_invar_geom\')

if AQYType == "Truncated"
    EA = E_ASW(1:5,1:420);

    t0_a(~ismember(t12_a.l,lrange),:) = [];
    t12_a(~ismember(t0_a.l,lrange),:) = [];
    r = 7; %in mm
    alph = t0_a.sASW;

    mtabNamesQY = {"histidine pos";...
        "glutamic acid pos";...
        "glutamine neg";...
        "tryptophan pos";...
        "kynurenine neg"};
    addpath("./quantum_invar_geom")
    AQYTable = table(mtabNamesQY,zeros(size(mtabNamesQY)),zeros(size(mtabNamesQY)),...
        'VariableNames',{'Metabolite','AQY','MaxLambda'});
    for jj = 1:length(mtabNamesQY)
        ii = find(mtabNames==mtabNamesQY{jj});
        if mtabNamesQY{jj}=="glutamic acid pos"
            continue
        end
        eps_quant_jj = eps_quant{:,jj+1};
        iBadEps = isnan(eps_quant_jj)|eps_quant_jj==0;
        eps_quant_jj(iBadEps) = [];
        EA_jj = EA(:,~iBadEps);
        lrange_jj = lrange(:,~iBadEps);
        LOQ_nM = LOQ_pM./1000;
        alph_jj = alph(~iBadEps);

        G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
        G2 = findgroups(sInfo.matrix(iTimesVSW), sInfo.timePoint(iTimesVSW), sInfo.sType(iTimesVSW));
        m1 = splitapply(@means,mtabData_exp(ii,iTimesASW)',G1)./1000;
        s1 = splitapply(@stds,mtabData_exp(ii,iTimesASW)',G1)./1000;
        m1(isnan(m1)|m1==0,:) = LOQ_nM(ii); s1(isnan(s1)|s1==0,:) = 0.05*LOQ_nM(ii);
        m2 = splitapply(@means,mtabData_exp(ii,iTimesVSW)',G2)./1000;
        s2 = splitapply(@stds,mtabData_exp(ii,iTimesVSW)',G2)./1000;
        m2(isnan(m2)|m2==0,:) = LOQ_nM(ii); s2(isnan(s2)|s2==0,:) = 0.05*LOQ_nM(ii);


        Q = 0.0001; %mean(eps_quant_jj./(2*max(eps_quant_jj))); % First Guess
        Qmax = 1;
        Qmin = 1e-7*Q;

        chiA = @(Q) chisq_quantum_invariant(m1,s1,EA_jj,eps_quant_jj, alph_jj,Q,Durations,lrange_jj,r);
        opts = optimset("PlotFcns", "optimplotfval", "MaxFunEvals", 3e7, "TolX", 1e-6);
        QAest = fmincon(chiA,Q, [],[],[],[],Qmin,Qmax,[],opts);

        AQYTable.AQY(jj) = QAest(1);
        AQYTable.MaxLambda(jj) = max(lrange_jj);

        figure
        plot(lrange_jj, QAest*ones(size(lrange_jj)), "LineWidth",2)
        xlabel("\lambda, nm"); ylabel("\Phi_{ASW}");
        xlim([280 max(lrange_jj)])
        set(gca, "YScale", "log")
        yyaxis right
        plot(lrange_jj, eps_quant_jj, "LineWidth",2)
        ylabel("\epsilon, M^{-1} cm^{-1}")
        set(gca, "YColor", "k", "YScale", "log")
    end

elseif AQYType == "Broad"
    % This section only calculates AQY using the full range of what was measured by
    % the radiometer. (280-699 nm)
    eps_quant = eps(ismember(eps.l, lrange),:);
    eps_quant = flip(eps_quant,1);
    t0_a(~ismember(t12_a.l,lrange),:) = [];
    t12_a(~ismember(t0_a.l,lrange),:) = [];
    r = 7; %in mm
    alph = t0_a.sASW;
    EA = E_ASW(1:5,1:420);
    EV = E_VSW(1:5,1:420);

    mtabNamesQY = {"histidine pos";...
        "glutamic acid pos";...
        "glutamine neg";...
        "tryptophan pos";...
        "kynurenine neg"};
    addpath("./quantum_invar_geom")
    AQYTable = table(mtabNamesQY,zeros(size(mtabNamesQY)),zeros(size(mtabNamesQY)),...
        'VariableNames',{'Metabolite','AQY','MaxLambda'});
    for jj = 1:length(mtabNamesQY)
        ii = find(mtabNames==mtabNamesQY{jj});
        eps_quant_jj = eps_quant{:,jj+1};
        iBadEps = isnan(eps_quant_jj)|eps_quant_jj==0;
        eps_quant_jj(iBadEps) = 0; %Set NaNs to zero
        EA_jj = EA(:,:); %revise these lines so that it's the full range
        lrange_jj = lrange(:,:);
        LOQ_nM = LOQ_pM./1000;

        alph_jj = alph;

        G1 = findgroups(sInfo.matrix(iTimesASW), sInfo.timePoint(iTimesASW), sInfo.sType(iTimesASW));
        G2 = findgroups(sInfo.matrix(iTimesVSW), sInfo.timePoint(iTimesVSW), sInfo.sType(iTimesVSW));
        m1 = splitapply(@means,mtabData_exp(ii,iTimesASW)',G1)./1000;
        s1 = splitapply(@stds,mtabData_exp(ii,iTimesASW)',G1)./1000;
        m1(isnan(m1)|m1==0,:) = LOQ_nM(ii); s1(isnan(s1)|s1==0,:) = 0.05*LOQ_nM(ii);
        m2 = splitapply(@means,mtabData_exp(ii,iTimesVSW)',G2)./1000;
        s2 = splitapply(@stds,mtabData_exp(ii,iTimesVSW)',G2)./1000;
        m2(isnan(m2)|m2==0,:) = LOQ_nM(ii); s2(isnan(s2)|s2==0,:) = 0.05*LOQ_nM(ii);


        Q = mean(eps_quant_jj./(2*max(eps_quant_jj))); % First Guess
        Qmax = 1;
        Qmin = 1e-7*Q;

        chiA = @(Q) chisq_quantum_invariant(m1,s1,EA_jj,eps_quant_jj, alph_jj,Q,Durations,lrange_jj,r);
        opts = optimset("PlotFcns", "optimplotfval", "MaxFunEvals", 3e7, "TolX", 1e-6);
        QAest = fmincon(chiA,Q, [],[],[],[],Qmin,Qmax,[],opts);

        AQYTable.AQY(jj) = QAest(1);
        AQYTable.MaxLambda(jj) = max(lrange_jj);

        figure
        plot(lrange_jj, QAest*ones(size(lrange_jj)), "LineWidth",2)
        xlabel("\lambda, nm"); ylabel("\Phi_{ASW}");
        xlim([280 max(lrange_jj)])
        set(gca, "YScale", "log")
        yyaxis right
        plot(lrange_jj, eps_quant_jj, "LineWidth",2)
        ylabel("\epsilon, M^{-1} cm^{-1}")
        set(gca, "YColor", "k", "YScale", "log")
    end
end


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
    ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", col{3});
    if ii~=1
        ub(ii).HandleVisibility = "off";
    end
    hold on
end
ii = find(mtabNames=="tryptophan pos");
scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 50, "filled",...
    'Color', col{1});
scBA = scatter(repelem(Durations(1),1,3), ba, 50, col{2},...
    "filled", "Marker","v");
scCA = scatter(repelem(Durations(6),1,3), ca,50,col{4},...
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
    "High Standard (" +string(round(MaxStd_pM(ii))) + " pM)"], "Location","northeast")

% subplot(1,2,2)
% for ii = 1:5
%     dCdt = @(t,Ct) -Ct.*int_solver(EV(ii,:)',eps_quant,QVest,lrange);
%     tc = 0:1/3600:Durations(ii+1);
%     [ty, yhat] = ode45(dCdt,tc,m2(1));
%     ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", col{3});
%     if ii~=1
%         ub(ii).HandleVisibility = "off";
%     end
%     hold on
% end
% ii = find(mtabNames=="tryptophan pos");
% scV = scatter(durVSW,mtabData_exp_nM(ii,iTimesVSW), 50, "filled",...
%     'Color', col{1});
% scBV = scatter(repelem(Durations(1),1,6), bv, 50,col{2},...
%     "filled", "Marker","v");
% scCV = scatter(repelem(Durations(6),1,3), cv,50,col{4},...
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
%     "High Standard (" +string(round(MaxStd_pM(ii))) + " pM")], "Location","northeast")

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
    ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", col{3});
    if ii~=1
        ub(ii).HandleVisibility = "off";
    end
    hold on
end
ii = find(mtabNames=="tryptophan pos");
scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 50, "filled",...
    'Color', col{1});
scBA = scatter(repelem(Durations(1),1,3), ba, 50, col{2},...
    "filled", "Marker","v");
scCA = scatter(repelem(Durations(6),1,3), ca,50,col{4},...
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
    "High Standard (" +string(round(MaxStd_pM(ii))) + " pM)"], "Location","northeast")

saveas(f, "../graphs/trpPlot/trpQY.png", "png")
%close(f)

%% AQY WIDEBAND CALCULATION
addpath("./quantum_invar_geom")
Q = mean(eps_quant./(2*max(eps_quant))); % First Guess 
Qmax = 1;
Qmin = 1e-7*Q;

chiA = @(Q) chisq_quantum_invariant(m1,s1,EA,eps_quant,alph,Q,Durations,lrange,r);
opts = optimset("PlotFcns", "optimplotfval", "MaxFunEvals", 3e7, "TolX", 1e-6);
QAest = fmincon(chiA,Q, [],[],[],[],Qmin,Qmax,[],opts);

figure
plot(lrange, QAest*ones(size(lrange)), "LineWidth",2)
xlabel("\lambda, nm"); ylabel("\Phi_{ASW}");
xlim([280 407])
set(gca, "YScale", "log")
yyaxis right
plot(lrange, eps_quant, "LineWidth",2)
ylabel("\epsilon, M^{-1} cm^{-1}")
set(gca, "YColor", "k", "YScale", "log")

%% Decay Estimation with AQY
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
    ub(ii) = plot(ty, yhat, "--", "LineWidth", 2, "Color", col{3});
    if ii~=1
        ub(ii).HandleVisibility = "off";
    end
    hold on
end
ii = find(mtabNames=="tryptophan pos");
scA = scatter(durASW,mtabData_exp_nM(ii,iTimesASW), 50, "filled",...
    'Color', col{1});
scBA = scatter(repelem(Durations(1),1,3), ba, 50, col{2},...
    "filled", "Marker","v");
scCA = scatter(repelem(Durations(6),1,3), ca,50,col{4},...
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
    "High Standard (" +string(round(MaxStd_pM(ii))) + " pM)"], "Location","northeast")

saveas(f, "../graphs/trpPlot/trpQY_invariant.png", "png")
%close(f)