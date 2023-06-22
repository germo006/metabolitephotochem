% This is to make a graph that compares the irradiance in the SunTest with
% simulated values from SMARTS at the BATS site.

Irrad = readtable("../datasets/SunTest_v_SMARTS_Jul11_BDA.xlsx", "Sheet", ...
    "Comparison");
STstd = std(Irrad{:,2:13}, [], 2);
STmean = Irrad.mean_ST;
SMARTS = Irrad.SMARTS_norm_dir;
l = Irrad.lambda;

load("AlbumMaps.mat", "chainsaw")
setDefaultFigs

Cf = figure;
A1 = area(l, [STmean-STstd,STmean, STmean+STstd]);
A1(1).FaceColor = "none";
A1(2).FaceColor = chainsaw{4};
A1(3).FaceColor = chainsaw{4};
A1(2).FaceAlpha = 0.3;
A1(3).FaceAlpha = 0.3;
A1(1).EdgeColor = chainsaw{4};
A1(2).EdgeColor = chainsaw{4};
A1(3).EdgeColor = chainsaw{4};
A1(3).LineStyle = "--";
A1(1).LineStyle = "--";
A1(1).HandleVisibility = "off";
A1(2).HandleVisibility = "off";

hold on

L1 = plot(l,SMARTS,"LineWidth",1.5,"Color",chainsaw{2});

ax = gca;
ax.YLabel.String = "Direct Normal Irradiance, W m^{-2} nm^{-1}";
ax.XLabel.String = "\lambda, nm";
ax.XLim = [300, 700];

legend({"SunTest (measured +/- \sigma)", "SMARTS (simulated, 11 July at BATS)"},"Location","northwest")

saveas(Cf, "../graphs/STvSMARTS.png", "png")