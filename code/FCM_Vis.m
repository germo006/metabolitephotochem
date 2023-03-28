function f = FCM_Vis
% Since Excel sucks for this, I'm coding this to read in the processed flow
% cytometry data and graph it.
load('AlbumMaps.mat')
FCM = readtable("../datasets/FCMData.csv");
FCM(isnan(FCM.t_h),:) = [];
FCM(FCM.Matrix=="MQ" | FCM.Matrix == "MQs", :) = [];

f = figure("Name","Flow Cytometry Results");
ASW_dark = FCM.Matrix == "ASW" & FCM.ctrl;
sASW = FCM.Matrix == "ASWs" & ~FCM.ctrl;
sASW_dark = FCM.Matrix == "ASWs" & FCM.ctrl;

VSW_dark = FCM.Matrix == "VSW" & FCM.ctrl;
sVSW = FCM.Matrix == "VSWs" & ~FCM.ctrl;
sVSW_dark = FCM.Matrix == "VSWs" & FCM.ctrl;

p1 = plot(FCM.t_h(sASW),FCM.subBlank(sASW), "Color", chainsaw{1}, "LineWidth",2);
hold on
p2 = plot(FCM.t_h(sVSW),FCM.subBlank(sVSW), "Color", chainsaw{3}, "LineWidth",2);
s1 = scatter(FCM.t_h(sASW_dark),FCM.subBlank(sASW_dark), 55, chainsaw{1}, "filled","s");
s2 = scatter(FCM.t_h(sVSW_dark),FCM.subBlank(sVSW_dark), 55, chainsaw{3}, "filled","s");
s3 = scatter(FCM.t_h(ASW_dark),FCM.subBlank(ASW_dark), 45, chainsaw{1}, "filled","v");
s4 = scatter(FCM.t_h(VSW_dark),FCM.subBlank(VSW_dark), 45, chainsaw{3}, "filled","^");

s1.MarkerEdgeColor = "k"; s2.MarkerEdgeColor = "k";
s3.MarkerEdgeColor = "k"; s4.MarkerEdgeColor = "k";


xlabel("time, h")
ylabel("cell-like objects mL^{-1}")
legend({"sASW","sVSW", "sASW dark", "sVSW dark","ASW dark","VSW dark"}, "location", "northwest")
xlim([-0.5,12.5])