setDefaultFigs
load('F:/Noah Germolus/Documents/MATLAB/NoahMaps/AlbumMaps.mat', "soft")
x = 0:0.1:24;
y = 1 + sin(x.*(pi/12) + pi);
subplot(4,1,1)
area(x,y,"FaceColor",soft{2}, "FaceAlpha",0.5,"EdgeColor","none")
ax = gca;

ax.YLim = [0,2.5]; ax.XLim = [0,24];
ylabel("[Metabolite], nM", "Interpreter", "latex") 
ax.TickLabelInterpreter = "latex";
yr = cos(x.*(pi/12) + pi).*(pi/12);
subplot(4,1,2)
area(x,yr,"FaceColor",soft{2}, "FaceAlpha",0.5,"EdgeColor","none")
ax = gca;
ax.YLim = [-1.5,1.5]; ax.XLim = [0,24];
ax.TickLabelInterpreter = "latex";
ylabel("$\frac{d[Metabolite]}{dt}$, nM h$^{-1}$", "Interpreter","latex")
legend("k(t), total", "Interpreter", "latex")
xp = [zeros(1,60),0.1:0.1:14,zeros(1,41)];
yrp = -0.5*sin(xp.*(pi/14));
subplot(4,1,3)
area(x,yrp,"FaceColor",soft{4}, "FaceAlpha",0.5,"EdgeColor","none")
hold on
xz = [0:0.1:3, 3.*ones(1,180), 3:-0.1:0.1];
yrz = 0.75*cos(xz.*(pi/6));
area(x,yrz,"FaceColor",soft{3}, "FaceAlpha",0.5,"EdgeColor","none")
yrm = 0.2.*cos(x.*(pi/12));
area(x,yrm,"FaceColor",soft{1}, "FaceAlpha",0.5,"EdgeColor","none")
legend(["k$_{photo}$","k$_{zoop}$","k$_{mix}$"], "Location", "southeast", "Interpreter", "latex");
yrb = yr - yrp - yrz - yrm;
ax = gca;
ax.TickLabelInterpreter = "latex";
ax.YLim = [-1,1]; ax.XLim = [0,24];
ylabel("$\frac{d[Metabolite]}{dt}$, nM h$^{-1}$", "Interpreter","latex")
%ax.YLabel.String = "^{d[Metabolite]}/_{dt}, nM h^{-1}";
subplot(4,1,4)
area(x',yrb', "FaceColor",soft{5},"FaceAlpha",0.5,"EdgeColor","none")
%colororder([chainsaw{5};chainsaw{2};chainsaw{3};chainsaw{1}]);
ax = gca;
ax.TickLabelInterpreter = "latex";
ax.YLim = [-2.5,1.5]; ax.XLim = [0,24];
ylabel("$\frac{d[Metabolite]}{dt}$, nM h$^{-1}$", "Interpreter","latex")
%ax.YLabel.String = "^{d[Metabolite]}/_{dt}, nM h^{-1}"; 
xlabel(ax, "time, h", "Interpreter","latex")
legend("k$_{bio}$ = k(t) - k$_{photo}$ - k$_{zoop}$ - k$_{mix}$", "Location","south", "Interpreter", "latex")