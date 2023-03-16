%% Dealing with Radiometry
clear 
clc
[l_SunTest,W,E] = loadradiometry;
load("AlbumMaps.mat", "chainsaw")
figure
for ii=1:3
    for jj = 1:4
        subplot(3,4,(ii-1)*4+jj)
        plot(l_SunTest, reshape(E(ii,jj,:),420,1), "LineWidth",2,...
            "Color",chainsaw{3})

        if jj==1
            ylabel("\mumol s^{-1} m^{-2}")
        end
        if ii==3
            xlabel("\lambda, nm")
        end
        xlim([280,699]); ylim([0,15]);
    end
end
clear ii jj
%% Well, let's see if we can interpolate the spectra across the vertical 
% swaths where our samples were. 
close all

E0_ASW = reshape(E(:,1,:),3,length(l_SunTest));
W0_ASW = reshape(E(:,1,:),3,length(l_SunTest));
E0_VSW = reshape(E(:,3,:),3,length(l_SunTest));
W0_VSW = reshape(E(:,3,:),3,length(l_SunTest));

E_ASW = zeros(7,length(l_SunTest));
W_ASW = E_ASW;
E_VSW = E_ASW;
W_VSW = E_ASW;

sloc = [1:2:13].'./14;
lloc = [1:2:5].'./6;

%Linear interpolation/extrapolation
E_ASW(1:3,:) = E0_ASW(2,:) - ((E0_ASW(2,:) - E0_ASW(1,:))./(lloc(2)-lloc(1))).*(lloc(2) - sloc(1:3));
E_ASW(4,:) = E0_ASW(2,:);
E_ASW(5:7,:) = E0_ASW(2,:) + ((E0_ASW(3,:) - E0_ASW(2,:))./(lloc(3)-lloc(2))).*(sloc(5:7) - lloc(2));

W_ASW(1:3,:) = W0_ASW(2,:) - ((W0_ASW(2,:) - W0_ASW(1,:))./(lloc(2)-lloc(1))).*(lloc(2) - sloc(1:3));
W_ASW(4,:) = W0_ASW(2,:);
W_ASW(5:7,:) = W0_ASW(2,:) + ((W0_ASW(3,:) - W0_ASW(2,:))./(lloc(3)-lloc(2))).*(sloc(5:7) - lloc(2));

E_VSW(1:3,:) = E0_VSW(2,:) - ((E0_VSW(2,:) - E0_VSW(1,:))./(lloc(2)-lloc(1))).*(lloc(2) - sloc(1:3));
E_VSW(4,:) = E0_VSW(2,:);
E_VSW(5:7,:) = E0_VSW(2,:) + ((E0_VSW(3,:) - E0_VSW(2,:))./(lloc(3)-lloc(2))).*(sloc(5:7) - lloc(2));

W_VSW(1:3,:) = W0_VSW(2,:) - ((W0_VSW(2,:) - W0_VSW(1,:))./(lloc(2)-lloc(1))).*(lloc(2) - sloc(1:3));
W_VSW(4,:) = W0_VSW(2,:);
W_VSW(5:7,:) = W0_VSW(2,:) + ((W0_VSW(3,:) - W0_VSW(2,:))./(lloc(3)-lloc(2))).*(sloc(5:7) - lloc(2));

UVB = 280:315; UVBi = find(ismember(l_SunTest, UVB'));
UVA = 315:400; UVAi = find(ismember(l_SunTest, UVA'));
Vis = 380:699; Visi = find(ismember(l_SunTest, Vis'));

W_ASW_UVB = trapz(W_ASW(:,UVBi),2);
W_ASW_UVA = trapz(W_ASW(:,UVAi),2);
W_ASW_Vis = trapz(W_ASW(:,Visi),2);

W_VSW_UVB = trapz(W_VSW(:,UVBi),2);
W_VSW_UVA = trapz(W_VSW(:,UVAi),2);
W_VSW_Vis = trapz(W_VSW(:,Visi),2);

save("irradiance.mat")