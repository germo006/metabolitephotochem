% This is meant to calculate an approximate terrigeneous contribution to
% the VSW DOM, based on the method of Fichot and Benner (2012)...sort of.

load('../datasets/absorbance.mat')

means_BLsub_a = a_means_BLsub; % Copy the variable of abs coefficients.
wavelengthRange = 275:295';
i_275_295 = ismember(means_BLsub_a.l, wavelengthRange);
a275_295 = means_BLsub_a(i_275_295,:);
a275_295 = flip(a275_295);
a0 = min(a275_295(:,1:end-1));
a275_295_subBL = a275_295;
a275_295_subBL(:,1:end-1) =  a275_295_subBL(:,1:end-1)-a0;
a275_295_subBL_log = a275_295_subBL;
a275_295_subBL_log(:,1:end-1) = log(a275_295_subBL_log(:,1:end-1));
a275_295_subBL_log_dM = a275_295_subBL_log{:,1:12};
a275_295_subBL_log_dM(isinf(a275_295_subBL_log_dM)) = NaN;
a275_295_subBL_log{:,1:12} = a275_295_subBL_log_dM;

clear i_275_295 a275_295_subBL_log_dM

% for ii=1:12
%     plot(a275_295_subBL_log.l, a275_295_subBL_log.(ii))
%     hold on
% end
%%

coefs = array2table(zeros(2,12),'VariableNames',a275_295_subBL_log.Properties.VariableNames(1:12));
x = wavelengthRange';
for ii = 1:12
    y = a275_295_subBL_log.(ii);
    ir = isnan(y);
    y(ir,:) = []; 
    xii = x(~ir);
    xii = [ones(size(xii)),xii];
    coefs{:,ii} = xii\y;
end
S_275_295 = coefs(2,:);
S_275_295{:,:} = -S_275_295{:,:};

TDLP9C = @(S) exp(3.172-(267.566.*S)) + exp(0.228.*S) -0.953*exp(S);
TDLP9C_est = S_275_295;
TDLP9C_est{:,:} = TDLP9C(TDLP9C_est{:,:});