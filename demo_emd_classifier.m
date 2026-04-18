% DEMO_EMD_CLASSIFIER  Visualise the classifier on synthetic curves.
 
rng(42);
shifts = (-20:20)';
 
% Successful match: parabola with a deep central dip
emd_success = 480 + 0.45*shifts.^2 - 120*exp(-shifts.^2/30) + randn(size(shifts))*8;
 
% Failed match: gentle bowl, no central dip
emd_fail = 600 + 40*exp(-shifts.^2/800) + 0.15*shifts.^2 + randn(size(shifts))*8;
 
[m1, r1] = helper.emd_match_classifier(emd_success);
[m2, r2] = helper.emd_match_classifier(emd_fail);
 
figure('Name','EMD Match Classifier Demo','Color','w','Position',[100 100 900 380]);
 
subplot(1,2,1);
plot(shifts, emd_success, 'Color',[0.11 0.62 0.46], 'LineWidth',2);
title(sprintf('Curve 1 — ratio=%.3f → %s', r1, result_str(m1)), ...
    'FontSize',13, 'FontWeight','normal');
xlabel('Time shift'); ylabel('EMD');
xline(0,'--','Color',[0.5 0.5 0.5],'Alpha',0.6,'LineWidth',1.2);
grid on; box off;
 
subplot(1,2,2);
plot(shifts, emd_fail, 'Color',[0.85 0.35 0.19], 'LineWidth',2);
title(sprintf('Curve 2 — ratio=%.3f → %s', r2, result_str(m2)), ...
    'FontSize',13, 'FontWeight','normal');
xlabel('Time shift'); ylabel('EMD');
xline(0,'--','Color',[0.5 0.5 0.5],'Alpha',0.6,'LineWidth',1.2);
grid on; box off;
 
sgtitle('EMD curve shape classifier', 'FontSize',14);
 
function s = result_str(is_match)
    if is_match, s = 'MATCH'; else, s = 'NO MATCH'; end
end
 