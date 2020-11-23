%% Quick Histogram of Muscles active at lowest threshold

load('D:\DATA\UH3 testing\LSP05\data_gen\threshes\1Hz-500us_thresholds_srtd.mat')
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};
RatThresh = zeros(size(threshes));
for i = 1:length(threshes)
    globaMinIndexes = find(threshes(i,:) == min(threshes(i,:)));
    RatThresh(i, [globaMinIndexes]) = 1;
end

Z = RatThresh(:,9:end);
bar(sum(Z)/26)
ylim([0,1])
xticklabels(musclesL)
xlabel('Residual Limb');
ylabel('Response Fraction');
title('Percent of Time Muscles Respond at Lowest Threshold');