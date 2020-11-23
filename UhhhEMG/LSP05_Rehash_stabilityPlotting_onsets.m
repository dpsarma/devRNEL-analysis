%% Stability and Onset plotting for 9 and 15
load('D:\FigRescources\UH3\LSP05\rehash\response_stability_e9.mat');
response9 = response([1 2]);
load('D:\FigRescources\UH3\LSP05\rehash\response_stability_e16.mat');
response16 = response([1 2]);
clear response;
musclesL = {'VM', 'VL', 'RF', 'BF', 'ST', 'TA', 'MG', 'LG',};
msubset = [9 12 14 16];
mlab = [musclesL(msubset-8)];

%% Make Onset and 

elecs = unique([response9(:).elec; response16(:).elec]);
response = [response9 response16];
for iD = 1:length(response)
    for iM = 1 : length(msubset)
        xinfo = response(iD).muscle(msubset(iM)).p2pResponse;
        ampinfo = [ 0 response(iD).StimAmps/1000];
        aucData(iD,iM) = sum(xinfo.*diff(ampinfo));
        for iter = 1:length(response(iD).StimAmps)-1
            xes = response(iD).muscle(msubset(iM)).p2pResponse(1:iter:end);
            amps = [ 0 response(iD).StimAmps(1:iter:end)/1000];
            aucds(iter) = sum(xes.*diff(amps));
        end
        aucerr(iD,iM) = std( aucds ) / sqrt( length(aucds) );
    end
end

for iD = 1:length(response)
    for iM = 1 : length(msubset)
        onsets(iD,iM) = mean(response(iD).muscle(msubset(iM)).onset(response(iD).muscle(msubset(iM)).threshIdx:end))*1000;
        tmpOns = response(iD).muscle(msubset(iM)).onset(response(iD).muscle(msubset(iM)).threshIdx:end)*1000;
        Onserror(iD,iM) = std( tmpOns ) / sqrt( length(tmpOns) );
        threshes(iD,iM) = (response(iD).muscle(msubset(iM)).threshold)/1000;
    end
end

% % hS6 = figure;maximize;
% % tiledlayout(1,2);
% % nexttile
% % bar([1:length(msubset)],aucData);ylabel('Cum EMG (uV*mA)');
% % legend({response.setDescrpt}, 'Location', 'southoutside', 'Orientation', 'horizontal','Interpreter','none');
% % title('Gross Recruitment');
% % xticklabels(mlab);
% % 
% % nexttile
% % bar([1:length(msubset)],onsets);ylabel('Response Time (ms)');
% % xticklabels(mlab);
% % legend({response.setDescrpt}, 'Location', 'southoutside', 'Orientation', 'horizontal','Interpreter','none');
% % title('Response Onset');
% % hold on;
% % XData = repmat([1:length(msubset)],6,1);
% % er = errorbar(XData,onsets,Onserror,Onserror);
% % er.Color = [0 0 0];                            
% % er.LineStyle = 'none';

% % nexttile
% % boxplot(onsets);
% % ylabel('Response Time (ms)')
% % xticklabels(musclesL);
% % title('Response Onset');

% % saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Stability\AUC_Onset' num2str(elecs(1)) '.png'])
% % saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Stability\AUC_Onset' num2str(elecs(1)) '.svg'])

%% Plot Superbar
C = linspecer(4);
hS6 = figure;maximize;
tiledlayout(1,2);
nexttile
% % superbar(aucData','E', 0.5*aucerr','BarFaceColor',permute(C, [3 1 2]));
% % ylabel('Cum EMG (uV*mA)');
superbar(threshes','BarFaceColor',permute(C, [3 1 2]));
ylabel('Threshold Stimulation Amplitude (mA)'); ylim([0 6]);
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:4]);
xticklabels(mlab);
legend({response.setDescrpt}, 'Location', 'southoutside', 'Orientation', 'horizontal','Interpreter','none');
title('Threshold Across Days');
hold on;
XData = repmat([1:length(msubset)],6,1);

nexttile;
superbar(onsets','E', Onserror','BarFaceColor',permute(C, [3 1 2]));
ylabel('Response Time (ms)');
% tmp = [' ', mlab, ' '];
% xticklabels(tmp);
xticks([1:4]);
xticklabels(mlab);
legend({response.setDescrpt}, 'Location', 'southoutside', 'Orientation', 'horizontal','Interpreter','none');
title('Response Onset');
hold on;
XData = repmat([1:length(msubset)],6,1);

% % saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Stability\AUC_Onset' num2str(elecs(1)) '-' num2str(elecs(end)) '.png'])
% % saveas(hS6,['D:\FigRescources\UH3\LSP05\rehash\Stability\AUC_Onset' num2str(elecs(1)) '-' num2str(elecs(end)) '.svg'])
% % 
