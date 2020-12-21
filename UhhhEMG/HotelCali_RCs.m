%% Plot RCs for LSP02b and LSP05
%VL-RF: '14274d' BF-ST: 'ea6c20'  LG-MG: '62799d' TA: 'f39d73'}
subjectName = 'LSP05';
switch subjectName
    case 'LSP05'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP05\LSP05_response_meta.mat')
% %         elecs = [1 2 3 4 7 8 9 10 11 12 13 14 15 16 18 19 20 23 24 25 26 27 28 29 30 31];
        elecs = [4 7 12 20 29];
        response5 = response;
        muscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'SO/TA', 'MG', 'LG',};

        
    
    case 'LSP02b'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM','Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left Ham', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP02b\LSP02b_response_meta.mat');
% %          elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34];
         elecs = [3 9 12 33];
% %         msubset = [1:16];
        response2b = response;
        muscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'Ham/TA', 'MG', 'LG',};
end
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];




%% Plot Comparison RC
% % msubset = [11 12 16];
% % elecs = [1 8 16 18 23 29];
% % for i = 1:length(msubset)
% %     m = msubset(i);
% %     h(i) = subplot(1,length(msubset),i);
% %     for r = 1:length(elecs)
% %         e = find([response.elec] == elecs(r));
% %         tmp(e) = max(response(e).muscle(m).p2pResponse);
% %     end
% %     maxR(m) = max(tmp);
% % end


hF = figure; maximize;
tiledlayout(length(elecs), 2, 'TileSpacing','compact');


msubset = [1 2 3 7 8]; % Flexors
C = linspecer(length(msubset));

for r = 1:length(elecs)
    h(r) = nexttile(1+2*(r-1));
    e = find([response.elec] == elecs(r));
    l = 0;
    for i = 1:length(msubset)
        m = msubset(i);
%         h(i) = subplot(1,length(msubset),i);
        plot(response(e).StimAmps/1000, response(e).muscle(m+8).p2pResponse,'Color',C(i,:), 'LineWidth',2);
        hold on;
        if strcmp(subjectName, 'LSP05') && m <8
            plot(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse,'Color',C(i,:), 'LineWidth',2, 'LineStyle', '--');
        end
        l = 1+2*(i-1);
        leg_labs(l) = strcat(muscles(m),'- Residual');
        leg_labs(l+1) = strcat(muscles(m),'- Intact');
    end
    title(['Electrode ' num2str(elecs(r))]);
    ylabel('EMG (uV)');
    xlabel('StimAmp (mA)');
        
    legend(leg_labs,'Location','westoutside');
    box off;
end
linkaxes([h(:)],'xy');

clear msubset C;
msubset = [4 5 6]; % Flexors
C = linspecer(length(msubset));

for r = 1:length(elecs)
    h(r) = nexttile(2+2*(r-1));
    e = find([response.elec] == elecs(r));
    l = 0;
    for i = 1:length(msubset)
        m = msubset(i);
%         h(i) = subplot(1,length(msubset),i);
        plot(response(e).StimAmps/1000, response(e).muscle(m+8).p2pResponse,'Color',C(i,:), 'LineWidth',2);
        hold on;
        plot(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse,'Color',C(i,:), 'LineWidth',2, 'LineStyle', '--');
        l = 1+2*(i-1);
        leg_labs(l) = strcat(muscles(m),'- Residual');
        leg_labs(l+1) = strcat(muscles(m),'- Intact');
    end
    title(['Electrode ' num2str(elecs(r))]);
    ylabel('EMG (uV)');
    xlabel('StimAmp (mA)');
        
    legend(leg_labs,'Location','eastoutside');
    box off;
end
linkaxes([h(:)],'xy');
tit = sgtitle(['Recruitment Curves ' subjectName]);
% % tit = sgtitle(['Recruitment Curves by Muscle for Electrodes ' num2str(elecs(1)) '-' num2str(elecs(end))]);
% % saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end)) '.png'])
% % savefig(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\P2P_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end))])
% % 
% % close all;