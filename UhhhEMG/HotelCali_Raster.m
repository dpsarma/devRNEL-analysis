%% Plot Comparison RC Stacked
% % clear all;
% % 
% % load('D:\FigRescources\UH3\LSP05\rehash\responses_500us_1Hz_allElecsRC.mat')
% % response = response2; clear response2;

%% Plot RCs for LSP02b and LSP05
%VL-RF: '14274d' BF-ST: 'ea6c20'  LG-MG: '62799d' TA: 'f39d73'}
subjectName = 'LSP02b';
switch subjectName
    case 'LSP05'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM', 'Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left TA', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP05\LSP05_response_meta.mat')
% %         elecs = [1 2 3 4 7 8 9 10 11 12 13 14 15 16 18 19 20 23 24 25 26 27 28 29 30 31];
        elecs = [1 4 7 9 11 18 13 20 16 23 26 29]; %[4 7 12 20 29];
        response5 = response;
        labelmuscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'SO/TA', 'MG', 'LG',};

        
    
    case 'LSP02b'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM','Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left Ham', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP02b\LSP02b_response_meta.mat');
% %          elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34];
         elecs = [1 3 6 9 12 33 34];
% %         msubset = [1:16];
        response2b = response;
        labelmuscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'Ham/TA', 'MG', 'LG',};
end
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];

% % 
% % labelmuscles = { 'VM', 'VL', 'RF', 'BF', 'ST', 'Ham', 'MG', 'LG'}; %Ham vs TA
msubset = [11 12 16 14]; %[9 11 10 12 13 14 15 16]; %[11 12 16 14];[3 4 8 6];
% % etmp= [response.elec];
% % elecs = etmp(1:4:end);
% % elecs = [response(15:end).elec];%[response.elec]; %[1 4 7 9 11 18 13 20 16 23 26 29]; %
C = linspecer(length(elecs));

for i = 1:length(msubset)
    m = msubset(i);
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        tmp(e) = max(response(e).muscle(m).p2pResponse);
    end
    musc_map(m,:) = tmp;
    maxR(m) = max(tmp);
end

hF = figure;% maximize;
for i = 1:length(msubset)
    m = msubset(i);
    
    for r = 1:length(elecs)
        h(i,r) = subplot(length(elecs),length(msubset),length(msubset)*(r-1)+i);
        e = find([response.elec] == elecs(r));
        p = area(response(e).StimAmps/1000, response(e).muscle(m).p2pResponse); %/maxR(m)
        p.FaceColor = C(r,:);
        hold on;
        if i==1
            ylabel('EMG (uV)');
        end
        if r==length(elecs)
            xlabel('StimAmp (mA)');
        end
        if r==1
            title(response(e).muscle(m).Muscle_ID);
        end
%         legend(['E' num2str(elecs(r))],'Location','northwest');
        box off
        grid off
        set(gca,'Visible','off')
        axis off;
        h(i,r).Color = 'none';
    end
    
     
    
end
% % for i = 1:length(msubset)
% % linkaxes([h(i,:)],'xy');
% % end
linkaxes([h],'xy');
tit = sgtitle(['Recruitment Curves by Muscle for Electrodes ' num2str(elecs(1)) '-' num2str(elecs(end))]);
set(0,'defaultAxesFontSize',12)
% % saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\Stacked_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end)) '.png'])
% % savefig(hF,['D:\FigRescources\UH3\LSP05\rehash\PAD\Stacked_RCs_e' num2str(elecs(1)) '-' num2str(elecs(end))])
% % 
% % close all;


