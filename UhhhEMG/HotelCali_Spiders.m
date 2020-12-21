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
        elecs = [7 11 18 16 26];
        elecsRast = [1 4 7 9 11 18 13 20 16 23 26 29];%[4 7 12 20 29];
        response5 = response;
        musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
        musclesR = {'VM', 'RF', 'VL', 'BF', 'ST', 'SO', 'MG', 'LG',};

        
    
    case 'LSP02b'
        Muscle_ID = {'Right VM', 'Right RF', 'Right VL', 'Right BF', 'Right ST', ...
        'Right TA', 'Right SO', 'Right LG','Left VM','Left RF','Left VL',...
        'Left BF', 'Left ST', 'Left Ham', 'Left MG', 'Left LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP02b\LSP02b_response_meta.mat');
% %          elecs = [1 2 3 4 5 6 7 8 9 10 11 12 33 34];
         elecs = [3 9 12 34];
         elecsRast = [1 3 6 9 12 33 34];
% %         msubset = [1:16];
        response2b = response;
        musclesL = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
        musclesR = {'VM', 'RF', 'VL', 'BF', 'ST', 'Ham', 'MG', 'LG',};
end
cmap = [0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.0784    0.1529    0.3020; 0.9176    0.4235    0.1255; 0.9176    0.4235    0.1255; 0.9529    0.6157    0.4510; 0.3843    0.4745    0.6157; 0.3843    0.4745    0.6157];

%% Get Max for Norm
msubset = [1:16];
for i = 1:length(msubset)
    m = msubset(i);
    h(i) = subplot(1,length(msubset),i);
    for r = 1:length(elecs)
        e = find([response.elec] == elecs(r));
        tmp(e) = max(response(e).muscle(m).p2pResponse);
    end
    maxR(m) = max(tmp);
end
%%
div = 1000;
msubset = [1:16];
for r = 1:length(response)
    tIdx = min([response(r).muscle.threshIdx])+1;
    for m = 1:length(msubset)
        if isnan(response(r).muscle(msubset(m)).threshIdx)
            threshold(r,m) = nan;
        else
            threshold(r,m) = response(r).muscle(msubset(m)).p2pResponse(tIdx);
        end
    end
end
% % threshold(find(isnan(threshold))) = 0;
% % if strcmp(subjectName, 'LSP05')
% %     threshold(:,8) = 0;
% % end

% % hS1 = figure; maximize;
spdrfill = 'on';
% % tiledlayout(length(elecs),2, 'TileSpacing','compact');
spdrmax(1:8) = round(max(max(threshold)));
spdrmin(1:8) = 0;%;0.01*max(max(threshold)); 
spdriter = 4; 

C = linspecer(length(elecsRast));

for i = 1:length(elecs)
    figure; %maximize; tiledlayout('flow', 'TileSpacing','compact');
    e = find([response.elec] == elecs(i));
    c = find([elecsRast] == elecs(i));
    
    spdrmax(1:8) = max(threshold(e,9:16));
    spdrmin(1:8) = .05*spdrmax(1:8);
    nexttile;
    spider_plot(threshold(e,9:16),...
        'AxesLimits', [spdrmin; spdrmax ],... %maxR(9:16)
        'AxesLabels', musclesL,'AxesPrecision', 0,...
        'FillOption', spdrfill,'FillTransparency', 0.5, ...
        'AxesInterval', spdriter, 'Direction', 'clockwise',...
        'Color',  C(c,:));
% %     
% %     nexttile;
% %     nexttile(2+2*(i-1));

% %     spider_plot(threshold(e,1:8),...
% %         'AxesLimits', [spdrmin; maxR(9:16)],...
% %         'AxesLabels', musclesR,'AxesPrecision', 0,...
% %         'FillOption', spdrfill,'FillTransparency', 0.5, ...
% %         'AxesInterval', spdriter, 'Direction', 'counterclockwise',...
% %         'Color', C(c,:));
    


%     title('Residual Limb (left)');
    % Legend properties
%     legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
    tmptit = sgtitle([subjectName ' - E' num2str(elecs(i))]);
    set(tmptit, 'Interpreter', 'none');
end

% % for i = 1:length(elecs)
% %     nexttile(1+2*(i-1));
% %     e = find([response.elec] == elecs(i));
% %     c = find([elecsRast] == elecs(i));
% %     spider_plot(threshold(e,9:16),...
% %         'AxesLimits', [spdrmin; spdrmax],...
% %         'AxesLabels', musclesL,'AxesPrecision', 0,...
% %         'FillOption', spdrfill,'FillTransparency', 0.2, ...
% %         'AxesInterval', spdriter, 'Direction', 'clockwise',...
% %         'Color',  C(c,:));
% %     %         'AxesLimits', [spdrmin(1:length(msubset)); spdrmax(1:length(msubset))],...
% % %     title('Residual Limb (left)');
% %     % Legend properties
% % %     legend(arrayfun(@num2str,elecs,'UniformOutput',false), 'Location', 'westoutside', 'Orientation', 'vertical');
% %      
% % end
% %     tmptit = sgtitle(['Thresholds for ' num2str(elecs(1)) '-' num2str(elecs(end)) ' (1Hz, 500us)']);
% %     set(tmptit, 'Interpreter', 'none');