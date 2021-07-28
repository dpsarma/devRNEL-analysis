%% Prevalence heatmaps 
subjectName = 'LSP02b';
switch subjectName
    case 'LSP05'
        Muscle_ID = {'Intact VM', 'Intact RF', 'Intact VL', 'Intact BF', 'Intact ST', ...
        'Intact TA', 'Intact SO', 'Intact LG','Resid. VM', 'Resid. RF','Resid. VL',...
        'Resid. BF', 'Resid. ST', 'Resid. TA', 'Resid. MG', 'Resid. LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP05\LSP05_response_meta.mat')
        elecs = [1:4, 7:11, 18 12 19 13 20 14 15 16 23 24 25 26 27 28 29 30 31];
        response5 = response;
        muscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG',};
        sensethresh = [4000 4000 3500 3000 3000 2000 2000 2500 2500 2000 ...
    2000 2000 2000 2000 3000 3000 2000 2000 2000 2000 1000 1000 2000 1000 1000 2000]; %LSP05     
    
    case 'LSP02b'
        Muscle_ID = {'Intact VM', 'Intact RF', 'Intact VL', 'Intact BF', 'Intact ST', ...
        'Intact TA', 'Intact SO', 'Intact LG','Resid. VM','Resid. RF','Resid. VL',...
        'Resid. BF', 'Resid. ST', 'Resid. Ham', 'Resid. MG', 'Resid. LG'};
        load('D:\FigRescources\UH3\HotelCali\LSP02b\LSP02b_response_meta.mat');
        elecs = [response.elec];
        response2b = response;
        muscles = {'VM', 'RF', 'VL', 'BF', 'ST', 'SM', 'SO', 'LG',};
end

% % %% Extract Hit Rate
% % % % elecs = [response.elec];
% % 
% % for i = 1:length(elecs)
% %     e =  find(([response.elec] == elecs(i)));
% %     lowthresh = min([response(e).muscle(:).threshIdx]);
% %     for m = 1:16
% %         if isnan(response(e).muscle(m).threshIdx)
% %             spinalmap_full(i,m) = 0;
% %             spinalmap_lowthresh(i,m) = 0;
% %         else
% %             spinalmap_full(i,m) = (length(response(e).StimAmps) - response(e).muscle(m).threshIdx+1)/length(response(e).StimAmps);
% %             spinalmap_lowthresh(i,m) = (length(response(e).StimAmps) - response(e).muscle(m).threshIdx+1)/(length(response(e).StimAmps)-lowthresh+1);
% %         end
% %     end
% % end
% % 
% % %% Figure
% % figure; 
% % cmap = magma;
% % colormap(cmap);
% % % % nexttile;
% % % % Z = normalize(spinalmap_full(:,[9:13,15:16]),'range');
% % Z = spinalmap_full(:,9:16);
% % % Z(:, [6 7]) = Z(:, [7 6]);
% % imagesc(Z);%(:,4:7));
% % caxis([0 1]);
% % colorbar
% % yticks([1:length(elecs)]);
% % yticklabels({elecs});
% % ylabel('electrodes');
% % xticks([1:8]);
% % xticklabels(muscles);
% % xtickangle(45)
% % set(gca,'xaxisLocation','top')
% % box off
% % grid off
% % set(gca,'TickLength',[0 0])
% % title([subjectName ' - Raw']);
% % 
% % % % %% Figure with lowthresh
% % % % nexttile;
% % % % % Z = normalize(spinalmap_lowthresh,'range');
% % % % Z = spinalmap_lowthresh(:,9:16);
% % % % imagesc(Z);%(:,4:7));
% % % % colorbar
% % % % caxis([0 1]);
% % % % yticks([1:length(elecs)]);
% % % % yticklabels({elecs});
% % % % xticks([1:7]);
% % % % xticklabels(muscles([1:5, 7:8]));
% % % % xtickangle(45)
% % % % set(gca,'xaxisLocation','top')
% % % % box off
% % % % grid off
% % % % set(gca,'TickLength',[0 0])
% % % % title([subjectName ' - Adjust to Lowest Thresh']);
% % % % ylabel('electrodes');

% % %% Grouped outcomes
% % if strcmp(subjectName,'LSP02b')
% %     g(1).elecs = [1:4]
% %     g(2).elecs = [5:8]
% %     g(3).elecs = [9:12]
% %     g(4).elecs = [33:35]
% % elseif strcmp(subjectName,'LSP05')
% %     g(1).elecs = [1:4];
% %     g(2).elecs = [7:9];
% %     g(3).elecs = [10:11 18 12 19 13 20];
% %     g(4).elecs = [14 15 23 16];
% %     g(5).elecs = [24:27];
% %     g(6).elecs = [29:31];
% % else
% %     return;
% % end
% % 
% % 
% % for m = 1:16   
% %     for i = 1:length(g)
% %         for f = length(g(i).elecs)
% %             e =  find(([response.elec] == g(i).elecs(f)));
% %             lowthresh = min([response(e).muscle(:).threshIdx]);
% %             if isnan(response(e).muscle(m).threshIdx)
% %                 spinalhit(f) = 0;
% %                 spinalcount(f) = 0;
% %             else
% %                 spinalhit(f) = (length(response(e).StimAmps) - response(e).muscle(m).threshIdx+1);
% %                 spinalcount(f) = length(response(e).StimAmps);%-lowthresh+1;
% %                 spinalHR(f) = spinalhit(f)/spinalcount(f);
% %             end
% %         end
% %         spinalmap_group(i,m) = sum(spinalhit)/sum(spinalcount);
% %         spinalmap_groupmn(i,m) = mean(spinalHR);
% %         leg_string{i} = num2str([g(i).elecs]);
% %     end
% % end
% % 
% % figure; 
% % cmap = magma;
% % colormap(cmap);
% % % % nexttile;
% % % % Z = normalize(spinalmap_full(:,[9:13,15:16]),'range');
% % Z = spinalmap_group(:,9:16);
% % % Z(:, [6 7]) = Z(:, [7 6]);
% % imagesc(Z);%(:,4:7));
% % caxis([0 1]);
% % colorbar
% % yticks([1:length(g)]);
% % yticklabels(leg_string);
% % ylabel('electrodes');
% % xticks([1:8]);
% % xticklabels(muscles);
% % xtickangle(45)
% % set(gca,'xaxisLocation','top')
% % box off
% % grid off
% % set(gca,'TickLength',[0 0])
% % title([subjectName ' - Grouped']);


% % %% Grouped muscles
% % if strcmp(subjectName,'LSP02b')
% %     ms(1).muscles = [1:3]+8
% %     ms(2).muscles = [4]+8
% %     ms(3).muscles = [8]+8
% %     muscles = {'Quadriceps', 'Biceps Femoris', 'Gastrocnemius'};
% % elseif strcmp(subjectName,'LSP05')
% %     ms(1).muscles = [1:3]+8
% %     ms(2).muscles = [4]+8
% %     ms(3).muscles = [8]+8
% %     ms(4).muscles = [6]+8
% %     muscles = {'Quadriceps', 'Biceps Femoris', 'Gastrocnemius', 'Tibialis Anterior'};
% % else
% %     return;
% % end
% % 
% % if strcmp(subjectName,'LSP02b')
% %     g(1).elecs = [1:3]
% %     g(2).elecs = [5:11]
% %     g(3).elecs = [12]
% %    
% % elseif strcmp(subjectName,'LSP05')
% %     g(1).elecs = [1:2];
% %     g(2).elecs = [3:4 7:10];
% %     g(3).elecs = [11 18 12 19 13 20 14:16 23:25];
% %     g(4).elecs = [26:27 29:31];
% % else
% %     return;
% % end
% % 
% % % % figure; tiledlayout(length(g),1)
% % % % cmap = magma;
% % % % colormap(cmap);
% % for i = 1:length(g)
% % % %     nexttile;
% %     figure;
% %     cmap = magma;
% %     colormap(cmap);
% %     for t = 1:length(g(i).elecs)
% %             e =  find(([response.elec] == g(i).elecs(t)));
% % % % for i = 1:length(elecs)
% % % %     e =  find([response.elec] == elecs(i));
% %             lowthresh = min([response(e).muscle(:).threshIdx]);  
% %             for f = 1:length(ms)
% %                 for n = 1:length(ms(f).muscles)
% %                     m = ms(f).muscles(n)
% % 
% %                     if isnan(response(e).muscle(m).threshIdx)
% %                         spinalhit(n) = 0;
% %                         spinalcount(n) = 0;
% %                         spinalHR(n) = 0
% %                     else
% %                         spinalhit(n) = (length(response(e).StimAmps) - response(e).muscle(m).threshIdx+1)
% %                         spinalcount(n) = length(response(e).StimAmps)-lowthresh+1;
% %                         spinalHR(n) = spinalhit(n)/spinalcount(n)
% %                     end
% %                 end
% %                 spinalmap_group(t,f) = sum(spinalhit)/sum(spinalcount)
% %                 spinalmap_groupmn(t,f) = mean(spinalHR)
% %                 clear spinalhit spinalcount spinalHR
% %         %         leg_string{i} = num2str([ms(i).muscles]);
% %             end
% %     end
% %     Z = spinalmap_groupmn;
% %     imagesc(Z);%(:,4:7));
% %     caxis([0 1]);
% % % %     colorbar
% %     yticks([1:length(g(i).elecs)]);
% %     yticklabels({g(i).elecs});
% % % %     ylabel('electrodes');
% %     xticks([1:length(ms)]);
% %     xticklabels(muscles);
% %     xtickangle(45)
% %     set(gca,'xaxisLocation','top')
% %     box off
% %     grid off
% %     set(gca,'TickLength',[0 0])
% %     clear Z spinalmap_groupmn spinalmap_group spinalhit spinalcount
% % end
% % % % sgtitle([subjectName ' - Grouped']);
% % % % % % yticks([1:length(g)]);
% % % % % % yticklabels(leg_string);
% % % % 
% % % % 
% % % % 
% % % % % % nexttile;
% % % % % % Z = normalize(spinalmap_full(:,[9:13,15:16]),'range');
% % % % Z = spinalmap_groupmn;
% % % % % Z(:, [6 7]) = Z(:, [7 6]);
% % % % imagesc(Z);%(:,4:7));
% % % % caxis([0 1]);
% % % % colorbar
% % % % yticks([1:length(elecs)]);
% % % % yticklabels({elecs});
% % % % ylabel('electrodes');
% % % % xticks([1:4]);
% % % % xticklabels(muscles);
% % % % xtickangle(45)
% % % % set(gca,'xaxisLocation','top')
% % % % box off
% % % % grid off
% % % % set(gca,'TickLength',[0 0])
% % % % title([subjectName ' - Grouped']);

%% Grouped muscles
if strcmp(subjectName,'LSP02b')
    ms(1).muscles = [1:3]+8
    ms(2).muscles = [4]+8
    ms(3).muscles = [8]+8
    muscles = {'Quadriceps', 'Biceps Femoris', 'Gastrocnemius'};
elseif strcmp(subjectName,'LSP05')
    ms(1).muscles = [1:3]+8
    ms(2).muscles = [4]+8
    ms(3).muscles = [8]+8
    ms(4).muscles = [6]+8
    muscles = {'Quadriceps', 'Biceps Femoris', 'Gastrocnemius', 'Tibialis Anterior'};
else
    return;
end

if strcmp(subjectName,'LSP02b')
    g(1).elecs = [1:3]
    g(2).elecs = [5:11]
    g(3).elecs = [12]
   
elseif strcmp(subjectName,'LSP05')
    g(1).elecs = [1:2];
    g(2).elecs = [3:4 7:10];
    g(3).elecs = [11 18 12 19 13 20 14:16 23:25];
    g(4).elecs = [26:27 29:31];
else
    return;
end

figure; tiledlayout(length(g),1)
cmap = magma;
colormap(cmap);
for i = 1:length(g)
    nexttile;
% %     figure;
% %     cmap = magma;
% %     colormap(cmap);
    for t = 1:length(g(i).elecs)
            e =  find(([response.elec] == g(i).elecs(t)));
% % for i = 1:length(elecs)
% %     e =  find([response.elec] == elecs(i));
            lowthresh = min([response(e).muscle(:).threshIdx]);  
            for f = 1:length(ms)
                for n = 1:length(ms(f).muscles)
                    m = ms(f).muscles(n)

                    if isnan(response(e).muscle(m).threshIdx)
                        spinalhit(n) = 0;
                        spinalcount(n) = 0;
                        spinalHR(n) = 0
                    else
                        spinalhit(n) = (response(e).muscle(m).threshIdx - lowthresh)%(length(response(e).StimAmps) - response(e).muscle(m).threshIdx+1)
                        spinalcount(n) = length(response(e).StimAmps)-lowthresh+1;
                        spinalHR(n) = 1-spinalhit(n)/spinalcount(n)
                    end
                end
                spinalmap_group(t,f) = 1-sum(spinalhit)/sum(spinalcount);
                spinalmap_groupmn(t,f) = mean(spinalHR)
                clear spinalhit spinalcount spinalHR
        %         leg_string{i} = num2str([ms(i).muscles]);
            end
    end
    Z = spinalmap_groupmn;
    imagesc(Z);%(:,4:7));
    caxis([0 1]);
% %     colorbar
    yticks([1:length(g(i).elecs)]);
    yticklabels({g(i).elecs});
% %     ylabel('electrodes');
    xticks([1:length(ms)]);
    if i == 1
        xticklabels(muscles);
        xtickangle(45)
        set(gca,'xaxisLocation','top')
    end
    box off
    grid off
    set(gca,'TickLength',[0 0])
    clear Z spinalmap_groupmn spinalmap_group spinalhit spinalcount
% %     clear spinalmap_groupmn spinalmap_group spinalhit spinalcount
end
sgtitle([subjectName ' - Grouped']);
% % yticks([1:length(g)]);
% % yticklabels(leg_string);



% % % % nexttile;
% % % % Z = normalize(spinalmap_full(:,[9:13,15:16]),'range');
% % Z = spinalmap_group;
% % % Z(:, [6 7]) = Z(:, [7 6]);
% % imagesc(Z);%(:,4:7));
% % caxis([0 1]);
% % colorbar
% % yticks([1:length(elecs)]);
% % yticklabels({elecs});
% % ylabel('electrodes');
% % xticks([1:4]);
% % xticklabels(muscles);
% % xtickangle(45)
% % set(gca,'xaxisLocation','top')
% % box off
% % grid off
% % set(gca,'TickLength',[0 0])
% % title([subjectName ' - Grouped']);
