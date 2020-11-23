% % load('D:\FigRescources\UH3\LSP05\rehash\responses_500us_1Hz_allElecsRC.mat')
% % response = response2; clear response2;
load('D:\FigRescources\UH3\HotelCali\LSP02b\LSP02b_response_meta.mat');
labelmuscles = { 'rVM', 'rRF', 'rVL', 'rBF', 'rST', 'rTA', 'rSO', 'rLG', 'lVM', 'lRF', 'lVL', 'lBF', 'lST', 'lHAM', 'lMG', 'lLG'};
msubset = [1:16] %[9 11 10 12 13 14 15 16]; %[11 12 16 14];


elecs = [response.elec];%[response.elec]; %[1 4 7 9 11 18 13 20 16 23 26 29]; %
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

for iE = 1:length(elecs)
    for iM = 1 : length(labelmuscles)
        xinfo = response(iE).muscle(msubset(iM)).p2pResponse;
        ampinfo = [ 0 response(iE).StimAmps];
        aucData(iE,iM) = sum(xinfo.*diff(ampinfo));
    end
end
figure;
tiledlayout('flow');

Z = musc_map(1:16,:); %aucData'; %musc_map(9:16,:);
ZZ = normalize(Z','range');
A = aucData'; %musc_map(9:16,:);
AA = normalize(A','range');

% % nexttile
% % imagesc(ZZ(1:12,:))
% % colormap(magma)
% % yticks(1:26)
% % yticklabels({response(1:12).elec})
% % xticks(1:length(msubset));
% % xticklabels(labelmuscles)
% % set(gca,'xaxisLocation','top')
% % box off
% % grid off
% % set(gca,'TickLength',[0 0])
% % colorbar
% % title('Lead 1 - P2P')
% % 
nexttile
imagesc(AA(1:12,:))
colormap(magma)
yticks(1:26)
yticklabels({response(1:12).elec})
xticks(1:length(msubset));
xticklabels(labelmuscles)
set(gca,'xaxisLocation','top')
box off
grid off
set(gca,'TickLength',[0 0])
colorbar
title('Lead 1 - AUC')

% % nexttile
% % imagesc(ZZ(13:end,:))
% % colormap(magma)
% % yticks(1:26)
% % yticklabels({response(13:end).elec})
% % xticks(1:length(msubset));
% % xticklabels(labelmuscles)
% % set(gca,'xaxisLocation','top')
% % box off
% % grid off
% % set(gca,'TickLength',[0 0])
% % colorbar
% % title('Lead 2 - P2P')

nexttile
imagesc(AA(13:end,:))
colormap(magma)
yticks(1:26)
yticklabels({response(13:end).elec})
xticks(1:length(msubset));
xticklabels(labelmuscles)
set(gca,'xaxisLocation','top')
box off
grid off
set(gca,'TickLength',[0 0])
colorbar
title('Lead 2 - AUC')
