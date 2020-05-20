%% PCA of Recruitement AUC for all electrodes:
clear all;

load('D:\FigRescources\UH3\LSP05\rehash\responses_500us_1Hz_allElecsRC.mat')
response = response2; clear response2;

imuscles = [9:16];
labelmuscles = { 'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG'};
 
%% PCA for AUC
for iE = 1:length(response)
    for iM = 1 : length(imuscles)
        xinfo = response(iE).muscle(imuscles(iM)).p2pResponse;
        ampinfo = [ 0 response(iE).StimAmps];
        aucData(iE,iM) = sum(xinfo.*diff(ampinfo));
    end
end

[coeff,score,latent,tsquared,explained,mu] = pca(aucData');
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

%% NNMF for AUC
for iE = 1:length(response)
    for iM = 1 : length(imuscles)
        xinfo = response(iE).muscle(imuscles(iM)).p2pResponse;
        ampinfo = [ 0 response(iE).StimAmps];
        aucData(iE,iM) = sum(xinfo.*diff(ampinfo));
    end
end
figure;maximize;
subplot(121)
% C = linspecer(8,'qualitative');
[W,H] = nnmf(aucData,3);
biplot(H','VarLabels',labelmuscles)
axis equal
title('NNMF AUC by Muscle');
    xlim([0 1]);
    ylim([0 1]);
    zlim([0 1]);
    
subplot(122)
% C = linspecer(8,'qualitative');
[W2,H2] = nnmf(aucData',3);
biplot(H2','VarLabels',arrayfun(@num2str,[response(:).elec],'UniformOutput',false))
axis equal
title('NNMF AUC by Electrode');
    xlim([0 1]);
    ylim([0 1]);
    zlim([0 1]);

    %% NNMF for by subgroup
stimAmps = [2250];
elecs = [16];
C = linspecer(6,'qualitative');
% % stimAmps = [1000,3000,6000];
clear resp h;
midx = [1 3 2 4 5 6 8 7];

hF = figure; maximize;
    
for ttmp = 1:6
% %     clearvars -except response ttmp hF resp labelmuscles stimAmps elecs C imuscles
    switch ttmp
    case 1
        idelecs = [7:8]
        titstring = 'Electrodes 7-8';
    case 2
        idelecs = [10:12] 
        titstring = 'Electrodes 10-12';
    case 3
        idelecs = [9 14:15]
        titstring = 'Electrodes 9, 14-15';
    case 4
        idelecs = [16]
        titstring = 'Electrode 16';
    case 5
        idelecs = [18 19]
        titstring = 'Electrodes 18-19';
    case 6
        idelecs = [20 23 24]
        titstring = 'Electrodes 20, 23-24';
    otherwise
        idelecs = [1 4 7 8];
    end
        
    for iS = 1:length(stimAmps) 
        for iE = 1:length(idelecs)
            eI = find([response(:).elec] == idelecs(iE));
            idx = find(response(eI).StimAmps == stimAmps(iS));
            if isempty(idx)
                continue;
            end
            for iM = 1 : length(imuscles)
                resp(ttmp).Matrix(iE,iM) =  response(eI).muscle(imuscles(iM)).p2pResponse(idx);
            end
        end
    end
            resp(ttmp).elecs = idelecs;
            subplot(2,3,ttmp)
        if length(idelecs) <= 2
            k = 1   
        elseif  length(idelecs) > 2
            k = 3
        end
        
        [resp(ttmp).W,resp(ttmp).H] = nnmf(resp(ttmp).Matrix,k);
        if isempty(resp(ttmp).H)
            continue;
        elseif size(resp(ttmp).H,1) <2
            h = scatter([1:8],resp(ttmp).H(midx)',[],C(ttmp,:),'filled');
%             h = biplot([resp(ttmp).H'  zeros(size(resp(ttmp).H))'],...
%                 'VarLabels',labelmuscles,'Color',C(ttmp,:));
        ylabel('Component 1')
        xticklabels(labelmuscles(midx));
        ylim([0 1]);
        else
            h = biplot(resp(ttmp).H','VarLabels',labelmuscles,'Color',C(ttmp,:));
            xlim([0 1]);
            ylim([0 1]);
            zlim([0 1]);
        end

        title(['NNMF for ' titstring]);
       
 
end
        
% %         xticklabels

%% PCA for all stimAmps
stimAmps = unique([response(:).StimAmps]);

for iE = 1:length(response)
    for iM = 1 : length(imuscles)
        for iS = 1:length(response(iE).StimAmps)
            idx = find(stimAmps==response(iE).StimAmps(iS));
        respMatrix(iM,idx,iE) =  response(iE).muscle(imuscles(iM)).p2pResponse(iS);
        end
    end
end

% % [coeff,score,latent,tsquared,explained,mu] = pca(respMatrix);
% % scatter3(score(:,1),score(:,2),score(:,3))
% % axis equal
% % xlabel('1st Principal Component')
% % ylabel('2nd Principal Component')
% % zlabel('3rd Principal Component')


[rows, columns, stimlevels] = size(respMatrix);
listOfResponses = double(reshape(respMatrix, rows * columns, stimlevels));
[coeff,score,latent,tsquared,explained,mu] = pca(listOfResponses);

transformedResponseList = listOfResponses * coeff;
% Column 1 is the values of principal component #1, column 2 is the PC2, and column 3 is PC3.
% Extract each column and reshape back into a rectangular image the same size as the original image.
pca1Image = reshape(transformedResponseList(:,1), rows, columns);
pca2Image = reshape(transformedResponseList(:,2), rows, columns);
pca3Image = reshape(transformedResponseList(:,3), rows, columns);
scatter3(score(:,1),score(:,2),score(:,3))
axis equal
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
zlabel('3rd Principal Component')

   
%% PCA for 1mA, 2mA, 3mA
stimAmps = [750,1500,2250];
% % stimAmps = [1000,3000,6000];
clear resp;
figure; maximize;
for iS = 1:length(stimAmps) 
    for iE = 1:length(response)-1
        idx = find(response(iE).StimAmps == stimAmps(iS));
    for iM = 1 : length(imuscles)
        resp(iS).Matrix(iE,iM) =  response(iE).muscle(imuscles(iM)).p2pResponse(idx);
        end
    end
    
    [resp(iS).coeff,resp(iS).score,resp(iS).latent,resp(iS).tsquared,resp(iS).explained,resp(iS).mu] = pca(resp(iS).Matrix');
    h = scatter3(resp(iS).score(:,1),resp(iS).score(:,2),resp(iS).score(:,3),'filled');
    axis equal
    hold on;
% %     h = scatter3(resp(iS).score(:,4),resp(iS).score(:,5),resp(iS).score(:,6),'filled');
end
    legend('750mA','1500mA','2250mA');
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
%% PCA by electrode across stimulations
clear resp;
figure; maximize;
for iE = 1:length(response)
    leglabs{iE} = ['E' num2str(response(iE).elec)];
    for iS = 1:length(response(iE).StimAmps) 
       for iM = 1 : length(imuscles)
            resp(iE).Matrix(iS,iM) =  response(iE).muscle(imuscles(iM)).p2pResponse(iS);
       end
    end
    
    [resp(iE).coeff,resp(iE).score,resp(iE).latent,resp(iE).tsquared,resp(iE).explained,resp(iE).mu] = pca(resp(iE).Matrix');
    h = scatter3(resp(iE).score(:,1),resp(iE).score(:,2),resp(iE).score(:,3),'filled');
    axis equal
    hold on;
% %     h = scatter3(resp(iS).score(:,4),resp(iS).score(:,5),resp(iS).score(:,6),'filled');
end
    legend(leglabs);
        xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
    
%% PCA by muscles
clear resp;
figure; maximize;

for iM = 1 : length(imuscles)
    for iE = 1:length(response)
        for iS = 1:length(response(iE).StimAmps) 
            resp(iM).Matrix(iE,iS) =  response(iE).muscle(imuscles(iM)).p2pResponse(iS);
       end
    end
    
    [resp(iM).coeff,resp(iM).score,resp(iM).latent,resp(iM).tsquared,resp(iM).explained,resp(iM).mu] = pca(resp(iM).Matrix');
    h = scatter3(resp(iM).score(:,1),resp(iM).score(:,2),resp(iM).score(:,3),'filled');
    axis equal
    hold on;
% %     h = scatter3(resp(iS).score(:,4),resp(iS).score(:,5),resp(iS).score(:,6),'filled');
end
legend(labelmuscles);
    xlabel('1st Principal Component')
    ylabel('2nd Principal Component')
    zlabel('3rd Principal Component')
    
    
    %% Factor Analysis for 1mA, 2mA, 3mA
stimAmps = [750,1500,2250,3000,4000,6000];
C = linspecer(length(stimAmps),'qualitative');
% % stimAmps = [1000,3000,6000];
clear resp h;
    figure; maximize;
for iS = 1:length(stimAmps) 
    for iE = 1:length(response)
        idx = find(response(iE).StimAmps == stimAmps(iS));
        if isempty(idx)
            continue;
        end
        for iM = 1 : length(imuscles)
            resp(iS).Matrix(iE,iM) =  response(iE).muscle(imuscles(iM)).p2pResponse(idx);
        end
    end
    subplot(2,length(stimAmps)/2,iS)
    [resp(iS).lamda,resp(iS).psi,resp(iS).T,resp(iS).stats,resp(iS).F] = factoran(resp(iS).Matrix,3);
    h = biplot(resp(iS).lamda,'VarLabels',labelmuscles,'Color',C(iS,:));
    axis equal
    title(['FA for ' num2str(stimAmps(iS)/1000) 'mA']);
        zlim([0 1]);
            ylim([0 1]);
end

%% NNMF for 1mA, 2mA, 3mA
stimAmps = [750,1500,2250,3000,4000,6000];
% % stimAmps = [1000,3000,6000];
C = linspecer(length(stimAmps),'qualitative');
clear resp h;
    figure; maximize;
for iS = 1:length(stimAmps) 
    for iE = 1:length(response)
        idx = find(response(iE).StimAmps == stimAmps(iS));
        if isempty(idx)
            continue;
        end
        for iM = 1 : length(imuscles)
            resp(iS).Matrix(iE,iM) =  response(iE).muscle(imuscles(iM)).p2pResponse(idx);
        end
    end
    subplot(2,length(stimAmps)/2,iS)
    [resp(iS).W,resp(iS).H] = nnmf(resp(iS).Matrix,3);
    h = biplot(resp(iS).H','VarLabels',labelmuscles,'Color',C(iS,:));%,'scores',resp(iS).W
    axis equal
    title(['NNMF for ' num2str(stimAmps(iS)/1000) 'mA']);
    xlim([0 1]);
    ylim([0 1]);
    zlim([0 1]);
end

