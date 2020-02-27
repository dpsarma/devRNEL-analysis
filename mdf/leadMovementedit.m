baseDir = 'D:\FigRescources\UH3\XRays';
subjectName = 'LSP02b';
subjDir = fullfile(baseDir,subjectName);
week0Img = dir(fullfile(subjDir,'*week 0*'));
week1Img = dir(fullfile(subjDir,'*week 1*'));
week2Img = dir(fullfile(subjDir,'*week 2*'));
week3Img = dir(fullfile(subjDir,'*week 3*'));
week4Img = dir(fullfile(subjDir,'*week 4*'));
%%
for icompare = 1:4
    switch icompare
        case 1
            img1 = imread(fullfile(week0Img.folder,week0Img.name));
            img2 = imread(fullfile(week1Img.folder,week1Img.name));
        case 2
            img1 = imread(fullfile(week2Img.folder,week2Img.name));
            img2 = imread(fullfile(week1Img.folder,week1Img.name));
        case 3
            img1 = imread(fullfile(week3Img.folder,week3Img.name));
            img2 = imread(fullfile(week1Img.folder,week2Img.name));
        case 4
            img1 = imread(fullfile(week4Img.folder,week4Img.name));
            img2 = imread(fullfile(week1Img.folder,week3Img.name));
    end            
    [movingPts,fixedPts] = cpselect(img1,img2,'WAIT',true);
    %%
    tform = fitgeotrans(movingPts,fixedPts,'affine');
    warpedImg1 = imwarp(img1,tform,'OutputView',imref2d(size(img2)));
    %%
    figure, hold on    
    C = imfuse(img2,warpedImg1,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
    imshow(C)     
    keyboard
    %%
    [x,y] = ginput(2);
    Contactlenpx = sqrt(diff(x)^2 + diff(y)^2); % 3 mm
    for ilead = 1:3
        [x,y] = ginput(2);
        VertMvmnt(icompare,ilead) = diff(y) * (3/Contactlenpx); % in mm
    end
    close all
end
%%
clear all
close all
clc
baseDir = 'R:\analysis\human\sensory_stim\X-rays';
subjects = {'USP02';'USP01b';'USP03'};
cmap = [0 0.447 0.741;0.85 0.325 0.098;0.929 0.694 0.125];
figure, hold on
for isubj = 1:3
    subjectName = subjects{isubj};
    subjDir = fullfile(baseDir,subjectName);
    load(fullfile(subjDir,'mvmnt.mat'));
%     temp(:,:,isubj) = VertMvmnt;
    isubj
    VertMvmnt
%     plot(VertMvmnt,'LineWidth',2,'Color',[0.7 0.7 0.7])
    if isubj == 2
        mean(VertMvmnt(:,[1 3]),2)
    else
        mean(VertMvmnt,2)
    end
    median(VertMvmnt,2)
%         plot((1:size(VertMvmnt,1))+(0.05*(isubj-2)),median(VertMvmnt,2),'LineWidth',2,'Color',cmap(isubj,:))
        plot((1:size(VertMvmnt,1)),median(VertMvmnt,2),'LineWidth',2,'Color',cmap(isubj,:))

    for ii = 1:size(VertMvmnt,1)        
%         scatter(repmat(ii+(0.05*(isubj-2)),1,size(VertMvmnt,2)),VertMvmnt(ii,:),120,cmap(isubj,:),'filled','MarkerEdgeColor','k')
        scatter(repmat(ii,1,size(VertMvmnt,2)),VertMvmnt(ii,:),120,cmap(isubj,:),'filled','MarkerEdgeColor','k')
    end
%     errorbar(1:size(VertMvmnt,1),mean(VertMvmnt,2),std(VertMvmnt,[],2),'LineWidth',2)
    clear VertMvmnt
end

% breakyaxis([45 70])

ax = gca;
ax.FontSize = 18;
ax.XLim = [0.5 4.5];
ax.XTick = 1:4;
ax.XLabel.String = 'Week of Testing';
ax.YLabel.String = 'Vertical Migration (mm)';
box off