%% Plotting Vertical Raster
clearvars -except response
close all

% 
% ttmp = input('Which Set? (1-6)');
% %% Params to be set from the user
%response = loadthestructure; 
% Put here the electrode number that you want for comparative analysis.
hF = figure; maximize;
for ttmp = 1:6
clearvars -except response ttmp hF
    switch ttmp
    case 1
        idelecs = [7:8]    %1:4    
    case 2
        idelecs = [10:12] %7:8
    case 3
        idelecs = [9 14:15] %10-15
    case 4
        idelecs = [16]  
    case 5
        idelecs = [18 19]
    case 6
        idelecs = [20 23:24]
    otherwise
        idelecs = [1 4 7 8];
end

%% Params to be set from the user
%response = loadthestructure; 
% Put here the electrode number that you want for comparative analysis. 
% idelecs = [1:4]  %[1 8 12 15 23 26 28];
elecs = [response(:).elec];
% nice comparing 7 or 8 with anything between 20 and 24 is great
% Select the muscles you want to look at (here all of them). Make sure
% labels correspond, in order
imuscles = [9:16];
labelmuscles = { 'VM', 'RF', 'VL', 'BF', 'ST', 'TA', 'MG', 'LG'};

% Motor pools position taken from a figure (numbers are pixel)
% Changing this will change the disposition of muscles on the y scale
MN_musclesYcoords = [
    140;... %VM -- invented
    165; ... % RF
    154; ...% VL
    334; ...% BF
    273; ...% ST
    206; ...%TA 
    283; ...%MG -- invented
    293; ...%LG
    ];
% Muscle position on the figure
musclesYcoords = [2 : 2: 16];
% Colors for plotting. Here there are 4 I like, in the following order
% dark blue
% dark cyan
% yellow
% orange
% put whatever you need
% colors = [[18 40 76];[2 144 159]; [255 213 0];[236 107 33]]./255;
% % for count = 1:length(elecs)
% %     idelecs = elecs(count);
C = linspecer(length(elecs),'qualitative');
for cc = 1:length(idelecs)
 colors(cc,:) = C(find(elecs == idelecs(cc)),:);
end


%% Retrieving params from the data structures
% Listing all electrodes in the data structure loaded
elecs_id_pos = [response(:).elec]';
% Convert electrode number in positional index in the structure
for ie = 1 : length(idelecs)
    [ielecs(ie)] = find(elecs_id_pos(:,1) == idelecs(ie));
    
end
% Making sure sorting of muscles is corrects
[sortedmuscles,sortedmuscleorder]  = sort([  max(MN_musclesYcoords) - MN_musclesYcoords]);

%% Plot 
   

    dimfig = 1.1*max(musclesYcoords);

    [sortedticks,sortedticksidx]  = sort([  dimfig - musclesYcoords]);

    % Just dots
    hS(ttmp) = subplot(6,1,ttmp); hold on;
    % This creates a bif of shift on the y axis between electrodes
    yshift =  linspace(-0.2,  0.2, length(ielecs));
    % Starting to plot stuff
    iEE = 0;
    for iE = ielecs
        iEE = iEE +1;
        for iM = 1 : length(imuscles)
            xinfo = response(iE).muscle(imuscles(iM)).p2pResponse; 
            ampinfo = [ 0 response(iE).StimAmps]; 


            yinfo = dimfig - ones(length(xinfo), 1)*musclesYcoords(sortedmuscleorder(iM)) + yshift(iEE);

            plot(sum(xinfo.*diff(ampinfo)), yinfo(1), 'o', ...
                'Markersize', 12,...
                'MarkerFaceColor', colors(iEE, :), ...
                'MarkerEdgeColor', colors(iEE, :))
    %         scatter((xinfo), yinfo, 'o', ...
    %             'MarkerFaceColor', colors(iEE, :), ...
    %             'MarkerEdgeColor', colors(iEE, :), ...
    %             'MarkerFaceAlpha', 0.5);

        end
% %          w = waitforbuttonpress;
    end

    sortedxticklabel = [labelmuscles];
    sortedxticklabel = sortedxticklabel(sortedmuscleorder);
    set(gca, 'ytick', sortedticks,...
        'yticklabel',(sortedxticklabel))
    xlabel('Recruitment AUC (mV*mA)')
    ylabel('Rostrocaudal Spinal Organization')
    title(['Electrodes: ' num2str(idelecs)]);
    xlim([0,4.5e6]);

    % NB: legend is difficult to make, check the colors order to know what
    % you're plotting
end
% % linkaxes([hS(:)],'x');
saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\RostrocaudalRC_AUC.png'])
saveas(hF,['D:\FigRescources\UH3\LSP05\rehash\RostrocaudalRC_AUC.svg'])