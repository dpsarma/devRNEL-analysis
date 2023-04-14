load('C:\data\Gait_BarPlotData3.mat')

% % XX = weave(allPeak',allPeak2');
% % 
% % 
% % function out = weave(varargin)
% % % Example: weave(jet,cool,autumn,spring);
% % N = numel(varargin);
% % out = nan(N * max(cellfun(@(s) (size(s, 1)), varargin)), size(varargin{1}, 2));
% % for ii = 1:N
% %   out(ii:N:size(varargin{ii}, 1) * N - (N - ii), :) = varargin{ii};
% % end
% % while all(isnan(out(end, :)))
% %   out(end, :) = [];
% % end
% % end


% % figure; nexttile; boxplot(allPeak'); ylim([0 3.5]);
% % nexttile; boxplot(allPeak2'); ylim([0 3.5]);
% % ttest2(allPeak(1,:),allPeak2(1,:));
%% Boxplot for Activation
% % allPeak  = cell2mat(rz);
% % allPeak2 = cell2mat(rz2);

XX = zeros(size(allPeak2));
    XX(1,:) = randsample(allPeak(1,:),length(allPeak2));
    XX(2,:) = allPeak2(1,:);
    XX(4,:) = allPeak2(2,:);
    XX(3,:) = randsample(allPeak(2,:),length(allPeak2));
    XX(6,:) = allPeak2(3,:);
    XX(5,:) = randsample(allPeak(3,:),length(allPeak2));
    XX(8,:) = allPeak2(4,:);
    XX(7,:) = randsample(allPeak(4,:),length(allPeak2));
    XX(10,:) = allPeak2(5,:);
    XX(9,:) = randsample(allPeak(5,:),length(allPeak2));
    XX(12,:) = allPeak2(6,:);
    XX(11,:) = randsample(allPeak(6,:),length(allPeak2));
    XX(14,:) = allPeak2(7,:);
    XX(13,:) = randsample(allPeak(7,:),length(allPeak2));
    XX(16,:) = allPeak2(8,:);
    XX(15,:) = randsample(allPeak(8,:),length(allPeak2));
   
    figure; boxplot(XX')
    pval = [];
    for i = 1:2:16
        pval = [pval ttest2(XX(i,:),XX(i+1,:),'Alpha', 0.001,'vartype','unequal')];
    end
    figure; bar(pval);
    %% Figure out Gait differences

   cycDur1 = []; cycDur2 = [];
        for f = 1:length(EP)
            for count = 1:size(EP(f).m(2).i,2)
            cycDur1 = [cycDur1 length(EP(f).m(2).i(count).etim)];  
            end
        end
  
        for f = 1:length(EP2)
            for count = 1:size(EP2(f).m(2).i,2)
            cycDur2 = [cycDur2 length(EP2(f).m(2).i(count).etim)];
            end
        end

        for m = mlist
           for f = 1:length(EP)
             for count = 1:size(EP(f).m(m).i,2)
                cycDur1 = [cycDur1 length(EP(f).m(2).i(count).etim)];  
             end
           end
  
            for f = 1:length(EP2)
                for count = 1:size(EP2(f).m(m).i,2)
                    cycDur2 = [cycDur2 length(EP2(f).m(m).i(count).etim)];
                end
            end
        end

    cycDur1 = randsample(cycDur1, length(cycDur2))/2000;
    cycDur2 = cycDur2/2000;
    figure;violinplot([cycDur2; cycDur1;]');
disp('step size ttest');
ttest2(cycDur1, cycDur2, 'Vartype','unequal', 'Alpha', 0.001)
% % celEP = {};
% %     %% Time Normalized Epochs
% %     for m = mlist
% %         for f = 1:length(EP)
% %             for count = 1:size(EP(f).m(m).i,2)
% %                 celEP{end+1} = EP(f).m(m).i(count).epoch;
% %             end
% %         end
% %     end
% %     celEP = celEP';
% %     maxlen=max(cellfun('prodofsize',celEP));
% %     for n=1:numel(celEP)
% %         current_elem=numel(celEP{n});
% %         if current_elem<maxlen
% %          celEP{n}((current_elem+1):maxlen)=NaN;
% %         end
% %     end
% % 
% %     extradim=1+ndims(celEP{1});
% %     bigarray=cat(extradim,celEP{:});
% %     tmp = mean(bigarray,extradim,'omitnan');
% %     e_tim = linspace(-20,130,length(tmp));


    %% All Epochs

    figure; tiledlayout(7,1); maximize;
    col_map = [0 20 89; 82 130 192; 245 77 0; 250 138 115]/255;
    mlist = [8 14 6 10 2 12 4];
    for m = mlist
        nexttile;
        if m == 10 || m == 2
            i = 1;
            pc = 'w:';
        elseif m == 12 || m == 4
            i = 3;
            pc = 'k:';
        elseif m == 14 || m == 6
            i = 4;
            pc = 'k:';
        elseif m == 8
            i = 2;
            pc = 'k:';
        else
            disp('muscle not found');
            continue;
        end

        % Stim
        celEP1 = {};
        for f = 1:length(EP)
            for count = 1:size(EP(f).m(m).i,2)
              celEP1{end+1} = EP(f).m(m).i(count).epoch;  
            plot(EP(f).m(m).i(count).etim, EP(f).m(m).i(count).epoch, 'Color', col_map(i,:));  
            hold on;
            end
        end
        xlim([0 100]); xticklabels([]); box off; ylim([0 inf]);
        title(mLabels(chan_remap(m)));
% %         ylabel({'Activation', '(z-scored EMG)'});

        celEP1 = celEP1';
        maxlen=max(cellfun('prodofsize',celEP1));
        for n=1:numel(celEP1)
            current_elem=numel(celEP1{n});
            if current_elem<maxlen
             celEP1{n}((current_elem+1):maxlen)=NaN;
            end
        end
        extradim=1+ndims(celEP1{1});
        bigarray=cat(extradim,celEP1{:});
        tmp = mean(bigarray,extradim,'omitnan');
        tmp_tim = linspace(0,100,length(tmp));
% %         plot(tmp_tim, tmp, pc, 'LineWidth',2);
        if m == mlist(end)
            xticklabels({0,'',20,'',40,'',60,'',80,'',100});
            xlabel('% of Gait Cycle')
        end
    end

    sgtitle({'Stimulation On','%Gait Phase', 'Right Heel Strike:Heel Strike'},'interpreter', 'none');

            % No Stim
figure; tiledlayout(7,1); maximize;
     for m = mlist
        nexttile;
        if m == 10 || m == 2
            i = 1;
            pc = 'w:';
        elseif m == 12 || m == 4
            i = 3;
            pc = 'k:';
        elseif m == 14 || m == 6
            i = 4;
            pc = 'k:';
        elseif m == 8
            i = 2;
            pc = 'k:';
        else
            disp('muscle not found');
            continue;
        end

        celEP2 = {};
        for f = 1:length(EP2)
            for count = 1:size(EP2(f).m(m).i,2)
              celEP2{end+1} = EP2(f).m(m).i(count).epoch;  
            plot(EP2(f).m(m).i(count).etim, EP2(f).m(m).i(count).epoch, 'Color', col_map(i,:));  
            hold on;
            end
        end
        xlim([0 100]); xticklabels([]); box off; ylim([0 inf]);
        title(mLabels(chan_remap(m)));
% %         ylabel({'Activation', '(z-scored EMG)'});

        celEP2 = celEP2';
        maxlen=max(cellfun('prodofsize',celEP2));
        for n=1:numel(celEP2)
            current_elem=numel(celEP2{n});
            if current_elem<maxlen
             celEP2{n}((current_elem+1):maxlen)=NaN;
            end
        end
        extradim2=1+ndims(celEP2{1});
        bigarray2=cat(extradim2,celEP2{:});
        tmp2 = mean(bigarray2,extradim2,'omitnan');
        tmp_tim2 = linspace(-20,120,length(tmp2));
% %         plot(tmp_tim2, tmp2, pc, 'LineWidth',2);

        if m == mlist(end)
            xticklabels({0,'',20,'',40,'',60,'',80,'',100});
            xlabel('% of Gait Cycle')
        end
     end
sgtitle({'No Stimulation', '%Gait Phase', 'Right Heel Strike:Heel Strike'},'interpreter', 'none');