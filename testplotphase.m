%% Plot EMG w/INSOLEs    
    disp('plotting')
% %     nexttile;
% %     figH = figure; maximize;
% %     tiledlayout(8,2)
% %     i = 0;
    c = [.875 .875 .875];
    
    emg_time = linspace(0,size(data.time,2)/fs,length(data.time)); %% EMG Time based on insole time length?


    stepcycle(f) = length(starts_R) - 1
    triallength(f) = round(emg_time(end))
    
figure; tiledlayout(4,1); maximize;
    for m = [3 4 6 8]%[9 1 11 3 12 4 14 6 16 8]
        z = emg(chan_remap,:)*1000*1000;
        mx_cop = max(max(z((m),:)));
        nexttile;
        plot(emg_time, z(m,:))
        hold on;
% %         plot(t,envelope(e,150,'rms'),'LineWidth',2);
% % %         plot(t,envelope(e,1000,'peak'),'LineWidth',2);
% %         hold on;
        e = emg(chan_remap(m),:); t = emg_time;
        e = e - mean(e);
%         v = find(islocalmin(envelope(e,500,'peak'))); 
        plot(t,abs(e));
        hold on;
        vline(t(idx), 'r-');
% %         mx_cop = max(abs(emg(chan_remap(m),:)))*1000;
% %         if m >=9
% %             yyaxis right; 
% %             plot(t_left, cop_L); ylabel('CoP');
% %             yyaxis left; hold on;
% %             for g =  2:length(starts_L)
% %                 x = [t_left(starts_L(g-1)) t_left(stops_L(g)) t_left(stops_L(g)) t_left(starts_L(g-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
% %                  h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
% %                 if g == length(starts_L)
% %                     continue;
% %                 else
% %                      x2 = [t_left(stops_L(g)) t_left(starts_L(g)) t_left(starts_L(g)) t_left(stops_L(g))]; 
% %                      h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
% %                 end
% %             end
% %         else
% %             yyaxis right; 
% %             plot(t_right, cop_R); ylabel('CoP');
% %             yyaxis left; hold on;
% % 
% %             for g =  2:length(starts_R)
% %                 x = [t_right(starts_R(g-1)) t_right(stops_R(g)) t_right(stops_R(g)) t_right(starts_R(g-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
% %                 h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
% %                 if g == length(starts_R)
% %                     continue;
% %                 else
% %                     x2 = [t_right(stops_R(g)) t_right(starts_R(g)) t_right(starts_R(g)) t_right(stops_R(g))]; 
% %                     h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
% %                 end
% %             end
% %         end
        box off;
    title(mLabels(chan_remap(m)));
    end
    tit = sgtitle({cell2mat(trialtype), cell2mat(filename)},'interpreter', 'none');