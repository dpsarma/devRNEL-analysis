t = emg_time; e = emg(chan_remap(),:);
e = e - mean(e);

 
figure;
nexttile;
plot(t,e);
nexttile;
yyaxis right;
plot(t_left, cop_L); ylabel('CoP');
yyaxis left; hold on;
for g =  2:length(starts_L)
    x = [t_left(starts_L(g-1)) t_left(stops_L(g)) t_left(stops_L(g)) t_left(starts_L(g-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
    h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
    if g == length(starts_L)
        continue;
    else
       x2 = [t_left(stops_L(g)) t_left(starts_L(g)) t_left(starts_L(g)) t_left(stops_L(g))]; 
       h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
    end
end
 if stimStatus == 1
% %             vline([cell2mat(stimEvts)-stimEvts{1}(1)],'r:');
                scatter([cell2mat(stimEvts)], 0, 0.5, 'r', '+');
                disp(['Plotting' mLabels(m) '- tile:0' num2str(i)]);     
 end
 ylabel([mLabels(chan_remap(m)) ' (mV)']);
 box off
 sgtitle(mLabels(chan_remap(12)));

figure;
plot(t,abs(e));
hold on;
plot(t,envelope(e,150,'rms'),'LineWidth',2);
plot(t,envelope(e,1000,'peak'),'LineWidth',2);

figure;
v = find(islocalmin(envelope(e,1000,'peak')));
plot(t,abs(e));
hold on;
vline(t(v), 'r-');
sgtitle(mLabels(chan_remap(12)));
