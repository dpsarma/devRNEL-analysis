clear all;
close all;

% Opening JSON file
f = py.open('\\share.files.pitt.edu\RnelShare\data_raw\human\uh3_stim\LNP02\data\Closed Loop\Insole_data\LNP02_CL_Ssn069_Set001_Blk001_Trl006.json');
 
% returns JSON object as
% a dictionary
data = py.json.load(f);


py.scipy.io.savemat('C:/data/tmp_insole.mat', mdict=struct('data', data));
    load('C:\data\tmp_insole.mat');

% % insole = 
 
% % % Iterating through the json
% % % list
% % for i in data['emp_details']:
% %     print(i)
% % end
% % % Closing file
% % f.close()

for i = 1:length(data)
    t_left(i) = data{i}.time_L;
    left(:,i) = data{i}.pressure_L;
    cop_L(i) = data{i}.cop_L;
    t_right(i) = data{i}.time_R;
    right(:,i) = data{i}.pressure_R;
    cop_R(i) = data{i}.cop_R;
end
t_left = t_left - t_left(1);
t_right = t_right - t_right(1);
%Find phases Left
mask = logical(cop_L(:).');    %(:).' to force row vector
stops_L = strfind([false, mask], [0 1]);
starts_L = strfind([mask, false], [1 0]);
centers_L = mean([starts_L;stops_L]);

%Find phases Right
mask = logical(cop_R(:).');    %(:).' to force row vector
stops_R = strfind([false, mask], [0 1]);
starts_R = strfind([mask, false], [1 0]);
centers_R = mean([starts_R;stops_R]);


figure;
nexttile;
stackedplot(t_left, left([1 15],:));
title('Left Insole')
nexttile;
stackedplot(t_right, right([1 15],:));
title('Right Insole')
nexttile;
plot(t_left, cop_L);
hold on; vline(t_left(starts_L),'r:','swing'); vline(t_left(stops_L), 'k:', 'stance');
% title('CoP: Left');
hold on; yyaxis right
plot(t_right, cop_R);
% title('CoP: Right');
sgtitle('LNP02_CL_Ssn069_Set001_Blk001_Trl006','interpreter', 'none');


figure;
plot(t_left, left([1 15],:)); hold on;
plot(t_right, right([1 15],:));
yyaxis right
plot(t_left, cop_L);
plot(t_right, cop_R);
vline(t_left(starts_L),'r:','Sw'); vline(t_left(stops_L), 'k:', 'St');
% % vline(t_right(starts_R),'r-','Sw'); vline(t_left(stops_R), 'k-', 'St');



% % figure; nexttile; plot(t_left, cop_L); hold on; vline(t_left(starts_L),'r:'); vline(t_left(stops_L), 'k:'); xlabel('Center of Pressure'); ylabel('time (sec)'); title('Phases based on Left CoP')
% % nexttile; plot(t_right, cop_R,'HandleVisibility', 'off'); hold on; vline(t_right(starts_R),'r:'); vline(t_left(stops_R), 'k:'); xlabel('Center of Pressure'); ylabel('time (sec)'); title('Phases based on Right CoP'); 
% % % legend('Swing', 'Stance');
% % 
% % for f =  2:length(starts_L)
% %     x = [t_left(starts_L(f-1)) t_left(stops_L(f)) t_left(stops_L(f)) t_left(starts_L(f-1))]; y = [-1 -1 1 1];
% %     h1 = fill(x, y, 'cyan','FaceAlpha',0.3);
% %     hf.FaceColor = [.875 .875 .875];
% % end

figure; nexttile; 
% Phases for Left side
plot(t_left, cop_L); hold on; 
% % vline(t_left(starts_L),'r:'); 
% % vline(t_left(stops_L), 'k:'); 

xlabel('Center of Pressure'); ylabel('time (sec)'); title('Phases based on Left CoP')
mx_cop = max(abs(cop_L));
c = [.875 .875 .875];

for f =  2:length(starts_L)
    x = [t_left(starts_L(f-1)) t_left(stops_L(f)) t_left(stops_L(f)) t_left(starts_L(f-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
    h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
    if f == length(starts_L)
        continue;
    else
        x2 = [t_left(stops_L(f)) t_left(starts_L(f)) t_left(starts_L(f)) t_left(stops_L(f))]; 
        h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
    end
end

nexttile;
plot(t_right, cop_R); hold on; 
% % vline(t_right(starts_R),'r:'); 
% % vline(t_right(stops_R), 'k:'); 

xlabel('Center of Pressure'); ylabel('time (sec)'); title('Phases based on right CoP')
mx_cop = max(abs(cop_R));
c = [.875 .875 .875];

for f =  2:length(starts_R)
    x = [t_right(starts_R(f-1)) t_right(stops_R(f)) t_right(stops_R(f)) t_right(starts_R(f-1))]; y = [-mx_cop -mx_cop mx_cop mx_cop];
    h1 = fill(x, y, c,'FaceAlpha',0.3, 'EdgeColor','none');
    if f == length(starts_R)
        continue;
    else
        x2 = [t_right(stops_R(f)) t_right(starts_R(f)) t_right(starts_R(f)) t_right(stops_R(f))]; 
        h2 = fill(x2, y, 'magenta', 'FaceAlpha',0.05,  'EdgeColor','none');
    end
end
