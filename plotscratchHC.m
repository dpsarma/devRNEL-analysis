% % i = 11;
% % figure;tiledlayout(4,1)
% %     for m = [11 12 13]
% %         nexttile;
% %         if meta(i).muscle(m).exist == 1
% %             plot(meta(i).muscle(m).first2timVec,meta(i).muscle(m).first2trace);
% %             title([meta(i).muscle(m).Muscle_ID '-'  num2str(meta(i).Amp)])
% %         end
% %     end


clear all
load('C:\figs\HotelCali\LNP02\LNP02_PAD_Elec4_20Hz.mat');

t = 11;
fs = 7.5e3;
begin = round(1.15*fs); %sec*fs
stop = round(1.25*fs);
figure;tiledlayout(3,1)
for m = [11 12 13] %VM BF ST
    ax(t) = nexttile;
    x = trial(t).timeVec(begin:stop);
    plot(x, trial(t).emgF(m,begin:stop));
    
    stim = trial(t).stims(1)/fs;
    hold on;
    vline([stim stim+0.05])
    ylim = [-500 500];
    

end
linkaxes([ax(:)],'xy');

clear all
load('C:\figs\HotelCali\LSP02b\LSP02b_PAD_Elec1_2Hz.mat');

t = 16;
fs = 30e3;
begin = round(1.13*fs); %sec*fs
stop = round(1.28*fs);

begin1 = trial(t).stims(1) - round(0.05*fs); %sec*fs
stop1 = trial(t).stims(1) + floor(0.05*fs);
begin2 = trial(t).stims(2);
stop2 = trial(t).stims(2) + floor(0.005*fs);
begin3 = trial(t).stims(1) + floor(0.05*fs) + floor(0.005*fs);
stop3 = trial(t).stims(1) + floor(0.05*fs) + floor(0.005*fs) + floor(0.045*fs);
figure;tiledlayout(4,1)
for m = [9 10 12 13] %VM RF BF ST
    ax(t) = nexttile;
    y = [trial(t).emgF(m,begin1:stop1) trial(t).emgF(m,begin2:stop2-1) trial(t).emgF(m,begin3:stop3-1)];
    x = trial(t).timeVec(begin:stop);
    plot(x, y);
    stim = trial(t).stims(1)/fs;
    hold on;
    vline([stim stim+0.05])
end
linkaxes([ax(:)],'xy');
    
clear all
load('C:\figs\HotelCali\LSP05\LSP05_PAD_Elec9_20Hz.mat');

t = 108;
fs = 30e3;
begin = trial(t).stims(1) - round(0.05*fs);; %sec*fs
stop = trial(t).stims(1) + floor(0.1*fs);
figure;tiledlayout(4,1)
for m = [9 12 16 14] %VM BF ST
    ax(t) = nexttile;
    x = trial(t).timeVec(begin:stop);
    plot(x, trial(t).emgF(m,begin:stop));
    ylim([-3e3 3e3])
end
linkaxes([ax(:)],'xy');




% % load('C:\figs\HotelCali\LSP05\LSP05_metaB_e9_1Hz.mat')
% % load('C:\figs\HotelCali\LSP05\LSP05_metaB_e9_2Hz.mat')
% % load('C:\figs\HotelCali\LSP05\LSP05_metaB_e9_5Hz.mat')
% % load('C:\figs\HotelCali\LSP05\LSP05_metaB_e9_10Hz.mat')

% t=108;
% figure; plot(trial(t).timeVec, trial(t).emgF(10,:));
% hold on;
% vline(trial(t).stims(1)/fs);
% 
% 
% t = 60; m = 16;
% figure; plot(meta(t).muscle(10).first2timVec,meta(t).muscle(m).first2trace)