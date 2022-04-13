%% Quick Plot for Delsys EMG for GRC Slides

data = py.pickle.load( py.open('C:/data/LL_UH3/LNP02_CL_Ssn056_Set001_Blk001_Trl015.pkl', 'rb'));
% data = py.pickle.load( py.open('C:/data/LL_UH3/LNP02_CL_Ssn055_Set001_Blk001_Trl002.pkl', 'rb')); %LNP02_CL_Ssn055_Set001_Blk001_Trl002 %LNP02_CL_Ssn053_Set001_Blk001_Trl007
py.scipy.io.savemat('C:/data/overground.mat', mdict=struct('data', data));
load('C:\data\overground.mat');
data.emg = double(data.emg);
% Filtering Settings:
fs = 2000;
lowCut = 75; %lowest frequency to pass
highCut = 750; %highest frequency to pass
Norder = 2;
Wp = [lowCut, highCut]/(.5*fs);
[b,a]=butter(Norder, Wp);

emg = filtfilt(b,a,data.emg')';

tend = length(data.time)/fs;
mLabels = {"Right TFL", "Right RF", "Right TA", "Right SO", "Right LG", "Right VL",...
    "Left TFL", "Left RF", "Left VL", "Right BF", "Left BF", "Left ST", "Left TA",...
    "Right ST", "Left SO", "Left LG"};
chan_remap = [1 2 6 10 14 3 4 5 7 8 9 11 12 13 15 16];
% stackedplot(data.time, data.emg(chan_remap,:)', "DisplayLabels", mLabels(chan_remap));
tiledlayout(4,1)
i = 0;
for m = [9, 10, 12, 16]
    i = i+1;
    ax(i) = nexttile
    plot(linspace(0,tend,length(data.time)), emg(chan_remap(m),:)*1000)
    hold on;
    
    ylabel([mLabels(chan_remap(m)) ' (mV)'])
    box off
end

linkaxes([ax(:)], 'x'); xlim([0,tend]); %ylim([-0.5,.5]); 
 
xlabel ('time (sec)');

%stackedplot(outdoors,'Title','Weather Data','DisplayLabels', mLabels)


% % trialInfo = load('C:/data/LL_UH3/LNP02_CL_Ssn056_Set001_Blk001_Trl015_x0085.nev');
% % nevStimCh = trialInfo.Stim_params(1).NSChannels{1,1}(1);

% % [stimEvts] = read_stimEvents('C:/data/LL_UH3/LNP02_CL_Ssn056_Set001_Blk001_Trl015_x0085');
% % 
% % [ns_status, hFile] = ns_OpenFile('C:/data/LL_UH3/LNP02_CL_Ssn056_Set001_Blk001_Trl015_x0085'); 
% %         [ns_RESULT, nsFileInfo] = ns_GetFileInfo(hFile);
% %         ns_CloseFile(hFile);