%% R01 EMG Intramuscular Data vs. Sleeve EMG - Part 1
%Load Data from a subject and plot raw data of all muscles...


%Example Data files: load('R:\data_raw\human\emg_decoding\EMGC14\20170424_s15\Exports\EMG\RNEL_Raw_30khz\RNEL_EMGC14_EMG_DATA_Raw_30kHz_trial_0008.mat')
disp('Starting Basic Plotting for HAPTIX Instramuscular EMG for R01 Resub\n');

%Load data folder
mainpath = 'R:\data_raw\human\emg_decoding';
disp('Please Select the Subject Directory (EXPORT folder)');
subjpath = uigetdir(mainpath);
oldfolder = cd(subjpath);
% emgpath = uigetdir(mainpath, 'Select EMG Data Folder');

disp('Please Select the EMG Data Folder (RAW 30K)');
[filenames, emgpathname] = uigetfile('*mat','Select all Files','MultiSelect', 'on');
disp('Select the Trials Metadata File');
[taskFile, taskpath] = uigetfile('*csv','Select Task Meta File','MultiSelect', 'off');

trials = readtable(fullfile(taskpath,taskFile),'ReadRowNames',true);

%For Subject 14, remove trial 0016 from trial table
trials(9,:) = [];

if size(filenames,2) ~= height(trials)
   disp('There is a mismatch between files and trials!!!!');
   return;
end

for fNum = 1:height(trials)
    load(fullfile(emgpathname,filenames{fNum}));
    data = emg_30khz';
    tRef = linspace(0, length(data)/30000, length(data));
    
    disp(['Plotting subplot raster for File ' num2str(fNum) ' - ' trials.sTrialType{fNum}]);
    figure; hold on;
    maximize;
    suplabel('Muscles (uV)','y');
    headTit = suptitle([trials.sSubject{fNum} ' - ' trials.sTrialType{fNum} ' (' trials.sFile{fNum} ')']);
    set(headTit,'Interpreter', 'none');
    for iChan = 1:size(data, 1)
        subplot(8,2,iChan); 
        plot(tRef,data(iChan,:)); 
%         legend(chanelsLabel{iChan},'Location','WestOutside');
        aY = ylabel([chanelsLabel{iChan}]);
        set(aY,'Interpreter', 'none');
        xlabel('Time (sec)');
    end

    
    disp(['Plotting stacked plot for File ' num2str(fNum) ' - ' trials.sTrialType{fNum}]);
    figure; maximize; 
    xlabel('Time (sec)'); 
    suplabel('Muscles (uV)','y');
      
    stackedplot(tRef, data, 'PlotLabels', chanelsLabel,'PlotSpans','equal', 'ShowBaseline',false); xlabel('Time (sec)'); suplabel('Muscles (uV)','y');
    headTit = title({[trials.sSubject{fNum} ' - ' trials.sTrialType{fNum}], ['(' trials.sFile{fNum} ')']});
    set(headTit,'Interpreter', 'none'); 
    
    saveallopenfigs(['C:\Users\dsarma\Box Sync\R-01 Resubmit Figs\1stPass\EMGC14-' num2str(fNum) '-']); %'D:\Figures\R01-Resub\1stPass\EMGC14-'
end


