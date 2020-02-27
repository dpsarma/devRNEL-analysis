function res = dataPrepSingle(uuid,params)
  %
  % in input we get a structure with subject, date, electrode and list of uuids of the trials objects
  % create one mdf object with snippets for each channel/muscles
  %
  % Input
  % - uuid (string): uuid of the trial that we need to extract the epochs
  %                  from
  % - params (struct): over write default values of parameters.
  %                    List and default vales is listed below:
  %
  %   * frequencyLow  = 10
  %   * frequencyHigh = 7500
  %   * preTime       = 50e-3
  %   * postTime      = 500e-3
  %   * rmsWindowTime = 50e-3
  %   * blankTime     = 5e-3
  %   * fileBase      = ''. It uses the file base from trial
  %

  res = 0;
  
  om = mdfManage.getInstance();

  % default values
  frequencyLow = 10;
  frequencyHigh = 7500;
  preTime = 50e-3;
  postTime = 500e-3;
  rmsWindowTime = 50e-3;
  blankTime = 5e-3;
  fileBase = '';

  if nargin > 1
     if isfield(params,frequencyLow)
       frequencyLow = params.frequencyLow;
     end
     if isfield(params,frequencyHigh)
       frequencyHigh = params.frequencyHigh;
     end
     if isfield(params,preTime)
       preTime = params.preTime;
     end
     if isfield(params,postTime)
       postTime = params.postTime;
     end
     if isfield(params,rmsWindowTime)
       rmsWindowTime = params.rmsWindowTime;
     end
     if isfield(params,blankTime)
       blankTime = params.blankTime;
     end
     if isfield(params,fileBase)
       fileBase = params.fileBase;
     end
  end

  % load trial object
  oTrial = mdf.load(uuid);

  % prepares the file name if needed
  if strlength(fileBase)==0
     files = oTrial.getFiles(false);
     temp = split(strrep(files.base,'\','/'),'/');
     temp{3} = 'epochs';
     temp{4} = [temp{4} '_Epc'];
     filebase = join(temp,'/');
     filebase = filebase{1};
  end

  % load stim object
  oStim = oTrial.stim;

  % Extract stim settings
  
  % mA
  stimAmp = oTrial.stimParams.cathAmp/1000;
  % Hz
  stimFreq = oTrial.stimParams.freq;
  % uS
  stimPulse = 2*oTrial.stimParams.width/1000000;


%  epochTime = linspace(-preTime, (stimPulse) + postTime, preSamples + postSamples + stimlength );


  % loop on all the emg channels
  for i = 1:oTrial.getLen('EMG')
     % retrieve individual emg channel
     oEMG = oTrial.EMG(i);

     % Select data
     emgRaw =  oEMG.data.wf;
     emgTime = oEMG.data.time;

     % extract sampling frequency
     fs = oEMG.fs;
     
     % indexes of the stim event
     stimOns = floor(oStim.stimT*fs);
     % Set Epoch Window for Stim-triggered Averaging
     preSamples = floor(preTime*fs);
     postSamples = floor(postTime*fs);
     stimLength = floor(2*stimPulse*fs)+2;
     epochTime = linspace(-preTime, (stimPulse) + postTime, preSamples + postSamples + stimLength );
  
     % Bandpass Data
     emgFiltered = bandpass(emgRaw,[frequencyLow,frequencyHigh],fs);

     %Rectify & Smooth Data
     emgRect = abs(emgFiltered);

     % Epoch Data into Segments
     emgFilteredEpochs = cell2mat( ...
       arrayfun( ...
         @(x) emgFiltered((stimOns(x)-preSamples):(stimOns(x)+(stimLength+postSamples-1))), ...
         1:length(stimOns), ...
         'UniformOutput',false)');
     emgRectEpochs = cell2mat( ...
       arrayfun( ...
         @(x) emgRect((stimOns(x)-preSamples):(stimOns(x)+(stimLength+postSamples-1))), ...
         1:length(stimOns), ... 
         'UniformOutput',false)');

     % Find Mean 
     % NOTES: we do not use it
     meanFilteredTrace = mean(emgFilteredEpochs);
     meanRectTrace = mean(emgRectEpochs);

     %%% To Save

     % % % Verify by plotting
     % % figure;plot(epochtime,emgepochs); hold on; plot(epochtime,meantrace,'k');
     % % figure;plot(epochtime,emgepochs_rect); hold on; plot(epochtime,meantrace_rect,'k');
     % % figure;
     % % for i = 1:size(emgepochs,1)
     % %     subplot(15,1,i);
     % %     plot(epochtime,emgepochs(i,:)); hold on
     % %     plot(epochtime,meantrace);
     % % end

     % Calculate Values for Recruitment Curve and Onsets
     rmsWindow = floor(rmsWindowTime * fs); %in ms
     blankSamples = floor(blankTime * fs); %in ms post stim
     windowSamples = preSamples + stimLength + blankSamples;

     rmsFilteredEpochs = rms(emgFilteredEpochs(:,windowSamples:windowSamples+rmsWindow)')';
     p2pFilteredEpochs = peak2peak(emgFilteredEpochs(:,windowSamples:windowSamples+rmsWindow)')';

     rmsMean = mean(rmsFilteredEpochs); 
     rmsError = std(rmsFilteredEpochs)/sqrt(length(rmsFilteredEpochs));
     p2pMean = mean(p2pFilteredEpochs); 
     p2pError = std(p2pFilteredEpochs)/sqrt(length(p2pFilteredEpochs));

     % create mdf object and populate it
     oEpochs = mdfObj('Epochs');
     % files
     oEpochs.setFiles([filebase oEMG.EMGlabel]);
     % data
     oEpochs.data.filteredEpochs = emgFilteredEpochs;
     oEpochs.data.rmsFilteredEpochs = rmsFilteredEpochs;
     oEpochs.data.p2pFilteredEpochs = p2pFilteredEpochs;
     oEpochs.data.rectEpochs = emgRectEpochs;
     oEpochs.data.time = epochTime;
     oEpochs.data.meanFilteredTrace = meanFilteredTrace;
     oEpochs.data.meanRectTrace = meanRectTrace;
     % metadata
     oEpochs.metadata.rmsMean = rmsMean;
     oEpochs.metadata.rmsError = rmsError;
     oEpochs.metadata.p2pMean = p2pMean;
     oEpochs.metadata.p2pError = p2pError;
     oEpochs.metadata.stimAmp = stimAmp;
     oEpochs.metadata.stimFreq = stimFreq;
     oEpochs.metadata.stimPulse = stimPulse;
     oEpochs.metadata.subject = oTrial.subject;
     oEpochs.metadata.time = oTrial.time;
     oEpochs.metadata.nevChannel = oEMG.NEVchan;
     oEpochs.metadata.muscle = oEMG.EMGlabel;

     oEpochs.save();

     % now establish parent child relationship with trial
     mdf.apcr(oTrial,oEpochs,'epochs');
     oEpochs.save();

     om.clear(oEpochs);

  end %for

  % save trial
  oTrial.save();
  om.clear(oTrial);     
  
  res = 1;

end %function

