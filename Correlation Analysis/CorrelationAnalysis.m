clear; clc; close all;
addpath('C:\Users\calvi\Documents\EEG_adaptation\eeglab_current\eeglab2023.0');
eeglab;

%% Some constants
CHANNELS = 32; REF3 = 13; FS = 1000;
MUSCLE_DENOISING_MODE = 'CCA'; % CCA or EEMD-CCA
CORR_MODE = 'CLEANED'; % CLEANED or RAW

%% Subjects to analyze
subjects = [3 1 5 6 ...
            7 9 10 11 13 ...
            15 16 17 18 19 ... 
            21 22 23 24 25]; 
subjects = [3,21];

%% Analyze each subject
for iSubject = 1:length(subjects)
ALLEEG=[];
subjectID = subjects(iSubject);
if subjectID < 10
   subjectIDstr = strcat('S0', num2str(subjectID));
else
   subjectIDstr = strcat('S', num2str(subjectID));
end
disp(['Analyzing ', subjectIDstr]);

% Get the path to files
subjectDataDir = append(strcat('./subjects/', subjectIDstr, '/eeg/*.vhdr'));
nameList = dir(subjectDataDir); nTrials = length(nameList);
path = strcat('./subjects/', subjectIDstr, '/eeg/'); 

for trial = 1:nTrials 
    % Load the EEG file
    fileName = nameList(trial).name
    EEG = pop_loadbv(path, fileName, [], []);
    [ALLEEG EEG CURRENTSET] = ...
        pop_newset(ALLEEG, EEG, 0,'setname',fileName,'gui','off'); 
    
    % Set channel locations
    EEG=pop_chanedit(EEG, 'rplurchanloc',CURRENTSET,'load',...
        {'ChannelLoc1208.ced','filetype','autodetect'});
    [ALLEEG EEG CURRENTSET] = ...
    pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
    
    % Re-reference
    ref = EEG.data(REF3, :);
    for channel = 1:CHANNELS
        EEG.data(channel, :) = EEG.data(channel, :) - ref;
    end
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    clear ref
    
    % Remove motion-artifact channels
    EEG = pop_select( EEG, 'rmchannel',{'TP9'});
    EEG = pop_select( EEG, 'rmchannel',{'TP10'});
    EEG = pop_select( EEG, 'rmchannel',{'SPL'});
    EEG = pop_select( EEG, 'rmchannel',{'SPR'});
    EEG = pop_select( EEG, 'rmchannel',{'FP1'});
    EEG = pop_select( EEG, 'rmchannel',{'FP2'});
    [ALLEEG EEG CURRENTSET] = ...
    pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
    
    % Extract triggers
    info = ALLEEG(CURRENTSET).event;
    info = struct2table(info);
    triggerIndicesRaw = info.latency;
    type = info.type; trig = [];
    for i = 1:length(type)
        if strcmp(type{i}, 'M  1')
           trig = [trig triggerIndicesRaw(i)]; 
        end
    end
    triggerIndicesRaw = trig';
    triggerIndicesNew = zeros(length(triggerIndicesRaw)/2, 2);
    for iTrigger = 1:length(triggerIndicesRaw)
        if mod(iTrigger,2)==1
            iRow = ceil(iTrigger/2);
            triggerIndicesNew(iRow,1) = triggerIndicesRaw(iTrigger);
        end
        if mod(iTrigger,2)==0
            iRow = iTrigger/2;
            triggerIndicesNew(iRow,2) = triggerIndicesRaw(iTrigger);
        end
    end
    
    % Slice the EEG data based on triggers 
    for iSegment = 1:size(triggerIndicesNew, 1)
        EEGSegment = pop_select(EEG, 'point', triggerIndicesNew(iSegment, :));
        % Band-pass 1-80 Hz
        EEGSegment = pop_eegfiltnew(EEGSegment,'locutoff',1,...
                          'hicutoff',80);
        eegData = double(EEGSegment.data);
        
        % Notch filter at 60 Hz
        eegData = notch_filter(eegData, 1000); 
        EEGSegment.data = eegData;
        
        if iSegment == 1
            OUTEEG = EEGSegment;
        end
        if iSegment > 1
            OUTEEG = pop_mergeset(OUTEEG, EEGSegment);
        end
        length(OUTEEG.times)
    end
    CURRENTSET
    [ALLEEG EEG CURRENTSET] = ...
        pop_newset(ALLEEG, OUTEEG, CURRENTSET,'overwrite','on','gui','off'); 
end

%% Merge files
EEG = pop_mergeset(ALLEEG,[1:nTrials]);
[ALLEEG EEG CURRENTSET] = pop_newset(...
    ALLEEG, EEG, 0,'setname','Merged Set','gui','off');
EEG = eeg_checkset(EEG); 

if strcmp(CORR_MODE, 'CLEANED')
    %% ICA & ICLabel
    EEG = pop_runica(ALLEEG(end), 'icatype', 'sobi'); %sobi  runica
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    EEG = eeg_checkset(EEG); EEG = pop_iclabel(EEG, 'default');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    cut = 0.5;
    EEG = pop_icflag(EEG,[NaN NaN;cut 1;cut 1; cut 1;cut 1;cut 1;NaN NaN]);
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset(EEG);

    EEG = pop_subcomp(EEG,'',0,0); 
    [ALLEEG EEG CURRENTSET] = ...
    pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');

    %% Muscular denoising
    if strcmp(MUSCLE_DENOISING_MODE,'EEMD-CCA')
        % EEMD-CCA
        numIMFs = 10; 
        [eemdEEG] = EEMD_CCA(EEG.data, numIMFs, FS);
        EEG.data = eemdEEG; % Downsampled to 500 Hz after this step
    elseif strcmp(MUSCLE_DENOISING_MODE,'CCA')
        % CCA
        eegData = double(EEG.data); gems = mean(eegData,2); % zero mean
        eegData = eegData-gems*ones(1,size(eegData,2)); 
        [y,w,r] = ccaqr(eegData,1);
        A = pinv(w'); nCCA=18; %<=27
        A(:,end-nCCA+1:end) = 0; B = A*y;
        ccaEEG = [B B(:,end-1+1:end)]+gems*ones(1,size(eegData,2));    
        EEG.data = ccaEEG;
    end
end

%% Correlation analysis
% Channels rearranged
rearranged = [7,14,18, ...
              9,13,15, ...
              11,12,16,17, ...
              8,19,20, ...
              5,21,22, ...
              4,24, ...
              1,2,3,25,26, ...
              6,23];
% rearranged = [4,25, ...
%               8,15,19, ...
%               10,14,16, ...
%               12,13,17,18, ...
%               9,20,21, ...
%               6,22,23, ...
%               5,26, ...
%               1,2,3,27,28];

eegData = double(EEG.data);
ts = 1/500:1/500:length(eegData(1,:))/500;

% Adaptation
indices = logical((ts>=240).*(ts<=1440));
eegCorr = [];
for i = 1:(CHANNELS-6)
   temp = eegData(i,:);
   eegCorr = [eegCorr; temp(indices)];
end

cor = corr(eegCorr'); 
corInterest = cor(rearranged,rearranged);
tbl = array2table(corInterest);
info = EEG.chanlocs;
info = struct2table(info);
lblInterest = info.labels(rearranged);
xValues = lblInterest; yValues = lblInterest;
h = heatmap(xValues,yValues,abs(corInterest));
exportgraphics(gcf,strcat(['./corr/clean/adapt/', ...
                            subjectIDstr,'_cca8.png']));
close all;

% Pre-quiet standing
indices = ts<=4; eegCorr = [];
for i = 1:(CHANNELS-6)
   temp = eegData(i,:);
   eegCorr = [eegCorr; temp(indices)];
end
cor = corr(eegCorr'); 
corInterest = cor(rearranged,rearranged);
h = heatmap(xValues,yValues,abs(corInterest));
exportgraphics(gcf,strcat(['./corr/clean/pre/', ...
                           subjectIDstr,'_cca8.png']));
close all;
end

disp(['Correlation Analysis Completed']);


