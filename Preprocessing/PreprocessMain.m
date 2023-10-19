clear; clc; close all;
addpath('C:\Users\calvi\Documents\EEG_adaptation\eeglab_current\eeglab2023.0');
eeglab;

%% Some constants
CHANNELS = 32; REF3 = 13; FS = 1000;
MUSCLE_DENOISING_MODE = 'CCA'; % CCA or EEMD-CCA

%% Select subjects
subjects = [3 1 5 6 ... %  
            7 9 10 11 13 ...
            15 16 17 18 19 ... 
            21 22 23 24 25]; % 4 is not usable 
% subjects = [16];

%% Each subject
for iSubject = 1:length(subjects)
ALLEEG=[];
subjectID = subjects(iSubject)
if subjectID < 10
   subjectIDstr = strcat('S0', num2str(subjectID));
else
   subjectIDstr = strcat('S', num2str(subjectID));
end
subjectDataDir = append(strcat('./subjects/', subjectIDstr, '/eeg/*.vhdr'));
nameList = dir(subjectDataDir); nTrials = length(nameList);
path = strcat('./subjects/', subjectIDstr, '/eeg/'); 

% Each trial
for trial = 1:nTrials 
    fileName = nameList(trial).name;
    EEG = pop_loadbv(path, fileName, [], []);
    [ALLEEG EEG CURRENTSET] = ...
        pop_newset(ALLEEG, EEG, 0,'setname',fileName,'gui','off'); 
    
    % Set channel locations
    EEG=pop_chanedit(EEG, 'rplurchanloc',CURRENTSET,'load',...
        {'ChannelLoc1208.ced','filetype','autodetect'});
    [ALLEEG EEG CURRENTSET] = ...
    pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
    
    % Re-referencing
    ref = EEG.data(REF3, :);
    for channel = 1:CHANNELS
        EEG.data(channel, :) = EEG.data(channel, :) - ref;
    end
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    clear ref
    
    % Remove motion-artifact contaminated channels
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
    trigger_indexes_raw = info.latency;
    type = info.type; trig = [];
    for i = 1:length(type)
        if strcmp(type{i}, 'M  1')
           trig = [trig trigger_indexes_raw(i)]; 
        end
    end
    trigger_indexes_raw = trig';
    trigger_indexes = zeros(length(trigger_indexes_raw)/2, 2);
    for i_trigger = 1:length(trigger_indexes_raw)
        if mod(i_trigger,2)==1
            i_row = ceil(i_trigger/2);
            trigger_indexes(i_row,1) = trigger_indexes_raw(i_trigger);
        end
        if mod(i_trigger,2)==0
            i_row = i_trigger/2;
            trigger_indexes(i_row,2) = trigger_indexes_raw(i_trigger);
        end
    end
    
    % Each segment
    for iSegment = 1:size(trigger_indexes, 1)
        
        EEGThis = pop_select(EEG, 'point', trigger_indexes(iSegment, :));
        
        % Band pass filtering
        EEGThis = pop_eegfiltnew(EEGThis,'locutoff',1,...
                                  'hicutoff',80);

        % Notch filtering
        eegData = double(EEGThis.data);
        eegData = notch_filter(eegData, 1000); 
        EEGThis.data = eegData;
       
        if iSegment == 1
            OUTEEG = EEGThis;
        end
        if iSegment > 1
            OUTEEG = pop_mergeset(OUTEEG, EEGThis);
        end
        length(OUTEEG.times)
    end
    CURRENTSET
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, OUTEEG, CURRENTSET,'overwrite','on','gui','off'); 
end

%% Merge all files
EEG = pop_mergeset(ALLEEG,[1:nTrials]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname','Merged Set','gui','off');
EEG = eeg_checkset(EEG);

%% ASR
% EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5,...
%                             'ChannelCriterion',0.8,...
%                             'LineNoiseCriterion',4, ...
%                             'Highpass','off', ...
%                             'BurstCriterion',20, ...
%                             'WindowCriterion',0.25, ...
%                             'BurstRejection','on', ...
%                             'Distance','Euclidian', ...
%                             'WindowCriterionTolerances',[-Inf 7] );
% 
% [ALLEEG EEG CURRENTSET] = ...
%     pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');

%% ICA & ICLabel
% Use the entire merged set. 
EEG = pop_runica(ALLEEG(end), 'icatype', 'sobi'); % sobi or runica
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% IC Label
EEG = eeg_checkset(EEG); EEG = pop_iclabel(EEG, 'default');
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

% ICA Flag
cut = 0.05;
EEG = pop_icflag(EEG,[NaN NaN;cut 1;cut 1; ...
                      cut 1;cut 1;cut 1;NaN NaN]);
[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset(EEG);

% Remove flagged components
EEG = pop_subcomp(EEG,'',0,0); % [] or ' means removing components flaged for rejection
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'gui','off'); 
[ALLEEG EEG CURRENTSET] = ...
pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');

%% Additional muscle denoising
if strcmp(MUSCLE_DENOISING_MODE, 'EEMD-CCA')
    % EEMD-CCA
    numIMFs = 10; 
    [clean_eeg] = EEMD_CCA(EEG.data, numIMFs, FS);
    EEG.data = clean_eeg;
elseif strcmp(MUSCLE_DENOISING_MODE, 'CCA')
    % CCA
    eegData = double(EEG.data); 
    gems = mean(eegData,2); % zero mean
    eegData = eegData-gems*ones(1,size(eegData,2)); 
    [y,w,r] = ccaqr(eegData,1);
    A = pinv(w'); nCCA = 14; % <=26
    A(:,end-nCCA+1:end)=0; B = A*y;
    ccaEEG = [B B(:,end-1+1:end)]+gems*ones(1,size(eegData,2));    
    EEG.data = ccaEEG;
end

%% Additional line noise removal
% EEG = pop_cleanline(EEG,'bandwidth',2, ...
%                      'chanlist',[1:26],...
%                      'compu tepower',1,...
%                      'linefreqs',60,'newversion',0,...
%                      'normSpectrum',0,'p',0.01,'pad',2,...
%                      'plotfigures',0,'scanforlines',0,...
%                      'sigtype','Channels','taperbandwidth',2,...
%                      'tau',100,'verb',1,'winsize',4,'winstep',1);

%% Saving final results
% Save EEG set
EEG = pop_saveset(ALLEEG(end),'filename', ...
    'ALL_fnl.set','filepath', strcat(path));

% Save EEG structure
save(strcat('./subjects/', subjectIDstr,...
    '/eeg/ALLEEG_cca.mat'), 'ALLEEG');

eeglab redraw

end
