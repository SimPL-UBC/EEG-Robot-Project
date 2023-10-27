%{
This script extracts raw data segments from the trials and 
connect them.
%}
%%
clear; clc; close all;
addpath('C:\Users\calvi\Documents\EEG_adaptation\eeglab_current\eeglab2023.0');
eeglab;

%% Some constants
CHANNELS = 32; REF = 13; FS = 1000; 
BANDPASS_ENABLED = true;

%% Select subjects
subjects = [1 3 5 6 ... %  
            7 9 10 11 13 ...
            15 16 17 18 19 ... 
            21 22 23 24 25]; % 4 is not usable 

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
    ref = EEG.data(REF, :);
    for channel = [1:4,6:26,28:CHANNELS]
        EEG.data(channel, :) = EEG.data(channel, :) - ref;
    end
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    clear ref
    
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
    
    if subjectID == 3
        if trial == 1
            trigger_indexes(1,:) = [];
        end
    end
    
    % Each segment
    for iSegment = 1:size(trigger_indexes, 1)
        
        EEGThis = pop_select(EEG, 'point', trigger_indexes(iSegment, :));
        
        if BANDPASS_ENABLED
        % Band pass filtering
        EEGThis = pop_eegfiltnew(EEGThis,'locutoff',1,...
                                  'hicutoff',50);
        end
        
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

%% Saving final results
% Save EEG set
EEG = pop_saveset(ALLEEG(end),'filename', ...
    'ALL_raw.set','filepath', strcat(path));

% Save EEG structure
save(strcat('./subjects/', subjectIDstr,...
    '/eeg/ALLEEG_raw.mat'), 'ALLEEG');

eeglab redraw
end



