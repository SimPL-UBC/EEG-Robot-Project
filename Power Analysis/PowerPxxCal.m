clear; clc; close all;

%% Some constants
FS = 1000; PREPROCESSED_MODE = 'EEMD-CCA';
if strcmp(PREPROCESSED_MODE, 'EEMD-CCA')
    FS = 500; 
end

CHANNELS = 1:26;

WINTIME = 1; % 1-sec window
DELTA=[1,4]; THETA=[4,8]; ALPHA=[8,14]; BETA=[14,30]; 
GAMMA=[30,50]; TTL=[1,50]; %vgamma = [80,160]; 
SIGNAL_LIMIT = 100; NFBINS = 129;

%% Select subjects
subjects = [1 3 5 6 ... %  
            7 9 10 11 13 ...
            15 16 17 18 19 ... % 17
            21 22 23 24 25]; % 4,17
% subjects = [3];

%%
powerM = containers.Map(); % Structure to save the powers

for iSubject = 1:length(subjects)
subjectID = subjects(iSubject);
if subjectID < 10
   subjectIDstr = strcat('S0', num2str(subjectID));
else
   subjectIDstr = strcat('S', num2str(subjectID));
end
powerM(subjectIDstr) = containers.Map();

subjectDataDir = append(strcat('./subjects/', subjectIDstr, '/eeg/*.vhdr'));
nameList = dir(subjectDataDir); nTrials = length(nameList);
path = strcat('./subjects/', subjectIDstr, '/eeg/'); 
load(strcat(path, 'ALLEEG_cca.mat'));

% Extract powers
overlapPercent = 0.5; timeGap = WINTIME/(1/overlapPercent);
timeVec = []; prevTrialEndTime = 0; 
triggers = []; prevTrialEndTrig = 0;
for iChannel = 1:26
    pwelchChannels{iChannel} = [];
end

% Each trial
for trial = 1:nTrials
    fprintf(strcat(subjectIDstr, " ", 'Trial', num2str(trial), '\n'));
    data = double(ALLEEG(trial).data); 
    
    windowNPts = WINTIME*FS;
    
    % Extract triggers
    latencies = ALLEEG(trial).event; latencies = struct2table(latencies);
    triggerIndexesRaw = latencies.latency; trigger_types = latencies.type;
    selectedIndexes = [];
    for trigger = 1:length(trigger_types)
        thisTrigType = trigger_types{trigger};
        if strcmp(thisTrigType, 'M  1')
            selectedIndexes = [selectedIndexes trigger];
        end
    end
    triggerIndexesRaw = triggerIndexesRaw(selectedIndexes);
    triggerIndexes = zeros(length(triggerIndexesRaw)/2, 2);
    for iTrigger = 1:length(triggerIndexesRaw)
        if mod(iTrigger,2)==1
            i_row = ceil(iTrigger/2);
            triggerIndexes(i_row,1) = triggerIndexesRaw(iTrigger);
        end
        if mod(iTrigger,2)==0
            i_row = iTrigger/2;
            triggerIndexes(i_row,2) = triggerIndexesRaw(iTrigger);
        end
    end
    triggers = [triggers; triggerIndexes+prevTrialEndTrig];
    
    % Each segment
    prevSegEndTime = 0;
    for segment = 1:size(triggerIndexes,1)
        segTriggers = triggerIndexes(segment, :);
        segData = data(:, segTriggers(1):segTriggers(2));
        segStartTime = (segTriggers(1)-1)/1000; 
        segNPts = segTriggers(2)-segTriggers(1)+1;
        
        % Process if the segment > 1 sec
        if segNPts >= 1000 
            alphaChannels = []; betaChannels = []; gammaChannels = []; 
            thetaChannels = []; deltaChannels = []; ttlChannels = [];
            for channel = 1:length(CHANNELS) 
                pwelchBufs = [];
                [dataWins,z] = buffer(segData(channel, :), windowNPts, ...
                               round(windowNPts*overlapPercent), 'nodelay');
                        
               for buf = 1:size(dataWins,2)
                    dataBuf = dataWins(:,buf);  
                    % Process if there is no large artifact > 100 uV
                    if sum(abs(dataBuf) > SIGNAL_LIMIT) == 0
                        if length(dataBuf) > 1
                            disp(['Computing powers ...']);
                            [pxx, f] = pwelch(dataBuf,[],[],[],FS); 
                            pwelchBufs = [pwelchBufs; pxx'];              
                        else
                            disp(['#### Buf insufficient ####', num2str(buf)])
                            pwelchBufs = [pwelchBufs; NaN(1,NFBINS)];
                        end
                    else
                        disp(['Artifact found > 100 uV ', ...
                              subjectIDstr, ...
                              ' Trl ', num2str(trial), ...
                              ' Seg ', num2str(segment), ...
                              ' Cnl ', num2str(channel)])
                        pwelchBufs = [pwelchBufs; NaN(1,NFBINS)];
                    end
                    % Build a time vector
                    if channel == 1
                        timeVec = [timeVec ...
                            (buf)*timeGap+segStartTime+prevTrialEndTime];
                    end
               end
               pwelchBufs = fillmissing(pwelchBufs, 'nearest');
               pwelchChannels{channel} = [pwelchChannels{channel};
                         pwelchBufs];
            end
            timeAdd = 0;
        else
            fprintf('# of pts < 1000 \n');
            timeAdd = length(segData(1, :))/1000;
        end
        fprintf(strcat(num2str(timeVec(end)), ' Sec', '\n'));
    end
    
    prevTrialEndTrig = triggerIndexes(end, 2)+prevTrialEndTrig;
    prevTrialEndTime = timeVec(end) + timeAdd;

end
subjectPowers = powerM(subjectIDstr);
subjectPowers('pxx') = pwelchChannels;
subjectPowers('TimeVec') = timeVec;
subjectPowers('Triggers') = triggers;
end

%% Save the power estimates
powerM('Overlap') = overlapPercent;
save(strcat('./', 'powers_cca_pxx'), 'powerM');

