clear; clc; close all;
addpath('C:\Users\calvi\Documents\EEG_adaptation\eeglab_current\eeglab2023.0');
[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;

%% Some constants
CHANNELS = 32; REF3 = 13; FS = 1000;
MUSCLE_DENOISING_MODE = 'EEMD-CCA'; % CCA or EEMD-CCA
CORR_MODE = 'CLEANED'; % CLEANED or RAW

%% Subjects to analyze
subjects = [1 3 5 6 7 ...
            9 10 11 13 ...
            15 16 17 18 19 ... 
            21 22 23 24 25]; 

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

%% Get the path to files
subjectDataDir = append(strcat('./subjects/', subjectIDstr, '/eeg/*.vhdr'));
nameList = dir(subjectDataDir); nTrials = length(nameList);
path = strcat('./subjects/', subjectIDstr, '/eeg/'); 

%% Load the EEG file
EEG = pop_loadset('filename','ALL_raw.set','filepath', path);
[ALLEEG, EEG, CURRENTSET] = eeg_store(ALLEEG, EEG, 0);

%% Remove motion-artifact contaminated channels
EEG = pop_select( EEG, 'rmchannel',{'TP9'});
EEG = pop_select( EEG, 'rmchannel',{'TP10'});
EEG = pop_select( EEG, 'rmchannel',{'SPL'});
EEG = pop_select( EEG, 'rmchannel',{'SPR'});
EEG = pop_select( EEG, 'rmchannel',{'FP1'});
EEG = pop_select( EEG, 'rmchannel',{'FP2'});
[ALLEEG EEG CURRENTSET] = ...
pop_newset(ALLEEG, EEG, CURRENTSET,'overwrite','on','gui','off');
   
%% Cleaning
if strcmp(CORR_MODE, 'CLEANED')
    %% ICA & ICLabel
    EEG = pop_runica(EEG, 'icatype', 'sobi'); %sobi  runica
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
        new = [];
        for c = 1:26
            temp = eemdEEG(c,:);
            temp = upsample(temp, 10);
            new = [new; temp];
        end
        EEG.data = eemdEEG;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    elseif strcmp(MUSCLE_DENOISING_MODE,'CCA')
        % CCA
        eegData = double(EEG.data); 
        eegData = (resample(eegData.', 1, 10)).';
        gems = mean(eegData,2); % zero mean
        eegData = eegData-gems*ones(1,size(eegData,2)); 
        [y,w,r] = ccaqr(eegData,1);
        A = pinv(w'); nCCA = 18;
        A(:,end-nCCA+1:end) = 0; B = A*y;
        ccaEEG = [B B(:,end-1+1:end)]+gems*ones(1,size(eegData,2));    
        EEG.data = ccaEEG;
        FS = 100;
        [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
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
eegData = EEG.data;
ts = 1/FS:1/FS:length(eegData(1,:))/FS;

% Extract channels from the raw data
rawALLEEG = load(strcat([path, 'ALLEEG_raw.mat']));
rawEEGData = double(rawALLEEG.ALLEEG(end).data);
SPL = rawEEGData(5, :); SPR = rawEEGData(27, :);

% Combine and save
eegData =[SPL; SPR; eegData];

% Adaptation
indices = logical((ts>=240).*(ts<=1440));
eegCorr = [];
for i = 1:(CHANNELS-4)
   temp = eegData(i,:);
   eegCorr = [eegCorr; temp(indices)];
end

info = EEG.chanlocs;
info = struct2table(info);
lblInterest = info.labels(rearranged);
lblInterest = [ {'SPL'}; {'SPR'}; lblInterest];
xValues = lblInterest; yValues = lblInterest;
rearranged = rearranged + 2;
rearranged = [1 2 rearranged];

cor = corr(eegCorr'); 
corInterest = cor(rearranged,rearranged);

figure;
h = heatmap(xValues,yValues,abs(corInterest));
exportgraphics(gcf,strcat(['./corr/clean/adapt/', ...
                           subjectIDstr,'_', ...
                           MUSCLE_DENOISING_MODE, ...
                           '.png']));
close all;

% Pre-quiet standing
indices = ts<=4; eegCorr = [];
for i = 1:(CHANNELS-4)
   temp = eegData(i,:);
   eegCorr = [eegCorr; temp(indices)];
end
cor = corr(eegCorr'); 
corInterest = cor(rearranged,rearranged);
figure;
h = heatmap(xValues,yValues,abs(corInterest));
exportgraphics(gcf,strcat(['./corr/clean/pre/', ...
                           subjectIDstr,'_', ...
                           MUSCLE_DENOISING_MODE, ...
                           '.png']));
close all;

%% Power spectrum of cerebellar channels
CBZcleaned = eegData(4,:); 
CBLcleaned = eegData(3,:);
CBRcleaned = eegData(5,:);
CBLraw = rawEEGData(10,:); 
CBZraw = rawEEGData(17,:);
CBRraw = rawEEGData(21,:); 

figure('Position',[50,50,1840,280]);
hold on; plot(CBZraw); plot(CBZcleaned); 
xlim([120000,130000]); legend('Raw', 'Cleaned')
exportgraphics(gcf,strcat(['./corr/clean/time/', ...
                           subjectIDstr,'_', ...
                           MUSCLE_DENOISING_MODE, ...
                           '.png']));
close all;

figure; hold on;  
[p,f] = pspectrum(CBZraw,FS); plot(f,p);
[p,f] = pspectrum(CBZcleaned,FS); plot(f,p);
[p,f] = pspectrum(CBLraw,FS); plot(f,p);
[p,f] = pspectrum(CBLcleaned,FS); plot(f,p);
[p,f] = pspectrum(CBRraw,FS); plot(f,p);
[p,f] = pspectrum(CBRcleaned,FS); plot(f,p);
[p,f] = pspectrum(SPL,1000); plot(f,p); 
[p,f] = pspectrum(SPR,1000); plot(f,p);
xlim([1,80]); %ylim([0,80]);
xlabel('Hz'); grid on; legend('Raw CBZ', 'Cleaned CBZ', ...
                              'Raw CBL', 'Cleaned CBL', ...
                              'Raw CBR', 'Cleaned CBR', ...
                              'SPL', 'SPR');
hold off;
exportgraphics(gcf,strcat(['./corr/clean/psd/', ...
                           subjectIDstr,'_', ...
                           MUSCLE_DENOISING_MODE, ...
                           '.png']));
close all;

clear SPL SPR CBZcleaned CBZraw ALLEEG ALLCOM
clear eegData eegStruct eemdEEG h rawEEGData
clear temp ts p f eegCorr

end

disp(['Correlation Analysis Completed']);


