clear; close all; clc;

%% Some constants
DELTA=[1,4]; THETA=[4,8]; ALPHA=[8,14]; BETA=[14,30]; 
GAMMA=[30,50]; TTL=[1,50]; NBANDS = 6;

CHANNELS = 1:26;

CERE = [7,14,18]; FRONT = [1,2,3,25,26]; 
MOTOR = [5,21,22]; PARIETAL = [11,12,25,26]; O = [9,13,15];

WINPTS = 120; LN_WID = 4; NMIN = 28;

%% Load the file
load('./powers_rm_pxx.mat'); 
load('./f.mat'); % Frequency Vector

%% Select subjects
subjects = [1 3 5 6 ... %  
            7 9 10 11 13 ...
            15 16 17 18 19 ... % 17
            21 22 23 24 25];
% subjects = [7 10 23];

%% Subject results
alpCereSubjects = []; betCereSubjects = []; delCereSubjects = [];
gamCereSubjects = []; theCereSubjects = []; ttlCereSubjects = [];

alpFronSubjects = []; betFronSubjects = []; delFronSubjects = [];
gamFronSubjects = []; theFronSubjects = []; ttlFronSubjects = [];

alpMotSubjects = []; betMotSubjects = []; delMotSubjects = [];
gamMotSubjects = []; theMotSubjects = []; ttlMotSubjects = [];

alpPSubjects = []; betPSubjects = []; delPSubjects = [];
gamPSubjects = []; thePSubjects = []; ttlPSubjects = [];

alpOSubjects = []; betOSubjects = []; delOSubjects = [];
gamOSubjects = []; theOSubjects = []; ttlOSubjects = [];
subjectDataPts = zeros(length(subjects), 1);

%% Power averaging
for iSubject = 1:length(subjects)
subjectID = subjects(iSubject)
if subjectID < 10
   subjectIDstr = strcat('S0', num2str(subjectID));
else
   subjectIDstr = strcat('S', num2str(subjectID));
end

subjectPowers = powerM(subjectIDstr); 
pxxs = subjectPowers('pxx');

% Structure to save powers for all channels
for iBand = 1:NBANDS
    powAllChnls{iBand} = [];
end

nWins = ceil(size(pxxs{1},1)/WINPTS);

% Each channel
for iChannel = 1:length(CHANNELS)
    disp(['Channel: ', num2str(iChannel)]);
    pxxChnl = pxxs{iChannel};

    powChnlDel = []; powChnlThe = [];  
    powChnlAlp = []; powChnlBet = [];
    powChnlGam = []; powChnlTtl = [];
    for iWin = 1:nWins
        % Average pxx based on windows
        if iWin == nWins
            meanPWelch = nanmean(pxxChnl((iWin-1)*WINPTS+1:end,:));
        else
            meanPWelch = nanmean(pxxChnl((iWin-1)*WINPTS+1:iWin*WINPTS,:));
        end

        if sum(isnan(meanPWelch)) == 129
            pDel = NaN; pThe = NaN; pAlp = NaN;
            pBet = NaN; pGam = NaN; pTtl = NaN;
        else
            pDel = bandpower(meanPWelch,f,DELTA,"psd");
            pThe = bandpower(meanPWelch,f,THETA,"psd");
            pAlp = bandpower(meanPWelch,f,ALPHA,"psd");
            pBet = bandpower(meanPWelch,f,BETA,"psd");
            pGam = bandpower(meanPWelch,f,GAMMA,"psd");
            pTtl = bandpower(meanPWelch,f,TTL,"psd");
        end
        powChnlDel = [powChnlDel pDel];
        powChnlThe = [powChnlThe pThe];
        powChnlAlp = [powChnlAlp pAlp];
        powChnlBet = [powChnlBet pBet];
        powChnlGam = [powChnlGam pGam];
        powChnlTtl = [powChnlTtl pTtl];
    end
    powAllChnls{1} = [powAllChnls{1}; powChnlDel];
    powAllChnls{2} = [powAllChnls{2}; powChnlThe];
    powAllChnls{3} = [powAllChnls{3}; powChnlAlp];
    powAllChnls{4} = [powAllChnls{4}; powChnlBet];
    powAllChnls{5} = [powAllChnls{5}; powChnlGam];
    powAllChnls{6} = [powAllChnls{6}; powChnlTtl];
end

% Get the average across cerebellar channels
delCereSubject = nanmean(powAllChnls{1}(CERE,:));
theCereSubject = nanmean(powAllChnls{2}(CERE,:));
alpCereSubject = nanmean(powAllChnls{3}(CERE,:));
betCereSubject = nanmean(powAllChnls{4}(CERE,:));
gamCereSubject = nanmean(powAllChnls{5}(CERE,:));
ttlCereSubject = nanmean(powAllChnls{6}(CERE,:));

% Adjust the time axis to be consistent for plotting
delCereSubject = reform(delCereSubject);
theCereSubject = reform(theCereSubject);
alpCereSubject = reform(alpCereSubject);
betCereSubject = reform(betCereSubject);
gamCereSubject = reform(gamCereSubject);
ttlCereSubject = reform(ttlCereSubject);

delCereSubjects = [delCereSubjects; delCereSubject];
theCereSubjects = [theCereSubjects; theCereSubject];
alpCereSubjects = [alpCereSubjects; alpCereSubject];
betCereSubjects = [betCereSubjects; betCereSubject];
gamCereSubjects = [gamCereSubjects; gamCereSubject];
ttlCereSubjects = [ttlCereSubjects; ttlCereSubject];

end

%% Visualization
% Power across subjects
hold on;
for iSubject = 1:size(theCereSubjects)
   plot(ttCereSubjects(iSubject,:)); 
end
grid on; xticks([0,4,24,28]);

% Subject average power plot
figure('Position',[50,50,1840,280]);
x = 1:NMIN; y = nanmean(ttlCereSubjects); 
sigma = nanstd(dataPlot); 
curve1 = y + sigma; curve2 = y - sigma;
x2 = [x, fliplr(x)]; inBetween = [curve1, fliplr(curve2)];
fill(x2, inBetween, 'k','FaceAlpha',0.2, ...
    'LineStyle','none'); hold on;
plot(x, y, 'k', 'LineWidth', LN_WID);
xticks([1:NMIN]); xlim([1,NMIN]); 
ylabel('mV^2/Hz');  
set(gca,'FontSize',22); grid on;

% % Cross band plot
% gray1 = [.8 .8 .8]; gray2 = [.6 .6 .6];
% figure('Position',[200,50,1600,800]); x = times'; 
% y = nanmean(ttlCereSubjects); 
% plot(x, y, 'k', 'LineWidth', 4); hold on;
% y = nanmean(ttlFronSubjects); 
% plot(x, y, 'LineWidth', 4, 'Color', gray1); hold on;
% y = nanmean(ttlMotSubjects); 
% plot(x, y, 'LineWidth', 4, 'Color', gray2); hold on;
% xticks([0, 4, 24, 28]); xlim([0,28]);
% xlabel('Time (minute)'); ylabel(yl); 
% set(gca,'FontSize',22); grid on;
% % title('Total power change across brain regions');
% legend('Cerebellar', 'Frontal', 'Motor');

% Subject average power at different bands
plot_band_powers(ttlCereSubjects, ...
                 delCereSubjects, ...
                 theCereSubjects, ...
                 alpCereSubjects, ...
                 betCereSubjects, ...
                 gamCereSubjects, 'mV^2/Hz');

% Subject-by-subject plot
figure('Position',[200,50,1440,840]); 
x = 1:length(subjects); 
for i = 1:length(subjects)
    subplot(5,4,i);
%     y = deltaCereSubjects(i,:); 
%     plot(x, y, 'LineWidth', lw); hold on;
%     y = deltaCereSubjects(i,:); 
%     plot(x, y, 'LineWidth', lw); hold on;
%     y = alphaCereSubjects(i,:); 
%     plot(x, y, 'LineWidth', lw); hold on;
%     y = betaCereSubjects(i,:); 
%     plot(x, y, 'LineWidth', lw); hold on;
%     y = gammaCereSubjects(i,:); 
%     plot(x, y, 'LineWidth', lw); hold on;
    y = ttlCereSubjects(i,:); 
    plot(x, y, 'LineWidth', LN_WID); hold on;
    xticks([1:NMIN]); xlim([1,NMIN]);
    xlabel('Time (minute)'); ylabel('mV^2/Hz'); 
    grid on; 
end

