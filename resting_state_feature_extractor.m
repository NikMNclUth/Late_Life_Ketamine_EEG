function [freqs,predslope,spectra,frontaltable,centraltable,posteriortable,temporaltable] = resting_state_feature_extractor(EEGfinal)
% resting pipeline. performs automated artifact removal of resting state
% EEG data.
% Contains the following externally developed functions:
% get_FAA.m
% Michael Tesar <tesarm@jcu.cz>
% Ceske Budejovice, 2016
%
% fitPowerLaw3steps.m
% Colombo, M., 2019 the spectral exponent of the resting EEG indexes the presence of consciousness during unresponsiveness induced by propofol, xenon, and ketamine
%
% Inputs:
%   EEG_final = pre-processed EEG structure variable
% Output:
%   freqs = frequencies estimated
%   predslope = exponential fit
%   spectra = relative power spectra
%   tables = roi data tables
%
% Routine composed by Dr. Nicholas Murphy, Baylor College of Medicine, 2021

%% step 1, determine regions of interest
F_ROI={'FP1';'FP2';'F7';'F3';'FZ';'F4';'F8';'CZ'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatef=find(ismember(chan_IDs,F_ROI));

C_ROI={'FC3';'FC1';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4';'CP1';'CP2';'CPZ'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatec=find(ismember(chan_IDs,C_ROI));

P_ROI = {'P7','P3','PZ','P4','P8','PO9','PO10',...
    'O1','OZ','O2','P1','P5','P2','P6','PO7','PO3','POZ','PO4','PO8'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatep=find(ismember(chan_IDs,P_ROI));

T_ROI = {'T7','T8','TP9','TP10','FT9','FT10','FT7','FT8','TP7','TP8'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatet=find(ismember(chan_IDs,T_ROI));

%% step 2, get subject code
naam = EEGfinal.filename(1:end-4);

%% step 3, detrending & demeaning
dataf = zeros(size(EEGfinal.data(:,:,:)));
for ii = 1:size(EEGfinal.data(:,:),3)
    dataf(:,:,ii) = EEGfinal.data(:,:,ii)-nanmean(EEGfinal.data(:,:,ii),2);
    temp = detrend(dataf(:,:,ii)');
    dataf(:,:,ii) = temp';
    clear temp
end

%% step 4, compute spectra
for u = 1:size(dataf,1)
    if u == 1
        [s,freqs] = pwelch(dataf(u,:),hanning(EEGfinal.srate),EEGfinal.srate/2,(EEGfinal.srate),EEGfinal.srate);
        spectra = zeros(size(dataf,1),length(s));
        spectra(u,:) = s;
    else
        [spectra(u,:),~] = pwelch(dataf(u,:),hanning(EEGfinal.srate),EEGfinal.srate/2,(EEGfinal.srate),EEGfinal.srate);
    end
end
% convert to relative power
spectra = spectra./nansum(spectra,2);

%% step 5, typical spectra (calculate peak alpha separately)
delta = find(freqs<4);
theta = find(freqs>=4 & freqs<7);
alpha = find(freqs>=7 & freqs<14);
beta = find(freqs>=25 & freqs<34);
gamma = find(freqs>=35 & freqs<45);

% Frontal
F.delta = nanmean(nanmean(spectra(candidatef,delta),2),1);
F.theta = nanmean(nanmean(spectra(candidatef,theta),2),1);
F.alpha= nanmean(nanmean(spectra(candidatef,alpha),2),1);
F.beta = nanmean(nanmean(spectra(candidatef,beta),2),1);
F.gamma = nanmean(nanmean(spectra(candidatef,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatef,alpha),1));
F.alpha_peak_power = y;
F.alpha_peak_frq = freqs(alpha(iF));

% Central
C.delta = nanmean(nanmean(spectra(candidatec,delta),2),1);
C.theta = nanmean(nanmean(spectra(candidatec,theta),2),1);
C.alpha= nanmean(nanmean(spectra(candidatec,alpha),2),1);
C.beta = nanmean(nanmean(spectra(candidatec,beta),2),1);
C.gamma = nanmean(nanmean(spectra(candidatec,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatec,alpha),1));
C.alpha_peak_power = y;
C.alpha_peak_frq = freqs(alpha(iF));

% Posterior
P.delta = nanmean(nanmean(spectra(candidatep,delta),2),1);
P.theta = nanmean(nanmean(spectra(candidatep,theta),2),1);
P.alpha= nanmean(nanmean(spectra(candidatep,alpha),2),1);
P.beta = nanmean(nanmean(spectra(candidatep,beta),2),1);
P.gamma = nanmean(nanmean(spectra(candidatep,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatep,alpha),1));
P.alpha_peak_power = y;
P.alpha_peak_frq = freqs(alpha(iF));

% Temporal
T.delta = nanmean(nanmean(spectra(candidatet,delta),2),1);
T.theta = nanmean(nanmean(spectra(candidatet,theta),2),1);
T.alpha= nanmean(nanmean(spectra(candidatet,alpha),2),1);
T.beta = nanmean(nanmean(spectra(candidatet,beta),2),1);
T.gamma = nanmean(nanmean(spectra(candidatet,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatet,alpha),1));
T.alpha_peak_power = y;
T.alpha_peak_frq = freqs(alpha(iF));

%% step 6, estimate slope
% spectra = nanmean(spectra,1);
% doPlot= 1;
% thisCol= [0 0 1];
% spectra = spectra./sum(spectra);
frBand=[3 45]; % frequencies to use
frBins= dsearchn(freqs, frBand(1) ):dsearchn(freqs,frBand(2));
XX= freqs(frBins); % cutdown frequencies% [global_power_vars_log_orig] = measure_power(global_orig_spectra.log_spectra,global_orig_spectra.freqs);
robRegMeth= 'ols';
for ii = 1:size(spectra,1)
    YY= spectra(ii,frBins); % spectrum cutdown by bins
    [intSlo, ~, Pows, ~,  ~, ~] = fitPowerLaw3steps(XX,YY, robRegMeth,  0, []);
    spectra2(ii,:) = Pows.obs;
    psd_slope(ii,1) = intSlo(2);% exponent
    psd_slope(ii,2) = intSlo(1);% intercept
    if ii == 1
        predslope = zeros(size(psd_slope,1),length(Pows.pred),4);
    end
    predslope(ii,:,1) = Pows.pred; % this is the background fit
    predslope(ii,:,2) = Pows.obs;
    predslope(ii,:,3) = Pows.res;
    predslope(ii,:,4) = Pows.frex;
end
F.psd_slope = nanmean(psd_slope(candidatef,:),1);
C.psd_slope = nanmean(psd_slope(candidatec,:),1);
P.psd_slope = nanmean(psd_slope(candidatep,:),1);
T.psd_slope = nanmean(psd_slope(candidatet,:),1);

%% step 7, frontal alpha asymmetry
[FAA,~,~,FAA_cat] = get_FAA(EEGfinal);

%% Export
frontaltable = table({naam},F.delta,F.theta,F.alpha,F.beta,...
    F.gamma,F.alpha_peak_frq,F.alpha_peak_power,F.psd_slope(1),FAA,FAA_cat);
frontaltable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent','alpha_asymm','alpha_asymm_direction'};
centraltable = table({naam},C.delta,C.theta,C.alpha,C.beta,...
    C.gamma,C.alpha_peak_frq,C.alpha_peak_power,C.psd_slope(1));
centraltable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent'};
posteriortable = table({naam},P.delta,P.theta,P.alpha,P.beta,...
    P.gamma,P.alpha_peak_frq,P.alpha_peak_power,P.psd_slope(1));
posteriortable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent'};
temporaltable = table({naam},T.delta,T.theta,T.alpha,T.beta,...
    T.gamma,T.alpha_peak_frq,T.alpha_peak_power,T.psd_slope(1));
temporaltable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent'};
end





function  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
% %%
% Colombo, M. A. et al. The spectral exponent of the resting EEG indexes the presence of consciousness during unresponsiveness induced by propofol, xenon, and ketamine. NeuroImage 189, 631–644 (2019).
% % fit a line on log-freq vs log-PSD, excluding peaks and their base.
% % USAGE EXAMPLE:
% % first compute the PSD
% epLen= 2* sRate; epShift= 1*srate;numFFT=[];
%  [myPSD,frex]= pwelch( myEEGch  , epLen, epShift,numFFT, sRate);
%  frBand=[1 40];
%  frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
% XX= frex(frBins);
% YY= myPSD(frBins);
% robRegMeth= 'ols';
% doPlot= 1;
% thisCol= [0 0 1];
%  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
%%%%%%%%% INPUT
% XX ---> are the frequency bins of interest (hz)
% YY ---> are the PSD relative to the bins of interest (mV^2 /Hz)
% robRegMeth --->is the method for robust regression fit (default: ols)
%                type help robustfit for more details
% doPlot  --->    0 (default) means no plot ;
%                 1 is plot XX and YY on log-log scale;
%                 2 is plot log(XX) and log(YY) on linear scale;
% thisCol --->    if plotting, is the rgb triplet
%%%%%%%%% OUTPUT
% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit  i.e. the spectral exponent

% stat is a structure with informations on the 2nd fit, (see 2nd argument of robustfit)

% Pows is a structure with
%           pred: 1*length(frBins) doubles, predicted PSD by the power-law fit
%            obs: 1*length(frBins) doubles, observed PSD    = YY
%            res: 1*length(frBins) doubles, residuals (pred-obs) PSD
%           frex: 1*length(frBins) doubles, observed frex   = XX
%       meanPred: 1 double, mean of the predicted PSD
%        meanObs: 1 double, mean of the observed PSD
%        meanRes: 1 double, mean of the residuals of the PSD
%     meanResPos: 1 double, mean of the residuals of the PSD (exclude negative values)
%       cenFrObs: 1 double, central Frequency of the observed PSD
%       cenFrRes: 1 double, central Frequency of the residual PSD

% Deviants is a structure with
%       res2: 1*length(frBins) doubles,
%        res: 1*length(frBins) doubles,
%        rej:  length(frBins)*1 logicals, 1 if bins are rejected, 0 if kept   for fitting the 2nd power-law
%         fr: min and max freqs considered
%       thre: threshold used for bins adjacent to the peaks
%     clusts: structure with information on the contiguous freq. bins

% stat0 is a structure with informations on the 1st fit, (see 2nd argument of robustfit)

% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOW DOES IT WORK?
% transform in log-log and upsample by 4 times
% fit a 1st line, find all residual >0
% small peaks (whose residuals are smaller than 1*mad(residuals)) do not count as peaks
% find clusters of rejected frequency bins
% consider only clusters of bins where there is a peak
% i.e. reject large enough peaks and their base (adjacent residuals >0)
% fit a second line, on the remaining residuals (those closer to a power-law)

if exist('robRegMeth','var')
else
    robRegMeth= 'ols';
end
%% LOG TRANSFORM
%%%%%%%%%%%%%%%%%
% vectorize frex and psd to avoid dimension mismatch later
XX= XX(:)';
YY= YY(:)';

X= log10(XX);
Y= log10(YY);

%%% INTERPOLATE IN LOG-LOG: X, Y --> Xi, Yi
XXi= logspace( X(1), X(end),  length(X)*4); % 2^10*2;--> auto-chosen
Xi= log10(XXi);
Yi= interpn( X, Y, Xi   );% ,'spline'
YYi= 10.^(Yi);

%% STEP 1, FIT 1st LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL 1: minimize Y residuals. OLS +variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ intSlo0, stat0] = robustfit(Xi,Yi,robRegMeth);
YRes0= stat0.resid;
YPred0=  ( intSlo0(1)+intSlo0(2)*(Xi))';

%% FIND DEVIANT RESIDUALS
%%%%%%%%%%%%%%%%%%%%%%%%%
threRes=0;
boolYdev= YRes0 > threRes ;%
%% ONLY CONSIDER DEVIANTS WHERE ALSO POWER-PEAKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      %(I.E. if power is high, but it decays so that it does not peak (e.g.peak is lower than earliest freq bin) --> not deviant
% peaks are searched on log-power .... (why not on residuals? avoid small bumps on res)
[pks0,locs0] = findpeaks( Yi ); %, 'MINPEAKDISTANCE', 8, 'NPEAKS', 10, 'SORTSTR','descend') ;
% REJECT SMALL PEAKS
% peaks are rejected if their residuals (detrended log-power) are low ... (why not on log-power? there is no meaningful threshold on log-pow)
pks=pks0; locs=locs0;

threResPks= 1* mad(YRes0,1); %    median(YRes0) + 1 *mad(myres,1);
threResPks= max(threResPks, .1);

rejectPks= (YRes0(locs) < threResPks )';
pks(rejectPks)=[];
locs(rejectPks)=[];

% now exclude deviant-residuals that do not contain any peak
Clusts = bwconncomp(boolYdev);
includeClust= false( 1, Clusts.NumObjects);
for cl= 1: Clusts.NumObjects
    thisIdx=Clusts.PixelIdxList{cl};
    thisclust= Xi( thisIdx );
    thisloc= (Xi(locs));
    includeClust(cl)=   sum( ismember ( thisclust , thisloc));
    boolYdev( thisIdx )=includeClust(cl);
end

%% STEP 2 ,2nd line after excluding peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MODEL 1: minimize Y residuals. OLS +variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try % try an OLS regression first
    [ intSlo, stat] = robustfit(Xi(~boolYdev),Yi(~boolYdev),robRegMeth);
catch % if not enough points for regression,(too small freq.band or deviations with too much band-width) just do not exclude any peak ans do ols
    [ intSlo, stat] = robustfit(Xi,Yi,'ols');
end
%%%% GET RESIDUALS ALSO AT EXCLUDED POINTS
YPred0= (Xi.*intSlo(2) +intSlo(1));
YRes2= Yi - YPred0;    %     YRes2= stat.resid;


%% STORE
stat.frBins= [ XXi(1)   XXi(end)  ];
stat.myName= [ num2str(round( stat.frBins(1) )) '-'  num2str(round(  stat.frBins(2)) ) ' Hz: '   num2str(intSlo(2)) ]; %num2str(intSlo')
stat.robRegMeth=robRegMeth;

Clusts.includeClust=includeClust;
Clusts.pks=pks;
Clusts.locs=locs;
Clusts.threResPks=threResPks;
Clusts.pks0=pks0;
Clusts.locs0=locs0;
Clusts.rejectPks0=rejectPks;

Deviants.res2= YRes2;
Deviants.res= YRes0;
Deviants.rej= boolYdev;
Deviants.fr= stat.frBins;
Deviants.thre= threRes;
Deviants.clusts=Clusts;

YYpred= 10.^( intSlo(1)+intSlo(2)*(Xi))';
resPow= YYi- YYpred' ; % figure; plot(XXi, resPow)

Pows.pred= YYpred';
Pows.obs= YYi;
Pows.res= resPow;
Pows.frex= XXi;

Pows.meanPred= nanmean( YYpred);
Pows.meanObs= nanmean( YYi);
Pows.meanRes= mean( resPow  ); %YYpred

idxPos= ( resPow) >0;%YYpred
Pows.meanResPos= mean( 10.^(Yi(idxPos))- YYpred(idxPos)'); %YYpred(idxPos)'
Pows.cenFrObs= sum(XXi.*YYi)./sum(YYi);
Pows.cenFrRes= sum(XXi.*resPow)./sum(resPow);

%% PLOTS
%%%%%%%%%
%   figure
if ~exist('doPlot','var')
    doPlot=0;
end
if doPlot
    %     figure;
    if ~exist('thisCol','var')
        thisCol= 'k';
    end
    myLW=2;
end
if doPlot ==1 % PLOT ON LOG-LOG
    
    loglog(XXi, YYi, ':','color', thisCol); hold on;
    
    xPlot= XXi; xPlot(boolYdev)= nan;
    loglog(xPlot, 10.^(Yi ), '-',   'LineWidth',myLW, 'color', thisCol); hold on;
    %%% PLOT FIRST SLOPE
    YYpred0= 10.^( intSlo0(1)+intSlo0(2)*(Xi))';
    
    %%% PLOT SECOND (FINAL) SLOPE
    hp(1)= plot( XXi([1 end]), YYpred([1 end]),'-','LineWidth',myLW,'color', thisCol );
    legend(hp,stat.myName);
    
    hp(2)= plot( XXi([1 end]), YYpred0([1 end]),'--','LineWidth',myLW,'color', [.5 .5 .5] );
    
    %%%% DOTS ON PEAKS
    %     xPlot= XXi; xPlot(~boolYdev)= nan;
    %     loglog( xPlot, 10.^(Yi ), 'o','markersize',4, 'color', thisCol); hold on;
    %     loglog(xPlot(Clusts.locs), 10.^(Yi(Clusts.locs) ) ,'.', 'color', thisCol); %.k
    
elseif doPlot==2 % PLOT ON LINEAR SCALES (but use log(X) log(Y) )
    
    plot( (Xi),  (Yi), ':','color', thisCol); hold on;
    
    xPlot=  (Xi); xPlot(boolYdev)= nan;
    plot(xPlot,  (Yi ), '-', 'linewidth',myLW,'color', thisCol);   %'-'
    xPlot=  (Xi); xPlot(~boolYdev)= nan;
    plot( xPlot,  (Yi ), ':', 'linewidth',2, 'color', thisCol ); % 'o' 'markersize',4,
    
    
    Ypred= ( intSlo(1)+intSlo(2)*(Xi));
    hp= plot(  (Xi), Ypred,'-','LineWidth',myLW,'color', thisCol );
    legend(hp,stat.myName);
    
    %%% PLOT FIRST SLOPE
    Ypred0=  ( intSlo0(1)+intSlo0(2)*(Xi))';
    hp(1)= plot( Xi([1 end]), Ypred0([1 end]),'--','LineWidth',myLW,'color', [.5 .5 .5] );
    
    %%%% DOTS ON PEAKS
    %  plot(xPlot(Clusts.locs), Yi(Clusts.locs) ,'.','color',thisCol);
    %     axis tight;
    
end
%%%%%%%%%%%%%%
%% debug PLOTS
%%%%%%%%%%%%%%
debugPlot=0;
if debugPlot ==1
    figure;
    sp(1)=   subplot(121);  % loglog(XXi  , Yi ,':k')
    loglog(XX  , YY ,':k');hold on;      loglog(XXi  , 10.^(YPred0) ,':b')
    plot (XXi(~boolYdev) ,  (YYi(~boolYdev) ),'.b');
    YYthrePks=   10.^( intSlo0(1)+threResPks+intSlo0(2)*(Xi))';
    loglog(XXi  ,YYthrePks ,':k');
    loglog(XXi  , (YYpred) ,':r')
    
    sp(2)=  subplot(122);     plot(XXi, YRes0,':k'); hold on;
    plot (XXi(~boolYdev) ,  (YRes0(~boolYdev) ),'.b'); axis tight;
    plot(xlim,[threRes threRes],':b');  plot(xlim,[threResPks threResPks],':k');
    plot(XXi  , YRes2 ,':r')
    
    try
        subplot(121);
        plot(XXi(locs0), YYi(locs0),'.','color',[.5 .5 .5],'markersize',20 );
        subplot(122);
        plot(XXi(locs0), YRes0(locs0),'.','color',[.5 .5 .5],'markersize',20);
        
        subplot(121);
        plot(XXi(locs), YYi(locs),'.k' ); plot(XXi(locs(1)),YYi(locs(1)),'.r','markersize',20),axis tight;
        subplot(122);
        plot(XXi(locs), YRes0(locs),'.k'); plot(XXi(locs(1)), YRes0(locs(1)),'.r','markersize',20)
        linkaxes(sp,'x');
    catch
    end
end
end


function [FAA,POW_L,POW_R,FAA_cat] = get_FAA(EEG)
% simple function of FAA plugin from EEGLab see details below:
% Compute FAA index from given data
% ==================================
% Computes frontal alpha asymmetry from specified data - channels,
% frequency and latency.
%
% Michael Tesar <tesarm@jcu.cz>
% Ceske Budejovice, 2016
%
% Output is log transformed power difference Right - Left 
% FAA_cat = if sign == 1 right>left. if sign == -1 right<left. if sign == 0
% hemispheres are equal.
%%
% Channel_labels={'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F5' 'F6' 'T3' 'T4' 'P5' 'P6'};
% channels not in smaller chanlist: [1 2 11 12 15 16]
% large channel list left = FP1, right  = FP2
% small channel list left = F3, right  = F4
indLeft = 1;
indRight = 2;
tmpData = mean(EEG.data, 3);
L = tmpData(indLeft, :);
R = tmpData(indRight, :);
%% Compute power spectrum for left channel
WIND = hamming(floor(length(L))/1.5);   % Get 1.5 sec time windows
OVER = floor((length(L))/1.5/2);        % 50% overlap
SIGN = L';                              % Get signal
[s, freqs, t, power] = spectrogram(SIGN, WIND, OVER, [], EEG.srate);
indFreqs = find(freqs>=8 & freqs<=14);
POW_L = power(indFreqs);
%% Compute power spectrum for right channel
WIND = hamming(floor(length(R))/1.5);   % Get 1.5 sec time windows
OVER = floor((length(R))/1.5/2);        % 50% overlap
SIGN = R';                              % Get signal
[s, freqs, t, power] = spectrogram(SIGN, WIND, OVER, [], EEG.srate);
indFreqs = find(freqs>=8 & freqs<=14);
POW_R = power(indFreqs);
%% Compute whole FAA
FAA = mean(abs(log(POW_R)-log(POW_L)));function [freqs,predslope,spectra,frontaltable,centraltable,posteriortable,temporaltable] = resting_state_feature_extractor(EEGfinal)
% resting pipeline. performs automated artifact removal of resting state
% EEG data.
% Contains the following externally developed functions:
% get_FAA.m
% Michael Tesar <tesarm@jcu.cz>
% Ceske Budejovice, 2016
%
% fitPowerLaw3steps.m
% Colombo, M., 2019 the spectral exponent of the resting EEG indexes the presence of consciousness during unresponsiveness induced by propofol, xenon, and ketamine
%
% Inputs:
%   EEG_final = pre-processed EEG structure variable
% Output:
%   freqs = frequencies estimated
%   predslope = exponential fit
%   spectra = relative power spectra
%   tables = roi data tables
%
% Routine composed by Dr. Nicholas Murphy, Baylor College of Medicine, 2021

%% step 1, determine regions of interest
F_ROI={'FP1';'FP2';'F7';'F3';'FZ';'F4';'F8';'CZ'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatef=find(ismember(chan_IDs,F_ROI));

C_ROI={'FC3';'FC1';'FC2';'FC4';'C3';'C1';'CZ';'C2';'C4';'CP1';'CP2';'CPZ'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatec=find(ismember(chan_IDs,C_ROI));

P_ROI = {'P7','P3','PZ','P4','P8','PO9','PO10',...
    'O1','OZ','O2','P1','P5','P2','P6','PO7','PO3','POZ','PO4','PO8'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatep=find(ismember(chan_IDs,P_ROI));

T_ROI = {'T7','T8','TP9','TP10','FT9','FT10','FT7','FT8','TP7','TP8'};
chan_IDs={EEGfinal.chanlocs.labels};
candidatet=find(ismember(chan_IDs,T_ROI));

%% step 2, get subject code
naam = EEGfinal.filename(1:end-4);

%% step 3, detrending & demeaning
dataf = zeros(size(EEGfinal.data(:,:,:)));
for ii = 1:size(EEGfinal.data(:,:),3)
    dataf(:,:,ii) = EEGfinal.data(:,:,ii)-nanmean(EEGfinal.data(:,:,ii),2);
    temp = detrend(dataf(:,:,ii)');
    dataf(:,:,ii) = temp';
    clear temp
end

%% step 4, compute spectra
for u = 1:size(dataf,1)
    if u == 1
        [s,freqs] = pwelch(dataf(u,:),hanning(EEGfinal.srate),EEGfinal.srate/2,(EEGfinal.srate),EEGfinal.srate);
        spectra = zeros(size(dataf,1),length(s));
        spectra(u,:) = s;
    else
        [spectra(u,:),~] = pwelch(dataf(u,:),hanning(EEGfinal.srate),EEGfinal.srate/2,(EEGfinal.srate),EEGfinal.srate);
    end
end
% convert to relative power
spectra = spectra./nansum(spectra,2);

%% step 5, typical spectra (calculate peak alpha separately)
delta = find(freqs<4);
theta = find(freqs>=4 & freqs<7);
alpha = find(freqs>=7 & freqs<14);
beta = find(freqs>=25 & freqs<34);
gamma = find(freqs>=35 & freqs<45);

% Frontal
F.delta = nanmean(nanmean(spectra(candidatef,delta),2),1);
F.theta = nanmean(nanmean(spectra(candidatef,theta),2),1);
F.alpha= nanmean(nanmean(spectra(candidatef,alpha),2),1);
F.beta = nanmean(nanmean(spectra(candidatef,beta),2),1);
F.gamma = nanmean(nanmean(spectra(candidatef,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatef,alpha),1));
F.alpha_peak_power = y;
F.alpha_peak_frq = freqs(alpha(iF));

% Central
C.delta = nanmean(nanmean(spectra(candidatec,delta),2),1);
C.theta = nanmean(nanmean(spectra(candidatec,theta),2),1);
C.alpha= nanmean(nanmean(spectra(candidatec,alpha),2),1);
C.beta = nanmean(nanmean(spectra(candidatec,beta),2),1);
C.gamma = nanmean(nanmean(spectra(candidatec,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatec,alpha),1));
C.alpha_peak_power = y;
C.alpha_peak_frq = freqs(alpha(iF));

% Posterior
P.delta = nanmean(nanmean(spectra(candidatep,delta),2),1);
P.theta = nanmean(nanmean(spectra(candidatep,theta),2),1);
P.alpha= nanmean(nanmean(spectra(candidatep,alpha),2),1);
P.beta = nanmean(nanmean(spectra(candidatep,beta),2),1);
P.gamma = nanmean(nanmean(spectra(candidatep,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatep,alpha),1));
P.alpha_peak_power = y;
P.alpha_peak_frq = freqs(alpha(iF));

% Temporal
T.delta = nanmean(nanmean(spectra(candidatet,delta),2),1);
T.theta = nanmean(nanmean(spectra(candidatet,theta),2),1);
T.alpha= nanmean(nanmean(spectra(candidatet,alpha),2),1);
T.beta = nanmean(nanmean(spectra(candidatet,beta),2),1);
T.gamma = nanmean(nanmean(spectra(candidatet,gamma),2),1);
[y,iF]=max(nanmean(spectra(candidatet,alpha),1));
T.alpha_peak_power = y;
T.alpha_peak_frq = freqs(alpha(iF));

%% step 6, estimate slope
% spectra = nanmean(spectra,1);
% doPlot= 1;
% thisCol= [0 0 1];
% spectra = spectra./sum(spectra);
frBand=[3 45]; % frequencies to use
frBins= dsearchn(freqs, frBand(1) ):dsearchn(freqs,frBand(2));
XX= freqs(frBins); % cutdown frequencies% [global_power_vars_log_orig] = measure_power(global_orig_spectra.log_spectra,global_orig_spectra.freqs);
robRegMeth= 'ols';
for ii = 1:size(spectra,1)
    YY= spectra(ii,frBins); % spectrum cutdown by bins
    [intSlo, ~, Pows, ~,  ~, ~] = fitPowerLaw3steps(XX,YY, robRegMeth,  0, []);
    spectra2(ii,:) = Pows.obs;
    psd_slope(ii,1) = intSlo(2);% exponent
    psd_slope(ii,2) = intSlo(1);% intercept
    if ii == 1
        predslope = zeros(size(psd_slope,1),length(Pows.pred),4);
    end
    predslope(ii,:,1) = Pows.pred; % this is the background fit
    predslope(ii,:,2) = Pows.obs;
    predslope(ii,:,3) = Pows.res;
    predslope(ii,:,4) = Pows.frex;
end
F.psd_slope = nanmean(psd_slope(candidatef,:),1);
C.psd_slope = nanmean(psd_slope(candidatec,:),1);
P.psd_slope = nanmean(psd_slope(candidatep,:),1);
T.psd_slope = nanmean(psd_slope(candidatet,:),1);

%% step 7, frontal alpha asymmetry
[FAA,~,~,FAA_cat] = get_FAA(EEGfinal);

%% Export
frontaltable = table({naam},F.delta,F.theta,F.alpha,F.beta,...
    F.gamma,F.alpha_peak_frq,F.alpha_peak_power,F.psd_slope(1),FAA,FAA_cat);
frontaltable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent','alpha_asymm','alpha_asymm_direction'};
centraltable = table({naam},C.delta,C.theta,C.alpha,C.beta,...
    C.gamma,C.alpha_peak_frq,C.alpha_peak_power,C.psd_slope(1));
centraltable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent'};
posteriortable = table({naam},P.delta,P.theta,P.alpha,P.beta,...
    P.gamma,P.alpha_peak_frq,P.alpha_peak_power,P.psd_slope(1));
posteriortable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent'};
temporaltable = table({naam},T.delta,T.theta,T.alpha,T.beta,...
    T.gamma,T.alpha_peak_frq,T.alpha_peak_power,T.psd_slope(1));
temporaltable.Properties.VariableNames = {'subj_code','delta_power','theta_power'...
    ,'alpha_power','beta_power','gamma_power','alpha_peak_frq','alpha_peak_power',...
    'psd_exponent'};
end





function  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
% %%
% Colombo, M. A. et al. The spectral exponent of the resting EEG indexes the presence of consciousness during unresponsiveness induced by propofol, xenon, and ketamine. NeuroImage 189, 631–644 (2019).
% % fit a line on log-freq vs log-PSD, excluding peaks and their base.
% % USAGE EXAMPLE:
% % first compute the PSD
% epLen= 2* sRate; epShift= 1*srate;numFFT=[];
%  [myPSD,frex]= pwelch( myEEGch  , epLen, epShift,numFFT, sRate);
%  frBand=[1 40];
%  frBins= dsearchn( frex, frBand(1) ):  dsearchn( frex, frBand(2));
% XX= frex(frBins);
% YY= myPSD(frBins);
% robRegMeth= 'ols';
% doPlot= 1;
% thisCol= [0 0 1];
%  [intSlo, stat, Pows, Deviants,  stat0, intSlo0] = fitPowerLaw3steps(XX,YY, robRegMeth,  doPlot, thisCol)
%%%%%%%%% INPUT
% XX ---> are the frequency bins of interest (hz)
% YY ---> are the PSD relative to the bins of interest (mV^2 /Hz)
% robRegMeth --->is the method for robust regression fit (default: ols)
%                type help robustfit for more details
% doPlot  --->    0 (default) means no plot ;
%                 1 is plot XX and YY on log-log scale;
%                 2 is plot log(XX) and log(YY) on linear scale;
% thisCol --->    if plotting, is the rgb triplet
%%%%%%%%% OUTPUT
% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit  i.e. the spectral exponent

% stat is a structure with informations on the 2nd fit, (see 2nd argument of robustfit)

% Pows is a structure with
%           pred: 1*length(frBins) doubles, predicted PSD by the power-law fit
%            obs: 1*length(frBins) doubles, observed PSD    = YY
%            res: 1*length(frBins) doubles, residuals (pred-obs) PSD
%           frex: 1*length(frBins) doubles, observed frex   = XX
%       meanPred: 1 double, mean of the predicted PSD
%        meanObs: 1 double, mean of the observed PSD
%        meanRes: 1 double, mean of the residuals of the PSD
%     meanResPos: 1 double, mean of the residuals of the PSD (exclude negative values)
%       cenFrObs: 1 double, central Frequency of the observed PSD
%       cenFrRes: 1 double, central Frequency of the residual PSD

% Deviants is a structure with
%       res2: 1*length(frBins) doubles,
%        res: 1*length(frBins) doubles,
%        rej:  length(frBins)*1 logicals, 1 if bins are rejected, 0 if kept   for fitting the 2nd power-law
%         fr: min and max freqs considered
%       thre: threshold used for bins adjacent to the peaks
%     clusts: structure with information on the contiguous freq. bins

% stat0 is a structure with informations on the 1st fit, (see 2nd argument of robustfit)

% intSlo(1) intercept of 2nd (final) powerLaw Fit
% intSlo(2) slope of 2nd (final) powerLaw Fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HOW DOES IT WORK?
% transform in log-log and upsample by 4 times
% fit a 1st line, find all residual >0
% small peaks (whose residuals are smaller than 1*mad(residuals)) do not count as peaks
% find clusters of rejected frequency bins
% consider only clusters of bins where there is a peak
% i.e. reject large enough peaks and their base (adjacent residuals >0)
% fit a second line, on the remaining residuals (those closer to a power-law)

if exist('robRegMeth','var')
else
    robRegMeth= 'ols';
end
%% LOG TRANSFORM
%%%%%%%%%%%%%%%%%
% vectorize frex and psd to avoid dimension mismatch later
XX= XX(:)';
YY= YY(:)';

X= log10(XX);
Y= log10(YY);

%%% INTERPOLATE IN LOG-LOG: X, Y --> Xi, Yi
XXi= logspace( X(1), X(end),  length(X)*4); % 2^10*2;--> auto-chosen
Xi= log10(XXi);
Yi= interpn( X, Y, Xi   );% ,'spline'
YYi= 10.^(Yi);

%% STEP 1, FIT 1st LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL 1: minimize Y residuals. OLS +variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ intSlo0, stat0] = robustfit(Xi,Yi,robRegMeth);
YRes0= stat0.resid;
YPred0=  ( intSlo0(1)+intSlo0(2)*(Xi))';

%% FIND DEVIANT RESIDUALS
%%%%%%%%%%%%%%%%%%%%%%%%%
threRes=0;
boolYdev= YRes0 > threRes ;%
%% ONLY CONSIDER DEVIANTS WHERE ALSO POWER-PEAKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%      %(I.E. if power is high, but it decays so that it does not peak (e.g.peak is lower than earliest freq bin) --> not deviant
% peaks are searched on log-power .... (why not on residuals? avoid small bumps on res)
[pks0,locs0] = findpeaks( Yi ); %, 'MINPEAKDISTANCE', 8, 'NPEAKS', 10, 'SORTSTR','descend') ;
% REJECT SMALL PEAKS
% peaks are rejected if their residuals (detrended log-power) are low ... (why not on log-power? there is no meaningful threshold on log-pow)
pks=pks0; locs=locs0;

threResPks= 1* mad(YRes0,1); %    median(YRes0) + 1 *mad(myres,1);
threResPks= max(threResPks, .1);

rejectPks= (YRes0(locs) < threResPks )';
pks(rejectPks)=[];
locs(rejectPks)=[];

% now exclude deviant-residuals that do not contain any peak
Clusts = bwconncomp(boolYdev);
includeClust= false( 1, Clusts.NumObjects);
for cl= 1: Clusts.NumObjects
    thisIdx=Clusts.PixelIdxList{cl};
    thisclust= Xi( thisIdx );
    thisloc= (Xi(locs));
    includeClust(cl)=   sum( ismember ( thisclust , thisloc));
    boolYdev( thisIdx )=includeClust(cl);
end

%% STEP 2 ,2nd line after excluding peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MODEL 1: minimize Y residuals. OLS +variants
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try % try an OLS regression first
    [ intSlo, stat] = robustfit(Xi(~boolYdev),Yi(~boolYdev),robRegMeth);
catch % if not enough points for regression,(too small freq.band or deviations with too much band-width) just do not exclude any peak ans do ols
    [ intSlo, stat] = robustfit(Xi,Yi,'ols');
end
%%%% GET RESIDUALS ALSO AT EXCLUDED POINTS
YPred0= (Xi.*intSlo(2) +intSlo(1));
YRes2= Yi - YPred0;    %     YRes2= stat.resid;


%% STORE
stat.frBins= [ XXi(1)   XXi(end)  ];
stat.myName= [ num2str(round( stat.frBins(1) )) '-'  num2str(round(  stat.frBins(2)) ) ' Hz: '   num2str(intSlo(2)) ]; %num2str(intSlo')
stat.robRegMeth=robRegMeth;

Clusts.includeClust=includeClust;
Clusts.pks=pks;
Clusts.locs=locs;
Clusts.threResPks=threResPks;
Clusts.pks0=pks0;
Clusts.locs0=locs0;
Clusts.rejectPks0=rejectPks;

Deviants.res2= YRes2;
Deviants.res= YRes0;
Deviants.rej= boolYdev;
Deviants.fr= stat.frBins;
Deviants.thre= threRes;
Deviants.clusts=Clusts;

YYpred= 10.^( intSlo(1)+intSlo(2)*(Xi))';
resPow= YYi- YYpred' ; % figure; plot(XXi, resPow)

Pows.pred= YYpred';
Pows.obs= YYi;
Pows.res= resPow;
Pows.frex= XXi;

Pows.meanPred= nanmean( YYpred);
Pows.meanObs= nanmean( YYi);
Pows.meanRes= mean( resPow  ); %YYpred

idxPos= ( resPow) >0;%YYpred
Pows.meanResPos= mean( 10.^(Yi(idxPos))- YYpred(idxPos)'); %YYpred(idxPos)'
Pows.cenFrObs= sum(XXi.*YYi)./sum(YYi);
Pows.cenFrRes= sum(XXi.*resPow)./sum(resPow);

%% PLOTS
%%%%%%%%%
%   figure
if ~exist('doPlot','var')
    doPlot=0;
end
if doPlot
    %     figure;
    if ~exist('thisCol','var')
        thisCol= 'k';
    end
    myLW=2;
end
if doPlot ==1 % PLOT ON LOG-LOG
    
    loglog(XXi, YYi, ':','color', thisCol); hold on;
    
    xPlot= XXi; xPlot(boolYdev)= nan;
    loglog(xPlot, 10.^(Yi ), '-',   'LineWidth',myLW, 'color', thisCol); hold on;
    %%% PLOT FIRST SLOPE
    YYpred0= 10.^( intSlo0(1)+intSlo0(2)*(Xi))';
    
    %%% PLOT SECOND (FINAL) SLOPE
    hp(1)= plot( XXi([1 end]), YYpred([1 end]),'-','LineWidth',myLW,'color', thisCol );
    legend(hp,stat.myName);
    
    hp(2)= plot( XXi([1 end]), YYpred0([1 end]),'--','LineWidth',myLW,'color', [.5 .5 .5] );
    
    %%%% DOTS ON PEAKS
    %     xPlot= XXi; xPlot(~boolYdev)= nan;
    %     loglog( xPlot, 10.^(Yi ), 'o','markersize',4, 'color', thisCol); hold on;
    %     loglog(xPlot(Clusts.locs), 10.^(Yi(Clusts.locs) ) ,'.', 'color', thisCol); %.k
    
elseif doPlot==2 % PLOT ON LINEAR SCALES (but use log(X) log(Y) )
    
    plot( (Xi),  (Yi), ':','color', thisCol); hold on;
    
    xPlot=  (Xi); xPlot(boolYdev)= nan;
    plot(xPlot,  (Yi ), '-', 'linewidth',myLW,'color', thisCol);   %'-'
    xPlot=  (Xi); xPlot(~boolYdev)= nan;
    plot( xPlot,  (Yi ), ':', 'linewidth',2, 'color', thisCol ); % 'o' 'markersize',4,
    
    
    Ypred= ( intSlo(1)+intSlo(2)*(Xi));
    hp= plot(  (Xi), Ypred,'-','LineWidth',myLW,'color', thisCol );
    legend(hp,stat.myName);
    
    %%% PLOT FIRST SLOPE
    Ypred0=  ( intSlo0(1)+intSlo0(2)*(Xi))';
    hp(1)= plot( Xi([1 end]), Ypred0([1 end]),'--','LineWidth',myLW,'color', [.5 .5 .5] );
    
    %%%% DOTS ON PEAKS
    %  plot(xPlot(Clusts.locs), Yi(Clusts.locs) ,'.','color',thisCol);
    %     axis tight;
    
end
%%%%%%%%%%%%%%
%% debug PLOTS
%%%%%%%%%%%%%%
debugPlot=0;
if debugPlot ==1
    figure;
    sp(1)=   subplot(121);  % loglog(XXi  , Yi ,':k')
    loglog(XX  , YY ,':k');hold on;      loglog(XXi  , 10.^(YPred0) ,':b')
    plot (XXi(~boolYdev) ,  (YYi(~boolYdev) ),'.b');
    YYthrePks=   10.^( intSlo0(1)+threResPks+intSlo0(2)*(Xi))';
    loglog(XXi  ,YYthrePks ,':k');
    loglog(XXi  , (YYpred) ,':r')
    
    sp(2)=  subplot(122);     plot(XXi, YRes0,':k'); hold on;
    plot (XXi(~boolYdev) ,  (YRes0(~boolYdev) ),'.b'); axis tight;
    plot(xlim,[threRes threRes],':b');  plot(xlim,[threResPks threResPks],':k');
    plot(XXi  , YRes2 ,':r')
    
    try
        subplot(121);
        plot(XXi(locs0), YYi(locs0),'.','color',[.5 .5 .5],'markersize',20 );
        subplot(122);
        plot(XXi(locs0), YRes0(locs0),'.','color',[.5 .5 .5],'markersize',20);
        
        subplot(121);
        plot(XXi(locs), YYi(locs),'.k' ); plot(XXi(locs(1)),YYi(locs(1)),'.r','markersize',20),axis tight;
        subplot(122);
        plot(XXi(locs), YRes0(locs),'.k'); plot(XXi(locs(1)), YRes0(locs(1)),'.r','markersize',20)
        linkaxes(sp,'x');
    catch
    end
end
end


function [FAA,POW_L,POW_R,FAA_cat] = get_FAA(EEG)
% simple function of FAA plugin from EEGLab see details below:
% Compute FAA index from given data
% ==================================
% Computes frontal alpha asymmetry from specified data - channels,
% frequency and latency.
%
% Michael Tesar <tesarm@jcu.cz>
% Ceske Budejovice, 2016
%
% Output is log transformed power difference Right - Left 
% FAA_cat = if sign == 1 right>left. if sign == -1 right<left. if sign == 0
% hemispheres are equal.
%%
% Channel_labels={'Fp1' 'Fp2' 'F3' 'F4' 'C3' 'C4' 'P3' 'P4' 'O1' 'O2' 'F5' 'F6' 'T3' 'T4' 'P5' 'P6'};
% channels not in smaller chanlist: [1 2 11 12 15 16]
% large channel list left = FP1, right  = FP2
% small channel list left = F3, right  = F4
indLeft = 1;
indRight = 2;
tmpData = mean(EEG.data, 3);
L = tmpData(indLeft, :);
R = tmpData(indRight, :);
%% Compute power spectrum for left channel
WIND = hamming(floor(length(L))/1.5);   % Get 1.5 sec time windows
OVER = floor((length(L))/1.5/2);        % 50% overlap
SIGN = L';                              % Get signal
[s, freqs, t, power] = spectrogram(SIGN, WIND, OVER, [], EEG.srate);
indFreqs = find(freqs>=8 & freqs<=14);
POW_L = power(indFreqs);
%% Compute power spectrum for right channel
WIND = hamming(floor(length(R))/1.5);   % Get 1.5 sec time windows
OVER = floor((length(R))/1.5/2);        % 50% overlap
SIGN = R';                              % Get signal
[s, freqs, t, power] = spectrogram(SIGN, WIND, OVER, [], EEG.srate);
indFreqs = find(freqs>=8 & freqs<=14);
POW_R = power(indFreqs);
%% Compute whole FAA
FAA = mean(abs(log(POW_R)-log(POW_L)));
%% Categorical Output
FAA_cat = sign(FAA);
    





end
%% Categorical Output
FAA_cat = sign(FAA);
    





end
