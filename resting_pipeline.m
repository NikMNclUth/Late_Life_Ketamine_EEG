function [EEGfinal] = resting_pipeline(current_dir,datafile)
% resting pipeline. performs automated artifact removal of resting state
% EEG data.
% Requirements:
%   EEGLab must be in your path. 
%   The MARA ICA artifact removal plugin must be installed and in the path.
%
% Inputs:
%   current_dir = directory of files to load (string)
%   datafile = name of datafile to load, will be appended to the directory
%   string (string)
% Output:
%   EEGfinal = cleaned EEG data
%
% Routine composed by Dr. Nicholas Murphy, Baylor College of Medicine, 2021
% References:
% EEGLab: A Delorme & S Makeig. EEGLAB: an open source toolbox for analysis of single-trial EEG dynamics, Journal of Neuroscience Methods 134:9-21 (2004)
% MARA: https://irenne.github.io/artifacts/, Irene Winkler, Stefan Haufe and Michael Tangermann. Automatic Classification of Artifactual ICA-Components for Artifact Removal in EEG Signals. Behavioral and Brain Functions, 7:30, 2011. 
%   Irene Winkler, Stephanie Brandl, Franziska Horn, Eric Waldburger, Carsten Allefeld and Michael Tangermann. Robust artifactual independent component classification for BCI practitioners. Journal of Neural Engineering, 11 035013, 2014. 
% faster: Nolan, H., Whelan, R.,&  Reilly, R.B. (2010). FASTER: Fully Automated Statistical Thresholding for EEG artifact Rejection./Journal of Neuroscience Methods, 192/, 152-162.



%% step 1, load data using appropriate routine for file extension
try
    EEG = loadcurry([current_dir, datafile, '/', name],'CurryLocations','True'); % keep true to load channel locations
catch
    EEG = pop_biosig([current_dir, datafile, '/', name]);
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup','/data/rcho/TOOLS/eeglab14_1_1b/functions/resources/Standard-10-5-Cap385_witheog.elp','lookup','/data/rcho/TOOLS/eeglab14_1_1b/plugins/dipfit2.3/standard_BESA/standard-10-5-cap385.elp');
    EEG = eeg_checkset( EEG );
end
chan_IDs_old={EEG.chanlocs.labels};
if sum(ismember(chan_IDs_old,'CB1') | ismember(chan_IDs_old, 'CB2'))== 2;
    %load black cap
    load('/data/smathew/AV-101/data/black_cap_adj2.mat');
    eloc = elocB2;
    % re-label the mastoid electrodes to match brainvision 10-20 labels
    for ii = 1:length(eloc)
        if strmatch('M1',eloc(ii).labels)
            eloc(ii).labels = 'TP9';
        end
        if strmatch('M2',eloc(ii).labels)
            eloc(ii).labels = 'TP10';
        end
    end
else
    eloc = readlocs('/data/smathew/AV-101/data/easycap channel configuration.ced');
end
chanlocs = eloc;
EEG.chanlocs=eloc;
chan_IDs={chanlocs.labels};
EEG = eeg_checkset( EEG );


%% step 2, filter the data, downsample, and remove line noise
lcut = 1;hcut = 50;
% Low pass filter
EEG = pop_eegfiltnew(EEG, [],hcut);
EEG = pop_resample(EEG, 250);
%High-pass filter the data at 1-Hz. Note that EEGLAB uses pass-band edge, therefore 1/2 = 0.5 Hz.
EEG = pop_eegfiltnew(EEG, lcut, []);
EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',[1:EEG.nbchan],'computepower',0,'linefreqs',...
    [60 120 180 240] ,'normSpectrum',0,'p',0.01,'pad',2,'plotfigures',0,'scanforlines',1,'sigtype',...
    'Channels','tau',100,'verb',0,'winsize',4,'winstep',1, 'ComputeSpectralPower','False');
EEG = eeg_checkset(EEG);

%% step 3, remove artifacts using ASR
EEGoriginal = EEG; % you need to create a separate EEG structure for interpolation purposes later on.
EEG = clean_rawdata(EEG, 5, [0.25 0.75], 0.7, 4, 40, 0.3); % recommended = 30
numbc = size(EEGoriginal.data,1)-size(EEG.data,1); % bad channels

%% step 4, Interpolate all the removed channels
EEG = pop_interp(EEG, EEGoriginal.chanlocs, 'spherical');
EEG = eeg_checkset( EEG ); % change this for r21 - dev notes.


%% step 5, Re-reference the data to average
EEG.nbchan = EEG.nbchan+1;
EEG.data(end+1,:) = zeros(1, EEG.pnts);
EEG.chanlocs(1,EEG.nbchan).labels = 'initialReference';
EEG = pop_reref(EEG, []);
EEG = pop_select( EEG,'nochannel',{'initialReference'});
EEG = eeg_checkset( EEG ); % change this for r21 - dev notes.

%% step 6, segment and reject poor epochs
EEG = eeg_regepochs(EEG,'recurrence',1); % use epochs of length 1 second for identifying artifacts
EEG = eeg_checkset( EEG ); % change this for r21 - dev notes.
EEGtemp = EEG;
[badtrials1] = find_bad_trials(EEGtemp,[]);
[~, badtrials2, ~, ~] = eegthresh( EEGtemp.data, EEGtemp.pnts, ...
    1:length(EEGtemp.chanlocs), -500, 500, [EEGtemp.xmin,EEGtemp.xmax], EEGtemp.xmin,EEGtemp.xmax);
[OUTEEG, ~, ~, ~] = ...
    pop_jointprob( EEGtemp, 1, 1:length(EEGtemp.chanlocs), ...
    6, 2, 0, 0, 0);
badtrials3 = find(OUTEEG.reject.rejjp == 1);
badtrials = sort(unique([badtrials1,badtrials2,badtrials3]));
EEG = pop_select(EEG, 'notrial', badtrials);
EEG = eeg_checkset( EEG ); % change this for r21 - dev notes.

%% step 7, Run ICA
if isfield(EEG.etc, 'clean_channel_mask')
    dataRank = min([rank(double(EEG.data(:,:)')) sum(EEG.etc.clean_channel_mask)]);
else
    dataRank = rank(double(EEG.data(:,:)'));
end
EEG = pop_runica(EEG, 'icatype','runica','extended',1,'pca',dataRank);

EEG = eeg_checkset( EEG ); % change this for r21 - dev notes.

%% step 8, Run MARA tp reject bad components
[~,temp,~]=processMARA ( EEG,EEG,EEG, [0, 0, 0, 0 , 0] );
% test version with plotting
% [~,temp,~]=processMARA ( EEG,EEG,EEG, [0, 0, 1,1 , 0] );
mbc = zeros(size(temp.reject.gcompreject));
mbc(temp.reject.MARAinfo.posterior_artefactprob > 0.5) = 1;
ibc = [];
badcomps = sort(unique([bad_comps',find(mbc),ibc]));
EEG = pop_subcomp( EEG, badcomps, 0);
EEGfinal = EEG;
end

function [badtrials] = find_bad_trials(EEGepoch,bad_chans)

data = double(EEGepoch.data);
data(bad_chans,:) = NaN;
[badindsA] = bt_amp(data);
[badindsV] = bt_var(data);
[badindsE] = bt_emd(data);
badtrials = sort(unique([badindsA,badindsV,badindsE]));
end

function [badindsA] = bt_amp(data)
for t = 1:size(data,1)
    for u = 1:size(data,3)
        ampdiffs(t,u) = max(data(t,:,u)) - min(data(t,:,u));
    end
end
ad = nanmean(ampdiffs,1);
adz = (ad-nanmean(ad))./nanstd(ad);
badindsA = find(adz>=3);
end

function [badindsV] = bt_var(data)
ev = nanmean(squeeze(nanvar(data,0,2)));
evz = (ev-nanmean(ev))./nanstd(ev);
badindsV = find(evz>=3);
end

function [badindsE] = bt_emd(data)
emd = zeros(1,size(data,3));
means = squeeze(nanmean(data(:,:),2));
for u = 1:size(data,3)
	emd(u) = nanmean(abs(squeeze(nanmean(data(:,:,u),2)) - means));
end
emdz = (emd-nanmean(emd))./nanstd(emd);
badindsE = find(emdz>=3);
end
