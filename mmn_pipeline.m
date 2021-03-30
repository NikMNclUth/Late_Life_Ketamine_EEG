function [EEGfinal] = mmn_pipeline(current_dir,datafile)
% mmn pipeline. performs automated artifact removal of mmn
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

%% step 2, get event codes
mmn_evcodes = [1,2,99];
mmn_altcodes = [5,6,103];
alleventtypes = {EEG.event.type}';
evs=double(cell2mat(alleventtypes));
evs2 = unique(evs);
% edit alternative codes
% in this section we search for event codes that are equal to the
% alternative event codes and then replace them with their
% equivalent standard index.
if sum(ismember(evs2,mmn_altcodes))>0
    f1 = find(evs==mmn_altcodes(1));
    f2 = find(evs==mmn_altcodes(2));
    f3 = find(evs==mmn_altcodes(3));
    if ~isempty(f1)
        for ac1 = 1:length(f1);EEG.event(f1(ac1)).type = mmn_evcodes(1);end
    end
    if ~isempty(f2)
        for ac1 = 1:length(f2);EEG.event(f2(ac1)).type = mmn_evcodes(2);end
    end
    if ~isempty(f3)
        for ac1 = 1:length(f3);EEG.event(f3(ac1)).type = mmn_evcodes(3);end
    end
end
clear f1 f2 f3 ac1 evs evs2 alleventtypes iter f fff


%% step 3, filter the data, downsample, and remove line noise
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

%% step4, pre-ICA
% epoch the data
EEG = pop_epoch(EEG, {1,2,99}, [-0.1,0.4], 'verbose', 'no', 'epochinfo', 'yes');
EEG=eeg_checkset(EEG);
% create a backup
% EEGbu = EEG;
% find bad channels using kurtosis
[~, badinds2,~] = pop_rejchan(EEG, 'elec',[1:EEG.nbchan],'threshold',[-3 3],...
    'norm','on','measure','kurt');
tempchans = 1:EEG.nbchan;tempchans(badinds2) = [];
% if bad chans exist remove and find second wave of bad channels
if ~isempty(badinds2)
    [~, badinds3,~] = pop_rejchan(EEG, 'elec',tempchans,'threshold',[-3 3],...
        'norm','on','measure','kurt');
    badinds2 = [badinds2,tempchans(badinds3)];
end
% create index of bad channels for ICA exclusion
badchans = sort(badinds2);
bad_channels_removed= {EEG.chanlocs(badchans).labels};
icchans = 1:length(EEG.chanlocs);
icchans(badchans) = [];
EEG = pop_select(EEG,'channel', icchans);


%% step 5, ICA
% run extended infomax, repeatable random seed, pca to adjust dims to
% matrix rank
EEG = pop_runica(EEG, 'extended',1,'interupt','off','reset_randomseed','off','pca',rank(EEG.data(:,:)));
EEG = eeg_checkset( EEG );
% MARA bad component check
[~,temp,~]=processMARA( EEG,EEG,EEG, [0, 0, 0, 0 , 0] );
mbc = zeros(size(temp.reject.gcompreject));
mbc(temp.reject.MARAinfo.posterior_artefactprob > 0.5) = 1;
if ~isempty(bad_comps)
    bad_comps = sort(unique([bad_comps',find(mbc)]));
    EEGprerej = EEG;
    EEG = pop_subcomp( EEG, bad_comps, 0);
    artifact_ICs=bad_comps;
else
    badcompsmara=find(temp.reject.MARAinfo.posterior_artefactprob > 0.5);
    EEG = pop_subcomp( EEG, badcompsmara, 0);
    artifact_ICs=badcompsmara;
end
%% step 6, post-ICA
% clean up of dim reduced data using statistical thresholding
bch = logical(zeros(size(EEG.data,1),1));
try
    X = EEG.data;
    [X,bTR,~,~]=finaldatacleanup(X,EEG.srate,EEG.chanlocs,bch,3,3);
    badtrials0 = find(bTR==1);
    EEG = pop_select(EEG, 'notrial', badtrials0);
    EEG.data = X;
catch
    [badtrials] = find_bad_trials(EEG,[]);
    EEG = pop_select(EEG, 'notrial', badtrials);
    EEG = eeg_checkset(EEG);
end
EEG = eeg_checkset(EEG);

%% step 7, interpolate bad channels
EEG = pop_interp(EEG, full_selected_channels, 'spherical');
EEG = eeg_checkset(EEG );
EEGfinal = EEG;
EEGfinal.nbchan = EEGfinal.nbchan+1;
EEGfinal.data(end+1,:,:) = zeros(1, EEGfinal.pnts,size(EEGfinal.data,3));
EEGfinal.chanlocs(1,EEGfinal.nbchan).labels = 'initialReference';
EEGfinal = pop_reref(EEGfinal, []);
EEGfinal = pop_select(EEGfinal,'nochannel',{'initialReference'});

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
