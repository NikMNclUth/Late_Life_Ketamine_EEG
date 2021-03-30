function [s1avgall,s2avgall,mmnavgall,mmnamps,mmnlats] = mmn_feature_extractor(dirlist,ROI)
% MMN feature extraction pipeline. 
% EEG data.
% Requirements:
%
% Inputs:
%   dirlist = directory list of files to load (string)
%   ROI = name of the channel to load (string)
% Output:
%   s1avgall = standard group data
%   s2avgall = deviant group data
%   mmnavgall = mmn group data
%   mmnamps = amplitudes for mmn at group level
%   mmnlats = latencies for mmn at group level
% Routine composed by Dr. Nicholas Murphy, Baylor College of Medicine, 2021

%% step 1, process group signals from chosen ROI
for iter = 1:length(dirlist)
    try
        load(dirlist{iter})
        for i = 1:EEGfinal.trials
            EEGfinal.data(:,:,i) = detrend(EEGfinal.data(:,:,i)')';
        end
        lcut = 3;%hcut = 30;
        % Low pass filter
        EEGfinal = pop_eegfiltnew(EEGfinal, lcut,[]);
        cond={EEGfinal.event.type}';
        evs=double(cell2mat(cond));
        % deal with potential bad data labelling
        f=find(evs>2);
        evs(f)=[];
        X=EEGfinal.data;
        sr=EEGfinal.srate;
        times=EEGfinal.times;
        bl=(1:26);
        EEGbl=baseline(X,bl,2);
        chanlocs=EEGfinal.chanlocs;
        chan_IDs={chanlocs.labels};
        candidate=find(ismember(chan_IDs,ROI));
        Standard=find(evs==1);
        Deviant=find(evs==2);
        s1=nanmean(EEGbl(candidate,:,Standard),3);
        s2=nanmean(EEGbl(candidate,:,Deviant),3);
        mmn = s2-s1;
        s1avg=nanmean(s1,1);
        s2avg=nanmean(s2,1);
        mmnavg = nanmean(mmn,1);
        s1avgall(:,iter)=gaussmooth_signal(s1avg,10);
        s2avgall(:,iter)=gaussmooth_signal(s2avg,10);
        mmnavgall(:,iter) = gaussmooth_signal(mmnavg,10);
    catch
    end
end

[Mbl_1,~,~,SEbl_1]=groupav(s1avgall,[],2);
[Mbl_2,~,~,SEbl_2]=groupav(s2avgall,[],2);
[Mbl_M,~,~,SEbl_M]=groupav(mmnavgall,[],2);
gMbl_1=gaussmooth_signal(Mbl_1,10);
gMbl_2=gaussmooth_signal(Mbl_2,10);
gMbl_M=gaussmooth_signal(Mbl_M,10);
gSEbl_1=gaussmooth_signal(SEbl_1,10);
gSEbl_2=gaussmooth_signal(SEbl_2,10);
gSEbl_M=gaussmooth_signal(SEbl_M,10);

%% step 2, MMN feature extraction
wind = find(times>=100 & times<=250);
[mmnamps,mmnlats] = min(mmnavgall(wind,:));
mmnlats = times(wind(mmnlats));


end

function [ smoothed_signal ] = gaussmooth_signal( signal,window )
%Gaussian smoothing - simple gaussian smoothing of eeg signal (low pass
%filtering alternative)
%   signal = signal to filter
%   window = number of points in smoothing window
% output = smoothed_signal
g = gausswin(window); % define smoothing parameters
g = g/sum(g); % create kernel
smoothed_signal = conv(signal,g,'same'); % convolve with original signal


end
