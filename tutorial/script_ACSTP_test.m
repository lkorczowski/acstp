% A TUTORIAL TO ACSTP FOR MATLAB
%
% *** History: 2016-05-04
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2016
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA
% SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for
% ERP Analysis", 2016
% *** Contact: louis.korczowski@gmail.com
% *** Website: louis-korczowski.org
%
%
% see also : ACSTP, ACSTPshow

%%%%%%%%%%%%%%%%%% READ ME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Welcome to the tutorial for using the companion method for Event-Related
% Potential Analysis
%
%%%% INTRODUCTION :
% Features (effects):
% -Regression-based trials averaging (reduce the effect of overlapping ERP)
% -Automatic Subspace reduction (find the best subspace to reject the
% artifacts while keeping the maximum of useful information)
% -Automatic Jitter Correction at the trial level (ERP peaks are sharper)
% -Automatic Trial Weighting (Amplitude correction and rejection of high-energy
% trials often related to artefacts)
%
% Results :
% -less trials required to have a clean ERP estimation
% -removal of artifacts such as EMG, EOG, blink, ... that could have distorded the
% real shape of the ERP with the standard Arithmetic Ensemble Average (AEA).
% -estimation of linear spatio-temporal filter maximizing the signal noise
% ratio of the evoked potential(s) according to some user-defined masks
% (see below).
%
% Critical parameters required from the user and information :
% -a well-chosen temporal mask for the evoked potential of interest
% -a wel-chosen spatial mask for the evoked potential of interest
% The estimated ERP epochs Xhat are returned such as :
% for X the EEG epochs with corrected latencies and with weights W, the ACSTP filter, does:
% Xhat(:,:,k)=As*Bs'*W(k)*X(:,:,k)*Bt*At';
%
% Actual limitations and advices :
% -while the code was heavily tested by the authors, it is still an ongoing
% work and bugs will happen. Please contact us if you have troubles to make
% it works for your data. (louis.korczowski@gmail.com)
% -we recommand to downsample your signal to remove most of the
% non-informative frequencies (e.g. Fs=128Hz should be fine for most ERP
% studies).
% -temporal filter cut-off frequencies should be considered with care as it can distord the
% waves. Unfortunately, this method heavily relies on empirical covariance
% matrices and low-frequencies implies often ill-conditionned
% covariances matrices on short epochs. The ACSTP has shown succesful
% results on zero-phase distorsion band-pass filtered data between [0.5-30] Hz.
% -this code is working for 2 classes (Target versus Non-Target) data but
% the method is the same of any number of classes. (NOT TESTED, please report any
% issue at louis.korczowski@gmail.com)
% 

clc
close all
clear all

% we propose a toy data set. An oddball visual ERP paradigm with only 2
% classes : Target (TA) and Non-Target (NT).
load('EEG_ex.mat')
% You can load your own data and put it in the following format:
% EEG is a structure with
%              Fs: scalar (sample rate in Hz)
%         Trigger: [nb samples x1 ] Trigger channel of '0' with '1' at the start
%                   of each sweep. There are [nb epochs] '1'.
%      EpochClass: [nb epochs x1] class of the sweeps (0 for Non-TARGET, 1
%                   for TARGET).
%        Channels: [nb samples x nb channels] preprocessed continuous EEG recordings
%   NoiseTrigger*: the equivalent of Trigger but for the sweep of the noise.
%                   By default, it takes the same.
% ElectrodesName*: {1 x nb channels} the names of the electrodes (useful
%                  for figures)

%% %%%%%%%%%%%%%%%% 0 SIGNAL CHARACTERISTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N= size(EEG.Channels,2); % number of electrodes
T= size(EEG.Channels,1); % number of samples in the continuous EEG recording
% if your structure is correct, the following function should show the
% ensemble average for TA/NT (i.e. TARGET (1) and Non-TARGET (0)).
figure
ACSTPshow(EEG)
% NOTE ABOUT EEG_ex.mat, a toy dataset :
% OBSERVATION1 
% -a quick overview of the Ensemble Average for these data showns that there
% are some blinks distorting heavily the TA estimation in frontal around 550 ms.
% RECOMMANDATION1
% -do not include the frontal electrodes in the spatial mask

% OBSERVATION2 
% -one can see that the TA evoked potential is mostly in perieto-occipital
% with maximum amplitude between 100-500ms (300ms) in Cz-P7-P3-Pz-P4-P8-01-0Z-02 
% RECOMMANDATION2
% -the aforementioned electrodes should be included in the spatial mask and
% increasing slighly the scope of the temporal filter to 50-550ms is a safe
% bet.
%% %%%%%%%%%%%%%%%% 1 CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ACSTPoptions is a structure with
%      Epoch_size: scalar, the length of the epoch window (in samples)
% LatencyCorr_max: scalar, the maximum of samples for the latency
%                   correction. Set 0 to disable the Latency correction.
% Mask_Electrodes: vector of the selectionned electrodes. Useful for the
%                  automatic subspace selection (BestPz) and latency
%                  correction (Latency).
%       Mask_Time: vector of the selectionned sample. Useful for the
%                  automatic subspace selection (BestPz) and latency
%                  correction (Latency).
% MaxIterLatency*: scalar, the maximum iteration allowed to compute the
%                   latency correction. Default: 10.
%    SubspaceDim*: vector, containing all the subspace dimension(s) (<nb electrodes)
%                  to test in descent order.
%                   By default, it is equal to (nb_channels:-1:(nb_channels/2))
%computeClassLat*: vector, containing all the class tag in which you want
%                   to compute the latency correction. By default, it
%                   computes it for all classes but it could be long (for
%                   instance you can skip the non-target).
%        Weights*: Default: true.
%                  option1(given) [nb epochs x1] vector, containing the weights for each
%                   epoch if it is computed from an external function.
%                  option2 (true/false) boolean. If true (default) the ACSTP compute the
%                   weight for each epoch. If false, the weights are
%                   desactivated (i.e. set to 1).
%        overlap*: Default: true (1). If set to 0, it will consider that
%                   the epochs are not overlapping and it will compute the standard
%                   ensemble average. If the window of your epochs is big (>1024 samples), it could be
%                   also interesting to disable this function to speed up the process.
%
%   *optional
ACSTPoptions.Epoch_size=1*EEG.Fs; %(11) the trials' length will be 1s
ACSTPoptions.LatencyCorr_max=0; %(12) +/- nb shifted samples allowed for the jitter correction 

%%% CRITICAL PARAMETER : MASK
%user parameters for the ACSTP to improve converge :
ACSTPoptions.Mask_Electrodes=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
% exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16) for
% EEG_ex.mat (see section 0 SIGNAL CHARACTERISTS for explanation)
ACSTPoptions.Mask_Time=[floor((0.05)*EEG.Fs):ceil((0.55)*EEG.Fs)]; %(14)* time window (in sample) used for latency calculation and Pz selection
% exemple: 50ms to 550ms (see section 0 SIGNAL CHARACTERISTS for explanation)
%%%

ACSTPoptions.computeClassLat=1; %compute the latencies for class TARGET (1) only.
% It is adviced to skip latencies for Non-Target (0) if you can as this
% process is the slowest of the ACSTP function.

% Subspace Dimension, example 1 (integer arrary in descent order) :
ACSTPoptions.SubspaceDim=[N:-1:N-6]; % find the optimal subspace
% between the number of electrodes N and N-6 (it should always be
% sorted in descent order). Should not be greater than N.

% % Subspace Dimension, example 2 :
%     ACSTPoptions.SubspaceDim=N-3;


%% optionnal test for comparison between some Noise Estimation
% uncomment one of the following section to test the ACSTP with different setup
%
% TEST1 : with complete random epochs for the noise covariance matrices estimation (performance should
% be much worse).
%
%     EEG.NoiseTrigger=size(EEG.Channels,1);
%     while max(EEG.NoiseTrigger)+ACSTPoptions.Epoch_size>size(EEG.Channels,1)
%             EEG.NoiseTrigger=(sort(randi([1,size(EEG.Channels,1)],length(EEG.EpochClass),1)));
%     end

% TEST2: noise covariance matrices estimation from only the class Target
%     EEG.TriggerPos=find(EEG.Trigger);
%     EEG.NoiseTrigger=EEG.TriggerPos(EEG.EpochClass==1);
%
% TEST3: noise covariance matrices estimation from only the class Non-Target
%     EEG.TriggerPos=find(EEG.Trigger);
%     EEG.NoiseTrigger=EEG.TriggerPos(EEG.EpochClass==0);
%
% TEST4: Standard (performance should be the best)
%     EEG.NoiseTrigger=[];

% % optionnal test for input debug
% tmpclass=0
%      EEG.TriggerPos=find(EEG.Trigger);
%      EEG.Trigger=EEG.TriggerPos(EEG.EpochClass==tmpclass);
% EEG.EpochClass(randi(length(EEG.EpochClass),30,1))=3;

%% ACSTP Algorithm
[EEG.EpochsACSTP ACSTPfilter]=ACSTP(EEG,ACSTPoptions);
% show results
close all
% the figure below will show the "Standard" Ensemble Average as well as the
% ACSTP Ensemble Average. The ERP distorsion due to the blinks in frontal 
% is removed with ACSTP
figure
ACSTPshow(EEG,ACSTPoptions,ACSTPfilter)
% the figure below will show the distribution of the latencies and weights.
% please note that with the current parameters, the latencies for NT
% (non-Target) are not estimated, see ACSTPoptions.computeClassLat for
% details.
figure
ACSTPshow(EEG,ACSTPoptions,ACSTPfilter,2)

% eeglab required
if exist('readlocs')
    % if needed, you can find the missing .locs files in ./data/loc/
    % in the github repo and add it in the loc folder of eeglab.
    ACSTPshow(EEG,ACSTPoptions,ACSTPfilter,3)
end

%% TO GO FURTHER
% -try different temporal masks (e.g. 250 to 350ms, a really selective mask
% for P300, expect to better estimate the P300 but lose the N2-P2 complexes).
% -try different spatial masks 
% -increase the maximum sample for jitter correction (e.g. 15).
