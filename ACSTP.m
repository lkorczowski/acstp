function [Xhat ACSTPstruct]=ACSTP(EEG,ACSTPoptions)
%%%% Adaptive Common Spatio-Temporal Pattern
% [Xhat ACSTPstruct]=ACSTP(EEG,ACSTPoptions)
%
% Compute the Adaptive Common Spatio-Temporal Pattern Filter (ACSTP) from the
% continuous recording EEG.
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
% -while the code was heavely tested by the authors, it is still an ongoing
% work and bugs will happen. Please contact us if you have troubles to make
% it works for your data. (louis.korczowski@gmail.com)
% -we recommand to downsample your signal to remove most of the
% non-informative frequencies (e.g. Fs=128Hz should be fine for most ERP
% studies).
% -temporal filter cut-off frequencies should be considered with care as it can distord the
% waves. Unfortunately, this method heavely relies on empirical covariance
% matrices and low-frequencies implies often ill-conditionned
% covariances matrices on short epochs. The ACSTP has shown succesful
% results on zero-phase distorsion band-pass filtered data between [0.5-30] Hz.
%
%
%
%%%% INPUTS :
% ------
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
%
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
%
%%%% OUTPUT :
% ------
% Xhat with corrected latencies, weighted and filtered
% ACSTPstruct is a structure with
%              EA: the ensemble average before ACSTP
%          EAcstp: the ensemble average corrected with latencies, weighted
%                   and filtered + with the effect of overlapping
%                   correction
%     As Bs Bt At: such as Xhat(:,:,k)=As{indClass}*Bs{indClass}'*W(k)*X(:,:,k)*Bt{indClass}*At{indClass}'
%                   Each filter is a cell containing a matrix filter for
%                   each class
%           Class: The tag and the order in which the filter are sorted
%          BestPz: Orders of the best subspace for the ACSTP for the given
%                  Class
%         Weights: [nb epochs x1] the weights of each epoch
%         Latency: [nb epochs x1] the corrected offset of each epoch
%      Epoch_size: scalar, the length of the epoch window (in samples)
%
%
% *** History: 11-Nov-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA
% SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for
% ERP Analysis", 2016
% *** Version: 0.9
% *** Contact: louis.korczowski@gmail.com
% *** Website: louis-korczowski.org
%
% see also : ACSTPshow, script_ACSTP_test

Xhat=[];
ACSTPstruct=struct;
EEG.Channels=double(EEG.Channels);
if ~isfield(ACSTPoptions,'DISPLAY')
    DISPLAY=0; %put '1' to plot the results at the end of the computation or '0' to not.
else
    DISPLAY=ACSTPoptions.DISPLAY;
end
if ~isfield(ACSTPoptions,'overlap')
    ACSTPoptions.overlap=1;
end
DISPLAYtime=0;%set on for debug (optimization)
if( DISPLAYtime) tic; end
%%%%%%%%%%%%%%%%%% INPUT EXTRACTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%check if the row are the electrodes and the columns are the samples
if size(EEG.Channels,1)>size(EEG.Channels,2),EEG.Channels=EEG.Channels';end


if isfield(EEG,'NoiseTrigger'),NoiseTrigger=EEG.NoiseTrigger;else NoiseTrigger=[];end

if (size(EEG.EpochClass,1)==1),EEG.EpochClass=EEG.EpochClass';end %transpose vector if necessary
if (size(EEG.Trigger,1)==1),EEG.Trigger=EEG.Trigger';end %transpose vector if necessary


%Check if the EEG INPUTS are correct
if size(EEG.Channels,2)<size(EEG.Channels,1)
    error(['EEG INPUT ERROR1: the temporal dimension of EEG.Channels is smaller than the spatial dimension.' ...
        ' Please check that EEG.Channels'' size is [nb samples x nb channels].' ]);
end

% check EEG.Trigger Channel
if ~(size(EEG.Channels,2)==length(EEG.Trigger))
    warning('ACSTP WARNING2: Size of EEG.Trigger not equals to nb samples of EEG.Channels')
    if length(EEG.Trigger)==length(EEG.EpochClass)
        %we assume that EEG.Trigger gives the position of the EEG.EpochClass instead of
        %being a trigger channel at the sample rate EEG.Fs
        tmp=zeros(size(EEG.Channels,2),1);
        tmp(EEG.Trigger)=1;
        EEG.Trigger=tmp;
        clear tmp
        warning('EEG.Trigger has been assumed to be the EEG.EpochClass locations instead of a trigger channel at the sample rate EEG.Fs')
    else
        error('ACSTP ERROR3: EEG INPUT ERROR.')
        
    end
end

% check NoiseTrigger Channel (if exist)
if ~isempty(NoiseTrigger)
    if ~(size(EEG.Channels,2)==length(NoiseTrigger))
        warning('ACSTP WARNING2b: Size of EEG.NoiseTrigger not equals to nb samples of EEG.Channels')
        %we assume that NoiseTrigger gives the position of the EEG.EpochClass instead of
        %being a trigger channel at the sample rate EEG.Fs
        tmp=zeros(size(EEG.Channels,2),1);
        tmp(NoiseTrigger)=1;
        NoiseTrigger=tmp;
        clear tmp
        warning('EEG.NoiseTrigger has been assumed to be the EEG.EpochClass locations instead of a trigger channel at the sample rate EEG.Fs')
    end
end


%%%%%%%%%%%%%%%%%% CSTP PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<2
    Window=1*EEG.Fs; %(11) the sweeps window will be 1s
    Delays=4; %(12) +/- nb shifted samples allowed for the jitter correction
    winElec=[7,9,10,11,12,13,14,15,16]; %(13)* electrodes used for latency calculation and Pz selection
    % exemple: Cz(7),P7(9),P3(10),Pz(11),P4(12),P8(13),O1(14),Oz(15),O2(16)
    winTime=[floor((0.05)*128):ceil((0.550)*128)]; %(14)* time window (in sample) used for latency calculation and Pz selection
    % exemple: 50ms to 550ms
else
    Window=ACSTPoptions.Epoch_size;
    if isfield(ACSTPoptions,'LatencyCorr_max')
        Delays=ACSTPoptions.LatencyCorr_max;
        if isempty(Delays)
            Delays=0;
        end
    else
        Delays=0;
    end
    if ~isfield(ACSTPoptions,'Mask_Electrodes')
        winElec=1:size(EEG.Channels,1);%put all electrodes
    elseif isempty(ACSTPoptions.Mask_Electrodes)
        winElec=1:size(EEG.Channels,1);%put all electrodes
    else
        winElec=ACSTPoptions.Mask_Electrodes;
    end
    
    if ~isfield(ACSTPoptions,'Mask_Time')
        winTime=1:Window;%put all samples
    elseif isempty(ACSTPoptions.Mask_Time)
        winTime=1:Window;%put all samples
        
    else
        winTime=ACSTPoptions.Mask_Time;
    end
end
%Check if the ACSTP INPUTS are correct
if any(find(EEG.Trigger)+Window-1>size(EEG.Channels,2))
    error('ACSTP ERROR4: ACSTP INPUT ERROR. The size of the epochs choosen in ACSTPoptions.Epoch_size is too big for the last event in EEG.Trigger.');
end
if ~isempty(NoiseTrigger)
    if any(find(NoiseTrigger)+Window>size(EEG.Channels,2))
        error('ACSTP ERROR4b: ACSTP INPUT ERROR. The size of the epochs choosen in ACSTPoptions.Epoch_size is too big for the last event in EEG.NoiseTrigger.');
    end
end
if any(find(EEG.Trigger)+Window+Delays-1>size(EEG.Channels,2))
    error('ACSTP ERROR5: ACSTP INPUT ERROR. The delays choosen in ACSTPoptions.LatencyCorr_max is too big for the last event in EEG.Trigger.');
    
end
if any(find(EEG.Trigger)-Delays<1)
    error('ACSTP ERROR6: ACSTP INPUT ERROR. The delays choosen in ACSTPoptions.LatencyCorr_max is too big for the first event in EEG.Trigger.');
    
end
if any(max(winElec)>size(EEG.Channels,1))
    error('ACSTP ERROR7: ACSTP INPUT ERROR. ACSTPoptions.Mask_Electrodes incorrect');
    
end
if any(max(winTime)>Window)
    error('ACSTP ERROR8: ACSTP INPUT ERROR. ACSTPoptions.Mask_Time incorrect');
    
end
%%%%%%% ADDITIONNALS STRUCTURES FOR ACSTP (no user input needed)%%%%%%%%%%%
WindowB=Window+2*Delays;
offset=-Delays;
[Xall]=epoch_p300(EEG.Channels,EEG.Trigger,WindowB,offset); %prepare epoch for full window+latency
if ~isfield(ACSTPoptions,'MaxIterLatency')
    maxIter_Latcor=10;
else
    maxIter_Latcor=ACSTPoptions.MaxIterLatency
end
if ~isfield(ACSTPoptions,'SubspaceDim')
    Pzfs=[size(Xall,1):-1:size(Xall,1)/2];%choose the screening of the Pz
else
    if isempty(ACSTPoptions.SubspaceDim)
        Pzfs=[size(Xall,1):-1:size(Xall,1)/2];%choose the screening of the Pz
    elseif max(ACSTPoptions.SubspaceDim)>size(EEG.Channels,1)
        error('ACSTP ERROR9: ACSTP INPUT ERROR. The subspace dimension(s) choosen in ACSTPoptions.SubspaceDim is/are greater than the number of electrodes.');
    elseif min(ACSTPoptions.SubspaceDim)<1
        error('ACSTP ERROR9: ACSTP INPUT ERROR. The subspace dimension(s) choosen in ACSTPoptions.SubspaceDim is/are smaller than 1.');
    else
        Pzfs=ACSTPoptions.SubspaceDim;
    end
end

%{
%the constant offset has been desactivated. If you'll like to add an
offset, please consider modifying EEG.Trigger as well as Epoch_size and Mask_Time
accordingly in ACSTPoptions.

if ~isfield(ACSTPoptions,'EpochsOffset')
    EpochsOffset=0;%offset has been set to zero
else
    EpochsOffset=ACSTPoptions.EpochsOffset;
    disp(['ACSTP INFO: The epochs will be with an offset of: ' num2str(EpochsOffset) ' samples.']);
    if (find(EEG.Trigger,1)-Delays+EpochsOffset<1)
        disp('ACSTP ERROR9: ACSTP INPUT ERROR. ACSTPoptions.EpochsOffset is too big for the first event in EEG.Trigger.');
    end
end
%}

% NOT GOOD what if computeClassLat is empty ????
if Delays>0
    if ~isfield(ACSTPoptions,'computeClassLat')
        ClassLat=unique(EEG.EpochClass);%all the classes
        disp(['ACSTP INFO: Latencies will be computed for class(es): ' num2str(ClassLat') '.']);
        disp('ACSTP INFO: Please consider only classes that need Latency correction to speed up the procedure.');
        disp( 'ACSTP INFO: See setting ACSTPoptions.computeClassLat');
        
    elseif isfield(ACSTPoptions,'computeClassLat')
        ClassLat=ACSTPoptions.computeClassLat;
        disp(['ACSTP INFO: Latencies will be computed only for class(es): ' num2str(ClassLat)]);
        
    end
end

% CHECK the Weights' estimation (true/false/given).
if ~isfield(ACSTPoptions,'Weights') % true (default)
    ComputeWeights=1; %ACSTP will estimate the Weights
    disp(['ACSTP INFO: Weights will be computed: ' num2str(ComputeWeights) ' (1=true/0=false)']);
else
    if(length(ACSTPoptions.Weights)==length(EEG.EpochClass)) %given
        ComputeWeights=0; %ACSTP will NOT estimate the Weights
        Weights=ACSTPoptions.Weights;
        warning(['ACSTP WARNING10: ACSTP INPUT WARNING. ACSTPoptions.Weights has been given and won''t be computed.'...
            ' Please check carefully your weights to avoid a scaling issue in output of the CSTP.']);
    elseif(length(ACSTPoptions.Weights)==1) %set true/false
        ComputeWeights=ACSTPoptions.Weights;%ACSTP does(true)/doesn't(false) estimate the Weights
        disp(['ACSTP INFO: Weights will be computed: ' num2str(ComputeWeights) ' (1=true/0=false)']);
    else %wrong size
        error('ACSTP ERROR11: ACSTP INPUT ERROR. ACSTPoptions.Weights size incorrect');
    end
    
end
if (~exist('Weights','var') && ~ComputeWeights) %if Weights estimation false
    Weights=ones(length(EEG.EpochClass),1);
end
%% ACSTP loop Algorithm

disp(['ACSTP INFO: Latency estimation: ' num2str(~(Delays==0)) ' (1=activated/0=disabled)']);

% Heuristic criterion for the Latency convergence:
CriteriaConvLatencies=1/length(find(EEG.EpochClass))*Delays; %nb latencies correction allowed for convergence.


[ACSTPstruct.EA]=EnsembleAverage(EEG.Channels,EEG.EpochClass,EEG.Trigger,Window);

%EEG.EpochClass=ones(size(EEG.EpochClass));
[Xbarz Class]=meanOverlap(EEG.Channels,EEG.Trigger,Window,EEG.EpochClass,[],ACSTPoptions.overlap);%estimate the arithmetic mean (0.10)
Output{1}=Xbarz; %save the AEA for future visualization
if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
disp((['COMPUTE ACSTP... STATE20CHAR         ' dtime]))
printCurrentState(['START               ' dtime])
iter=1;%for each iteration, a Pz
for Pzf=Pzfs;
    %%
    
    
    %%(A.1) #################### DATA PREPARATION ####################
    X=epoch_p300(EEG.Channels,EEG.Trigger,Window); %get epochs
    
    if ~isempty(NoiseTrigger)
        Xnoise=epoch_p300(EEG.Channels,NoiseTrigger,Window);
    else
        Xnoise=X;
    end
    
    %if needed, we can estimate the noise covariance from differents epochs
    %(e.g. resting state or randomly selected epochs)
    
    
    %%(A.2) #################### LATENCY INITIALIZATION ####################
    
    LatencyCorrection(:,iter)=zeros(length(EEG.EpochClass),1); %initialization latency at zero
    
    %%(A.3) #################### WEIGHTS INITIALIZATION ####################
    if (ComputeWeights)
        Weights=WeightsEstimation(X,EEG.EpochClass);
    end
    
    %%(A.4) #################### Xbarz INITIALIZATION ####################
    
    [Xbarz Class]=meanOverlap(EEG.Channels,EEG.Trigger,Window,EEG.EpochClass,Weights,ACSTPoptions.overlap);
    
    if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
    %%(B) #################### ACSTP Initialization ####################
    printCurrentState(['Pzf ' num2str(Pzf) 'ACSTP          ' dtime])
    [Bs Bt As At Class eigV]=CSTP(X,EEG.EpochClass,Xbarz,Weights,winElec,winTime,Pzf,Xnoise);
    
    %%(3) #################### Filtered Ensemble Average Initialization ####################
    Zbarz{iter}=applyCSTP(Xbarz,Bs,Bt,As,At,Class);
    
    Output{2}{iter}=Zbarz{iter}; %save after CSTP (initialized)
    Output{3}{iter}={Bs Bt As At};%save CSTP (initialized)
    
    if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
    printCurrentState(['Pzf ' num2str(Pzf) 'WEIGHT         ' dtime])
    
    %%(3) #################### WEIGHTS ESTIMATION ####################
    if (ComputeWeights)
        Weights=WeightsEstimation(X,EEG.EpochClass,Bs,Bt,As,At);
    end
    
    %%(4) #################### LATENCY CORRECTION LOOP ####################
    if Delays>0 %ONLY IF THE DELAYS NEED TO BE COMPUTED
        for indCl=1:length(ClassLat)
            class_for_lat=ClassLat(indCl);
            classLat_indices=(EEG.EpochClass==class_for_lat);%just compute latencies for class target
            tempoLatency=[];
            STOPlat=1;
            iterLat=1;
            if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
            printCurrentState(['Pzf ' num2str(Pzf) 'LATEN' num2str(class_for_lat) '         ' dtime])
            
            LatencyCorrection(classLat_indices,iter)=zeros(size(LatencyCorrection(classLat_indices,iter))); %initialize latencies
            while STOPlat %converge criteria
                
                [tempoLatency Crit_Trace]=LatencyEstimation(Xall(:,:,classLat_indices),X(:,:,classLat_indices),EEG.EpochClass(classLat_indices),Weights(classLat_indices),Bs(end),Bt(end),As(end),Bt(end),winElec,winTime);
                
                Conv(iter,iterLat)=ConvergenceLatencies(tempoLatency,LatencyCorrection(classLat_indices,iter));
                if (iterLat>1 && Conv(iter,iterLat)<CriteriaConvLatencies) || iterLat>=maxIter_Latcor
                    STOPlat=0;
                end
                LatencyCorrection(classLat_indices,iter)=tempoLatency;
                TriggerCorrected=CorrectLatency(EEG.Trigger,LatencyCorrection(:,iter));
                X=epoch_p300(EEG.Channels,TriggerCorrected,Window);
                iterLat=iterLat+1;
                
            end
            
        end
        
        
        %%(4) #################### CORRECT LATENCIES ####################
        TriggerCorrected=CorrectLatency(EEG.Trigger,LatencyCorrection(:,iter));
        X=epoch_p300(EEG.Channels,TriggerCorrected,Window);
        
        [Xbarz Class]=meanOverlap(EEG.Channels,TriggerCorrected,Window,EEG.EpochClass,Weights,ACSTPoptions.overlap);
        
        % Xbarz is then the EAE averaged with weighted epochs and with
        % corrected lattencies
        
        if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
        printCurrentState(['Pzf ' num2str(Pzf) 'FINALI         ' dtime])
        
    else
        if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
        
        printCurrentState(['NOLAT_esti         ' dtime])
    end
    
    %%(3) #################### ACSTP FINAL ####################
    
    [Bs Bt As At Class eigV fullBs{iter} fullBt{iter} fullAs{iter} fullAt{iter}]=CSTP(X,EEG.EpochClass,Xbarz,Weights,winElec,winTime,Pzf,Xnoise);
    
    %%(3) #################### FINALIZATION ####################
    % for visualization after weights, latency correction and CSTP.
    Zbarz{iter}=applyCSTP(Xbarz,Bs,Bt,As,At,Class);
    % output will be usefull to find the best subspace-dimension
    Output{2}{iter}=Zbarz{iter};
    Output{3}{iter}={Bs Bt As At};
    Output{4}{iter}=Weights;
    %end
    iter=iter+1;
    
    
end
% find the best Pz for each class
ACSTPstruct.EAoverlap=Output{1};%standard arithmetic ensemble average
for indClass=1:size(Output{2}{1},3)
    [Pzind SNR]=best_Pz(Output{2},winElec,winTime,indClass);
    ACSTPstruct.Bs{indClass}=Output{3}{Pzind}{1}{indClass};
    ACSTPstruct.Bt{indClass}=Output{3}{Pzind}{2}{indClass};
    ACSTPstruct.As{indClass}=Output{3}{Pzind}{3}{indClass};
    ACSTPstruct.At{indClass}=Output{3}{Pzind}{4}{indClass};
    ACSTPstruct.EAcstp(:,:,indClass)=mean(Output{2}{Pzind}(:,:,indClass),3);
    ACSTPstruct.BestPz(indClass)=Pzfs(Pzind);
end
ACSTPstruct.Latency=LatencyCorrection(:,Pzind);
ACSTPstruct.Weights=Output{4}{Pzind};
ACSTPstruct.Class=Class;
ACSTPstruct.fullBs=fullBs{Pzind};
ACSTPstruct.fullBt=fullBt{Pzind};
ACSTPstruct.fullAs=fullAs{Pzind};
ACSTPstruct.fullAt=fullAt{Pzind};
ACSTPstruct.PzSNR=SNR;

% do the final filtering from initial data from the optimal
% SubspaceDimension
TriggerCorrected=CorrectLatency(EEG.Trigger,ACSTPstruct.Latency);
X=epoch_p300(EEG.Channels,TriggerCorrected,Window); %apply latency correction
Xw=applyWeights(X,Output{4}{Pzind}); %apply the weights
Xhat=applyCSTP(Xw,ACSTPstruct.Bs,ACSTPstruct.Bt,ACSTPstruct.As,ACSTPstruct.At,EEG.EpochClass); % apply the ACSTP

if(DISPLAYtime) dtime=[' t(' num2str(toc) ')']; else dtime=''; end
printCurrentState(['ACSTP has finished  ' dtime])

% ACSTPstruct is a structure with
% EnsembleAverage: the ensemble average corrected with latencies, weighted
%                   and filtered + with the effect of overlapping
%                   correction
%     As Bs Bt At: such as Xhat(:,:,k)=As{indClass}*Bs{indClass}'*W(k)*X(:,:,k)*Bt{indClass}*At{indClass}'
%                   Each filter is a cell containing a matrix filter for
%                   each class
%           Class: The tag and the order in which the filter are sorted
%         Weights: [nb epochs x1] the weights of each epoch
%         Latency: [nb epochs x1] the corrected offset of each epoch



%% find best Pz and PLOT
if DISPLAY
    disp('WARNING : display not available, please see ACSTPshow.m')
end


% include all needed functions
% *** History: 5-June-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% for the full commented versions please contact me :
% louis.korczowski [at] gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [Bs Bt As At Class eigV fullBs fullBt fullAs fullAt]=CSTP(X,Y,P,W,winElec,winTime,Pzf,Xnoise)
        if nargin<8
            Xnoise=X;
            
        end
        if nargin<5 || isempty(winElec) || isempty(winTime)
            winElec=1:size(X,1);
            winTime=1:size(X,2);
        end
        
        
        if nargin<4 || isempty(W);
            W=ones(1,size(X,3));
        end
        
        
        if nargin<3 || isempty(P) %check and verify the average P300
            P=[];Bool=1; %non overlapping method
        else
            Bool=0; %use already computed mean P300
            if ~iscell(P) %check if the ERP is with a cell structure
                tmpP=P;
                nbClasses=size(P,3);
                clear P
                for z=1:nbClasses
                    P{z}=tmpP(:,:,z);%one cell for each class
                end
            end
        end
        % Algo Parameters
        eigWhite1=1:size(X,1);
        eigWhite2=1:size(X,2);
        
        if nargin<2 || isempty(Y)
            Y=ones(size(X,3),1);
        end
        
        
        
        Class=unique(Y);
        nbClasses=length(Class);
        Nbe=size(X,1);
        
        if Bool==1
            for z=1:nbClasses
                %estimate the EAE if needed
                P{z}=EnsembleAverage(X(:,:,Y==Class(z))); %(0.10)
            end
        end
        
        %%%%%%%%%%%%%%%% START COMPUTING CSPT %%%%%%%%%%%%%%%%%%%%%%%%
        %compute Sample covariance matrices for each sweep
        for k=1:size(Xnoise,3)
            %             Cs_k(:,:,k)=Xnoise(:,:,k)*Xnoise(:,:,k)'/(size(Xnoise,2)-1); %scm
            Cs_k(:,:,k)=cov(Xnoise(:,:,k)');
            %             Ct_k(:,:,k)=Xnoise(:,:,k)'*Xnoise(:,:,k)/(size(Xnoise,1)-1); %scm
            Ct_k(:,:,k)=cov(Xnoise(:,:,k));
        end
        
        %mean covariance matrices
        Cs=mean(Cs_k,3); %(0.15)
        Ct=mean(Ct_k,3); %(0.15)
        
        %%%%%%%%%%%%%%%%%%%%% BILINEAR WHITENING %%%%%%%%%%%%%%%%%%%%%%%%
        [Ux Ws]=eig(Cs); %(0.18)
        [Vx Wt]=eig(Ct); %(0.18)
        
        [Ws,ind] = sort(diag(Ws),'descend');
        sumWs=cumsum(Ws);
        indWhite1tmp=find(sumWs>max(sumWs)*(1-1e-12),1);
        % suppress only the smallest spatial eigenvectors to keep (1-1e-12)% of the
        % power of the signal
        
        if max(eigWhite1)>indWhite1tmp
            eigWhite1=1:indWhite1tmp;
        end
        if any(Ws(eigWhite1)<=0); disp('Warning negative SPATIAL eigenvalue(s)'); end
        
        Ux = Ux(:,ind);%sort spatial eigenvector
        Ux = Ux(:,eigWhite1);%whitening reduction
        [Wt,ind] = sort(diag(Wt),'descend');
        sumWt=cumsum(Wt);
        indWhite2tmp=find(sumWt>max(sumWt)*(1-1e-12),1);
        % suppress only the smallest temporal eigenvectors to keep (1-1e-12)% of the
        % power of the signal
        
        if max(eigWhite2)>indWhite2tmp
            eigWhite2=1:indWhite2tmp;
        end
        if any(Wt(eigWhite2)<=0); disp('Warning negative TEMPORAL eigenvalue(s)'); end
        
        Vx = Vx(:,ind); %sort temporal eigenvector
        Vx = Vx(:,eigWhite2);%whitening reduction
        
        Fsp = (sqrt(pinv(diag(Ws(eigWhite1)))) * Ux')'; %(0.19)
        Ft =(sqrt(pinv(diag(Wt(eigWhite2)))) * Vx')'; %(0.19)
        Gs =Ux*sqrt((diag(Ws(eigWhite1)))) ; %(0.20)
        Gt =Vx*sqrt((diag(Wt(eigWhite2)))) ; %(0.20)
        %%%%%%%%%%%%%% WHITENING MATRICES COMPUTED %%%%%%%%%%%%%%%
        
        
        %% %%%%%%%%%%%% COMPUTE CSTP FILTERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for z=1:length(P) %for each class
            Z{z}=Fsp'*P{z}*Ft; %EAE whitening 
            [Uz{z},Wz{z},Vz{z}]=svd((Z{z})); 
                                    
            if exist('Pzf') %if Pzf has be given by the user (adviced)
                eigV=1:Pzf;
            else
                %%%%%%%%%%%%%% FIND OPTIMAL Pz (optional) %%%%%%%%%%%%%%
                Fdenum=norm(Z{z},'fro').^2;
                for Pz=Nbe*0.75:Nbe-1
                    
                    Uz_r{z}=Uz{z}(:,1:Pz);
                    Vz_r{z}=Vz{z}(:,1:Pz);
                    Zz=Uz_r{z}*Wz{z}(Pz,Pz)*Vz_r{z}';
                    Fnum=norm(Zz(winElec,winTime),'fro').^2;
                    indPz(Pz)=Fnum/Fdenum;
                    
                end
                [~, Pze]=max(indPz);
                eigV=1:Pze;
                %%%%%%%%%%%%%% FIND OPTIMAL Pz END (optional) %%%%%%%%%%%%%%
            end
            
            
            if isempty(eigV) %if no optimal Pz
                eigV=1;
            end
            if max(eigV)>size(Uz{z},1) | max(eigV)>size(Vz{z},1)
                error('Subspace dimension (Pz) too high, estimated covariance matrix is lower rank. Try reduce the SubspaceDim')
            end
            Uz_r{z}=Uz{z}(:,eigV); %subspace reduction
            Vz_r{z}=Vz{z}(:,eigV); %subspace reduction
            
            
            % Compute Bs, Bt, As and At with both whitening and subspace reduction
            Bs{z}=Fsp*Uz_r{z};%
            Bt{z}=Ft*Vz_r{z};%
            As{z}=Gs* Uz_r{z}; %
            At{z}=Gt* Vz_r{z}; %
            
            fullBs{z}=Fsp*Uz{z}; %only for display purpose
            fullBt{z}=Ft*Vz{z}; %only for display purpose
            fullAs{z}=Gs* Uz{z}; %only for display purpose
            fullAt{z}=Gt* Vz{z}; %only for display purpose
        end
    end

    function Xw=applyCSTP(X,Bs,Bt,As,At,Y)
        if nargin<6
            Y=ones(1,size(X,3));
        end
        Classes=unique(Y);
        for k=1:size(X,3)
            c=find(Y(k)==Classes);
            Xw(:,:,k)=As{c}*Bs{c}'*X(:,:,k)*Bt{c}*At{c}';
        end
    end

    function [P Class]=EnsembleAverage(E,Y,Flash,Window,W)
        if nargin>2
            if isempty(Y)
                Y=ones(length(find(Flash)),1);
            end
            if isempty(Flash) || isempty(Window) || length(find(Flash))~=length(Y) || nargin<4
                error('Invalid Inputs in EnsembleAverage please check Flash or Window')
            end
            Xep=epoch_p300(E,Flash,Window);
            
        end
        
        if nargin<3
            Window=size(E,2);
            Xep=E;
            if nargin<2
                Y=ones(1,size(Xep,3));
            end
        end
        
        if nargin>4
            Xep=applyWeights(Xep,W);
        end
        
        
        
        
        Class=unique(Y);
        P=zeros(size(E,1),Window,length(Class));
        for z=1:length(Class)
            c=find(Y==Class(z));
            P(:,:,z)=mean(Xep(:,:,c),3);
        end
        
    end

    function [Emean Class]=meanOverlap(E,Flash,Window,TAG,Weights,doregression)
        if nargin<6 || isempty(doregression)
            doregression=1;
        end
        if nargin<5 || isempty(Weights)
            Weights=ones(length(find(Flash)),1);
        end
        
        if nargin<4 || isempty(TAG)
            TAG=ones(length(find(Flash)),1);
        end
        N=size(E,2);
        Class=unique(TAG);
        
        indFlash=find(Flash);%all Flashes (target and non-target)
        if (Window*length(Class))>1024, warning(['number of classes (' num2str(length(Class)) ') or length epochs (' num2str(Window) ' samples) are too big, cancelled regression.']);doregression=0;end
        
        %check if the trials are overlapping
        if any(indFlash(1:end-1)+Window>indFlash(2:end)) & doregression  %LS regression
            %build the weighted toeplitz matrix
            Toep=[];
            for iz=1:length(Class)
                FlashZ=zeros(N,1);
                indFZ=indFlash(TAG==Class(iz));
                FlashZ(indFZ)=Weights(TAG==Class(iz)); %include Weights in the Toeplitz matrix
                Toep=[Toep;toeplitz(zeros(1,Window),FlashZ)];% Toeplitz matrix for each class
                WeigthsCalib(iz)=(sum(Weights(TAG==Class(iz)).^2)/sum(Weights(TAG==Class(iz))));
            end
            
            Emean=E*Toep'*pinv(Toep*Toep');
            
            Emean=reshape(Emean,[size(E,1),Window,length(Class)]);
            for iz=1:length(Class)
                Emean(:,:,iz)=Emean(:,:,iz)*WeigthsCalib(iz); %(0.19) Weight correction
            end
        else %arithmetic ensemble average
            %     disp('AEA')
            indFlash=find(Flash);
            for idf=1:length(indFlash)
                epoch(:,:,idf)=E(:,indFlash(idf):indFlash(idf)+Window-1)*Weights(idf);
            end
            for iz=1:length(Class)
                Emean(:,:,iz)=mean(epoch(:,:,TAG==Class(iz)),3); %(0.19) 1/sum(Weights(TAG==Class(iz)))*
            end
        end
    end

    function [X isBad] = epoch_p300(Signal,Target,window,offset,Artefacts)
        if (nargin<5 || isempty(Artefacts))
            Artefacts=zeros(1,length(Target));
        end
        if (nargin < 4 || isempty(offset))
            offset = 0;
        end
        ixTarget = find(Target);
        while (ixTarget(end)+window-1)>size(Signal,2)
            ixTarget(end) = [];
        end
        NTarget = length(ixTarget);
        X = zeros(size(Signal,1),window,NTarget);
        isBad=zeros(1,NTarget);
        for i=1:NTarget
            X(:,:,i) = Signal(:,ixTarget(i)+offset:ixTarget(i)+window-1+offset);
            isBad(i)= length(find(Artefacts(ixTarget(i)+offset:ixTarget(i)+window-1+offset)>0));
        end
    end

    function [W,Xhatk]=WeightsEstimation(X,Y,Bs,Bt,As,At)
        
        W=ones(1,size(X,3));
        if nargin<2
            Y=ones(size(X,3),1);
        end
        Classes=unique(Y); % find all the classes
        
        if nargin<3 %%no CSTP weights (initialization)
            for k=1:size(X,3) % for each epoch
                W(k)=1/norm(X(:,:,k),'fro'); %compute the initilized weights
            end
        else
            Xhatk=zeros(size(X));
            for k=1:size(X,3) % for each epoch
                %Xk=X(:,:,k);
                c=(Y(k)==Classes); % class of the current epoch
                Xhatk(:,:,k)=As{c}*Bs{c}'*X(:,:,k)*Bt{c}*At{c}'; % signal estimation
                W(k)=norm(Xhatk(:,:,k) ,'fro')/norm( X(:,:,k)-Xhatk(:,:,k),'fro'); % SNR estimation
            end
            
        end
        % Weights normalization
        for z=1:length(Classes) % for each class
            Wz=W(Y==Classes(z));
            W(Y==Classes(z))=length(Wz)*Wz/sum(Wz); % normalize each class weights
        end
    end

    function FlashCorrected=CorrectLatency(Flash,Latencies,RND)
        if nargin<3
            RND=1:length(find(Flash)); %no random seed
        end
        if size(Flash,1) == 1, Flash = Flash'; end;
        FlashCorrected=zeros(size(Flash));
        indF=find(Flash);
        if length(indF)~=Latencies
            indF(RND)=indF(RND)+Latencies;
            indFcorrected=indF;
        else
            indFcorrected=indF+Latencies;
        end
        FlashCorrected(indFcorrected)=1;
    end

    function Conv=ConvergenceLatencies(Lat1,Lat2)
        Conv=mean(abs(Lat1-Lat2));
    end

    function [Latency Crit_TraceN]=LatencyEstimation(Xall,X,Y,W,Bs,Bt,As,At,winElec,winTime)
        
        WindowSize=size(X,2);
        SizeMax=size(Xall,2);
        finalInd=SizeMax-WindowSize;
        delay0=floor(finalInd/2)+1;
        if nargin<3
            Y=ones(size(X,3),1);
        end
        if nargin<4
            W=ones(size(Y));
        end
        Xweighted=applyWeights(X,W);
        Xallw=applyWeights(Xall,W);
        Classes=unique(Y);
        
        if nargin<5 || isempty(Bs)
            for i=1:length(Classes)
                Bs{i}=eye(size(X,1));
                Bt{i}=eye(size(X,2));
            end
        end
        if nargin<7 || isempty(As)
            for i=1:length(Classes)
                
                As{i}=eye(size(X,1));
                At{i}=eye(size(X,1));
            end
        end
        if nargin<9
            winElec=1:size(X,1);
        end
        if nargin<10
            winTime=1:size(X,2);
        end
        
        for k=1:size(Xall,3)
            
            z=find(Classes==Y(k));
            Xother=Xweighted(:,:,[1:k-1 k+1:end]);
            Tag=Y([1:k-1 k+1:end]);
            Xother=Xother(:,:,Tag==Classes(z));
            Pother=mean(Xother,3); %P300 estimation
            Ybar=Bs{z}'*Pother*Bt{z};            
            %%%%%%% WARNING: SLOW, this section will be optimized (in Fourier) %%%
            %%%%%%% Find the maximum correlation to estimate the latencies 
            for lat=1:finalInd+1
                Xtmp=Xallw(:,lat:lat+WindowSize-1,k);
                Xwk=zeros(size(Xtmp));
                Xwk(winElec,winTime)=Xtmp(winElec,winTime);
                Yk=Bs{z}'*Xwk*Bt{z};                
                Criterion_Trace(k,lat)=trace(Ybar*Yk');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            Crit_TraceN(k,:)=Criterion_Trace(k,:)-min(Criterion_Trace(k,:));%normalize criteria
            [Pik Peaks]=findpeaks(Crit_TraceN(k,:)); %criteria
            
            if isempty(Peaks)
                Peaks=delay0;
            else
                indP=find(Pik>max(Pik)*0.66); %accept only peak of 66% of max (heuristic criteria)
                Peaks=Peaks(indP);
            end
            
            [trash IndM]=min(abs(Peaks-delay0)); %choose closest peak from allowed ones
            Latency(k)=Peaks(IndM)-delay0;
            
        end
        
        Latency=Latency';
        
    end

    function [Pzind SNR]=best_Pz(Zbarz,Electrodes,Window,indClass)
        % it is important that the Pz are sorted in descend order
        if nargin<4
            indClass=size(Zbarz{1},3); %take the "target" class that should be the last
        end
        for i=1:length(Zbarz)
            noise=Zbarz{i}(:,:,indClass);
            signal=noise(Electrodes,Window);
            SNR(i)=norm(signal,'fro')/norm(noise,'fro');
        end
        if(length(Zbarz)==1),Pzind=1;return;return;end
        
        SNR = smooth(SNR,3);% smooth SNR (moving average 3 points)
        SNR=SNR-min(SNR);
        if length(Zbarz)>3 & length(unique(SNR))>1
            [maxs INDtmp]=findpeaks(SNR);
            %  plot(SNR);hold on;plot(INDtmp,maxs,'ro')
            
        else
            INDtmp=[];
        end
        % remove the first and last for the computation of the best Pz 
        [tmpPz Pzind]=max(SNR(2:end-1));Pzind=Pzind+1;
        
        if ~isempty(INDtmp)
            IND=INDtmp(maxs>(tmpPz*0.66));
            if ~isempty(IND)
                Pzind=IND(1);%take the greatest valid Pz
            end
        end
        
        %check that you don't take the maximum or minimum Pz
        if length(Zbarz)>=3
            if Pzind==length(Zbarz),Pzind=Pzind-1;end
            if Pzind==1,Pzind=Pzind+1;end
        end
        
    end

    function printCurrentState(INFO,previousLine)
        if nargin<2
            previousLine=true;
        end
        str = [sprintf(INFO)];
        if previousLine
            str = [sprintf(char(repmat(double('\b'),1,length(str)+1))) str];
        end
        disp(str);
    end

    function Xw=applyWeights(X,W)
        
        for k=1:length(W)
            Xw(:,:,k)=X(:,:,k)*W(k);
        end
        
    end

end