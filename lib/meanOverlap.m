function [Emean Class]=meanOverlap(E,Flash,Window,TAG,Weights,doregression)
%MEANOVERLAP Ensemble average of epochs that may overlap (regression or arithmetic).
%   [Emean, Class] = meanOverlap(E, Flash, Window)
%   [Emean, Class] = meanOverlap(E, Flash, Window, TAG, Weights, doregression)
%
%   If any two consecutive epochs overlap (onset(k)+Window > onset(k+1)),
%   the averages of all classes are estimated jointly by least-squares
%   multivariate regression on the continuous signal (eq. 1.9 of Congedo et
%   al. 2016): Emean = E*T'*pinv(T*T'), where T is the stacked Toeplitz
%   design matrix of all classes, holding the epoch weights. Each class
%   average is then rescaled by j_z = sum(w.^2)/sum(w). Otherwise (or if
%   doregression is 0) it returns the arithmetic mean of the weighted epochs.
%
%   INPUTS :
%       E : continuous EEG, size [nb electrodes x nb samples]
%       Flash : trigger vector of length [nb samples], non-zero at the first
%               sample of each epoch. Only the positions are used (find(Flash));
%               the values are ignored, so pass epoch weights in Weights.
%       Window : scalar, epoch length in samples
%       TAG* : class label of each epoch, vector of length nnz(Flash).
%              Default: all ones (one class).
%       Weights* : weight of each epoch, vector of length nnz(Flash).
%              Default: all ones.
%       doregression* : 1 (default) to use the regression when epochs
%              overlap, 0 to always use the arithmetic mean. The regression
%              is skipped automatically (with a warning) when
%              Window * nb classes > 1024, to limit memory use.
%   *optional
%
%   OUTPUTS :
%       Emean : ensemble average of each class, size
%               [nb electrodes x Window x nb classes], 3rd dimension sorted
%               as in Class
%       Class : the sorted class labels, unique(TAG)
%
%   Example (two classes, 1 s epochs, EEG struct as in the ACSTP tutorial):
%       [ERP, Class] = meanOverlap(EEG.Channels', EEG.Trigger, EEG.Fs, EEG.EpochClass);
%
%   Note: ACSTP.m uses its own nested copy of this function; this file is
%   used by ACSTPshow(EEG) and can be called directly.
%
% *** History: 2015-10-28
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-temporal
% common pattern: A companion method for ERP analysis in the time domain", J. Neurosci. Methods, 2016
%
% see also: ACSTP, ACSTPshow
if nargin<6 || isempty(doregression)
    %by default, it will use the regression to estimate the EAE
    doregression=1;
end

if nargin<5 || isempty(Weights)
    %no weights
    Weights=ones(length(find(Flash)),1);
end

if nargin<4 || isempty(TAG)
    %only one class
    TAG=ones(length(find(Flash)),1);
end
N=size(E,2);
Class=unique(TAG);

indF=find(Flash);%all Flashes (target and non-target)

if (Window*length(Class))>1024, warning(['number of classes (' num2str(length(Class)) ') or length epochs (' num2str(Window) ' samples) are too big, cancelled regression.']);doregression=0;end
    
%check if the trials are overlapping
if any(indF(1:end-1)+Window>indF(2:end)) & doregression  %LS regression
    %build the weighted toeplitz matrix
    %     disp('Regression')
    Toep=[];
    for z=1:length(Class)
        FlashZ=zeros(N,1);
        indFZ=indF(TAG==Class(z));
        FlashZ(indFZ)=Weights(TAG==Class(z)); %include Weights in the Toeplitz matrix
        Toep=[Toep;toeplitz(zeros(1,Window),FlashZ)];% Toeplitz matrix for each class
        WeigthsCalib(z)=(sum(Weights(TAG==Class(z)).^2)/sum(Weights(TAG==Class(z))));
    end
    

    
    Emean=E*Toep'*pinv(Toep*Toep'); %seems faster
    %     Emean=(Toep'\E')'; %(0.19) %multiclass LS-regression for the ensemble average estimation
    
    Emean=reshape(Emean,[size(E,1),Window,length(Class)]);
    for z=1:length(Class)
        Emean(:,:,z)=Emean(:,:,z)*WeigthsCalib(z); %(0.19) Weight correction
    end
else %arithmetic ensemble average
    %     disp('AEA')
    indF=find(Flash);
    for i=1:length(indF)
        epoch(:,:,i)=E(:,indF(i):indF(i)+Window-1)*Weights(i);
    end
    for z=1:length(Class)
        Emean(:,:,z)=mean(epoch(:,:,TAG==Class(z)),3); %(0.19) 1/sum(Weights(TAG==Class(z)))*
    end
end