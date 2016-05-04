function [Emean Class]=meanOverlap(E,Flash,Window,TAG,Weights,doregression)
%[Emean Class]=meanOverlap(E,Flash,Window,TAG,Weights)
% Ensemble Averaging of the EEG data E with
% overlapping samples for z classes using multivariate
%
%   INPUTS :
%       E : the EEG data matrix,  size [nb electrodes x nb samples]
%       Flash : a zero vector with only 1 at the initial sample for each epoch, size [nb samples]
%               for weighted epochs, replace 1 by the weight of each epoch
%       Window : a scalar giving the length of the epochs
%       TAG* : in case of several classes, TAG is a vector of the length of the
%           number of 1 contained in Flash.
%       Weights* : a vector of the size of TAG for weightened epochs
%   *optional
%
%   Output :
%       Emean : the matrix containing the ensemble average for each class,
%           with size [nb electrodes x (Window) x(nb classes)] with the 3rth dim sorted in
%           order with the output Class
%       Class : this output give the order of each class in Emean
%
%
% *** History: 2015-10-28
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
%
% see also: script_meanOverlap, tuto_meanOverlap
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