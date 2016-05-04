function h=plotEEG(sig,Scale,fs,LABELS,FontSize,Color,Marker,LineStyle,LineWidth)
%h=plotEEG(sig,Scale,fs,LABELS,FontSize,Color,Marker,LineStyle)
% *** History: 19-Mar-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
% *** Related work: L. KORCZOWSKI, M. CONGEDO, C. JUTTEN Single-Trial Classification of Multi-User P300-Based 
% Brain-Computer Interface Using Riemannian Geometry" (IEEE EMBC, 2015)
if nargin<9 || isempty(LineWidth)
   LineWidth=1 ;
end
if nargin<8 || isempty(LineStyle)
   LineStyle='-' ;
end
if nargin<7 || isempty(Marker)
   Marker='none' ;
end
if nargin<6 || isempty(Color)
    Color=[.1 .1 .1];
end
if nargin<5 || isempty(FontSize)
    FontSize=24;
end
sig=sig*1; % *1 for positive on top, -1 for positive down
NbSamples=size(sig,2);
NbDim=size(sig,1);
if nargin<2 || isempty(Scale)
    Scale = abs(max(max(sig,[],2)));
    
end
mi=repmat(-Scale,[size(sig,1) 1]);
    ma=repmat(Scale,[size(sig,1) 1]);
    
if nargin<4 || isempty(LABELS)
    LABELS=(1:NbDim)';
    
end
if nargin<3 || isempty(fs)
    fs=128;
end
LABELS=flipud(LABELS);

t = (0:(NbSamples-1))./fs;
%sig = rand(NbDim,NbSamples);

% calculate shift

shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
shift = repmat(shift,1,NbSamples);

%plot 'eeg' data
h=plot(t,sig-shift,'LineWidth',1.2,'Color',Color,'Marker',Marker,'LineStyle',LineStyle,'linewidth',LineWidth);


% edit axes
set(gca,'ytick',sort(mean(sig-shift,2)),'yticklabel',LABELS,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic');
%xlhand = get(gca,'ylabel')
%set(xlhand,'fontsize',20)

grid on
ylim([ min(min(-shift-Scale*2)) max(max(-shift+Scale*2))])
%drawline