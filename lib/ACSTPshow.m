function ACSTPshow(EEG,ACSTPoptions,ACSTPfilter,Type)
% ACSTPshow(EEG,ACSTPoptions,ACSTPfilter)
%
% Show ensemble average Target(TA) Non-Target(NT)
% and their difference (TA-NT) before and after ACSTP.
%
% ACSTPshow(EEG,ACSTPoptions,ACSTPfilter,Type)
%
% With Type a scalar that defines the plot to do :
%      1: ensemble average Target(TA) Non-Target(NT) (default)
%      2: latencies and Weights
%      3: Filters components topoplot (WARNING : will create one figure of each class, eeglab REQUIRED)
%
% *** History: 29-Oct-2015
% *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
% *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
if isfield(EEG, 'ElectrodesName')
    Mapping=(EEG.ElectrodesName');
else
    Mapping=[];
end

if nargin<4
    Type=1;
end
if nargin<2 % just plot Enseble Average
    offset=0;
    FontSize=10;
    xticklabels = (offset/EEG.Fs):.25:(EEG.Fs+offset)/EEG.Fs;
    xticks = linspace(0, (1.01)*(EEG.Fs)/EEG.Fs, numel(xticklabels));
    set(gcf,'units','normalized','outerposition',[0.1 0.1 0.2 .9])
    set(gcf, 'InvertHardCopy', 'off');
    
    [ERP Class]=meanOverlap(EEG.Channels',EEG.Trigger,EEG.Fs,EEG.EpochClass);
    Scale=FindScaling(ERP(:,:,2));
    for IND=1:2
        if IND==1
            LineStype='-';
            LW=1;
            Color=[0.9,0.1,0.1];
        else
            LineStype='-';
            LW=2;
            Color=[0.1,0.1,0.9];
        end
        plotEEG(ERP(:,:,IND),Scale,EEG.Fs,Mapping,LW,Color,'none',LineStype);
        hold on
        set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
    end
    annotation('textbox', [0.3,0.01,0.5,.05],'units','normalized',...
        'String',['NT(' num2str(sum(EEG.EpochClass==0)) ') TA(' num2str(sum(EEG.EpochClass==1)) ') ; EEG scaling: ' ...
        num2str(Scale*2) 'µV' ' Blue (TA), Red (NT)'],'fontsize',FontSize,'fontname','times new roman');
else
    
    switch Type
        case 1
            %% Show Target/Non-Target + difference TA-NT
            
            offset=0;
            Scale=FindScaling(ACSTPfilter.EA(:,:,2));
            %close all
            ClassName={'NT','TA'};
            nbLINE=9;
            nbCOLUMN=8;
            FontSize=12;
            set(gcf, 'PaperPosition', [0 0 25*nbCOLUMN/5 60*nbLINE/9],'units','normalized','outerposition',[0.1 0.1 0.8 .9])
            
            for IND=1:4
                
                for classIND=1:2
                    
                    
                    
                    
                    t = (0:(ACSTPoptions.Epoch_size-1))./EEG.Fs;
                    
                    
                    switch IND
                        case 1
                            %PLOTX=mean(X(:,:,EpochClass==1),3);
                            PLOTX=ACSTPfilter.EA(:,:,classIND);
                            info=['AEA'];
                            if ACSTPfilter.Class(classIND)==0;
                                info=[info ''];
                            else
                                info=[info ''];
                            end
                        case 2
                            %PLOTX=Xcstp(:,:,EpochClass==1);
                            PLOTX=ACSTPfilter.EAcstp(:,:,classIND);
                            info=['ACSTP(' num2str(ACSTPfilter.BestPz(classIND)) ')'];
                        case 3
                            PLOTX=ACSTPfilter.EA(:,:,2)-ACSTPfilter.EA(:,:,1);
                            info=['AEA TA-NT'];
                        case 4
                            PLOTX=ACSTPfilter.EAcstp(:,:,2)-ACSTPfilter.EAcstp(:,:,1);
                            info=['ACSTP TA-NT'];
                            
                            
                    end
                    xticklabels = (offset/EEG.Fs):.25:(ACSTPoptions.Epoch_size+offset)/EEG.Fs;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT EEG %%%%%%%%%%%%%%%
                    TemporalCells=(1:nbCOLUMN:nbCOLUMN*(nbLINE-2));TemporalCells=sort([TemporalCells,TemporalCells+1]);
                    if classIND==1
                        LineStype='-';
                        LW=1;
                        Color=[0.9,0.1,0.1];
                    else
                        
                        LineStype='-';
                        LW=2;
                        Color=[0.1,0.1,0.9];
                        
                        
                        
                    end
                    if IND>2
                        Color=[0.1,0.1,0.1];
                    end
                    hold on;
                    subplot(nbLINE,nbCOLUMN,TemporalCells+2*(IND-1));h(:,IND)=plotEEG(PLOTX,Scale,EEG.Fs,Mapping,LW,Color,'none',LineStype); % Xbar TA
                    hold off;
                    
                    if(IND~=1),set(gca,'yticklabel',[]);end
                    
                    if mod(IND-1,2)>0
                        set(gca, 'color', [0.95 .95 .95])
                    end
                    set(gcf, 'color', [1 1 1])
                    %                 xlim([0 1-offset/EEG.Fs])
                    
                    title(info,'FontSize',FontSize)
                    %xlabel('Time (s)')
                    xticks = linspace(0, (1.01)*(ACSTPoptions.Epoch_size)/EEG.Fs, numel(xticklabels));
                    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
                    set(gcf, 'InvertHardCopy', 'off');
                    % global field power on average
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%% SUBPLOT GFP mean %%%%%%%%%%%%%%%
                    if IND>2
                        TemporalCells=(nbCOLUMN*(nbLINE-2)+1:nbCOLUMN:nbCOLUMN*(nbLINE-1+1));
                    else
                        TemporalCells=(nbCOLUMN*(nbLINE-classIND)+1:nbCOLUMN:nbCOLUMN*(nbLINE-classIND+1));
                    end
                    TemporalCells=sort([TemporalCells,TemporalCells+1]);
                    hold all;
                    subplot(nbLINE,nbCOLUMN,TemporalCells+2*(IND-1))
                    area(t,global_field_power(PLOTX'),'LineWidth',2,'FaceColor',[.5 .5 .5]); %mean GFP Xbar TA
                    hold off;
                    
                    if IND==1
                        %ylabel('GFP','FontSize',FontSize,'fontname','times new roman','FontAngle','italic')
                        SCALEgfp=Scale*0.99;%max(global_field_power(PLOTX'));
                    end
                    %                 xlim([0 1-offset/EEG.Fs])
                    ylim([0 SCALEgfp]);grid on;
                    if IND==1
                        ylabel(ClassName(classIND),'FontSize',FontSize,'fontname','times new roman','FontAngle','italic');
                    end
                    %set(gca,'YtickLabel',[])
                    %xticks = linspace(0, (1.01)*(ACSTPoptions.Epoch_size)/EEG.Fs, numel(xticklabels));
                    set(gca, 'XTick', xticks, 'XTickLabel', xticklabels,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
                    set(gca, 'color', [0.95 .95 .95])
                    
                    
                end
                annotation('textbox', [0.3,0.01,0.5,.05],'units','normalized',...
                    'String',['NT(' num2str(sum(EEG.EpochClass==0)) ') TA(' num2str(sum(EEG.EpochClass==1)) ') ; Top: EEG scaling: ' ...
                    num2str(Scale*2) 'µV' ' ; Bottom: GFP (µV^2) ; Blue (TA), Red (NT), Black (TA-NT)'],'fontsize',FontSize,'fontname','times new roman');
                
                
            end
            
            
            
            
        case 2
            set(gcf, 'PaperPosition', [0 0 25 40],'units','normalized','outerposition',[0.6 0.1 0.375 .9])
            FontSize=24;
            subplot(221)
            hist(ACSTPfilter.Latency(EEG.EpochClass==0),(-ACSTPoptions.LatencyCorr_max:+ACSTPoptions.LatencyCorr_max))
            
            xlabel('Latencies NT','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
            subplot(222)
            hist(ACSTPfilter.Latency(EEG.EpochClass==1),(-ACSTPoptions.LatencyCorr_max:+ACSTPoptions.LatencyCorr_max))
            xlabel('Latencies TA','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
            subplot(2,2,[3 4])
            hist(ACSTPfilter.Weights)
            xlabel('Weights','fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
            set(gcf, 'color', [1 1 1])
            
        case 3
            if ~exist('readlocs')
                error('ERROR ACSTPshow, READLOCS.m not found, please verify that eeglab is installed')
            else
                FontSize=24;
                
                for indC=1:length(ACSTPfilter.Class)
                    figure
                    topoplot_components(ACSTPfilter.fullAt(indC), ACSTPfilter.fullAs(indC),ACSTPfilter.Class(indC),EEG.Fs,4,5,EEG.ElectrodesName,[])
                    annotation('textbox', [0.3,0.01,0.5,.05],'units','normalized',...
                        'String',['Class' num2str(ACSTPfilter.Class(indC))],'fontsize',FontSize,'fontname','times new roman');
                end
            end
        otherwise
            error('ACSTPshow: Type incorrect')
    end
end




%**************** NESTED ***********************
    function [gfp]=global_field_power(P)
        for k=1:size(P,3);
            u=P(:,:,k)';
            
            gfp(:,k)=sqrt(mean(u.^2,1));
        end
    end

    function h=plotEEG(sig,Scale,fs,LABELS,FontSize,Color,Marker,LineStyle)
        %plotEEG(sig,Scale,fs,LABELS)
        % *** History: 19-Mar-2015
        % *** Author: Louis KORCZOWSKI (c), GIPSA-Lab, 2015
        % *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
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
        
        ts = (0:(NbSamples-1))./fs;
        %sig = rand(NbDim,NbSamples);
        
        % calculate shift
        
        shift = cumsum([0; abs(ma(1:end-1))+abs(mi(2:end))]);
        shift = repmat(shift,1,NbSamples);
        
        %plot 'eeg' data
        h=plot(ts,sig-shift,'LineWidth',2,'Color',Color,'Marker',Marker,'LineStyle',LineStyle);
        
        
        % edit axes
        set(gca,'ytick',sort(mean(sig-shift,2)),'yticklabel',LABELS,'fontsize',FontSize,'fontname','times new roman','FontAngle','italic')
        %xlhand = get(gca,'ylabel')
        %set(xlhand,'fontsize',20)
        
        grid on
        ylim([ min(min(-shift-Scale*2)) max(max(-shift+Scale*2))])
        %drawline
    end
    function topoplot_components(TemporalComponents, SpatialComponents,Class,Fs,nb_rows,nb_columns,localization,savePath)
        % TemporalComponents, TemporalComponents, cells [nbClasses x 1] are the
        % ACSTP filter (Bt and Bs respectively). The components are in column.
        
        % localization is either the localization structure from eeglab (see
        % readlocs.m) or the electrodes names (cells).
        %
        % *** History: 2015-04-19
        % *** Author: Louis KORCZOWSKI, GIPSA-Lab, 2015
        % *** Related work: M. CONGEDO, L. KORCZOWSKI, A. DELORME, F. LOPES DA SILVA "Spatio-Temporal Common Pattern a Reference Companion Method for ERP Analysis" (submitted)
        
        if nargin<8
            savePath=[]; %do not save the figure
        end
        
        if nargin<7
            error('Please consider to give the localization of electrodes or at least the electrodes names')
        elseif iscell(localization) %localization is a cell of electrodes names i.g. {'Fz','Cz','Pz'}
            [~, tmpo]=findElectrodes(localization, localization);
            clear localization
            localization=tmpo;
        end
        
        %% TOPOPLOT FOR THE COMPONENTS
        CompLABELS=Generate_Components_Numbers(1:size(SpatialComponents{1},2))';
        for classINDc=1:length(Class)
            FontSizec=24;
            [As]=(SpatialComponents{classINDc});
            As=normPOW(As');
            [At]=(TemporalComponents{classINDc});
            At=(normEEG(At(:,1:size(As,1))'));
            
            TemporalCellsc=(1:nb_columns:nb_columns*nb_rows);
            subplot(nb_rows,nb_columns,TemporalCellsc);plotEEG(At,4,Fs,CompLABELS);xlabel('Time (s)')
            for i=1:size(As,1)
                SpatialCells=(1:nb_columns*nb_rows);
                SpatialCells=SpatialCells(~ismember(SpatialCells,TemporalCellsc));
                subplot(nb_rows,nb_columns,[SpatialCells(i)]);
                topoplot(abs(As(i,:)),localization,'plotrad',0.55);title(CompLABELS{i},'FontSize',FontSizec);
                set(gcf, 'PaperPosition', [0 0 40 60]);
                caxis([0 1]);
            end
            hb=colorbar;
            set(hb,'Units','normalized', 'position', [0.92 0.3 0.02 0.6],'FontSize',FontSizec);
            colormap(gray);
            colormap(flipud(colormap));
            
            
            
            
            %spaceplots
            if nargin>6
                if ~isempty(savePath)
                    if ~exist(savePath)
                        mkdir(savePath);
                    end
                    saveas(htopo,[savePath 'class' num2str(classINDc) '.tiff'])
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NESTED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function AUsers=Generate_Components_Numbers(Numbers)
            for nb=Numbers
                if nb<10
                    AUsers{nb}=['c0' num2str(nb)];
                else
                    AUsers{nb}=['c' num2str(nb)];
                end
            end
            AUsers=AUsers(find(~cellfun(@isempty,AUsers)));
        end
        
        function Xnorm=normPOW(X)
            Xnorm=abs(X);
            Xnorm=Xnorm-repmat((min(Xnorm,[],2)), [1,size(Xnorm,2)]);
            Xnorm=Xnorm./repmat((max(Xnorm,[],2)), [1,size(Xnorm,2)]);
        end
        
        function [IndicesElectrodes LocsElectrodes]=findElectrodes(ElectrodesName, SelectedElectrodes)
            % From a file of electrodes in a randomized order ElectrodesName, find the
            % indices of the wanted electrodes SelectedElectrodes
            %
            % exemple
            % ElectrodesName={'Fp1','Fp2','Cz','O1','O2'}
            % SelectedElectrodes= {'Cz','O1'}
            if exist('readlocs')
                alllocs=readlocs('Standard-10-20-Cap81.locs');
            else
                error('ERROR FINDELECTRODES, READLOCS.m not found, please check your eeglab install')
            end
            %LocsElectrodes=struct;
            for ind=1:length(SelectedElectrodes)
                tmp(ind,:)=strcmpi(SelectedElectrodes{ind},ElectrodesName);
                if exist('alllocs')
                    tmp2=alllocs(strcmpi({alllocs.labels},SelectedElectrodes{ind}));
                    if ~isempty(tmp2)
                        LocsElectrodes(ind)=alllocs(strcmpi({alllocs.labels},SelectedElectrodes{ind})); %need eeglab installed
                        
                    end
                end
            end
            IndicesElectrodes=find(sum(tmp,1));
        end
        
        
    end

    function Scale=FindScaling(PLOTX,winElec,winTime)
        if nargin<2
            winElec=1:size(PLOTX,1);
        end
        if nargin<3
            winTime=1:size(PLOTX,2);
        end
        Scale=max(max(abs(mean(PLOTX(winElec,winTime,:),3))))/2;
    end

    function Xnorm=normEEG(X,Method,Param)
        if nargin<2
            Method='';
        end
        if strcmp(Method,'fro')
            
            for k=1:size(X,3)
                coef(k)=norm(X(:,:,k),'fro')/length(X(:,:,k));
                Xnorm(:,:,k)=X(:,:,k)/coef(k);
            end
            
        elseif strcmp(Method,'baseline')
            if nargin<2 || isempty(Param)
                winBaseline=96:128; % HARDCODED NOT GOOD !!!!!!!!!!!!! ! ! !
            else
                winBaseline=Param;
            end
            for k=1:size(X,3)
                coef(:,:,k)=repmat(mean(X(:,winBaseline,k),2),1,size(X,2)); %compute the baseline
            end
            Xnorm=X-coef;
            
        else
            Xnorm=X./repmat(sqrt(var(X,[],2)), [1,size(X,2)]);
            
            for k=1:size(Xnorm,3)
                Xnorm(:,:,k)=Xnorm(:,:,k)-repmat(mean(Xnorm(:,:,k),2),[1 size(Xnorm,2)]);
            end
            
        end
    end

% END OF THE FUNCTION
end