function [adj,rslts]= mjdemetra(data,varargin)

%% Performs seasonal adjustment on the data according to the options chosen
%
% Inputs: 
%
% [CLASS                                            ][ARGUMENT NAME]
%  -------------------------------------------------  -------------
% [TsData (JD+)                                     ] data
% [string (Matlab)  -maybe this is better           ] option
% [integer (Matlab)                                 ] forecastHorizon
%
% https://blogs.mathworks.com/loren/2009/05/05/nice-way-to-set-function-defaults/
%
% Outputs:
%
% [CLASS                                            ][ARGUMENT NAME]
%  -------------------------------------------------  -------------
% [timeseries (JDemetra+ object) structure          ] rslts
% [timeseries (Matlab]                              ] adj
%
% Note 1: Make sure Matlab uses the appropiate Java version
%        (type version -java to find out which version is used)
%         If Matlab is using a version of Java that is not compatible with
%         JDemetra+, then 
% Note 2: Edit the classpath.txt file (type which classpath.txt to find its location)
%         to make sure the paths containing your .jar libraries are listed. For example, 
%         my classpath.txt file includes the following path, containing the java compiled sofware:
%         L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB\demetra-tstoolkit-2.2.2.jar
%         If you don't want to modify the classpath.txt file, then just
%         uncomment this line:
%         javaclasspath('L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB\')
% 
%                                        Param1          Param2                       Param3
% Examples:                            __________   __________________      ______________________
%                                     /          \ /                  \    /                      \
%        [sa, rslts]= mjdemetra2(data,'horizon',20,'Method','TramoSeats','CalendarOption','RSAfull')
%        [sa, rslts]= mjdemetra2(data2,            'Method','X13'      );
%        [sa, rslts]= mjdemetra2(data,'horizon',20,'Method','TramoSeats','CalendarOption','RSA5')
%        [sa, rslts]= mjdemetra2(data,'horizon',20,'Method','X13'       ,'CalendarOption','RSA5c')
%        [sa, rslts]= mjdemetra2(data)
%        [sa, rslts]= mjdemetra2(data,                                  ,'CalendarOption','RSA0')
%        [sa, rslts]= mjdemetra2(data,                                                          , 'grafico',false)
%
%
% By default, the method 'TramoSeats' is used (with specification'RSAfull', unless otherwise stated)
% *Calendar options for TramoSeats:
%  RSAfull (default):
%  RSA0             :
%  RSA1             :
%  RSA2             :
%  RSA3             :
%  RSA4             :
%  RSA5             :
% *Calendar options for X13:
%  RSAX11 (default) :
%  RSA0             :
%  RSA1             :
%  RSA2c            :
%  RSA3             :
%  RSA4c            :
%  RSA5c            :
%
%
% Code written by David de Antonio Liedo
% Feel free to modify this function or create alternative outputs
% Dont' hesitate to contact me if you need assistance 
% Contact email: david.deantonioliedo@nbb.be
% Github       : https://github.com/Liedo


p = inputParser;


% Data matrix is a required input
addRequired(p,'data',@ismatrix);

% Seasonal adjusment method is option (TrameSeats is the default)
defaultSaMethod = 'TramoSeats';
validSaMethods = {'TramoSeats','X13'};
checkSaMethod = @(x) any(validatestring(x,validSaMethods));
addParameter(p,'Method',defaultSaMethod,checkSaMethod)

% Seasonal adjusment options 
defaultCalendarOption = 'RSAfull';
validCalendarOption = {'RSAfull','RSA0','RSA1','RSA2','RSA3','RSA4','RSA5','RSA4c','RSA5c','RSA2c','RSAX11'};
checkCalendarOption = @(x) any(validatestring(x,validCalendarOption));
addParameter(p,'CalendarOption',defaultCalendarOption,checkCalendarOption)


defaultHorizon=12 ;% 12 months, 12 years, 12 quarters...
addParameter(p,'horizon',defaultHorizon,@isnumeric);


defaultPlot=true ;% 12 months, 12 years, 12 quarters...
addParameter(p,'plot',defaultPlot,@islogical);

%defaultSaMethod='TramoSeats' % 12 months, 12 years, 12 quarters...
%addParameter(p,'Method',defaultSaMethod,@ischar);

%defaultCalendarOption='RSAfull' % 12 months, 12 years, 12 quarters...
%addParameter(p,'CalendarOption',defaultCalendarOption,@ischar);


p.KeepUnmatched = true;
 

%Pass the values of all of the function inputs.
parse(p,data,varargin{:});

horizon=p.Results.horizon; % parameter
saMethod=p.Results.Method; % optional
saOption=p.Results.CalendarOption; % optional 
grafico = p.Results.plot 

disp('...........................................................................')
disp(['Seasonal Adjustment Method: ',p.Results.Method])
disp(['Options for calendar effects and outlier-detection: ',p.Results.CalendarOption])
disp('...........................................................................')

if ~isempty(fieldnames(p.Unmatched))
   disp('Extra inputs:')
   disp(p.Unmatched)
end
if ~isempty(p.UsingDefaults)
   disp('Using defaults: ')
   disp(p.UsingDefaults)
end
 

if strcmp(saOption,'RSAX11') | strcmp(saOption,'RSA2c')| strcmp(saOption,'RSA4c')| strcmp(saOption,'RSA5c')
        saOptionJD = ec.satoolkit.x13.X13Specification.fromString(saOption);
% not available for x13???
%        saOptionJD.getSeatsSpecification().setPredictionLength(horizon)

        rslts = ec.satoolkit.algorithm.implementation. ...
            X13ProcessingFactory.process(data, saOptionJD);
        if  strcmp(saMethod,'TramoSeasts')  
            disp('The SA options used correspond to the X13 method, so your call for the TramoSeats method is being ignored')
        end

elseif  strcmp(saMethod,'X13')
        saOptionJD = ec.satoolkit.x13.X13Specification.fromString(saOption);
% not available for x13???
%        saOptionJD.getSeatsSpecification().setPredictionLength(horizon)
        rslts = ec.satoolkit.algorithm.implementation. ...
            X13ProcessingFactory.process(data, saOptionJD);
else %saMethod=='TramoSeats'
        %saOptionJD = ec.satoolkit.tramoseats.TramoSeatsSpecification.fromString('RSAfull')
        saOptionJD = ec.satoolkit.tramoseats.TramoSeatsSpecification.fromString(saOption);
        saOptionJD.getSeatsSpecification().setPredictionLength(horizon);
        rslts = ec.satoolkit.algorithm.implementation. ...
            TramoSeatsProcessingFactory.process(data, saOptionJD);
end


%%  Outliers Identification
%   This java object has the following properties: description,
%   coefficient, stdError, pValue
frecuencia = data.getFrequency().intValue() 
  
multiplicative=rslts.getData('log',java.lang.Boolean(1).getClass())
 


temp = java.lang.Integer(1); % generate a Java integer (output of the function next line)
nout=  rslts.getData('preprocessing.regression.nout', temp.getClass()) % number of outliers

temp=ec.tstoolkit.information.RegressionItem('descripcion',0,0,0)  % generate Regression item (e.g. outlier)
for i=1:nout % extract all outliers
eval(['outlier{i}=char(rslts.getData(''preprocessing.regression.out(',num2str(i),')'', temp.getClass()))'])
description{i}=outlier{i}(1:2)
idx1=find(outlier{i}=='(')
idx2=find(outlier{i}==')')
dateString{i}=  outlier{i}(idx1+1:idx2-1) 
    if frecuencia==4
        yearO   =outlier{i}(idx2-4:idx2-1) 
        quarterO=outlier{i}(idx1+1:idx2-5)    
        if  strcmp(quarterO,'I-')
        quarter='-Q1'
        elseif strcmp(quarterO,'II-')
        quarter='-Q2'
        elseif strcmp(quarterO,'III-')
        quarter='-Q3'
        elseif strcmp(quarterO,'IV-')
        quarter='-Q4'    
        end
        dateOutlier(i)=datenum([year,quarter],'YYYY-QQ')
    else
       % dateOutlierStr{i}=outlier{i}(idx1+1:idx2-1);
        dateOutlier(i)=datenum(dateString{i},'mm-YYYY');
    end
end

 

% Forecasts and estimation uncertainty
if  strcmp(saMethod,'TramoSeats')

    % Adjusted data
    adjusted = rslts.getData('sa', data.getClass());
    T=adjusted.getLength();


    % SA data F
    saTsF = rslts.getData('decomposition.sa_lin_f', data.getClass());
    saTsF_se = rslts.getData('decomposition.sa_lin_ef', data.getClass());
    saF    = nan(T+horizon,1);     %  initialization with NAN
    saF_se = nan(T+horizon,1);     %  initialization with NAN
    
    if multiplicative
         
        for i=T:(T+horizon-1)
        saF(i+1,1)= saTsF.get(i-T);
        temp      = saTsF_se.get(i-T);
        saF_se(i+1,1)=ec.tstoolkit.modelling.arima.LogForecasts.expStdev(temp,saF(i+1,1));
        end
    
    else
        for i=T:(T+horizon-1)
        saF(i+1,1)= saTsF.get(i-T);
       % temp      = saTsF_se.get(i-T);
       % saF_se(i+1,1)=ec.tstoolkit.modelling.arima.LogForecasts.expStdev(temp,saF(i+1,1));
        saF_se(i+1,1)=saTsF_se.get(i-T);
        end
        
    end

    % SA data
    saTs = rslts.getData('decomposition.sa_lin', data.getClass());
    saTs_se = rslts.getData('decomposition.sa_lin_e', data.getClass());    
    ycalTs = rslts.getData('ycal', data.getClass());
 

    ycal   = nan(T+horizon,1); 
    adj   = nan(T+horizon,1);     %  initialization with NAN
    adj2_F   = nan(T+horizon,1);     %  initialization with NAN
    sa    = nan(T+horizon,1);     %  initialization with NAN
    sa_se = nan(T+horizon,1);     %  initialization with NAN
%    saF_se = nan(T+horizon,1);     %  initialization with NAN

    if multiplicative
        
        for i=0:T-1
        adj(i+1,1)  =adjusted.get(i);
        sa(i+1,1)   =saTs.get(i);
        temp        =saTs_se.get(i);
        sa_se(i+1,1)=ec.tstoolkit.modelling.arima.LogForecasts.expStdev(temp, sa(i+1,1));
        ycal(i+1,1) = (ycalTs.get(i))
        end
    
        adj2  = exp(sa);
        adj2U = adj2 + 2*sa_se;
        adj2L = adj2 - 2*sa_se;

        adj2_F  = exp(saF);
        adj2U_F   = adj2_F + 2*saF_se;
        adj2L_F   = adj2_F - 2*saF_se;

    else
        for i=0:T-1
        adj(i+1,1)  =adjusted.get(i);
        sa(i+1,1)   =saTs.get(i);
        %temp        =saTs_se.get(i);
        %sa_se(i+1,1)=ec.tstoolkit.modelling.arima.LogForecasts.expStdev(temp, sa(i+1,1));
        sa_se(i+1,1)=saTs_se.get(i);
        ycal(i+1,1) =ycalTs.get(i)
        end
    
        %adj2  = exp(sa);
        adj2  = sa;
        adj2U = adj2 + 2*sa_se;
        adj2L = adj2 - 2*sa_se;

        adj2_F  = exp(saF);
        adj2U_F   = adj2_F + 2*saF_se;
        adj2L_F   = adj2_F - 2*saF_se;

    
    
    end
     
    
    % raw data
    rawData       =nan(T+horizon,1);
    for i=0:T-1
    rawData(i+1,1)   =data.get(i);
    end

else % x13
    
    % Adjusted data
    adjusted = rslts.getData('sa', data.getClass());
    ycalTs = rslts.getData('ycal', data.getClass());
    T=adjusted.getLength();

    % SA data
  %->  saTs = rslts.getData('decomposition.sa_lin', data.getClass())

    adj   = nan(T,1);     %  initialization with NAN
    ycal   = nan(T,1);     %  initialization with NAN
    
%->    sa    = nan(T+horizon,1);     %  initialization with NAN

    for i=0:T-1
    adj(i+1,1)   =adjusted.get(i);
    ycal(i+1,1)   =ycalTs.get(i);
%->    sa(i+1,1)    =saTs.get(i);
    end

%->    adj2  = exp(sa);
    
    % convert it into time series
    rawData       =nan(T,1);
    for i=0:T-1
    rawData(i+1,1)   =data.get(i);
    end

    
end






if grafico

   for j=1:nout
      oType{j}=['   ',description{j}] 
   end    
    
%sa = rslts.getData('decomposition.sa_lin', data.getClass())
%saF = rslts.getData('decomposition.sa_lin_f', data.getClass())
%saE = rslts.getData('decomposition.sa_lin_e', data.getClass())
%saEf = rslts.getData('decomposition.sa_lin_ef', data.getClass());
if  strcmp(saMethod,'TramoSeats')

    %% Create timeseries

    dominio=data.getDomain(); 
    dominioF= saTsF.getDomain(); 

    T=data.getLength();
    H= saTsF.getLength();
    for i=0:(T-1)
        fechas{i+1}=char(dominio.get(i).toString())
    end
    for i=T:(T+H-1)
        fechas{i+1}=char(dominioF.get(i-T).toString())
    end
   % tiempo0=datestr(datenum(fechas,'mm-yyyy'))
    
    if frecuencia==4
    for i=0:T+H-1
    yearF=fechas{i+1}(end-3:end)
    quarterF=fechas{i+1}(1:end-4)
    if  strcmp(quarterF,'I-')
    quarter='-Q1'
    elseif strcmp(quarterF,'II-')
    quarter='-Q2'
    elseif strcmp(quarterF,'III-')
    quarter='-Q3'
    elseif strcmp(quarterF,'IV-')
    quarter='-Q4'    
    end

    fechasQ{i+1}=[yearF,quarter]
    end
    
    timeFormat='YYYY-QQ';
    tiempo00=datenum(fechasQ,timeFormat)    
   
    else
    timeFormat='mm-yyyy'
    tiempo00= datenum(fechas,timeFormat)

    end
    
    % SA data and events (outliers)
    ts_adj=timeseries(adj,tiempo00)
    %[~,index] = max(ts_adj.Data(:,1));
    %[~,indexT] = min(ts_adj.Data(:,1));
    % eventTime= tiempo00(index,:)
    % eventTimeT= tiempo00(indexT,:)
     
    
    for j=1:nout
        index = find(dateOutlier(j)==ts_adj.Time)
        new_event = tsdata.event(description{j},ts_adj.Time(index))
        new_event.Units = 'seconds'
        ts_adj = addevent(ts_adj,new_event); 
        indice(j)=index;
    end
    

    saOutputSE=[ adj2U  adj2L  adj2U_F  adj2L_F];
    saOutput=[rawData adj  adj2 ycal];
    ts=timeseries(saOutput,tiempo00)
    ts.TimeInfo.Format=timeFormat
    tsCI=timeseries(saOutputSE,tiempo00)
    tsCI.TimeInfo.Format=timeFormat
    %figure,plot(ts)

    % event https://nl.mathworks.com/help/matlab/ref/timeseries.plot.html
    
 
     
   %%
   

   
   figure
   % subplot(2,1,1)
  %  plot([rawData adj  adj2 ]),legend('raw data', 'sa','exp(decomposition.sa-lin)') % add forecasts and interval for both forecasts and past
  %  hold on
  %  plot([adj2U  adj2L  adj2_F adj2U_F  adj2L_F],'k:')%,legend('raw data', 'sa','exp(decomposition.sa-lin)') % add forecasts and interval for both forecasts and past
  %  subplot(2,1,2)
   % plot(ts),legend('raw data', 'sa','exp(decomposition.sa-lin)')
   % grid on
   % hold on
   % plot(tsCI,':')
    
   % subplot(2,1,2)
   
    plot(ts_adj,'.-b'),  datetick('x', timeFormat), xlabel('time'),title(['JD+ adjustment with TRAMO-SEATS-',p.Results.CalendarOption]);
    hold on
    % too complicated: use tiempo00 instead
    plot(tiempo00,ts.Data(:,1),'k',... % datenum(fechas,'mm-yyyy'),ts.Data(:,2),'b.-',...        
         tiempo00,ts.Data(:,3),'c',...
         tiempo00,adj+rawData-ycal ,'r',... 
         tiempo00,tsCI.Data(:,1),'c:',...
         tiempo00,tsCI.Data(:,2),'c:',...
         tiempo00,tsCI.Data(:,3),'c:',...
         tiempo00,tsCI.Data(:,4),'c:',...
         tiempo00,mean(tsCI.Data(:,3:4),2),'c')
         legend('Tramo-Seats SA  & Cal','raw','exp(decomposition.sa-lin)','Tramo-Seats SA')
    datetick('x', timeFormat);
    hold on
     if nout>0       
      text( dateOutlier  ,ts_adj.Data(indice,1),oType,'Color','red')     
     end

    grid minor
%    dateOutlier
    %         tiempo00(3),ts.Data(3,1),'r-s',...%,'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6],...
    %         tiempo00(8),ts.Data(8,1),'r-s')%,'MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6]),...

else
    %% Create timeseries

    dominio=data.getDomain(); 
    %->dominioF= saTsF.getDomain(); 

    T=data.getLength();
    %->H= saTsF.getLength();
    for i=0:(T-1)
        fechas{i+1}=char(dominio.get(i).toString())
    end
    %->for i=T:(T+H-1)
    %->    fechas{i+1}=char(dominioF.get(i).toString())
    %->end
   
    % tiempo0=datestr(datenum(fechas,'mm-yyyy'))
    
    if frecuencia==4
        for i=0:T-1
        yearF=fechas{i+1}(end-3:end)
        quarterF=fechas{i+1}(1:end-4)
        if  strcmp(quarterF,'I-')
        quarter='-Q1'
        elseif strcmp(quarterF,'II-')
        quarter='-Q2'
        elseif strcmp(quarterF,'III-')
        quarter='-Q3'
        elseif strcmp(quarterF,'IV-')
        quarter='-Q4'    
        end

        fechasQ{i+1}=[yearF,quarter]
        end
        
        timeFormat='YYYY-QQ'
        tiempo00=datenum(fechasQ,timeFormat)
   
    else
        timeFormat='mm-yyyy'
        tiempo00= (datenum(fechas,timeFormat))

    end


    % SA data and events
    ts_adj=timeseries(adj,tiempo00)
 
    for j=1:nout
        index = find(dateOutlier(j)==ts_adj.Time)
        new_event = tsdata.event(description{j},ts_adj.Time(index))
        new_event.Units = 'seconds'
        ts_adj = addevent(ts_adj,new_event); 
        indice(j)=index;
    end
    
    
    
    saOutput=[rawData adj ycal];
    ts=timeseries(saOutput,tiempo00)
    ts.TimeInfo.Format=timeFormat
    %figure,plot(ts)

    %%
     
     figure
     plot(ts_adj,'.-b'),  datetick('x', timeFormat), xlabel('time'),title(['JD+ adjustment with X13-',p.Results.CalendarOption]);
     hold on
     plot( tiempo00 , rawData   ,'k') , datetick('x', timeFormat)
     hold on
     plot( tiempo00 , adj+rawData-ycal   ,'r') , datetick('x', timeFormat)
     if nout>0
         hold on %     
         text( dateOutlier  ,ts_adj.Data(indice,1),oType,'Color','red')     
     end
     legend('X13 SA  & Cal','raw','X13 SA')
     datetick('x', timeFormat);
     %plot( tiempo00(marcador) ,ts_adj.Data(marcador,1),'m-s','MarkerSize',10) 
     %hold on
end

end

%toc
  
%display(toc)