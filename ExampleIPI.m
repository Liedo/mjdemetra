% This is my working directory; choose your own
cd('L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB');
addpath('L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB');
clear all

% If needed, specify the location of the JDemetra+ compiled software: demetra-tstoolkit-2.2.2.jar
% javaclasspath('L:\DSXNPAPER\Project 2017\R model\models\JDinMATLAB\');

% Load the data
[A,B]=xlsread('IPI','sheet1');
y=A(:,1); % in this example,the data is stored in A  
timeCellArray=char(B(2:end,1));  % the dates are stored in B in a format that Matlab does not recognize
yearString=timeCellArray(1,1:4); % Let's specify the year corresponding to the first observation
monthString=timeCellArray(1,6:7);% Let's specify the month corresponding to the first observation

% Specify a time
firstPeriod = ec.tstoolkit.timeseries.simplets.TsPeriod(ec.tstoolkit.timeseries.simplets. ...
    TsFrequency.Monthly,...   % specify the frequency
    str2num(yearString), ...  % specify the year of the first observation (it's a string, so we convert it to number)
    str2num(monthString)-1);  % specify the month. NOTE that month 1 in Java is 0, so we use: str2num(monthString)-1
data = ec.tstoolkit.timeseries.simplets.TsData(firstPeriod , y, false); % Create the TsData object and name it 'data'

% 'data' is a TsData object. You can print it in Matlab
display(data) % you can print it in your command window but you cannot plot it, e.g. plot(data) will not work

%%  Seasonal Adjustment with JDemetra+
%   This function will take the TsData object as an input and produce the
%   seasonal adjusted data (sa) as an output that can be read in Matlab

[output, rslts]= mjdemetra(data); % It will run TramoSeats (with RSAfull option) and plot results
 
 
[output, rslts]= mjdemetra(data,'Method','X13'); % this will run X13 (with RSAfull option) and plot results
 
 
[output, rslts]= mjdemetra(data,'Method','X13','CalendarOption','RSA2C'); % this will run X13 (with RSA2C option) and plot results
 
tic
[output, rslts]= mjdemetra(data,'plot',false); % this doesn't make the plot, and runs the default method/option
toc
