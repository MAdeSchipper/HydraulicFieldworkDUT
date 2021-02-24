%% Example how to cross correlate two signals
% This example shows this for 2 sine waves  
% M.A. de Schipper, 2015

clear all
clc
close all

%% Make two example signals
dt_wave=1; % sampling frequency signal
t=1:dt_wave:500; % time vector

Twave=10; % Wave period waves in example
sig1=sin(6.28/Twave*t)+0.1*(rand(size(t))-0.5); % Signal 1 including random noise. Wave height =1 m. noise is 0.1 m

timelag=.4; % artificial time shift between signal 1 and signal 2
sig2=sin(6.28/Twave*(t+timelag))+0.1*(rand(size(t))-0.5); % Signal 1 including random noise. Wave height =1 m. noise is 0.1 m



%% Show  signal 1 in figure
figure; 
 plot(t,sig1,'x--k')
 xlabel('time')
 ylabel('surface elevation (m)')

 %% Show both signals in a figure
figure; 
 plot(t,sig1,'x--k')
 hold on
  plot(t,sig2,'x--b')
 xlabel('time')
 ylabel('surface elevation  (m)')
xlim([0 100])
 
 
%% Resample data to 10Hz (to be able to see the timeshift well)

dt=.1;         % new sampling rate
t_interp=t(1):dt:t(end); % new time vector
sig1_spline = interp1(t, sig1, t_interp, 'spline'); % Interpolated signal 
sig2_spline = interp1(t, sig2, t_interp, 'spline'); 

%% Input signals into cross correlation
a=sig1_spline-mean(sig1_spline); % first remove the mean of the timeseries
b=sig2_spline-mean(sig2_spline);

maxlag=10;      % Maximum time lag you expect in seconds
subwindow=5;    % number of subseries the time series is devided in

%This program is uses xcorr.m  %% Signal Processing toolbox!
%Jamie MacMahan 3/14/02

subwindow=subwindow+1;
index=round(linspace(1,length(a),subwindow));
maxlag=round((maxlag-1)./dt);
for i=1:length(index)-1;
         AA=detrend(a(index(i):index(i+1)));
         BB=detrend(b(index(i):index(i+1)));
         [C,LAGS]=xcorr(AA,BB,maxlag,'biased');
     if i==1
         C1=zeros(size(C));
     end
     C1=C./sqrt(var(AA).*var(BB))+C1;
	 %C1=C+C1;
end
c=C1./(subwindow-1);
lags=LAGS.*dt;

[temp,ind]=max(c);      %% find max correlation peak
Best_Time_Shift=lags(ind)
        
%% plot crosscorrelation
figure
plot(lags, c)
axlimits=axis;
hold on
plot([Best_Time_Shift Best_Time_Shift], [axlimits(3:4)],'r') % line at maximum correlation
grid on;
xlabel('Time shift in seconds')
ylabel('Correlation value')
title(['maximum correlation for a time shift of ' num2str(Best_Time_Shift) ' sec'])
