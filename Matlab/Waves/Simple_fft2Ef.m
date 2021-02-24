
function [f_spec,Ef] = Simple_fft2Ef(t,y)


%% Simple FFT Gives back the AMPLITUDE spectrum!! 
% for tutorial look below!

% Sampling frequency:
dt=t(end)-t(end-1);
Fs = 1/dt; % Hz			
nt=length(t);

f = Fs*(0:nt-1)/nt;	% Frequency vector 
f_spec=f(1:round(nt/2)); % pick only first half of the freq vector

Ffty=fft(y);
A_fft=abs(Ffty); % amplitudes
A_spec=1/nt*2*A_fft(1:round(nt/2)); % pick only first half of the freq vector

df=f_spec(2)-f_spec(1);
Ef=(1/df)*0.5*(abs(A_spec)).^2; % Variance density spectrum!

%% Figure 2 check
% figure
% plot(f_spec,A_spec)
% 



%% Tutorial

% %%%%%%%%%%%%%%%%%%_____________________%%%%%%%%%%%%%%%%%%%
% %% Simple example of FFT use
% 
% % Sampling frequency:
% dt=0.001;
% Fs = 1/dt; % Hz			
% 
% % Time vector;
% t = dt:dt:5;
% nt=length(t)
% 
% % Signal
% y=sin(2*pi*t*50) + 5*sin(2*pi*t*60);      
% 
% 	
% % Frequency vector
% % NFFT = length(y)
% f = Fs*(0:nt-1)/nt;	
% 
% A_fft=abs(fft(y));
% 
% %% 
% figure
% plot(f,A_fft)
% 
% f_spec=f(1:round(nt/2));
% A_spec=2/nt*A_fft(1:round(nt/2));
% 
% %%
% figure
% plot(f_spec,A_spec)