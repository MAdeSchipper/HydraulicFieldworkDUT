%% Tutorial

%%%%%%%%%%%%%%%%%%_____________________%%%%%%%%%%%%%%%%%%%
%% Simple example of FFT use

% Sampling frequency:
dt=0.001;
Fs = 1/dt; % Hz			

% Time vector;
t = dt:dt:5;
nt=length(t)

% Signal
y=sin(2*pi*t*50) + 5*sin(2*pi*t*60);      


figure
plot(t,y)
	
% Frequency vector
% NFFT = length(y)
f = Fs*(0:nt-1)/nt;	

A_fft=abs(fft(y));

%% 
figure
plot(f,A_fft)

f_spec=f(1:round(nt/2));
A_spec=2/nt*A_fft(1:round(nt/2));

%%
figure
plot(f_spec,A_spec)