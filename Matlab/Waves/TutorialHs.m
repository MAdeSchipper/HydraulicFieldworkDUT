%% Tutorial, different ways to calculate the HS out of a time signal

%% Load Signal

%%%-- synthetic data

dt=1/2; 
nt=1000;
t = dt:dt:nt*dt;
y=2.55*sin(2*pi*t*0.1) + 0.2*sin(2*pi*t*0.2-pi/2) + 0.1*sin(2*pi*t*0.3-pi/2);  



%% Important... Detrend the data to remove a tilt

y=detrend(y); 

%%
figure
plot (t,y)
%% By means of the variance of the time signal (Sierd's method)

disp ('Based on the variance of the signal ( var(signal) )' )
Hs_sierd=4*sqrt(var(y))


%% By means of the variance of the variance density spectrum 
disp ('Methods based on the variance density Spectrum' )

Fs=1/dt;                %Sampling rate

% Matt's vage scripts  
[f_spec,A_spec] = Simple_fft(t,y);
        %%% Simple_fft does the following:
        % A_fft=abs(fft(y));
        % f_spec=f(1:round(nt/2));
        % A_spec=2/nt*A_fft(1:round(nt/2));

df=f_spec(2)-f_spec(1);

M0=sum(0.5*A_spec.^2);
disp ('Based on Matthieus summation method of the spectrum' )
Hm0_matt=4*sqrt(M0)

%% Martijn's Trapeziodal method
E=(1/df)*0.5*(abs(A_spec)).^2; % convert amplitudes to variance(0.5*A^2), and variance to variance density (1/df)
m0=trapz(f_spec,E);
disp ('Based on Martijns summation method of the spectrum' )
Hm0_mart=4*sqrt(m0)



%% Plot signal


figure(1)
plot(t,y)
xlabel('time [s]')
ylabel('surface elevation [m]')
title('signal')
%xlim([0 100])

%% plot Spectrum

figure(2)
title('Amplitude Spectrum')
plot(f_spec,A_spec,'k')
xlabel('freq [Hz]')
ylabel('Amplitude [m]')

hold off
