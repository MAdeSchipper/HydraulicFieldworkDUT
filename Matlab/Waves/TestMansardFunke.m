%% SCRIPT OM DE WERKING VAN MANSARD & FUNKE SCHEIDING VAN INKOMENDE EN UITGAANDE GOLVEN TE TESTEN.
% In dit script wordt de scheidingsmethode van Mansard & funke getest.

clear 
close all
%% Parameters

% [zsi zsr] = mansard_funke(t,zs,xp,h)
% t = tijd in seconds (length burst)
% zs = tijdseries van laagfrequente signal (voor 3 drukdozen)
% x = afstand tussen de drukdozen. Vector met drie waarde [m] (1 = 0; 2= tov 1; 3 = tov 1) 
% h = waterdiepte [m] (die nemen wij gemiddeld over de drie drukdozen). 

h=2;

t=1:1:1800;
xp=[0 5 10];

f_lw=1/100; % long wave frequency
wx=2*pi*f_lw; % long wave angular 
kx = disper(wx,h,9.81); % long wave wavenumber (free long wave)

%% TEST 1, enkel een inkomende lange golf 
% Test om te kijken hoe goed enkel een inkomende golf wordt herkent.

%maak kunstmatige tijdseries per wavegauge
a=1; % amplitude inkomende golf
signal1=a*cos(kx*xp(1)-wx*t);
signal2=a*cos(kx*xp(2)-wx*t);
signal3=a*cos(kx*xp(3)-wx*t);

%%
% figuur om te testen hoe de kunstmatige tijdseries eruitzien

figure
plot(t,signal1,t,signal2)
xlabel('time [s]')
ylabel('surface elevation [m]')
legend('signal measured at station 1','signal measured at station 2')
xlim([0 200])

%%
% top.. het is dezelfde cosinus, maar enkel wat verschoven

% Nu eens kijken wat er gebeurd als we de inkomende en uitgaande golf
% proberen te vinden met de methode van Mansard&Funke.

%Combineer de signalen in een matrix
zs=[signal1; signal2; signal3];

%%bereken inkomend en uitgaand signaal
[zsi zsr] = mansard_funke(t,zs,xp,h);

% figuur tijdseries
figure
plot(t,zsi,t,zsr,t,signal1'-(zsr+zsi))
legend('zsi (incoming)','zsr (outgoing)','difference original and combined incoming&outgoing')
xlabel('time t [s]')
ylabel('surface elevation [m]')
title ('separated timeseries')

%% 
% Bulk statistieken van het ingevoerde signaal bij punt 1 en het
% de gescheide signalen.
% het inkomende signaal gheeft een amplitude $a$ gelijk aan 1. De variantie is dan:
var_based_on_amp=a^2/2

% op basis van de tijdserie komen we op een variantie van 
var_inputsignal=var(signal1)

% Dat komt mooi overeen.

%op basis van de tijdserie komen we op een variantie van 
var_incomingwaves=var(zsi)
var_outgoingwaves=var(zsr)

% difference beteen input variance and variance after processing is

Difference_variance=var_inputsignal - (var_incomingwaves+var_outgoingwaves)

%%
% We zien dus dat er een kleine toename is van de variantie na het scheiden.


%% TEST 2, een inkomende lange golf en een uitgaande lange golf
% Here we make an artificial timeseries with an incoming long wave with
% amplitude $a$ of 1 and a reflected wave with half the amplitude

a_inc=1;
a_out=0.5;

% this looks sort of like this
figure
plot(t,a_inc*cos(kx*2*xp(1)-wx*t),t,a_out*sin(-kx*xp(1)-wx*t+0.1*pi))

%maak kunstmatige tijdseries 
signal1=a_inc*cos(kx*xp(1)-wx*t)+a_out*sin(-kx*xp(1)-wx*t+0.1*pi);
signal2=a_inc*cos(kx*xp(2)-wx*t)+a_out*sin(-kx*xp(2)-wx*t+0.1*pi);
signal3=a_inc*cos(kx*xp(3)-wx*t)+a_out*sin(-kx*xp(3)-wx*t+0.1*pi);

 zs=[signal1; signal2; signal3];
%%
% figuur om te testen hoe de kunstmatige tijdseries eruitzien

figure
plot(t,signal1,t,signal2)
xlabel('time [s]')
ylabel('surface elevation [m]')
legend('signal measured at station 1','signal measured at station 2')
%xlim([0 800])

%%
% top.. het is dezelfde cosinus, maar enkel wat verschoven

% Nu eens kijken wat er gebeurd als we de inkomende en uitgaande golf
% proberen te vinden met de methode van Mansard&Funke.

%%bereken inkomend en uitgaand signaal
[zsi zsr] = mansard_funke(t,zs,xp,h);


figure
plot(t,zsi,t,zsr,t,signal1'-(zsr+zsi))
legend('zsi (incoming)','zsr (outgoing)','difference original and combined incoming&outgoing')
xlabel('time t [s]')
ylabel('surf elevation [m]')
title ('separated timeseries')

%% 
% Bulk statistieken van het ingevoerde signaal bij punt 1 en het
% de gescheide signalen.
% het inkomende signaal gheeft een amplitude $a$ gelijk aan 1. De variantie is dan:
var_inc_based_on_amp=(a_inc)^2/2
var_out_based_on_amp=(a_out)^2/2

% op basis van de tijdserie komen we op een variantie van 
var_inputsignal=var(signal1)

% Dat komt mooi overeen.

%op basis van de tijdserie komen we op een variantie van 
var_incomingwaves=var(zsi)
var_outgoingwaves=var(zsr)

% difference beteen input variance and variance after processing is

Difference_variance=var_inputsignal - (var_incomingwaves+var_outgoingwaves)

%%
% We zien dus dat er een kleine toename is van de variantie na het scheiden.



%% TEST 3, een inkomende lange golf en een uitgaande lange golf met ruis

%maak kunstmatige tijdseries 


signal1=cos(kx*xp(1)-wx*t)-0.5*cos(-kx/2*xp(1)-wx*t+0.25*pi)+(rand(1,length(t))-1)*0.2;
signal2=cos(kx*xp(2)-wx*t)-0.5*cos(-kx/2*xp(2)-wx*t+0.25*pi)+(rand(1,length(t))-1)*0.2;
signal3=cos(kx*xp(3)-wx*t)-0.5*cos(-kx/2*xp(3)-wx*t+0.25*pi)+(rand(1,length(t))-1)*0.2;

 zs=[signal1; signal2; signal3];

%%bereken inkomend en uitgaand signaal
[zsi zsr] = mansard_funke(t,zs,xp,h);

figure
plot(t,zsi,t,zsr,t,signal1'-(zsr+zsi))
legend('zsi (incoming)','zsr (outgoing)','difference original and combined incoming&outgoing')
xlabel('time t [s]')
ylabel('surf elevation [m]')
title ('separated timeseries')