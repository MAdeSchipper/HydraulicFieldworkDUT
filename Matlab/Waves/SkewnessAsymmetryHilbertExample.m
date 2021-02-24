%% Example for CD2

t=0:0.1:50;


T=8; % wave period
d_phi=pi/2 % pahse shift
a_1=1; % amplitude wave 1;
a_2=0.25; % amplitude wave higher harmonic;

eta=a_1*cos(t*2*pi/T)+a_2*cos(t*2*pi/(T/2)-d_phi);

% Sk volgens Brinkkemper
eta_m=eta-mean(eta);
sk=mean(eta_m.^3)/(mean(eta_m.^2))^1.5
eta_H=imag(hilbert(eta_m));
ass=mean(eta_H.^3)/(mean(eta_H.^2))^1.5

figure
plot(t,eta,'linewidth',2)
grid on 
hold on;
plot(t,eta_H,'linewidth',2)

xlim([15 35])
xlabel ('time [s]')
ylabel ('surf. elevation \eta [m]')
legend ('\eta=cos(2\pi/T x)+cos(2\pi/(T/2) x + \pi/2)','Hilbert(cos(\eta))')
set(gca,'fontsize',14)