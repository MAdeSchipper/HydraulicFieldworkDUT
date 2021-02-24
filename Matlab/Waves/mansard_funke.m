function [zsi zsr] = mansard_funke(t,zs,xp,h)
% Use three-gauge method to separate incident from reflected wave
%
% [zsi zsr] = mansard_funke(t,zs,xp,h)

% zero-base time
tz=t(:)-t(1);

% length of time vector
% lt = length(t);

% size of zs
sz = size(zs);
if min(sz)~=3
    error('you must supply three time series of zs');
end
if sz(1)==3
    zs=zs';  % convert to row based time
end

% detrend water level signal.. no?
zsd = detrend(zs);
zsrem = zs-zsd;

% number of measurement points
K=length(tz);

% frequency range
nf = floor(K/2);
df = 1/tz(end);
ff = df*[0:1:round(K/2) -1*nf+1:1:-1];
% fx = ff(1:nf);
wx = 2*pi*ff;
kx = disper(wx,h,9.81);
% need to flip this round to make similar to Mansard Funke calculations
kx=fliplr(kx);

% Fourier transform
Q = fft(zsd,[],1);
% B = abs(Q).*exp(1i.*angle(Q));
% B = abs(Q).*cos(-angle(Q)) - 1i.*abs(Q).*sin(-angle(Q));
% B = conj(Q);
B = Q;
% Bp = B(1:nf,:); % only positive side

% cut out hf noise
fp = ff(abs(Q(2:end,1))==max(abs(Q(2:end,1))))+ff(2);
fp = fp(fp>0);
ffreq = [5 10];
B(abs(ff)>ffreq(1)*fp,:) = B(abs(ff)>ffreq(1)*fp,:) - ...
               repmat(...
                       min(1,abs(ff(abs(ff)>ffreq(1)*fp))-ffreq(1)*fp./(ffreq(2)-ffreq(1)))',...
                       1,3) ...
               .*B(abs(ff)>ffreq(1)*fp,:);
% cut out lf noise
ffreq = [1/5 1/10];
B(abs(ff)<ffreq(1)*fp,:) = B(abs(ff)<ffreq(1)*fp,:) - ...
               repmat(...
                       min(1,ffreq(1)*fp-abs(ff(abs(ff)<ffreq(1)*fp))./(ffreq(1)-ffreq(2)))',...
                       1,3) ...
               .*B(abs(ff)<ffreq(1)*fp,:);



% set up output array
zsic = 0*B(:,1);
zsrc = 0*B(:,1);
for j=1:length(B)
    if kx(j)~=0
        % calculate beta and gamma from instrument locations
        beta = kx(j)*(xp(2)-xp(1));
        gamma = kx(j)*(xp(3)-xp(1));
        D = 2*((sin(beta))^2 + (sin(gamma))^2 + (sin(gamma-beta))^2);
        R1 = (sin(beta))^2 + (sin(gamma))^2;
        Q1 = sin(beta)*cos(beta) + sin(gamma)*cos(gamma);
        R2 = sin(gamma)*sin(gamma-beta);
        Q2 = sin(gamma)*cos(gamma-beta)-2*sin(beta);
        R3 = -sin(beta)*sin(gamma-beta);
        Q3 = sin(beta)*cos(gamma-beta)-2*sin(gamma);
        zsic(j) = 1/D*( B(j,1)*(R1+1i*Q1) + B(j,2)*(R2+1i*Q2) + B(j,3)*(R3+1i*Q3));
        zsrc(j) = 1/D*( B(j,1)*(R1-1i*Q1) + B(j,2)*(R2-1i*Q2) + B(j,3)*(R3-1i*Q3));
    end
    
end

zsi = real(ifft(zsic));
zsi = zsi+zsrem(:,1);
zsr = real(ifft(zsrc));
zsr = zsr+zsrem(:,1);