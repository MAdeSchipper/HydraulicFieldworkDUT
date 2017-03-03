%%%%Script to translate offshore wave heights and angles to nearshore values
%
% using linear shoaling and refraction, assuming parralel depth contours.
% 
% Wave angles > =- 90 at the buoy are assumed to move offshore and
% calculated wave heights inshore are therefore set to zero.

%%%% INPUT:
%%% offshore wave data:
% d_buoy  % in m ; depth at offshore station
% H0 % in m  ; wave height at offshore station
% Theta0 % Angle of wave incidence in degrees WITH RESPECT TO SHORENORMAL!!! If data is from bouy (in deg N) subtract shoreline orientation first
% T  % in s ; wave period at offshore station;
%%% parameters
% d_inshore % in m ; depth at inshore station.
%  (if d_inshore_in set to nan, then values at breakpoint are calculated)
% gamma_b % -  ; breaking criterium H/h
%  (if gamma_b is nan, then wave heights are not affected by breaking (and too large inshore))

%%%% EXAMPLE:
% 
%  Hrms0=[1 2 3 4]; % wave height
%  theta0=[0 10 5 20]; % wave angle wrt shorenormal
%  T=[10 7 8 9]; % wave period
%  d_buoy=32; % depth at the deep water location in m
%  d_inshore_in=5; % depth inshore
%  gamma_b=0.8; % breaking criterium h/h
% 
%  [h_inshore,d_inshore,theta_inshore]=Waves2Nearshore(d_buoy,Hrms0,theta0,T,d_inshore_in,gamma_b);

%%%%
% M.A. de Schipper 2010-2013


function [H_inshore,d_inshore,Theta_inshore]=Waves2Nearshore(d_buoy,H0,Theta0,T,d_inshore_in,gamma_b)

% deep water values
w=2*pi./T;
k0 = disper(w, d_buoy, 9.81); % k0
L0=2*pi./k0;
C0=L0./T;
n0=0.5*(1+2*k0.*87./sinh(2*k0.*d_buoy));
Cg0=n0.*C0;

%  values inshore
if ~isnan(d_inshore_in) % for finite values of d_inshore
    disp('---Nearshore values calculated ---')
    k = disper(w, d_inshore_in, 9.81);
    L=2*pi./k;
    C=L./T;
    n=0.5*(1+2*k.*d_inshore_in./sinh(2*k.*d_inshore_in));
    Cg=n.*C;
    Theta_inshore=asind(C./C0.*sind(Theta0));
    
    Ks=real(sqrt(Cg0./Cg));                    % shoaling coeff
    Kr=real(sqrt(cosd(Theta0)./cosd(Theta_inshore))); % ref coeff
    
    H_inshore=Kr.*Ks.*H0;                               % wave height inshore
    if ~isnan(gamma_b)
        H_inshore=min(Kr.*Ks.*H0,gamma_b*d_inshore_in);       %if waves are breaking at 5m waterdepth take the gamma h as waveheight
    end
    d_inshore=ones(size(H_inshore))*d_inshore_in;
    
else % for NaN value of d_inshore calculate values at breakpoint (solved iteratively).
    disp('--- values at Breakpoint calculated ---')
    h_table=0:0.05:max(H0)/gamma_b*3;  % range to search for breakerdepth
    
    for i_day =1:length(T)
       
        if abs(Theta0(i_day))<90 
            
        k_i = disper(w(i_day), h_table, 9.81);
        L_i=2*pi./k_i;
        C_i=L_i./T(i_day);
        n_i=0.5*(1+2*k_i.*h_table./sinh(2*k_i.*h_table));
        Cg_i=n_i.*C_i;
        Theta_i=asind(C_i./C0(i_day).*sind(Theta0(i_day)));
        Theta_i(h_table==0)=0; % avoid problems with for depths being equal to 0;
        
        Ks_i=sqrt(Cg0(i_day)./Cg_i);    % shoaling coefficient
        Kr_i=sqrt(cosd(Theta0(i_day))./cosd(Theta_i));% refraction coefficient
        
        Hs_i=Kr_i.*Ks_i.*H0(i_day);
        
        [dum,ind]=min(abs(Hs_i-gamma_b*h_table));  % search for the solution with the smallest error
        if ind==length(h_table)
            warning('max of table reached.. be careful')
        end
        
        
        % store best results
        H_inshore(i_day)=Hs_i(ind);
        d_inshore(i_day)=h_table(ind);
        Theta_inshore(i_day)=Theta_i(ind);
        
        else
         H_inshore(i_day)=NaN;
        d_inshore(i_day)=NaN;
        Theta_inshore(i_day)=NaN;
        end
    end
end

function k = disper(w, h, g)
%DISPER  Linear dispersion relation.
%
%   absolute error in k*h < 5.0e-16 for all k*h
%
%   Syntax:
%   k = disper(w, h, g)
%
%   Input:
%   w = 2*pi/T, were T is wave period
%   h = water depth
%   g = gravity constant
%
%   Output:
%   k = wave number
%
%   Example
%   k = disper(2*pi/5,5,9.81);

%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2010 <COMPANY>
%       G. Klopman, Delft Hydraulics, 6 Dec 1994


%%
if nargin < 3,
  g = 9.81;
end;

w2 = (w.^2) .* h ./ g;
q  = w2 ./ (1 - exp (-(w2.^(5/4)))) .^ (2/5);

for j=1:2,
  thq     = tanh(q);
  thq2    = 1 - thq.^2;
  a       = (1 - q .* thq) .* thq2;
  b       = thq + q .* thq2;
  c       = q .* thq - w2;
  arg     = (b.^2) - 4 .* a .* c;
  arg     = (-b + sqrt(arg)) ./ (2 * a);
  iq      = find (abs(a.*c) < 1.0e-8 * (b.^2));
  arg(iq) = - c(iq) ./ b(iq);
  q       = q + arg;
end;

k = sign(w) .* q ./ h;

ik    = isnan (k);
k(ik) = zeros(size(k(ik)));




