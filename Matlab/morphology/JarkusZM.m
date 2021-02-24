%%% Get jarkus for Zandmotor domain
clear all
close all
clc

loaded2009 =0
loaded2010 =0
loaded2011 =0
loaded2012 =0
%loaded2013 =1 % not in the repo yet

NCfile=     jarkus_url
%NCfile='http://dtvirt5.deltares.nl:8080/thredds/dodsC/opendap/rijkswaterstaat/jarkus/profiles/transect.nc';

transect_ids=nc_varget(NCfile,'id');
begin=10300
eind=11560

ind=find(transect_ids >= 9000000+begin & transect_ids <= 9000000+eind);
transects2use=transect_ids(ind)-9000000;

dates=nc_varget(NCfile,'time_bathy');

%%% 2009
if loaded2009==0
    disp (' *** 2009 *** ')
    % get jarkus data uit database
    Jarkus2009=[];
    for transectid=transects2use'; % ten noorden van de opgang
        try
            disp(['getting data  ' num2str(transectid)])
            d09 = jarkus_readTransectDataNetcdf (NCfile,'Delfland',transectid, 2009);
            Jarkus2009=[Jarkus2009; [d09.xRD d09.yRD d09.zi]];
        end
    end
    Jarkus2009(isnan(Jarkus2009(:,3)),:)=[];
       
    %%Interpolating to grid
    dy_cross=10; dx_alongs=10;
    ZMjarkus_2009.grids.dx_alongs=dx_alongs;    ZMjarkus_2009.grids.dy_cross=dy_cross;
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(7.05E4:dx_alongs:7.5E4,4.50E5:dy_cross:4.55E5);
    Zi=griddata(Jarkus2009(:,1),Jarkus2009(:,2),Jarkus2009(:,3),Xi,Yi) ;
    warning on
    ZMjarkus_2009.grids.Xi=Xi;
    ZMjarkus_2009.grids.Yi=Yi;
    ZMjarkus_2009.grids.Zi=Zi;
    toc
    disp(' done!...')
    
    % Cut off interpolation error
    if 0 % turn on if you want to redefine box
        figure
        th=1; % thinning if needed
        pcolor(ZMjarkus_2009.grids.Xi,ZMjarkus_2009.grids.Yi,ZMjarkus_2009.grids.Zi)
        shading interp
        hold on
        scatter(Jarkus2009(1:th:end,1),Jarkus2009(1:th:end,2),'.k')
        
        %ginput(12)
    end
    Con_data=1.0e+05 *[   0.7051    4.5223
    0.7049    4.4998
    0.7137    4.4999
    0.7286    4.5162
    0.7407    4.5300
    0.7502    4.5409
    0.7502    4.5498
    0.7257    4.5447
    0.7199    4.5397
    0.7146    4.5321
    0.7103    4.5281
    0.7079    4.5242]; % contour datapoints
    
    ind_data = inpolygon(ZMjarkus_2009.grids.Xi,ZMjarkus_2009.grids.Yi,Con_data(:,1),Con_data(:,2));
    ZMjarkus_2009.grids.int_mask=nan(size(ZMjarkus_2009.grids.Zi));
    ZMjarkus_2009.grids.int_mask(ind_data)=ind_data(ind_data==1); % Mask with 1's where there is data, NaN where there is not.
    
    
    %% Rotated
    disp(' Rotating!...')
    x0=72421.90; % opgang keet schelpenpad
    y0=451326.1;
    theta0 =-132+90;  %% gok
    ZMjarkus_2009.origin.x0=x0;
    ZMjarkus_2009.origin.y0=y0;
    ZMjarkus_2009.origin.theta0=theta0;
    
    [Xr,Yr]=RotateTopo(Jarkus2009(:,1),Jarkus2009(:,2),x0,y0,theta0);
    XYZ_r=[ Xr Yr Jarkus2009(:,3)];
    
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(-2500:dx_alongs:4500,-100:dy_cross:1800);
    Zi=griddata(XYZ_r(:,1),XYZ_r(:,2),XYZ_r(:,3),Xi,Yi) ;
    ZMjarkus_2009.grids.Xi_r=Xi;
    ZMjarkus_2009.grids.Yi_r=Yi;
    ZMjarkus_2009.grids.Zi_r=Zi;
    warning on
    toc
    
    % Cut off interpolation error
    if 0
        figure
        th=1
        pcolor( ZMjarkus_2009.grids.Xi_r, ZMjarkus_2009.grids.Yi_r,ZMjarkus_2009.grids.Zi_r)
        shading interp
        hold on
        scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
        %ginput(8)
    end
    
     Con_data=1.0e+03 *[  4.5   -0.0639
   -2.50   -0.0806
   -2.50    1.7861
    4.5     1.7861]; % contour datapoints
     
     ind_data = inpolygon( ZMjarkus_2009.grids.Xi_r, ZMjarkus_2009.grids.Yi_r,Con_data(:,1),Con_data(:,2));
     ZMjarkus_2009.grids.int_mask_r=nan(size(ZMjarkus_2009.grids.Zi_r));
     ZMjarkus_2009.grids.int_mask_r(ind_data)=ind_data(ind_data==1);
    
     %
     dates_yr=dates(45,ind)+datenum('01-01-1970');
     datesthisdata=datestr(round(nanmean(dates_yr)))
     ZMjarkus_2009.matlabdate=datestr(round(nanmean(dates_yr)));
    save ('ZMjarkus_2009.mat','ZMjarkus_2009','Jarkus2009')
else
    load ZMjarkus_2009.mat
end

%fig_org
figure
th=1; % thinning if needed
pcolor(ZMjarkus_2009.grids.Xi,ZMjarkus_2009.grids.Yi,ZMjarkus_2009.grids.Zi.*ZMjarkus_2009.grids.int_mask)
shading interp
hold on
scatter(Jarkus2009(1:th:end,1),Jarkus2009(1:th:end,2),'.k')
colorbar

%fig_rot
figure
th=1; %thinning
pcolor( ZMjarkus_2009.grids.Xi_r, ZMjarkus_2009.grids.Yi_r, ZMjarkus_2009.grids.Zi_r.*ZMjarkus_2009.grids.int_mask_r)
shading interp
hold on
scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
colorbar

%% 2010
if loaded2010==0
    disp (' *** 2010 *** ')
    % get jarkus data uit database
    Jarkus2010=[];
    for transectid=transects2use'; % ten noorden van de opgang
        try
            disp(['getting data  ' num2str(transectid)])
            d09 = jarkus_readTransectDataNetcdf (NCfile,'Delfland',transectid, 2010);
            Jarkus2010=[Jarkus2010; [d09.xRD d09.yRD d09.zi]];
        end
    end
    Jarkus2010(isnan(Jarkus2010(:,3)),:)=[];
    
    %%Interpolating to grid
    dy_cross=10; dx_alongs=10;
    ZMjarkus_2010.grids.dx_alongs=dx_alongs;    ZMjarkus_2010.grids.dy_cross=dy_cross;
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(7.05E4:dx_alongs:7.5E4,4.50E5:dy_cross:4.55E5);
    Zi=griddata(Jarkus2010(:,1),Jarkus2010(:,2),Jarkus2010(:,3),Xi,Yi) ;
    warning on
    ZMjarkus_2010.grids.Xi=Xi;
    ZMjarkus_2010.grids.Yi=Yi;
    ZMjarkus_2010.grids.Zi=Zi;
    toc
    disp(' done!...')
    
    % Cut off interpolation error
    if 0 % turn on if you want to redefine box
        figure
        th=1; % thinning if needed
        pcolor(ZMjarkus_2010.grids.Xi,ZMjarkus_2010.grids.Yi,ZMjarkus_2010.grids.Zi)
        shading interp
        hold on
        scatter(Jarkus2010(1:th:end,1),Jarkus2010(1:th:end,2),'.k')
        colorbar
               
        %ginput(12)
    end
    Con_data=1.0e+05 *[   0.7050    4.5000
    0.7140    4.5000
    0.7306    4.5175
    0.7499    4.5402
    0.7500    4.5498
    0.7303    4.5499
    0.7179    4.5359
    0.7051    4.5207]; % contour datapoints
    
    ind_data = inpolygon(ZMjarkus_2010.grids.Xi,ZMjarkus_2010.grids.Yi,Con_data(:,1),Con_data(:,2));
    ZMjarkus_2010.grids.int_mask=nan(size(ZMjarkus_2010.grids.Zi));
    ZMjarkus_2010.grids.int_mask(ind_data)=ind_data(ind_data==1); % Mask with 1's where there is data, NaN where there is not.
    
    
    %% Rotated
    disp(' Rotating!...')
    x0=72421.90; % opgang keet schelpenpad
    y0=451326.1;
    theta0 =-132+90;  %% gok
    ZMjarkus_2010.origin.x0=x0;
    ZMjarkus_2010.origin.y0=y0;
    ZMjarkus_2010.origin.theta0=theta0;
    
    [Xr,Yr]=RotateTopo(Jarkus2010(:,1),Jarkus2010(:,2),x0,y0,theta0);
    XYZ_r=[ Xr Yr Jarkus2010(:,3)];
    
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(-2500:dx_alongs:4500,-100:dy_cross:1800);
    Zi=griddata(XYZ_r(:,1),XYZ_r(:,2),XYZ_r(:,3),Xi,Yi) ;
    ZMjarkus_2010.grids.Xi_r=Xi;
    ZMjarkus_2010.grids.Yi_r=Yi;
    ZMjarkus_2010.grids.Zi_r=Zi;
    warning on
    toc
    
    % Cut off interpolation error
    if 0
        figure
        th=1
        pcolor( ZMjarkus_2010.grids.Xi_r, ZMjarkus_2010.grids.Yi_r,ZMjarkus_2010.grids.Zi_r)
        shading interp
        hold on
        scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
        %ginput(4)
    end
    
     Con_data=1.0e+03 *[  4.5   -0.0639
   -2.50   -0.0806
   -2.50    1.7861
    4.5     1.7861]; % contour datapoints
     
     ind_data = inpolygon( ZMjarkus_2010.grids.Xi_r, ZMjarkus_2010.grids.Yi_r,Con_data(:,1),Con_data(:,2));
     ZMjarkus_2010.grids.int_mask_r=nan(size(ZMjarkus_2010.grids.Zi_r));
     ZMjarkus_2010.grids.int_mask_r(ind_data)=ind_data(ind_data==1);
    
     %
     dates_yr=dates(46,ind)+datenum('01-01-1970');
     datesthisdata=datestr(round(nanmean(dates_yr)))
     ZMjarkus_2010.matlabdate=datestr(round(nanmean(dates_yr)));
     
        save ('ZMjarkus_2010.mat','ZMjarkus_2010','Jarkus2010')
else
    load ZMjarkus_2010.mat
end

%fig_org
figure
th=1; % thinning if needed
pcolor(ZMjarkus_2010.grids.Xi,ZMjarkus_2010.grids.Yi,ZMjarkus_2010.grids.Zi.*ZMjarkus_2010.grids.int_mask)
shading interp
hold on
scatter(Jarkus2010(1:th:end,1),Jarkus2010(1:th:end,2),'.k')
colorbar

%fig_rot
figure
th=1 %thinning
pcolor( ZMjarkus_2010.grids.Xi_r, ZMjarkus_2010.grids.Yi_r, ZMjarkus_2010.grids.Zi_r.*ZMjarkus_2010.grids.int_mask_r)
shading interp
hold on
scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
colorbar

pause(0.1)
%% 2011
if loaded2011==0
    disp (' *** 2011 *** ')
    % get jarkus data uit database
    Jarkus2011=[];
    for transectid=transects2use'; % ten noorden van de opgang
        try
            disp(['getting data  ' num2str(transectid)])
            d09 = jarkus_readTransectDataNetcdf (NCfile,'Delfland',transectid, 2011);
            Jarkus2011=[Jarkus2011; [d09.xRD d09.yRD d09.zi]];
        end
    end
    Jarkus2011(isnan(Jarkus2011(:,3)),:)=[];
    
    %%Interpolating to grid
    dy_cross=10; dx_alongs=10;
    ZMjarkus_2011.grids.dx_alongs=dx_alongs;    ZMjarkus_2011.grids.dy_cross=dy_cross;
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(7.05E4:dx_alongs:7.5E4,4.50E5:dy_cross:4.55E5);
    Zi=griddata(Jarkus2011(:,1),Jarkus2011(:,2),Jarkus2011(:,3),Xi,Yi) ;
    warning on
    ZMjarkus_2011.grids.Xi=Xi;
    ZMjarkus_2011.grids.Yi=Yi;
    ZMjarkus_2011.grids.Zi=Zi;
    toc
    disp(' done!...')
    
    % Cut off interpolation error
    if 0 % turn on if you want to redefine box
        figure
        th=1; % thinning if needed
        pcolor(ZMjarkus_2011.grids.Xi,ZMjarkus_2011.grids.Yi,ZMjarkus_2011.grids.Zi)
        shading interp
        hold on
        scatter(Jarkus2011(1:th:end,1),Jarkus2011(1:th:end,2),'.k')
        
        %ginput(12)
    end
    Con_data=1.0e+05 *[        0.7050    4.5000
    0.7140    4.5000
    0.7306    4.5175
    0.7499    4.5402
    0.7500    4.5498
    0.7303    4.5499
    0.7179    4.5359
    0.7051    4.5207   ]; % contour datapoints
    
    ind_data = inpolygon(ZMjarkus_2011.grids.Xi,ZMjarkus_2011.grids.Yi,Con_data(:,1),Con_data(:,2));
    ZMjarkus_2011.grids.int_mask=nan(size(ZMjarkus_2011.grids.Zi));
    ZMjarkus_2011.grids.int_mask(ind_data)=ind_data(ind_data==1); % Mask with 1's where there is data, NaN where there is not.
    
    
    %% Rotated
    disp(' Rotating!...')
    x0=72421.90; % opgang keet schelpenpad
    y0=451326.1;
    theta0 =-132+90;  %% gok
    ZMjarkus_2011.origin.x0=x0;
    ZMjarkus_2011.origin.y0=y0;
    ZMjarkus_2011.origin.theta0=theta0;
    
    [Xr,Yr]=RotateTopo(Jarkus2011(:,1),Jarkus2011(:,2),x0,y0,theta0);
    XYZ_r=[ Xr Yr Jarkus2011(:,3)];
    
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(-2500:dx_alongs:4500,-100:dy_cross:1800);
    Zi=griddata(XYZ_r(:,1),XYZ_r(:,2),XYZ_r(:,3),Xi,Yi) ;
    ZMjarkus_2011.grids.Xi_r=Xi;
    ZMjarkus_2011.grids.Yi_r=Yi;
    ZMjarkus_2011.grids.Zi_r=Zi;
    warning on
    toc
    
    % Cut off interpolation error
    if 0
        figure
        th=1
        pcolor( ZMjarkus_2011.grids.Xi_r, ZMjarkus_2011.grids.Yi_r,ZMjarkus_2011.grids.Zi_r)
        shading interp
        hold on
        scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
        %ginput(8)
    end
    
    Con_data=1.0e+03 *[  4.5   -0.0639
   -2.50   -0.0806
   -2.50    1.7861
    4.5     1.7861]; % contour datapoints
     
     ind_data = inpolygon( ZMjarkus_2011.grids.Xi_r, ZMjarkus_2011.grids.Yi_r,Con_data(:,1),Con_data(:,2));
     ZMjarkus_2011.grids.int_mask_r=nan(size(ZMjarkus_2011.grids.Zi_r));
     ZMjarkus_2011.grids.int_mask_r(ind_data)=ind_data(ind_data==1);
    
     %
     dates_yr=dates(47,ind)+datenum('01-01-1970');
     datesthisdata=datestr(round(nanmean(dates_yr)))
     ZMjarkus_2011.matlabdate=datestr(round(nanmean(dates_yr)));
     
      save ('ZMjarkus_2011.mat','ZMjarkus_2011','Jarkus2011')
else
    load ZMjarkus_2011.mat
end

%fig_org
figure
th=1; % thinning if needed
pcolor(ZMjarkus_2011.grids.Xi,ZMjarkus_2011.grids.Yi,ZMjarkus_2011.grids.Zi.*ZMjarkus_2011.grids.int_mask)
shading interp
hold on
scatter(Jarkus2011(1:th:end,1),Jarkus2011(1:th:end,2),'.k')
colorbar

%fig_rot
figure
th=1 %thinning
pcolor( ZMjarkus_2011.grids.Xi_r, ZMjarkus_2011.grids.Yi_r, ZMjarkus_2011.grids.Zi_r.*ZMjarkus_2011.grids.int_mask_r)
shading interp
hold on
scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
colorbar
%%

%% 2012
if loaded2012==0
    disp (' *** 2012 *** ')
    % get jarkus data uit database
    Jarkus2012=[];
    for transectid=transects2use'; % ten noorden van de opgang
        try
            disp(['getting data  ' num2str(transectid)])
            d09 = jarkus_readTransectDataNetcdf (NCfile,'Delfland',transectid, 2012);
            Jarkus2012=[Jarkus2012; [d09.xRD d09.yRD d09.zi]];
        end
    end
    Jarkus2012(isnan(Jarkus2012(:,3)),:)=[];
    
    %%Interpolating to grid
    dy_cross=10; dx_alongs=10;
    ZMjarkus_2012.grids.dx_alongs=dx_alongs;    ZMjarkus_2012.grids.dy_cross=dy_cross;
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(7.05E4:dx_alongs:7.5E4,4.50E5:dy_cross:4.55E5);
    Zi=griddata(Jarkus2012(:,1),Jarkus2012(:,2),Jarkus2012(:,3),Xi,Yi) ;
    warning on
    ZMjarkus_2012.grids.Xi=Xi;
    ZMjarkus_2012.grids.Yi=Yi;
    ZMjarkus_2012.grids.Zi=Zi;
    toc
    disp(' done!...')
    
    % Cut off interpolation error
    if 0 % turn on if you want to redefine box
        figure
        th=1; % thinning if needed
        pcolor(ZMjarkus_2012.grids.Xi,ZMjarkus_2012.grids.Yi,ZMjarkus_2012.grids.Zi)
        shading interp
        hold on
        scatter(Jarkus2012(1:th:end,1),Jarkus2012(1:th:end,2),'.k')
        
        %ginput(12)
    end
    Con_data=1.0e+05 *[    0.7050    4.5000
    0.7140    4.5000
    0.7306    4.5175
    0.7499    4.5402
    0.7500    4.5498
    0.7303    4.5499
    0.7179    4.5359
    0.7051    4.5207]; % contour datapoints
    
    ind_data = inpolygon(ZMjarkus_2012.grids.Xi,ZMjarkus_2012.grids.Yi,Con_data(:,1),Con_data(:,2));
    ZMjarkus_2012.grids.int_mask=nan(size(ZMjarkus_2012.grids.Zi));
    ZMjarkus_2012.grids.int_mask(ind_data)=ind_data(ind_data==1); % Mask with 1's where there is data, NaN where there is not.
    
    
    %% Rotated
    disp(' Rotating!...')
    x0=72421.90; % opgang keet schelpenpad
    y0=451326.1;
    theta0 =-132+90;  %% gok
    ZMjarkus_2012.origin.x0=x0;
    ZMjarkus_2012.origin.y0=y0;
    ZMjarkus_2012.origin.theta0=theta0;
    
    [Xr,Yr]=RotateTopo(Jarkus2012(:,1),Jarkus2012(:,2),x0,y0,theta0);
    XYZ_r=[ Xr Yr Jarkus2012(:,3)];
    
    disp(' gridding!...')
    tic
    warning off
    [Xi,Yi]=meshgrid(-2500:dx_alongs:4500,-100:dy_cross:1800);
    Zi=griddata(XYZ_r(:,1),XYZ_r(:,2),XYZ_r(:,3),Xi,Yi) ;
    ZMjarkus_2012.grids.Xi_r=Xi;
    ZMjarkus_2012.grids.Yi_r=Yi;
    ZMjarkus_2012.grids.Zi_r=Zi;
    warning on
    toc
    
    % Cut off interpolation error
    if 0
        figure
        th=1
        pcolor( ZMjarkus_2012.grids.Xi_r, ZMjarkus_2012.grids.Yi_r,ZMjarkus_2012.grids.Zi_r)
        shading interp
        hold on
        scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
        %ginput(8)
    end
    
    Con_data=1.0e+03 *[  4.5   -0.0639
   -2.50   -0.0806
   -2.50    1.7861
    4.5     1.7861]; % contour datapoints
     
     ind_data = inpolygon( ZMjarkus_2012.grids.Xi_r, ZMjarkus_2012.grids.Yi_r,Con_data(:,1),Con_data(:,2));
     ZMjarkus_2012.grids.int_mask_r=nan(size(ZMjarkus_2012.grids.Zi_r));
     ZMjarkus_2012.grids.int_mask_r(ind_data)=ind_data(ind_data==1);
    %
     dates_yr=dates(48,ind)+datenum('01-01-1970');
     datesthisdata=datestr(round(nanmean(dates_yr)))
     ZMjarkus_2012.matlabdate=datestr(round(nanmean(dates_yr)));
     
      save ('ZMjarkus_2012.mat','ZMjarkus_2012','Jarkus2012')
else
    load ZMjarkus_2012.mat
end

%fig_org
figure
th=1; % thinning if needed
pcolor(ZMjarkus_2012.grids.Xi,ZMjarkus_2012.grids.Yi,ZMjarkus_2012.grids.Zi.*ZMjarkus_2012.grids.int_mask)
shading interp
hold on
scatter(Jarkus2012(1:th:end,1),Jarkus2012(1:th:end,2),'.k')
colorbar

%fig_rot
figure
th=1 %thinning
pcolor( ZMjarkus_2012.grids.Xi_r, ZMjarkus_2012.grids.Yi_r, ZMjarkus_2012.grids.Zi_r.*ZMjarkus_2012.grids.int_mask_r)
shading interp
hold on
scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
colorbar

%% 2013
% if loaded2013==0
%     disp (' *** 2013 *** ')
%     % get jarkus data uit database
%     Jarkus2013=[];
%     for transectid=transects2use'; % ten noorden van de opgang
%         try
%             disp(['getting data  ' num2str(transectid)])
%             d09 = jarkus_readTransectDataNetcdf (NCfile,'Delfland',transectid, 2013);
%             Jarkus2013=[Jarkus2013; [d09.xRD d09.yRD d09.zi]];
%         end
%     end
%     Jarkus2013(isnan(Jarkus2013(:,3)),:)=[];
%     
%     %%Interpolating to grid
%     dy_cross=10; dx_alongs=10;
%     ZMjarkus_2013.grids.dx_alongs=dx_alongs;    ZMjarkus_2013.grids.dy_cross=dy_cross;
%     disp(' gridding!...')
%     tic
%     warning off
%     [Xi,Yi]=meshgrid(7.05E4:dx_alongs:7.5E4,4.50E5:dy_cross:4.55E5);
%     Zi=griddata(Jarkus2013(:,1),Jarkus2013(:,2),Jarkus2013(:,3),Xi,Yi) ;
%     warning on
%     ZMjarkus_2013.grids.Xi=Xi;
%     ZMjarkus_2013.grids.Yi=Yi;
%     ZMjarkus_2013.grids.Zi=Zi;
%     toc
%     disp(' done!...')
%     
%     % Cut off interpolation error
%     if 0 % turn on if you want to redefine box
%         figure
%         th=1; % thinning if needed
%         pcolor(ZMjarkus_2013.grids.Xi,ZMjarkus_2013.grids.Yi,ZMjarkus_2013.grids.Zi)
%         shading interp
%         hold on
%         scatter(Jarkus2013(1:th:end,1),Jarkus2013(1:th:end,2),'.k')
%         
%         %ginput(12)
%     end
%     Con_data=1.0e+05 *[       0.7050    4.5000
%     0.7140    4.5000
%     0.7306    4.5175
%     0.7499    4.5402
%     0.7500    4.5498
%     0.7303    4.5499
%     0.7179    4.5359
%     0.7051    4.5207]; % contour datapoints
%     
%     ind_data = inpolygon(ZMjarkus_2013.grids.Xi,ZMjarkus_2013.grids.Yi,Con_data(:,1),Con_data(:,2));
%     ZMjarkus_2013.grids.int_mask=nan(size(ZMjarkus_2013.grids.Zi));
%     ZMjarkus_2013.grids.int_mask(ind_data)=ind_data(ind_data==1); % Mask with 1's where there is data, NaN where there is not.
%     
%     
%     % Rotated
%     disp(' Rotating!...')
%     x0=72421.90; % opgang keet schelpenpad
%     y0=451326.1;
%     theta0 =-132+90;  %% gok
%     ZMjarkus_2013.origin.x0=x0;
%     ZMjarkus_2013.origin.y0=y0;
%     ZMjarkus_2013.origin.theta0=theta0;
%     
%     [Xr,Yr]=RotateTopo(Jarkus2013(:,1),Jarkus2013(:,2),x0,y0,theta0);
%     XYZ_r=[ Xr Yr Jarkus2013(:,3)];
%     
%     disp(' gridding!...')
%     tic
%     warning off
%     [Xi,Yi]=meshgrid(-2500:dx_alongs:4500,-100:dy_cross:1800);
%     Zi=griddata(XYZ_r(:,1),XYZ_r(:,2),XYZ_r(:,3),Xi,Yi) ;
%     ZMjarkus_2013.grids.Xi_r=Xi;
%     ZMjarkus_2013.grids.Yi_r=Yi;
%     ZMjarkus_2013.grids.Zi_r=Zi;
%     warning on
%     toc
%     
%     % Cut off interpolation error
%     if 0
%         figure
%         th=1
%         pcolor( ZMjarkus_2013.grids.Xi_r, ZMjarkus_2013.grids.Yi_r,ZMjarkus_2013.grids.Zi_r)
%         shading interp
%         hold on
%         scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
%         %ginput(8)
%     end
%     
%      Con_data=1.0e+03 *[  3.9758   -0.0639
%    -1.80   -0.0806
%    -1.80    1.7861
%     3.9758    1.7861]; % contour datapoints
%      
%      ind_data = inpolygon( ZMjarkus_2013.grids.Xi_r, ZMjarkus_2013.grids.Yi_r,Con_data(:,1),Con_data(:,2));
%      ZMjarkus_2013.grids.int_mask_r=nan(size(ZMjarkus_2013.grids.Zi_r));
%      ZMjarkus_2013.grids.int_mask_r(ind_data)=ind_data(ind_data==1);
%     
%       save ('ZMjarkus_2013.mat','ZMjarkus_2013','Jarkus2013')
% else
%     load ZMjarkus_2013.mat
% end
% 
% %fig_org
% figure
% th=1; % thinning if needed
% pcolor(ZMjarkus_2013.grids.Xi,ZMjarkus_2013.grids.Yi,ZMjarkus_2013.grids.Zi.*ZMjarkus_2013.grids.int_mask)
% shading interp
% hold on
% scatter(Jarkus2013(1:th:end,1),Jarkus2013(1:th:end,2),'.k')
% colorbar
% 
% %fig_rot
% figure
% th=1 %thinning
% pcolor( ZMjarkus_2013.grids.Xi_r, ZMjarkus_2013.grids.Yi_r, ZMjarkus_2013.grids.Zi_r.*ZMjarkus_2013.grids.int_mask_r)
% shading interp
% hold on
% scatter(XYZ_r(1:th:end,1),XYZ_r(1:th:end,2),'.k')
% colorbar
%%




