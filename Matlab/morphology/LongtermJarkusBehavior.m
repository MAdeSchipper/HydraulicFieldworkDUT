%% Simple file to download and show Jarkus profile over the years.
%Requires Open Earth


%NCfile='http://dtvirt5.deltares.nl:8080/thredds/dodsC/opendap/rijkswaterstaat/jarkus/profiles/transect.nc';
NCfile='http://opendap.deltares.nl/thredds/dodsC/opendap/rijkswaterstaat/jarkus/profiles/transect_r20171124.nc'

nc_dump(NCfile)

transectid=10883;

%% plot jaren
clear ztot xtot indyr
figure
hold all
set(gca,'XDir','reverse')
xlim([00 1400])
%ylim([400 800])
xlabel(' cross shore distance [m]' )
ylabel(' elevation [m NAP]' )
years=1965:2014
for i=1:length(years)
    iyr=years(i);
    try
        datayear = jarkus_readTransectDataNetcdf (NCfile,'Delfland',transectid, iyr);
        indyr=find(~isnan(datayear.zi)==1);
       plot(datayear.xi(indyr),datayear.zi(indyr))
        xtot(:,i)=datayear.xi;
        ztot(:,i)=datayear.zi;                
    end
   
end

%%
disp(' calculating slope...')
up=1;
down=-4;

dz=0.1;

b_b=up-0.5*dz;
b_up=up+0.5*dz;

b_b2=down-0.5*dz;
b_up2=down+0.5*dz;

B_offs_all= 1300;
B_ons_all = -60;
i_transect=1
  for i_date=1:length(xtot(1,:))
      try
            xi=xtot(:,i_date);
            zi=ztot(:,i_date);
            Vol_up(i_date,i_transect)=jarkus_getVolume(xi, zi, b_up, b_b ,B_ons_all,B_offs_all);
     x_up(i_date,i_transect)=Vol_up(i_date,i_transect)/(b_up-b_b);
     
                 Vol_down(i_date,i_transect)=jarkus_getVolume(xi, zi, b_up2, b_b2 ,B_ons_all,B_offs_all);
     x_down(i_date,i_transect)=Vol_down(i_date,i_transect)/(b_up2-b_b2);
%      if i_date<=23 && i_transect==11
%         x_up(i_date,i_transect)=NaN;
%      end
       slope(i_date,i_transect)=(x_down(i_date,i_transect)-x_up(i_date,i_transect))/(up - down);
      end
  end
disp('slope calculated!')

%%
figure
plot(years,slope)

mean(slope(1:21))

%%
disp(' calculating volume...')

b_b=-10;
b_up=3;


B_offs_all= 800;
B_ons_all = -100;
i_transect=1
  for i_date=1:length(xtot(1,:))
      try
            xi=xtot(:,i_date);
            zi=ztot(:,i_date);
            Vol(i_date,i_transect)=jarkus_getVolume(xi, zi, b_up, b_b ,B_ons_all,B_offs_all);
   end
  end
disp('vol calculated!')
%%
figure
plot(years, Vol(:,i_transect))