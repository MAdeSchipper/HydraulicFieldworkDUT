
% Turn the Leica HandHeld generated exported custom data into XYZ.mat in RDNAP
% Before running splitin separate files
% improved

% syntax [XYZ]= leica_out2xyz('outputfilename_via_handheld_customdata.txt',1.33,'quad')
clear all

filename_customdata='NCKsummer19.txt'

fid = fopen(filename_customdata);
pi = 0;
for t = 1:200000
    
    tline = fgetl(fid);
    if tline<0
        break
    end
    if ~strcmp(tline(1:4),'RTCM')   % get rid of RTK reference station LNR, for better plotting
        pi=pi+1;                    % counter for points
         [raw.point_id(pi) raw.surveydate(pi) raw.surveytime(pi) raw.x(pi) raw.y(pi) raw.h(pi) raw.q1(pi) raw.q2(pi) dum]...
             =  strread(tline,'%s %s %s %f %f %f %f %f %s','delimiter','	');
%                [raw.point_id(pi) raw.surveydate(pi) raw.x(pi) raw.y(pi) raw.h(pi) raw.q1(pi) raw.q2(pi)]...
 %                  =  strread(tline,'%s %s %f %f %f %f %f','delimiter','	');
    end
end

%% 
figure
scatter(raw.x,raw.y,10,raw.h,'filled')
xlabel('x [m]')
ylabel('y [m]')
grid on
colorbar

%% map extends
minx=min(raw.x)
maxx=max(raw.x)
miny=min(raw.y)
maxy=max(raw.y)

%%
[Xi,Yi]=meshgrid(minx:5:maxx,miny:5:maxy)
%%
Zi = griddata(raw.x,raw.y,raw.h,Xi,Yi)


%%
figure
pcolor(Xi,Yi,Zi)
shading interp
hold on
scatter(raw.x,raw.y,10,'w','filled')
xlabel('x [m]')
ylabel('y [m]')
grid on
colorbar


