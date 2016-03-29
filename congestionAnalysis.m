function out=congestionAnalysis(TripLat,TripLon,TripTime)
%congestionAnalysis takes lat/lon/time cell arrays for all nodes at a particular
%school, and returns the relationship between congestion and travel time
%for a particular road segment as defined by frequent travel.

%school pos indicates which position (start or end) of the lat/lon vector
%the school is. This may vary (for some reason).

%Input:
%TripLat <cellArray> : each cell is a <double array> of Latitude data for the morning trip
%TripLon <cellArray> : each cell is a <double array> of Longitude data for the morning trip
%TripLon <cellArray> : each cell is a <char array> of time data for the morning trip in the format YYYY-MM-DD HH:MM:SS.F

%Output:
%corridor <structArray>: of micro-trips within the  with fields:
% duration of microtrip (seconds)
% start time of microtrip (MATLAB datenum)
% endtime of microtrip (MATLAB datenum)
% distance of microtrip (km)
% mean velocity of microtrip (km/hr)

% sample usage:

%load('AnonDemoData.mat')
% corridorData=congestionAnalysis(RandTripLat{1},RandTripLon{1},RandTripTime{1}) for first school in sample data
% corridorData=congestionAnalysis(RandTripLat{1},RandTripLon{1},RandTripTime{1}) for second school in sample data

try

close all
clc
    
figure
ax1=gca;
hold on
dcm_obj = datacursormode(gcf);
set(dcm_obj,'enable','on')
set(ax1,'FontSize', 14)
set(dcm_obj,'UpdateFcn',@NewCallback_ds)
ylabel('Lat','FontSize', 14)
xlabel('Lon','FontSize', 14)

figure
ax2=gca;
hold on
dcm_obj = datacursormode(gcf);
set(dcm_obj,'enable','on')
set(dcm_obj,'UpdateFcn',@NewCallback_ds)
set(ax2,'FontSize', 14)
ylabel('Lat','FontSize', 14)

count_full=0;
header={'ID','Time','Lat','Lon'};
total=[];
mean_loc=[];
lats_vec=[];
lons_vec=[];
ind_vec=[];
time_vec=[];

inner_rad=0.005;
outer_rad=0.01;
delta_lat_thresh=0.004;
dist_thresh=0.005; %50 m distance threshold from street start/end to points

for i=1:length(TripLat)%loop through all sensors for a school

    %skip empty cells
    if isempty(TripLat{i})
        continue 
    end

    count_full=count_full+1;

    mat_time=datenum(TripTime{i});

    [sorted_time,sort_idx] = sort(mat_time);

    lats=TripLat{i}(sort_idx);
    lons=TripLon{i}(sort_idx);

    lats_vec=[lats_vec;lats];
    lons_vec=[lons_vec;lons];
    ind_vec=[ind_vec; i*ones(length(sorted_time),1)];
    time_vec=[time_vec;sorted_time];

    h1=plot(ax1,lons,lats,'.');
    plot(ax1,lons(end),lats(end),'*','MarkerSize',12);
    mean_loc(count_full,1)=lons(end);
    mean_loc(count_full,2)=lats(end);

    h2=plot(ax2,sorted_time,lats,'.');
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    datetick(ax2,'x','mm-dd HH:MM','keepticks')
    
end

%determine the school location from the endpoints of most of the trips
mean_school_lon=nanmedian(mean_loc(:,1));
mean_school_lat=nanmedian(mean_loc(:,2));
plot(ax1,mean_school_lon,mean_school_lat,'xk','MarkerSize',20,'LineWidth',3);
viscircles(ax1,[mean_school_lon,mean_school_lat],inner_rad); %inner circle - ~50 m at equator 
viscircles(ax1,[mean_school_lon,mean_school_lat],outer_rad);%outer circle, ~1km at equator

points_outer=pointsincircle({lons_vec lats_vec},outer_rad,[mean_school_lon mean_school_lat]);
points_inner=pointsincircle({lons_vec lats_vec},inner_rad,[mean_school_lon mean_school_lat]);

figure
plot(points_outer.in{1},points_outer.in{2},'.')
hold
plot(points_inner.in{1},points_inner.in{2},'.')

[latlon,row_outer,~]=intersect([points_outer.in{1},points_outer.in{2}],[points_inner.out{1},points_inner.out{2}],'rows');

ind_outer=points_outer.in_ind; %indices of points that are in the outer circle
ind_final=ind_outer(row_outer); %indicies of points that are in the outer circle but out of the inner circle

lat_points=latlon(:,2);
lon_points=latlon(:,1);
node_points=ind_vec(ind_final);
time_points=time_vec(ind_final);

% points within the radius
figure; plot(lon_points,lat_points,'.')

nodes_in=unique(node_points);

%plot individual node traces within radius
figure
hold
for ii=1:length(nodes_in)
    plot(lon_points(node_points==nodes_in(ii)),lat_points(node_points==nodes_in(ii)),'.')

end
set(gca,'FontSize', 14)
ylabel('Lat','FontSize', 14)
xlabel('Lon','FontSize', 14)

% first cluster all points into individual streets
X=[lon_points,lat_points];

Opt=evalclusters(X,'kmeans','CalinskiHarabasz','KList',[1:6]); %accomodate up to 6 roads

K=Opt.OptimalK;

opts = statset('Display','final');
[G,C] = kmeans(X,K,'Distance','cityblock',...
    'Replicates',5,'Options',opts);

clr = lines(K);

figure
scatter(X(:,1),X(:,2),12,clr(G,:),'Marker','.')
hold
scatter(C(:,1), C(:,2), 30, clr, 'Marker','x', 'LineWidth',3)
set(gca,'FontSize', 14)
ylabel('Lat','FontSize', 14)
xlabel('Lon','FontSize', 14)

%using Piecewise linear interpolation for 1-D interpolation (table lookup) http://www.mathworks.com/matlabcentral/fileexchange/40913-piecewise-linear-least-square-fit
figure
ax3=gca;
hold on
dcm_obj = datacursormode(gcf);
set(dcm_obj,'enable','on')
set(dcm_obj,'UpdateFcn',@NewCallback_ds)
set(ax3,'FontSize', 14)
ylabel('Lat','FontSize', 14)
xlabel('Lon','FontSize', 14)

street_start=zeros(K,2);
street_end=zeros(K,2);

for kk=1:K

    XI = linspace(min(lon_points(G==kk)),max(lon_points(G==kk)),30);
    YI = lsq_lut_piecewise( lon_points(G==kk), lat_points(G==kk), XI );
    
    %null events and lat's out of Singapore are removed
    XI(YI>1.8 | YI<1)=[];
    YI(YI>1.8 | YI<1)=[];
    
    %rate of change too high removed
    dlat=abs(diff(YI));
    
    ind_dlat_rem=find(dlat>delta_lat_thresh);
    ind_dlat_rem=ind_dlat_rem+1;
    
    XI(ind_dlat_rem)=[];
    YI(ind_dlat_rem)=[];
    

    plot(ax3, lon_points(G==kk),lat_points(G==kk),'.',XI,YI,'+-','LineWidth',2)
    
    %scenario where the start is defined by the beginning/end of the line.
    %This doesn't particularly handle situations well where the fitted
    %lines are irregular or the roads are circular around the school
    street_start(kk,1)=XI(1);
    street_start(kk,2)=YI(1);
    
    street_end(kk,1)=XI(end);
    street_end(kk,2)=YI(end);

end

% go through each of the nodes to find out their transit times along the routes
figure
hold
ylabel('Lat','FontSize',14);
xlabel('Lon','FontSize',14);
set(gca,'FontSize', 14)
success_count=1;

for jj=1:length(nodes_in)
    
    %initialize the parameters
    node_lons=lon_points(node_points==nodes_in(jj));
    node_lats=lat_points(node_points==nodes_in(jj));
    node_times=time_points(node_points==nodes_in(jj));
    
    node_lonlat=[node_lons,node_lats];
    if length(node_lonlat)> 6 %only run this code if sufficient points are available
        for jk=1:K
            [idx_start,d_start]=knnsearch(node_lonlat,street_start(jk,:),'K',1);
            [idx_end,d_end]=knnsearch(node_lonlat,street_end(jk,:),'K',1);
            if d_start<dist_thresh && d_end<dist_thresh %a successful match
                
                corridor{success_count}.duration=abs((node_times(idx_end)-node_times(idx_start))*3600*24); %in seconds
                
                if corridor{success_count}.duration < 3600 %implying that the trip was during the same day, and during the same trip
                    corridor{success_count}.start=node_times(idx_start);
                    corridor{success_count}.end=node_times(idx_end);

                    dist=lldistkm_vec([node_lats(idx_start),node_lons(idx_start)],[node_lats(idx_end),node_lons(idx_end)]);
                    corridor{success_count}.distance=dist;
                    corridor{success_count}.meanspeed=dist/(corridor{success_count}.duration/3600);
                    
                    if corridor{success_count}.meanspeed < 30 % threshold to ensure that meaningful mean speeds are obtained
                        % if not, overwrite previous success
                        success_count=success_count+1;
                        plot(node_lons(idx_start),node_lats(idx_start),'xr','MarkerSize',20);
                        plot(node_lons(idx_end),node_lats(idx_end),'xg','MarkerSize',20);
                    end
                end
            end
        end
    end
    
    plot(node_lons,node_lats,'.');
    

end

%plot the start times versus the duration of the trip
figure
ax4=gca;
ylabel('Absolute Duration (seconds)');
xlabel('Hour of travel');
hold

figure
ax5=gca;
ylabel('Distance (km)');
xlabel('Hour of travel');
hold

figure;
ax6=gca;
ylabel('Mean Speed (km/hr)');
xlabel('Hour of travel');
hold

figure;
ax7=gca;
ylabel('Mean Speed (km/hr)');
xlabel('Distance (km)');
hold

hrtm=zeros(success_count-1,1);
mn_vel=zeros(success_count-1,1);
dist=zeros(success_count-1,1);

for ss=1:success_count-1
    [~,~,~,H,MN,~] = datevec(corridor{ss}.start);

    plot(ax4,H+MN/60,corridor{ss}.duration,'.') %time spent in corridor
    plot(ax5,H+MN/60,corridor{ss}.distance,'.') %distance travelled in corridor
    plot(ax6,H+MN/60,corridor{ss}.meanspeed,'*') % mean speed in corridor
    plot(ax7,corridor{ss}.distance,corridor{ss}.meanspeed,'*') %distanc vs speed
    
    hrtm(ss)=H+MN/60;
    mn_vel(ss)=corridor{ss}.meanspeed;
    dist(ss)=corridor{ss}.distance;
end

idx_in=find(hrtm>6&hrtm<8);
meanspeed_in=mn_vel(idx_in);
dist_in=dist(idx_in);
hrtm_in=hrtm(idx_in);
[traveler_counts,traveler_edges,bin_idx] = histcounts(hrtm_in,24);
figure
centers=mean([traveler_edges(1:end-1);traveler_edges(2:end)]);
bar(centers,traveler_counts);
ylabel('Frequency','FontSize',14);
xlabel('Time of Travel (Fractional Hours from Midnight)','FontSize',14);
set(gca,'FontSize', 14)

for bs=1:length(centers)
    if sum(bin_idx==bs)~=0
        speed_tm(bs)=mean(meanspeed_in(bin_idx==bs));
        dist_tm(bs)=mean(dist_in(bin_idx==bs));
    end
end
figure
plot(centers,speed_tm,'.')
ylim([0,30])
ylabel('Mean Speed (km/hr)','FontSize',14);
xlabel('Time of Travel (Fractional Hours from Midnight)','FontSize',14);
set(gca,'FontSize', 14)
%get the coefficients and R^2 of the linear model
mdl=fitlm(centers,speed_tm)

figure
plot(dist_tm,speed_tm,'.')
ylim([0,30])
ylabel('Mean Speed (km/hr)','FontSize',14);
xlabel('Mean Segment Distance (km)','FontSize',14);
set(gca,'FontSize', 14)

out=corridor;

catch ME %debugging interface
    disp(ME.message)
    keyboard
end