function make_clm5_input_stations
%{
Function to create all clm5 met data input files required for PTCLM
simulations, based on OSHD COSMO1E data (2015/10-2020/07). At the end of
the script we create dummy data for one additional season (just use for
spin up).

This is the final script used to create PTCLM5 input files. It's run for
most stations  (575), albeit we later only use about 50 of them.

created by JM / June 2021

output:
monthly files including all met variables needed for simulations.
%}

clearvars
close all

%% basefolder output netcdf
bf_out='/media/malle/LaCie/CLM5_input/input_stations_OSHD_v2';
%% define basefolder of OSHD data
%station info file
station_info=load('/media/malle/LaCie/OSHD_info/STAT_LIST.mat');
%get info for each
loc_x=station_info.statlist.x';
loc_y=station_info.statlist.y';
elev=station_info.statlist.z';

%transform ch1903 to lat/long
xd = (loc_x - 600000)/1000000;
yd = (loc_y - 200000)/1000000;
lon = (2.6779094 + 4.728982 .* xd + 0.791484 .* xd .* yd + 0.1306 .* xd .* yd.^2 - 0.0436 .* xd.^3) .* 100 ./ 36;
lat = (16.9023892 + 3.238272 .* yd - 0.270978 .* xd.^2 - 0.002528 .* yd.^2 - 0.0447 .* xd.^2 .* yd - 0.0140 .* yd.^3) .* 100 ./ 36;

edge_n=ceil(lat); % north edge of grid cell within which point of interest resides
edge_s=floor(lat); % south edge of grid cell within which point of interest resides
edge_w=floor(lon); % west edge of grid cell within which point of interest resides
edge_e=ceil(lon); % east edge of grid cell within which point of interest resides
% bf_oshd='smb://phoenix/oshd_archive/DATA_COSMO_1E/PROCESSED_STAT_ANALYSIS/COSMO1';
bf_oshd='/media/malle/LaCie/OSHD_data_phoenix/OSHD_station_data';
dir_data_time=dir(bf_oshd);%what files are in there?
tf = ismember( {dir_data_time.name}, {'.', '..'});
dir_data_time(tf) = []; %remove current and parent directory

year_months={dir_data_time.name}';

for tii = 1:size(year_months,1)

    year_month=year_months{tii};
    time_stamp=[year_month(1:4) '-' year_month(6:7)]; %for file naming

    bf_data=fullfile(bf_oshd,year_month);
    dir_data=dir(bf_data);%what files are in there?
    tf = ismember( {dir_data.name}, {'.', '..'});
    dir_data(tf) = []; %remove current and parent directory

    %% get file names
    names_data={dir_data.name}';

    idx=(contains(names_data,'prec'));
    names_precip=names_data(idx);
    idx_wind=(contains(names_data,'wcor')); %corrected wind speed, m/s I think
    names_wcor=names_data(idx_wind);
    idx_temp=(contains(names_data,'tcor')); %corrected air temp, deg C
    names_temp=names_data(idx_temp);
    idx_rh=(contains(names_data,'rhum')); %relative humidity, %
    names_rh=names_data(idx_rh);
    idx_lwr=(contains(names_data,'ilwc')); %LWR, W/m2
    names_lwr=names_data(idx_lwr);
    idx_pair=(contains(names_data,'pair')); %pressure, hPa I think
    names_pair=names_data(idx_pair);
    idx_dif=(contains(names_data,'idfc'));
    names_dif=names_data(idx_dif);
    idx_dir=(contains(names_data,'idrc'));
    names_dir=names_data(idx_dir);

    if strcmp(time_stamp,'2016-02')||strcmp(time_stamp,'2020-02') %get rid of extra day of leap year...
        names_wcor(673:end)=[];
        names_temp(673:end)=[];
        names_rh(673:end)=[];
        names_lwr(673:end)=[];
        names_pair(673:end)=[];
        names_dif(673:end)=[];
        names_dir(673:end)=[];
        names_precip(673:end)=[];
    end

    %time stamps (does not need to be repeated for remaining variables)
    time_stamps=datetime(cellfun(@(x) x(6:17),names_precip,'UniformOutput',false),'InputFormat','yyyyMMddHHmm');
    length_time=size(time_stamps,1);
    time_days=days(time_stamps-time_stamps(1));%need time in days since first of month

    pair_oshd=nan(length(loc_x),length_time);
    lwr_oshd=nan(length(loc_x),length_time);
    rh_oshd=nan(length(loc_x),length_time);
    tair_oshd=nan(length(loc_x),length_time);
    wind_oshd=nan(length(loc_x),length_time);
    prec_oshd=nan(length(loc_x),length_time);
    swr_oshd=nan(length(loc_x),length_time);

    for mix=1:size(names_pair,1)
        load_pair=load(fullfile(bf_oshd,year_month,cell2mat(names_pair(mix))));
        load_lwr=load(fullfile(bf_oshd,year_month,cell2mat(names_lwr(mix))));
        load_rh=load(fullfile(bf_oshd,year_month,cell2mat(names_rh(mix))));
        load_tair=load(fullfile(bf_oshd,year_month,cell2mat(names_temp(mix))));
        load_wind=load(fullfile(bf_oshd,year_month,cell2mat(names_wcor(mix))));
        load_dir=load(fullfile(bf_oshd,year_month,cell2mat(names_dir(mix)))); %units w/m2
        load_dif=load(fullfile(bf_oshd,year_month,cell2mat(names_dif(mix)))); %units w/m2
        load_precip=load(fullfile(bf_oshd,year_month,cell2mat(names_precip(mix)))); %units are mm pro zeitschritt = mm/h => /3600 to obtain mm/s

        iswr=(load_dir.stat.data)'+(load_dif.stat.data)';%need to transpose to fit surface data
        precip_mm_s=(load_precip.stat.data/3600)';%need to transpose to fit surface data

        pair=(load_pair.stat.data)'; %might have to convert... check units (No, this is OK!)
        lwr=(load_lwr.stat.data)';
        rh=(load_rh.stat.data)';
        tair=(load_tair.stat.data+273.15)'; %convert to Kelvin here
        wind=(load_wind.stat.data)';

        pair_oshd(:,mix)=pair;
        lwr_oshd(:,mix)=lwr;
        rh_oshd(:,mix)=rh;
        tair_oshd(:,mix)=tair;
        wind_oshd(:,mix)=wind;
        swr_oshd(:,mix)=iswr;
        prec_oshd(:,mix)=precip_mm_s;
    end

    for miix=1:length(loc_x) %now need to make seperate file for each station
        lat_in=lat(miix); %latitude of point of interest
        long_in=lon(miix)+180; %longitude of point of interest.. put on 0 to 360 grid
        elev_in=elev(miix); %elevation of point of interest

        edge_n_in=edge_n(miix); % north edge of grid cell within which point of interest resides
        edge_s_in=edge_s(miix); % south edge of grid cell within which point of interest resides
        edge_w_in=edge_w(miix); % west edge of grid cell within which point of interest resides
        edge_e_in=edge_e(miix); % east edge of grid cell within which point of interest resides

        stat_name=station_info.statlist.acro{1,miix};
        stat_long_name=station_info.statlist.name{1,miix};

        if contains(stat_name,'*') %for some stations have several occurances, need to select which ones we are interested in
            stat_name=[stat_name(2:end) '_JM'];
        end

        %contains(stat_name,'#') ||  contains(stat_name,'*') || % have some
        %stations with this so redid it...

        if  contains(stat_name,'$') || contains(stat_name,'+')
            %disp(['station ' stat_long_name ' ignored for now... check later!'])
        else

            if ~exist(fullfile(bf_out,stat_name),'dir')
                mkdir(fullfile(bf_out,stat_name))
            end
            nc_file=fullfile(bf_out,stat_name,['ptclm5_' time_stamp '.nc']);

            if exist(nc_file,'file') == 2
                delete(nc_file)
            else

                nccreate(nc_file,'time','Dimensions',{'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'time', 'long_name', 'observation time');
                ncwriteatt(nc_file, 'time', 'units', ['days since ' datestr(time_stamps(1),31)]);
                ncwriteatt(nc_file, 'time', 'calendar', 'noleap');
                ncwrite(nc_file,'time',time_days);

                %long+lat
                nccreate(nc_file,'LONGXY','Dimensions',{'lon',1,'lat',1},'Format','classic','Datatype','double')
                ncwriteatt(nc_file, 'LONGXY', 'long_name', 'longitude');
                ncwriteatt(nc_file, 'LONGXY', 'units', 'degrees_east');
                ncwriteatt(nc_file, 'LONGXY', 'mode', 'time-invariant');
                ncwrite(nc_file,'LONGXY',long_in);

                nccreate(nc_file,'LATIXY','Dimensions',{'lon',1,'lat',1},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'LATIXY', 'long_name', 'latitude');
                ncwriteatt(nc_file, 'LATIXY', 'units', 'degrees_north');
                ncwriteatt(nc_file, 'LATIXY', 'mode', 'time-invariant');
                ncwrite(nc_file,'LATIXY',lat_in);

                %4 edges of grid
                nccreate(nc_file,'EDGEN','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'EDGEN', 'long_name', 'northern edge in atmospheric data');
                ncwriteatt(nc_file, 'EDGEN', 'units', 'degrees_north');
                ncwriteatt(nc_file, 'EDGEN', 'mode', 'time-invariant');
                ncwrite(nc_file,'EDGEN',edge_n_in);

                nccreate(nc_file,'EDGES','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'EDGES', 'long_name', 'southern edge in atmospheric data');
                ncwriteatt(nc_file, 'EDGES', 'units', 'degrees_north');
                ncwriteatt(nc_file, 'EDGES', 'mode', 'time-invariant');
                ncwrite(nc_file,'EDGES',edge_s_in);

                nccreate(nc_file,'EDGEE','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'EDGEE', 'long_name', 'eastern edge in atmospheric data');
                ncwriteatt(nc_file, 'EDGEE', 'units', 'degrees_east');
                ncwriteatt(nc_file, 'EDGEE', 'mode', 'time-invariant');
                ncwrite(nc_file,'EDGEE',edge_e_in);

                nccreate(nc_file,'EDGEW','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'EDGEW', 'long_name', 'western edge in atmospheric data');
                ncwriteatt(nc_file, 'EDGEW', 'units', 'degrees_east');
                ncwriteatt(nc_file, 'EDGEW', 'mode', 'time-invariant');
                ncwrite(nc_file,'EDGEW',edge_w_in);

                %% now individual data -- precip
                nccreate(nc_file,'PRECTmms','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'PRECTmms', 'long_name', 'PRECTmms total precipitation');
                ncwriteatt(nc_file, 'PRECTmms', 'units', 'mm H2O / sec');
                ncwriteatt(nc_file, 'PRECTmms', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'PRECTmms', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'PRECTmms', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'PRECTmms',prec_oshd(miix,k),[1 1 k]);
                end
                %% now individual data -- solar
                nccreate(nc_file,'FSDS','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'FSDS', 'long_name', 'total incident solar radiation');
                ncwriteatt(nc_file, 'FSDS', 'units', 'W/m**2');
                ncwriteatt(nc_file, 'FSDS', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'FSDS', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'FSDS', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'FSDS',swr_oshd(miix,k),[1 1 k]);
                end
                %% now individual data -- tair_rh
                nccreate(nc_file,'TBOT','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'TBOT', 'long_name', 'temperature at the lowest atm level');
                ncwriteatt(nc_file, 'TBOT', 'units', 'K');
                ncwriteatt(nc_file, 'TBOT', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'TBOT', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'TBOT', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'TBOT',tair_oshd(miix,k),[1 1 k]);
                end

                nccreate(nc_file,'RH','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'RH', 'long_name', 'relative humidity at the lowest atm level (RH)');
                ncwriteatt(nc_file, 'RH', 'units', '%');
                ncwriteatt(nc_file, 'RH', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'RH', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'RH', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'RH',rh_oshd(miix,k),[1 1 k]);
                end
                %% now individual data -- wind_zbot
                nccreate(nc_file,'WIND','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'WIND', 'long_name', 'wind at the lowest atm level');
                ncwriteatt(nc_file, 'WIND', 'units', 'm/s');
                ncwriteatt(nc_file, 'WIND', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'WIND', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'WIND', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'WIND',wind_oshd(miix,k),[1 1 k]);
                end

                z_meas=repmat(10,length_time,1); %m, observational height
                nccreate(nc_file,'ZBOT','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'ZBOT', 'long_name', 'observational height');
                ncwriteatt(nc_file, 'ZBOT', 'units', 'm');
                ncwriteatt(nc_file, 'ZBOT', 'mode', 'time-independent');
                for k=1:length_time
                    ncwrite(nc_file,'ZBOT',z_meas(k),[1 1 k]);
                end
                %% now individual data -- lwr_pa
                nccreate(nc_file,'FLDS','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'FLDS', 'long_name', 'incident longwave radiation');
                ncwriteatt(nc_file, 'FLDS', 'units', 'W/m**2');
                ncwriteatt(nc_file, 'FLDS', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'FLDS', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'FLDS', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'FLDS',lwr_oshd(miix,k),[1 1 k]);
                end

                nccreate(nc_file,'PSRF','Dimensions',{'lon',1,'lat',1,'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(nc_file, 'PSRF', 'long_name', 'surface pressure at the lowest atm level');
                ncwriteatt(nc_file, 'PSRF', 'units', 'Pa');
                ncwriteatt(nc_file, 'PSRF', 'mode', 'time-dependent');
                ncwriteatt(nc_file, 'PSRF', '_FillValue',9.999999616903162e+35);
                ncwriteatt(nc_file, 'PSRF', 'missing_value',9.999999616903162e+35);
                for k=1:length_time
                    ncwrite(nc_file,'PSRF',pair_oshd(miix,k),[1 1 k]);
                end

                %global attributes for entire .nc file
                ncwriteatt(nc_file,'/','institution','WSL Institute for snow and avalanche research SLF, snow hydrology group');
                ncwriteatt(nc_file,'/','history',['File Origin - This file was created by J.Malle on ' datestr(now)]);
                ncwriteatt(nc_file,'/','temp. resolution','1-hourly, based on COSMO1E analysis data');
                ncwriteatt(nc_file,'/','site_location',['Latitude: ' num2str(lat_in) ' Longitude: '...
                    num2str(long_in) ' Elevation (masl): ' num2str(elev_in) ' Station Name: ' stat_long_name]);
            end
        end
    end
end

%% if interested in dummy data, also run this part:
bf_all='/media/malle/LaCie/CLM5_input/input_stations_OSHD_v2';

run_data=dir(bf_all);
tf=ismember({run_data.name},{'.','..'});
run_data(tf)=[];
files={run_data.name};
old_year=[repmat(2017,4,1);repmat(2018,9,1)];
new_year=[repmat(2014,4,1);repmat(2015,9,1)];
old_month=[9:12,1:9]';

for mii=1:size(run_data,1)
    bf=fullfile(bf_all,files{mii});
    mkdir(fullfile(bf,'dummy_data'));
    for sii=1:size(old_year)
        old_file=fullfile(bf,['ptclm5_' num2str(old_year(sii)) '-' num2str(old_month(sii),'%02.f') '.nc']);
        new_file=fullfile(bf,'dummy_data',['ptclm5_' num2str(new_year(sii)) '-' num2str(old_month(sii),'%02.f') '.nc']);
        copyfile(old_file,new_file);
        new_date=datetime(new_year(sii),old_month(sii),1);
        ncwriteatt(new_file, 'time', 'units', ['days since ' datestr(new_date,31)]);
    end
end

end