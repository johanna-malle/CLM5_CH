function make_clm5_input_v2
%{
Function to create all clm5 met data input files required for swiss wide
simulations

created by JM / April 2021

26.7.2021: had to update surface dataset because of 0-360 grid

output:
precipitation file for each month
solar radiation file for each month
remaining met variable file for each month
%}
clearvars
close all

%% basefolder output netcdf
bf_out='/media/malle/LaCie/CLM5_input/OSHD_CH1km_v3';
%% define basefolder of OSHD data
bf_oshd='/media/malle/LaCie/OSHD_data_phoenix/OSHD_data/';
dir_data_time=dir(bf_oshd);%what files are in there?
tf = ismember( {dir_data_time.name}, {'.', '..'});
dir_data_time(tf) = []; %remove current and parent directory

year_months={dir_data_time.name}';

%reference surface dataset
%for lat long, get surface dataset
surf_data='/media/malle/LaCie/CLM5_input/domain.lnd.CH_1km_v4_navy.210727.nc';
surf_LONGXY=ncread(surf_data,'xc');
surf_LATIXY=ncread(surf_data,'yc');
edge_w=min(min(surf_LONGXY));
edge_e=max(max(surf_LONGXY));
edge_s=min(min(surf_LATIXY));
edge_n=max(max(surf_LATIXY));

%already make mask here
idfc_load=load(fullfile(bf_oshd,'2018.04/idfc_201804010800_cosmo1_A_201804020000.mat'));
data=idfc_load.grid.data';
idx_data=(data==-9999);
first_noNaN=nan(size(idx_data,1),1);
for siix=1:size(idx_data,1)
    first_noNaN(siix)=find(idx_data(siix,:)==0,1);
end


for tii = 41:size(year_months,1)

    year_month=year_months{tii};

    time_stamp=[year_month(1:4) '-' year_month(6:7)]; %for file naming

    bf_data=fullfile(bf_oshd,year_month);
    dir_data=dir(bf_data);%what files are in there?
    tf = ismember( {dir_data.name}, {'.', '..'});
    dir_data(tf) = []; %remove current and parent directory

    %% precipitation
    %name of netcdf file
    name_gen='clmforc.OSHD1km.Prec.';
    nc_precip=fullfile(bf_out,'precip',[name_gen time_stamp '.nc']);

    %get file names
    names_data={dir_data.name}';
    idx=(contains(names_data,'prec'));
    names_precip=names_data(idx);
    if strcmp(time_stamp,'2016-02')||strcmp(time_stamp,'2020-02') %get rid of extra day of leap year...
     names_precip(673:end)=[];
    end

    %time stamps (does not need to be repeated for remaining variables)
    time_stamps=datetime(cellfun(@(x) x(6:17),names_precip,'UniformOutput',false),'InputFormat','yyyyMMddHHmm');
    length_time=size(time_stamps,1);
    time_days=days(time_stamps-time_stamps(1));%need time in days since first of month

    %get precip data - open all precip .mat files and create matrix
    prec_oshd=nan(365,272,length_time);
    for mix=1:size(names_precip,1)
        name_mix=names_precip(mix);
        load_precip=load(fullfile(bf_oshd,year_month,cell2mat(name_mix))); %units are mm pro zeitschritt = mm/h => /3600 to obtain mm/s
        precip_mm_s=(load_precip.grid.data/3600)';%need to transpose to fit surface data

        %overwrite nans: take average of non-nan value for nan fields...
        for ii=1:size(precip_mm_s,1)
            cut_line=precip_mm_s(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            precip_mm_s(ii,:)=cut_line;
        end
        prec_oshd(:,:,mix)=precip_mm_s;
    end
    %% solar - dif and dir stored seperately
    %name of netcdf file
    name_gen_solar='clmforc.OSHD1km.Solr.';
    nc_solar=fullfile(bf_out,'solar',[name_gen_solar time_stamp '.nc']);

    %get file names
    idx_dif=(contains(names_data,'idfc'));
    names_dif=names_data(idx_dif);
    idx_dir=(contains(names_data,'idrc'));
    names_dir=names_data(idx_dir);

    if strcmp(time_stamp,'2016-02')||strcmp(time_stamp,'2020-02') %get rid of extra day of leap year...
        names_dif(673:end)=[];
        names_dir(673:end)=[];
    end

    %get dir+dif radiation data
    swr_oshd=nan(365,272,length_time);
    for mix=1:size(names_dif,1)
        name_mix_dir=names_dir(mix);
        name_mix_dif=names_dif(mix);

        load_dir=load(fullfile(bf_oshd,year_month,cell2mat(name_mix_dir))); %units w/m2
        load_dif=load(fullfile(bf_oshd,year_month,cell2mat(name_mix_dif))); %units w/m2

        iswr=(load_dir.grid.data)'+(load_dif.grid.data)';%need to transpose to fit surface data
        for ii=1:size(iswr,1)
            cut_line=iswr(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            iswr(ii,:)=cut_line;
        end
        swr_oshd(:,:,mix)=iswr;
    end
    %% remaining file names
    %tair+rh
    name_gen_tq='clmforc.OSHD1km.TQ.';
    nc_tq=fullfile(bf_out,'tair_rh',[name_gen_tq time_stamp '.nc']);

    %wind+zbot
    name_gen_wz='clmforc.OSHD1km.WZ.';
    nc_wz=fullfile(bf_out,'wind_zbot',[name_gen_wz time_stamp '.nc']);

    %lwr+pa
    name_gen_lp='clmforc.OSHD1km.LP.';
    nc_lp=fullfile(bf_out,'lwr_pa',[name_gen_lp time_stamp '.nc']);

    %% get file names
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

    if strcmp(time_stamp,'2016-02')||strcmp(time_stamp,'2020-02') %get rid of extra day of leap year...
        names_wcor(673:end)=[];
        names_temp(673:end)=[];
        names_rh(673:end)=[];
        names_lwr(673:end)=[];
        names_pair(673:end)=[];
    end

    pair_oshd=nan(365,272,length_time);
    lwr_oshd=nan(365,272,length_time);
    rh_oshd=nan(365,272,length_time);
    tair_oshd=nan(365,272,length_time);
    wind_oshd=nan(365,272,length_time);

    for mix=1:size(names_pair,1)
        load_pair=load(fullfile(bf_oshd,year_month,cell2mat(names_pair(mix))));
        load_lwr=load(fullfile(bf_oshd,year_month,cell2mat(names_lwr(mix))));
        load_rh=load(fullfile(bf_oshd,year_month,cell2mat(names_rh(mix))));
        load_tair=load(fullfile(bf_oshd,year_month,cell2mat(names_temp(mix))));
        load_wind=load(fullfile(bf_oshd,year_month,cell2mat(names_wcor(mix))));

        pair=(load_pair.grid.data)'; %might have to convert... check units
        lwr=(load_lwr.grid.data)';
        rh=(load_rh.grid.data)';
        tair=(load_tair.grid.data+273.15)'; %convert to Kelvin here
        wind=(load_wind.grid.data)';

        %overwrite nans: take average of non-nan value for nan fields...
        %pair:
        for ii=1:size(pair,1)
            cut_line=pair(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            pair(ii,:)=cut_line;
        end
        pair_oshd(:,:,mix)=pair;

        %lwr
        for ii=1:size(lwr,1)
            cut_line=lwr(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            lwr(ii,:)=cut_line;
        end
        lwr_oshd(:,:,mix)=lwr;

        %rh
        for ii=1:size(rh,1)
            cut_line=rh(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            rh(ii,:)=cut_line;
        end
        rh_oshd(:,:,mix)=rh;

        %tair
        for ii=1:size(tair,1)
            cut_line=tair(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            tair(ii,:)=cut_line;
        end
        tair_oshd(:,:,mix)=tair;

        %wind
        for ii=1:size(wind,1)
            cut_line=wind(ii,:);
            id_nan=find(idx_data(ii,:)==1);
            cut_line(id_nan)=cut_line(first_noNaN(ii));
            wind(ii,:)=cut_line;
        end
        wind_oshd(:,:,mix)=wind;
    end
    %% write netcdf files - time, lat, long, edges same for all three => loop
    nc_combo={nc_precip,nc_solar,nc_tq,nc_wz,nc_lp};

    for dix=1:5 %loop through precip, solar, rest file
        nc_file=cell2mat(nc_combo(dix));
        %time
        nccreate(nc_file,'time','Dimensions',{'time',length_time},'Format','classic','Datatype','double');
        ncwriteatt(nc_file, 'time', 'long_name', 'observation time');
        ncwriteatt(nc_file, 'time', 'units', ['days since ' datestr(time_stamps(1),31)]);
        ncwriteatt(nc_file, 'time', 'calendar', 'noleap');
        ncwrite(nc_file,'time',time_days);

        %long+lat
        nccreate(nc_file,'LONGXY','Dimensions',{'lon',365,'lat',272},'Format','classic','Datatype','double')
        ncwriteatt(nc_file, 'LONGXY', 'long_name', 'longitude');
        ncwriteatt(nc_file, 'LONGXY', 'units', 'degrees_east');
        ncwriteatt(nc_file, 'LONGXY', 'mode', 'time-invariant');
        ncwrite(nc_file,'LONGXY',surf_LONGXY);

        nccreate(nc_file,'LATIXY','Dimensions',{'lon',365,'lat',272},'Format','classic','Datatype','double');
        ncwriteatt(nc_file, 'LATIXY', 'long_name', 'latitude');
        ncwriteatt(nc_file, 'LATIXY', 'units', 'degrees_north');
        ncwriteatt(nc_file, 'LATIXY', 'mode', 'time-invariant');
        ncwrite(nc_file,'LATIXY',surf_LATIXY);

        %4 edges of grid
        nccreate(nc_file,'EDGEN','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
        ncwriteatt(nc_file, 'EDGEN', 'long_name', 'northern edge in atmospheric data');
        ncwriteatt(nc_file, 'EDGEN', 'units', 'degrees_north');
        ncwriteatt(nc_file, 'EDGEN', 'mode', 'time-invariant');
        ncwrite(nc_file,'EDGEN',edge_n);

        nccreate(nc_file,'EDGES','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
        ncwriteatt(nc_file, 'EDGES', 'long_name', 'southern edge in atmospheric data');
        ncwriteatt(nc_file, 'EDGES', 'units', 'degrees_north');
        ncwriteatt(nc_file, 'EDGES', 'mode', 'time-invariant');
        ncwrite(nc_file,'EDGES',edge_s);

        nccreate(nc_file,'EDGEE','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
        ncwriteatt(nc_file, 'EDGEE', 'long_name', 'eastern edge in atmospheric data');
        ncwriteatt(nc_file, 'EDGEE', 'units', 'degrees_east');
        ncwriteatt(nc_file, 'EDGEE', 'mode', 'time-invariant');
        ncwrite(nc_file,'EDGEE',edge_e);

        nccreate(nc_file,'EDGEW','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
        ncwriteatt(nc_file, 'EDGEW', 'long_name', 'western edge in atmospheric data');
        ncwriteatt(nc_file, 'EDGEW', 'units', 'degrees_east');
        ncwriteatt(nc_file, 'EDGEW', 'mode', 'time-invariant');
        ncwrite(nc_file,'EDGEW',edge_w);

        %global attributes for entire .nc file
        ncwriteatt(nc_file,'/','institution','WSL Institute for snow and avalanche research SLF, snow hydrology group');
        ncwriteatt(nc_file,'/','history',['File Origin - This file was created by J.Malle on ' datestr(now)]);
        ncwriteatt(nc_file,'/','grid location','Switzerland, based on 1km OSHD grid (BAFU_DEM_2020_1000.txt)');
        ncwriteatt(nc_file,'/','temp. resolution','1-hourly, based on COSMO1E analysis data');
    end
    %% now individual data -- precip
    nccreate(nc_precip,'PRECTmms','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_precip, 'PRECTmms', 'long_name', 'PRECTmms total precipitation');
    ncwriteatt(nc_precip, 'PRECTmms', 'units', 'mm H2O / sec');
    ncwriteatt(nc_precip, 'PRECTmms', 'mode', 'time-dependent');
    ncwriteatt(nc_precip, 'PRECTmms', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_precip, 'PRECTmms', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_precip,'PRECTmms',prec_oshd);
    %% now individual data -- solar
    nccreate(nc_solar,'FSDS','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_solar, 'FSDS', 'long_name', 'total incident solar radiation');
    ncwriteatt(nc_solar, 'FSDS', 'units', 'W/m**2');
    ncwriteatt(nc_solar, 'FSDS', 'mode', 'time-dependent');
    ncwriteatt(nc_solar, 'FSDS', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_solar, 'FSDS', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_solar,'FSDS',swr_oshd);
    %% now individual data -- tair_rh
    nccreate(nc_tq,'TBOT','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_tq, 'TBOT', 'long_name', 'temperature at the lowest atm level');
    ncwriteatt(nc_tq, 'TBOT', 'units', 'K');
    ncwriteatt(nc_tq, 'TBOT', 'mode', 'time-dependent');
    ncwriteatt(nc_tq, 'TBOT', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_tq, 'TBOT', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_tq,'TBOT',tair_oshd);

    nccreate(nc_tq,'RH','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_tq, 'RH', 'long_name', 'relative humidity at the lowest atm level (RH)');
    ncwriteatt(nc_tq, 'RH', 'units', '%');
    ncwriteatt(nc_tq, 'RH', 'mode', 'time-dependent');
    ncwriteatt(nc_tq, 'RH', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_tq, 'RH', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_tq,'RH',rh_oshd);

    %% now individual data -- wind_zbot
    nccreate(nc_wz,'WIND','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_wz, 'WIND', 'long_name', 'wind at the lowest atm level');
    ncwriteatt(nc_wz, 'WIND', 'units', 'm/s');
    ncwriteatt(nc_wz, 'WIND', 'mode', 'time-dependent');
    ncwriteatt(nc_wz, 'WIND', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_wz, 'WIND', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_wz,'WIND',wind_oshd);

    z_meas=repmat(10,365,272); %m, observational height
    nccreate(nc_wz,'ZBOT','Dimensions',{'lon',365,'lat',272},'Format','classic','Datatype','double');
    ncwriteatt(nc_wz, 'ZBOT', 'long_name', 'observational height');
    ncwriteatt(nc_wz, 'ZBOT', 'units', 'm');
    ncwriteatt(nc_wz, 'ZBOT', 'mode', 'time-independent');
    ncwrite(nc_wz,'ZBOT',z_meas);
    %% now individual data -- lwr_pa
    nccreate(nc_lp,'FLDS','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_lp, 'FLDS', 'long_name', 'incident longwave radiation');
    ncwriteatt(nc_lp, 'FLDS', 'units', 'W/m**2');
    ncwriteatt(nc_lp, 'FLDS', 'mode', 'time-dependent');
    ncwriteatt(nc_lp, 'FLDS', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_lp, 'FLDS', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_lp,'FLDS',lwr_oshd);

    nccreate(nc_lp,'PSRF','Dimensions',{'lon',365,'lat',272,'time',length_time},'Format','classic','Datatype','double');
    ncwriteatt(nc_lp, 'PSRF', 'long_name', 'surface pressure at the lowest atm level');
    ncwriteatt(nc_lp, 'PSRF', 'units', 'Pa');
    ncwriteatt(nc_lp, 'PSRF', 'mode', 'time-dependent');
    ncwriteatt(nc_lp, 'PSRF', '_FillValue',9.999999616903162e+35);
    ncwriteatt(nc_lp, 'PSRF', 'missing_value',9.999999616903162e+35);
    ncwrite(nc_lp,'PSRF',pair_oshd);

end

end
