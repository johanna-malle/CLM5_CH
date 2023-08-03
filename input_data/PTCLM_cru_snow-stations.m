function regrid_crujura_PTCLM_split
%{
Script to regrid precip data from global crujura V2.2. dataset for PTCLM5
simulations. Write out data with and without lapse rate.

need to specify swiss domain, and where input output data is. Loops through
all files, all time steps and creates new CLM5 input data sets.


%}


bf_out_raw='/media/malle/LaCie/CLM5_input/PTCLM/met_input/';
list={'lwr_pa','precip','solar','tair_rh','wind_zbot','tair_rh_lapse'};
% list={'tair_rh_lapse'};

%bf='/media/malle/LaCie/CLM5_input/surfdata_new/PTCLM5/domain';
bf='/media/malle/LaCie/CLM5_input/surfdata_new/PTCLM5/domain_new';
stats=dir(bf);
tf=ismember({stats.name},{'.','..','.3'});
stats(tf)=[];
tf1=contains({stats.name},{'.ocn.'});
stats(tf1)=[];
domains_stations={stats.name}';

%%
%% topo forcing to calculate temperature at sea level
%for input data
topo_forcing=('/media/malle/LaCie/CLM5_input/cru7/topodata_360x720_c060528.nc');
TOPO=ncread(topo_forcing,'TOPO');
topo_lat=ncread(topo_forcing,'LATIXY');
topo_long=ncread(topo_forcing,'LONGXY');

%need topo meshgrid
% shift both old and new grid from center of cell to lower left corner... needed for interp2
topo_long_shift=topo_long-0.25; %shift because we want lower left corner id for interp2
topo_lat_shift=topo_lat-0.25; %shift because we want lower left corner id for interp2
% make meshgrids for interpolation
[X_llc_topo,Y_llc_topo]=meshgrid(topo_long_shift(:,1),topo_lat_shift(1,:)');

stat_list=load('/media/malle/LaCie/OSHD_info/STAT_LIST.mat');
acros= stat_list.statlist.acro;
elev=stat_list.statlist.z; %grid cells outside of domain => just leave it at the sealevel temperature
%% loop through stations
for domix=1:size(domains_stations,1)
    dom_data=domains_stations{domix};
    if strcmp(dom_data(15),'_')==1
        loc=dom_data(12:14);
    elseif strcmp(dom_data(16),'_')==1
        loc=dom_data(12:15);
    else
        error('station does not exist')
    end
    %get elev
    id_stat=strcmp(loc,acros);
    if sum(id_stat)==0 && strcmp(loc,'DIS')==1
        id_stat=strcmp(['#' loc],acros);
    elseif sum(id_stat)==0
        id_stat=strcmp(['*' loc],acros); %apart from DIS, all other selected stations have a * infront... needed here to avoid confusion between stations at different elevations for various locations
    end

    if sum(id_stat)~=0
        elev_stat=elev(id_stat);
        %new grid (vertices)
        LONGXY_new=ncread(fullfile(bf,dom_data),'xc');
        LATIXY_new=ncread(fullfile(bf,dom_data),'yc');
        LONGXY_vert_new=ncread(fullfile(bf,dom_data),'xv');
        LATIXY_vert_new=ncread(fullfile(bf,dom_data),'yv');
        edge_w=min(min(min(LONGXY_vert_new)));
        edge_e=max(max(max(LONGXY_vert_new)));
        edge_s=min(min(min(LATIXY_vert_new)));
        edge_n=max(max(max(LATIXY_vert_new)));

        %need lower left of new grid set
        long_new_llc=squeeze(LONGXY_vert_new(1,:,:));
        lat_new_llc=squeeze(LATIXY_vert_new(1,:,:));

        %% paths
        out_folder1=fullfile(bf_out_raw,loc);
        if exist(out_folder1,'dir')~=7 %exist always returns 7 if it is a folder
            mkdir(out_folder1)
        end

        out_folder2=fullfile(bf_out_raw,loc,'cru_jra');
        if exist(out_folder2,'dir')~=7 %exist always returns 7 if it is a folder
            mkdir(out_folder2)
        end

        out_folder_overall=fullfile(bf_out_raw,loc,'cru_jra','split');
        if exist(out_folder_overall,'dir')~=7 %exist always returns 7 if it is a folder
            mkdir(out_folder_overall)
        end
        %%
        for oix=1:size(list,2)
            var=list{oix};
            %% paths
            if strcmp(var,'tair_rh_lapse')==1
                bf_in=fullfile('/media/malle/LaCie/CLM5_input/cru_jra_CLM5/','tair_rh');
            else
                bf_in=fullfile('/media/malle/LaCie/CLM5_input/cru_jra_CLM5/',var);
            end
            bf_out=fullfile(out_folder_overall,var);
            if exist(bf_out,'dir')~=7 %exist always returns 7 if it is a folder
                mkdir(bf_out)
            end
            %% loop through all files
            % old domain
            dir_in=dir(bf_in);%what files are in there?
            tf = ismember( {dir_in.name}, {'.', '..'});
            dir_in(tf) = []; %remove current and parent directory
            names_data={dir_in.name}';
            for six=1:size(names_data,1)
                name_six=cell2mat(names_data(six));
                file_world=fullfile(bf_in,(name_six));
                file_ch_all=fullfile(bf_out,(name_six));
                if exist(file_ch_all,'file')==2
                    delete(file_ch_all)
                end
                %% get needed variables for writing new file
                time=ncread(file_world,'time'); %for writing new file later
                time_att=ncreadatt(file_world,'time','units'); %for writing new file later
                length_time=size(time,1);
                %% grids
                %old grid
                LONGXY=ncread(file_world,'LONGXY');
                LATIXY=ncread(file_world,'LATIXY');
                %% shift both old and new grid from center of cell to lower left corner... needed for interp2
                LONGXY_shift=LONGXY-0.25; %shift because we want lower left corner id for interp2
                LATIXY_shift=LATIXY-0.25; %shift because we want lower left corner id for interp2
                %% make meshgrids for interpolation
                [X_llc,Y_llc]=meshgrid(LONGXY_shift(:,1),LATIXY_shift(1,:)');
                [X_new_llc,Y_new_llc]=meshgrid(long_new_llc(:,1),lat_new_llc(1,:)');
                %% write dataset - things that are the same
                nccreate(file_ch_all,'time','Dimensions',{'time',length_time},'Format','classic','Datatype','double');
                ncwriteatt(file_ch_all, 'time', 'long_name', 'observation time');
                ncwriteatt(file_ch_all, 'time', 'units',time_att);
                ncwriteatt(file_ch_all, 'time', 'calendar', 'noleap');
                ncwrite(file_ch_all,'time',time);

                %long+lat
                nccreate(file_ch_all,'LONGXY','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2)},'Format','classic','Datatype','double')
                ncwriteatt(file_ch_all, 'LONGXY', 'long_name', 'longitude');
                ncwriteatt(file_ch_all, 'LONGXY', 'units', 'degrees_east');
                ncwriteatt(file_ch_all, 'LONGXY', 'mode', 'time-invariant');
                ncwrite(file_ch_all,'LONGXY',LONGXY_new);

                nccreate(file_ch_all,'LATIXY','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2)},'Format','classic','Datatype','double');
                ncwriteatt(file_ch_all, 'LATIXY', 'long_name', 'latitude');
                ncwriteatt(file_ch_all, 'LATIXY', 'units', 'degrees_north');
                ncwriteatt(file_ch_all, 'LATIXY', 'mode', 'time-invariant');
                ncwrite(file_ch_all,'LATIXY',LATIXY_new);

                %4 edges of grid
                nccreate(file_ch_all,'EDGEN','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(file_ch_all, 'EDGEN', 'long_name', 'northern edge in atmospheric data');
                ncwriteatt(file_ch_all, 'EDGEN', 'units', 'degrees_north');
                ncwriteatt(file_ch_all, 'EDGEN', 'mode', 'time-invariant');
                ncwrite(file_ch_all,'EDGEN',edge_n);

                nccreate(file_ch_all,'EDGES','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(file_ch_all, 'EDGES', 'long_name', 'southern edge in atmospheric data');
                ncwriteatt(file_ch_all, 'EDGES', 'units', 'degrees_north');
                ncwriteatt(file_ch_all, 'EDGES', 'mode', 'time-invariant');
                ncwrite(file_ch_all,'EDGES',edge_s);

                nccreate(file_ch_all,'EDGEE','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(file_ch_all, 'EDGEE', 'long_name', 'eastern edge in atmospheric data');
                ncwriteatt(file_ch_all, 'EDGEE', 'units', 'degrees_east');
                ncwriteatt(file_ch_all, 'EDGEE', 'mode', 'time-invariant');
                ncwrite(file_ch_all,'EDGEE',edge_e);

                nccreate(file_ch_all,'EDGEW','Dimensions',{'scalar',1},'Format','classic','Datatype','double');
                ncwriteatt(file_ch_all, 'EDGEW', 'long_name', 'western edge in atmospheric data');
                ncwriteatt(file_ch_all, 'EDGEW', 'units', 'degrees_east');
                ncwriteatt(file_ch_all, 'EDGEW', 'mode', 'time-invariant');
                ncwrite(file_ch_all,'EDGEW',edge_w);

                %global attributes for entire .nc file
                ncwriteatt(file_ch_all,'/','institution','WSL Institute for snow and avalanche research SLF, snow hydrology group');
                ncwriteatt(file_ch_all,'/','history',['File Origin - This file was created by J.Malle on ' datestr(now)]);
                ncwriteatt(file_ch_all,'/','grid location','Switzerland, interpolated from global cru jura V2.2 data');
                ncwriteatt(file_ch_all,'/','temp. resolution','6-hourly, based on 0.5deg cru jura V2.2 data');

                %% now actually write specific data
                if strcmp(var,'lwr_pa')==1
                    PSRF=(ncread(file_world,'PSRF'));
                    FLDS=(ncread(file_world,'FLDS'));
                    sz_t=size(PSRF,3); %how many time samples?

                    PSRF_regrid=nan([size(LONGXY_new) sz_t]); %output file
                    FLDS_regrid=nan([size(LONGXY_new) sz_t]); %output file
                    % interpolate variables\
                    for miix=1:sz_t
                        PSRF_inter=interp2(X_llc,Y_llc,(PSRF(:,:,miix))',X_new_llc,Y_new_llc);
                        FLDS_inter=interp2(X_llc,Y_llc,(FLDS(:,:,miix))',X_new_llc,Y_new_llc);
                        PSRF_regrid(:,:,miix)=PSRF_inter';
                        FLDS_regrid(:,:,miix)=FLDS_inter';
                    end
                    nccreate(file_ch_all,'PSRF','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all,'PSRF', 'long_name',ncreadatt(file_world,'PSRF','long_name'));
                    ncwriteatt(file_ch_all,'PSRF', 'units', ncreadatt(file_world,'PSRF','units'));
                    ncwriteatt(file_ch_all,'PSRF', 'mode',ncreadatt(file_world,'PSRF','mode'));
                    ncwriteatt(file_ch_all,'PSRF', '_FillValue',ncreadatt(file_world,'PSRF','_FillValue'));
                    ncwriteatt(file_ch_all,'PSRF', 'missing_value',ncreadatt(file_world,'PSRF','missing_value'));
                    ncwrite(file_ch_all,'PSRF',PSRF_regrid);

                    nccreate(file_ch_all,'FLDS','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all,'FLDS', 'long_name',ncreadatt(file_world,'FLDS','long_name'));
                    ncwriteatt(file_ch_all,'FLDS', 'units', ncreadatt(file_world,'FLDS','units'));
                    ncwriteatt(file_ch_all,'FLDS', 'mode',ncreadatt(file_world,'FLDS','mode'));
                    ncwriteatt(file_ch_all,'FLDS', '_FillValue',ncreadatt(file_world,'FLDS','_FillValue'));
                    ncwriteatt(file_ch_all,'FLDS', 'missing_value',ncreadatt(file_world,'FLDS','missing_value'));
                    ncwrite(file_ch_all,'FLDS',FLDS_regrid);

                    %%
                elseif strcmp(var,'precip')==1
                    PRECIP=(ncread(file_world,'PRECTmms'));
                    PRECIP_regrid=nan([size(LONGXY_new) length_time]); %output file
                    for miix=1:size(PRECIP,3)
                        PRECIP_inter=interp2(X_llc,Y_llc,(PRECIP(:,:,miix))',X_new_llc,Y_new_llc); %potentially use 'nearest' for precipitation?
                        PRECIP_regrid(:,:,miix)=PRECIP_inter';
                    end

                    nccreate(file_ch_all,'PRECTmms','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all, 'PRECTmms', 'long_name', 'PRECTmms total precipitation');
                    ncwriteatt(file_ch_all, 'PRECTmms', 'units', 'mm H2O / sec');
                    ncwriteatt(file_ch_all, 'PRECTmms', 'mode', 'time-dependent');
                    ncwriteatt(file_ch_all, 'PRECTmms', '_FillValue',9.999999616903162e+35);
                    ncwriteatt(file_ch_all, 'PRECTmms', 'missing_value',9.999999616903162e+35);
                    ncwrite(file_ch_all,'PRECTmms',PRECIP_regrid);

                    %%
                elseif strcmp(var,'solar')==1
                    FSDS=(ncread(file_world,'FSDS'));
                    FSDS_regrid=nan([size(LONGXY_new) length_time]); %output file
                    for miix=1:size(FSDS,3)
                        FSDS_inter=interp2(X_llc,Y_llc,(FSDS(:,:,miix))',X_new_llc,Y_new_llc); %potentially use 'nearest' for precipitation?
                        FSDS_regrid(:,:,miix)=FSDS_inter';
                    end

                    nccreate(file_ch_all,'FSDS','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all, 'FSDS', 'long_name', 'total incident solar radiation');
                    ncwriteatt(file_ch_all, 'FSDS', 'units', 'W/m**2');
                    ncwriteatt(file_ch_all, 'FSDS', 'mode', 'time-dependent');
                    ncwriteatt(file_ch_all, 'FSDS', '_FillValue',9.999999616903162e+35);
                    ncwriteatt(file_ch_all, 'FSDS', 'missing_value',9.999999616903162e+35);
                    ncwrite(file_ch_all,'FSDS',FSDS_regrid);
                    %%
                elseif strcmp(var,'tair_rh')==1 || strcmp(var,'tair_rh_lapse')==1
                    QBOT=(ncread(file_world,'QBOT'));
                    TBOT=(ncread(file_world,'TBOT'));
                    sz_t=size(TBOT,3); %how many time samples?
                    TBOT_regrid=nan([size(LONGXY_new) sz_t]); %output file
                    TBOT_regrid_lapse=nan([size(LONGXY_new) sz_t]); %output file
                    QBOT_regrid=nan([size(LONGXY_new) sz_t]); %output file
                    for miix=1:sz_t
                        TBOT_inter=interp2(X_llc,Y_llc,(TBOT(:,:,miix))',X_new_llc,Y_new_llc);
                        TBOT_regrid(:,:,miix)=TBOT_inter';
                        %deal with lapse rate
                        TOPO_cru7=interp2(X_llc_topo,Y_llc_topo,TOPO',X_new_llc,Y_new_llc); %interpolate topo data to station
                        temp_interp=interp2(X_llc,Y_llc,TBOT(:,:,miix)',X_new_llc,Y_new_llc); %interpolate temp data to station
                        t_sealevel=temp_interp+(6.5/1000)*flip(TOPO_cru7); %global topo to bring temperature down to sealevel
                        T_lapse_rate=t_sealevel-(6.5/1000)*elev_stat;         %use DEM of Switzerland to "relapse"
                        TBOT_regrid_lapse(:,:,miix)=T_lapse_rate';

                        QBOT_inter=interp2(X_llc,Y_llc,(QBOT(:,:,miix))',X_new_llc,Y_new_llc);
                        QBOT_regrid(:,:,miix)=QBOT_inter';
                    end


                    nccreate(file_ch_all,'QBOT','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all, 'QBOT', 'long_name', 'specific humidity at the lowest atm level');
                    ncwriteatt(file_ch_all, 'QBOT', 'units', 'kg/kg');
                    ncwriteatt(file_ch_all, 'QBOT', 'mode', 'time-dependent');
                    ncwriteatt(file_ch_all, 'QBOT', '_FillValue',9.999999616903162e+35);
                    ncwriteatt(file_ch_all, 'QBOT', 'missing_value',9.999999616903162e+35);
                    ncwrite(file_ch_all,'QBOT',QBOT_regrid);

                    nccreate(file_ch_all,'TBOT','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all, 'TBOT', 'long_name', 'temperature at the lowest atm level');
                    ncwriteatt(file_ch_all, 'TBOT', 'units', 'K');
                    ncwriteatt(file_ch_all, 'TBOT', 'mode', 'time-dependent');
                    ncwriteatt(file_ch_all, 'TBOT', '_FillValue',9.999999616903162e+35);
                    ncwriteatt(file_ch_all, 'TBOT', 'missing_value',9.999999616903162e+35);
                    if strcmp(var,'tair_rh')==1
                        ncwrite(file_ch_all,'TBOT',TBOT_regrid);
                    else
                        ncwrite(file_ch_all,'TBOT',TBOT_regrid_lapse);
                    end

                    %%
                elseif strcmp(var,'wind_zbot')==1
                    Wind=(ncread(file_world,'WIND'));
                    sz_t=size(Wind,3); %how many time samples?
                    Wind_regrid=nan([size(LONGXY_new) sz_t]); %output file
                    for miix=1:sz_t
                        Wind_inter=interp2(X_llc,Y_llc,(Wind(:,:,miix))',X_new_llc,Y_new_llc);
                        Wind_regrid(:,:,miix)=Wind_inter';
                    end

                    nccreate(file_ch_all,'WIND','Dimensions',{'lon',size(LONGXY_new,1),'lat',size(LONGXY_new,2),'time',length_time},'Format','classic','Datatype','double');
                    ncwriteatt(file_ch_all, 'WIND', 'long_name', 'wind at the lowest atm level');
                    ncwriteatt(file_ch_all, 'WIND', 'units', 'm/s');
                    ncwriteatt(file_ch_all, 'WIND', 'mode', 'time-dependent');
                    ncwriteatt(file_ch_all, 'WIND', '_FillValue',9.999999616903162e+35);
                    ncwriteatt(file_ch_all, 'WIND', 'missing_value',9.999999616903162e+35);
                    ncwrite(file_ch_all,'WIND',Wind_regrid);

                end
            end
        end
    end
end

end
