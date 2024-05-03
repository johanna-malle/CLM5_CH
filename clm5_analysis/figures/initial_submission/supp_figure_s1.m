function make_wiggle_plot_new_data_v2
close all
clearvars
%{
function to aggregate data from PTCLM5 simulations for wiggle plot and
create wiggle plot for a) crujra based CLM5 runs b) OSHD based CLM5 runs c)
JIM-based simulations.
%}

addpath(genpath('/home/malle/Documents/MATLAB/'));

%% define paths/already load some data
bf_out='/home/lud11/malle/CLM5_CH'; %where to save plots/mat files
%% folders of clm5 runs:
bf_clm5='/media/malle/LaCie/PTCLM5_analysis_final/';
clm5_oshd_newsurf=fullfile(bf_clm5,'PTCLM_all_OSHD_newSurf');
clm5_oshd_origsurf=fullfile(bf_clm5,'PTCLM_all_OSHD_origSurf');

clm5_crujraL_newsurf=fullfile(bf_clm5,'PTCLM_all_CRUJRA_lapse_newSurf');
clm5_crujraL_origsurf=fullfile(bf_clm5,'PTCLM_all_CRUJRA_lapse_origSurf');

clm5_crujraNoL_newsurf=fullfile(bf_clm5,'PTCLM_all_CRUJRA_Nolapse_newSurf');
clm5_crujraNoL_origsurf=fullfile(bf_clm5,'PTCLM_all_CRUJRA_Nolapse_origSurf');
%%
bf_points=('/home/lud11/malle/CLM5_CH/FSM_new/analysed_points')
dir_in=dir(bf_points);
tf=ismember({dir_in.name},{'.','..'});
dir_in(tf)=[];
files_name={dir_in.name}';
locs1=cellfun(@(x) strsplit(x,'_'),files_name, 'UniformOutput',false)'; %all locations
locs=cellfun(@(x) x{2},locs1, 'UniformOutput',false); %all locations
elev=cellfun(@(x) x{4}(1:end-4),locs1, 'UniformOutput',false); %all locations
sz_locs=size(locs,2);
%%
bf_meas=('/home/lud11/malle/CLM5_CH/dvd_oshd');
dir_in=dir(bf_meas);
tf=ismember({dir_in.name},{'.','..','HS_'});
dir_in(tf)=[];
files_name_meas={dir_in.name}';
%% initialize matrices and loop through all locations of interest
mat_all_oshd_newSurf=nan(365,sz_locs,6);
mat_all_crujra_newSurf=nan(365,sz_locs,11);
mat_all_crujra_noLapse_newSurf=nan(365,sz_locs,11);
mat_all_oshd_origSurf=nan(365,sz_locs,6);
mat_all_crujra_origSurf=nan(365,sz_locs,11);
mat_all_crujra_noLapse_origSurf=nan(365,sz_locs,11);
mat_all_jim=nan(365,sz_locs,6);
mat_all_measured=nan(365,sz_locs,3);

elev_inter=nan(sz_locs,1);
time_int=cell(sz_locs,1);
comp_HS_all=cell(sz_locs,1);

for tix=1:sz_locs %loop through all locations of interest
    name_sta=locs{tix}; %location of interest
    display(name_sta)
    elev_inter(tix)=str2num(elev{tix});
    jim_in=table2timetable(readtable(fullfile(bf_points,files_name{tix})));
    meas_in=table2timetable(readtable(fullfile(bf_meas,files_name_meas{find(contains(files_name_meas,name_sta))})));
    jim_in_shift=retime(jim_in,'daily','previous');

    oshd_newsurf=load(fullfile(clm5_oshd_newsurf,[name_sta '_output_2014_1_1.mat']));
    crujra_newsurf=load(fullfile(clm5_crujraL_newsurf,[name_sta '_output_2010_1_1.mat']));
    crujra_nolapse_newsurf=load(fullfile(clm5_crujraNoL_newsurf,[name_sta '_output_2010_1_1.mat']));

    oshd_origsurf=load(fullfile(clm5_oshd_origsurf,[name_sta '_output_2014_1_1.mat']));
    crujra_origsurf=load(fullfile(clm5_crujraL_origsurf,[name_sta '_output_2010_1_1.mat']));
    crujra_nolapse_origsurf=load(fullfile(clm5_crujraNoL_origsurf,[name_sta '_output_2010_1_1.mat']));

    %% compare to CLM5 output
    %timetable(dateshift(jim_in.time_stamps_jim,'start','day'),jim_in{:,1:3})

    id_hs=strcmp(oshd_origsurf.var_combo,'SNOW_DEPTH');
    id_hs_crujra=strcmp(crujra_origsurf.var_combo,'SNOW_DEPTH');

    id_swe=strcmp(oshd_origsurf.var_combo,'H2OSNO');
    id_swe_crujra=strcmp(crujra_origsurf.var_combo,'H2OSNO');

    id_scf=strcmp(oshd_origsurf.var_combo,'FSNO');
    id_scf_crujra=strcmp(crujra_origsurf.var_combo,'FSNO');

    timeyears_round=dateshift(datetime(datevec(oshd_origsurf.timeyears)), 'start', 'minute', 'nearest');
    TT_ptclm1=timetable(timeyears_round,oshd_origsurf.data_combo{id_hs}',...
        oshd_origsurf.data_combo{id_swe}',oshd_origsurf.data_combo{id_scf}',...
        oshd_newsurf.data_combo{id_hs}',oshd_newsurf.data_combo{id_swe}',oshd_newsurf.data_combo{id_scf}',...
        'VariableNames',{'HS_CLM5_OSHD_origSurf','SWE_CLM5_OSHD_origSurf','SCF_CLM5_OSHD_origSurf',...
        'HS_CLM5_OSHD_newSurf','SWE_CLM5_OSHD_newSurf','SCF_CLM5_OSHD_newSurf'});

    timeyears_round_crujra=dateshift(datetime(datevec(crujra_origsurf.timeyears)), 'start', 'minute', 'nearest');
    TT_ptclm1_crujra=timetable(timeyears_round_crujra,crujra_origsurf.data_combo{id_hs_crujra}',...
        crujra_origsurf.data_combo{id_swe_crujra}',crujra_origsurf.data_combo{id_scf_crujra}',...
        crujra_newsurf.data_combo{id_hs_crujra}',crujra_newsurf.data_combo{id_swe_crujra}',crujra_newsurf.data_combo{id_scf_crujra}',...
        'VariableNames',{'HS_CLM5_crujra_origSurf','SWE_CLM5_crujra_origSurf','SCF_CLM5_crujra_origSurf',...
        'HS_CLM5_crujra_newSurf','SWE_CLM5_crujra_newSurf','SCF_CLM5_crujra_newSurf'});

    TT_ptclm1_crujra_nolapse=timetable(timeyears_round_crujra,crujra_nolapse_origsurf.data_combo{id_hs_crujra}',...
        crujra_nolapse_origsurf.data_combo{id_swe_crujra}',crujra_nolapse_origsurf.data_combo{id_scf_crujra}',...
        crujra_nolapse_newsurf.data_combo{id_hs_crujra}',...
        crujra_nolapse_newsurf.data_combo{id_swe_crujra}',crujra_nolapse_newsurf.data_combo{id_scf_crujra}',...
        'VariableNames',{'HS_CLM5_crujra_nolapse_origSurf','SWE_CLM5_crujra_nolapse_origSurf','SCF_CLM5_crujra_nolapse_origSurf',...
        'HS_CLM5_crujra_nolapse_newSurf','SWE_CLM5_crujra_nolapse_newSurf','SCF_CLM5_crujra_nolapse_newSurf'});

    tt_all=synchronize(meas_in,TT_ptclm1,jim_in_shift);
    tt_all.diff_all_newSurf=(tt_all.HS_meas-tt_all.HS_CLM5_OSHD_newSurf);
    tt_all.diff_all_origSurf=(tt_all.HS_meas-tt_all.HS_CLM5_OSHD_origSurf);
    tt_all.diff_jim=(tt_all.HS_meas-tt_all.HS_jim);
    daily_all=retime(tt_all,'daily','firstvalue');

    tt_all_crujra=synchronize(meas_in,TT_ptclm1_crujra,TT_ptclm1_crujra_nolapse);
    tt_all_crujra.diff_all_newSurf=(tt_all_crujra.HS_meas-tt_all_crujra.HS_CLM5_crujra_newSurf);
    tt_all_crujra.diff_all_origSurf=(tt_all_crujra.HS_meas-tt_all_crujra.HS_CLM5_crujra_origSurf);
    tt_all_crujra.diff_all_noLapse_newSurf=(tt_all_crujra.HS_meas-tt_all_crujra.HS_CLM5_crujra_nolapse_newSurf);
    tt_all_crujra.diff_all_noLapse_origSurf=(tt_all_crujra.HS_meas-tt_all_crujra.HS_CLM5_crujra_nolapse_origSurf);
    daily_all_crujra=retime(tt_all_crujra,'daily','firstvalue');

    years=categorical(year(daily_all.time_HS_meas));
    time_overall=daily_all.time_HS_meas;
    diff_overall=daily_all.diff_all_newSurf;
    diff_overall_orig=daily_all.diff_all_origSurf;
    diff_jim=daily_all.diff_jim;
    years_u=unique(years);
    HS_meas_overall=daily_all.HS_meas;

    %% diff plot
    time_all=tt_all.time_HS_meas;
    time_all_cru=tt_all_crujra.time_HS_meas;

    fh=figure('Position',[963 339 521 619]);hold on;grid on
    ah1=subplot(2,1,1);hold on;grid minor;grid on
    ah2=subplot(2,1,2);hold on;grid minor;grid on
    plot(ah1,time_all,tt_all.HS_jim-tt_all.HS_meas,'Color',rgb('dark gray'),...
        'LineStyle','-','LineWidth',1.5,'DisplayName','FSM2')
    plot(ah1,time_all,tt_all.HS_CLM5_OSHD_newSurf-tt_all.HS_meas,'Color',rgb('orange'),...
        'LineWidth',1.5,'DisplayName','Clim_{OSHD 1km}LU_{HR 1km}')
    plot(ah1,time_all,tt_all.HS_CLM5_OSHD_origSurf-tt_all.HS_meas,'Color',rgb('orange'),...
        'LineWidth',1.5,'DisplayName','Clim_{OSHD 1km}LU_{Gl 1km}','LineStyle',':')
    plot(ah1,time_all_cru,tt_all_crujra.HS_CLM5_crujra_newSurf-tt_all_crujra.HS_meas,'Color',rgb('grass'),...
        'LineWidth',1.5,'DisplayName','Clim_{CRU* 1km}LU_{HR 1km}')
    plot(ah1,time_all_cru,tt_all_crujra.HS_CLM5_crujra_origSurf-tt_all_crujra.HS_meas,'Color',rgb('grass'),...
        'LineWidth',1.5,'DisplayName','Clim_{CRU* 1km}LU_{Gl 1km}','LineStyle',':')
    plot(ah1,time_all_cru,tt_all_crujra.HS_CLM5_crujra_nolapse_newSurf-tt_all_crujra.HS_meas,'Color',rgb('cornflower'),...
        'LineWidth',1.5,'DisplayName','Clim_{CRU 1km}LU_{HR 1km}')
    plot(ah1,time_all_cru,tt_all_crujra.HS_CLM5_crujra_nolapse_origSurf-tt_all_crujra.HS_meas,'Color',rgb('cornflower'),...
        'LineWidth',1.5,'DisplayName','Clim_{CRU 1km}LU_{Gl 1km}','LineStyle',':')
    plot(ah1,time_all(1),0,'Color',rgb('scarlet'),...
        'LineStyle','-','LineWidth',1.6,'DisplayName','OBSERVED','LineStyle','-.')
    title(ah1,[name_sta ' - ' num2str(elev_inter(tix)) 'm'],'FontWeight','Normal')
    xlim(ah1,[datetime(2017,10,1) datetime(2018,7,1)])
    ylabel(ah1,'\Delta HS [m]')

    plot(ah2,time_all,tt_all.HS_meas,'Color',rgb('scarlet'),...
        'LineStyle','-','LineWidth',1.6,'DisplayName','OBSERVED','LineStyle','-.')
    plot(ah2,time_all,tt_all.HS_jim,'Color',rgb('dark gray'),...
        'LineStyle','-','LineWidth',1.5,'DisplayName','FSM2')
    plot(ah2,time_all,tt_all.HS_CLM5_OSHD_newSurf,'Color',rgb('orange'),...
        'LineWidth',1.5,'DisplayName','OSHD HighRes')
    plot(ah2,time_all_cru,tt_all_crujra.HS_CLM5_crujra_newSurf,'Color',rgb('grass'),...
        'LineWidth',1.5,'DisplayName','CRUJRA+ HighRes')
    plot(ah2,time_all_cru,tt_all_crujra.HS_CLM5_crujra_nolapse_newSurf,'Color',rgb('cornflower'),...
        'LineWidth',1.5,'DisplayName','CRUJRA HighRes')

    pos_2=ah2.Position;
    ah2.Position=[pos_2(1) pos_2(2)-0.05 pos_2(3:4)];

    le2=legend(ah2,'Location','southoutside','Orientation','horizontal');
    set(le2,'visible','off')
    xlim(ah2,[datetime(2017,10,1) datetime(2018,7,1)])
    ylabel(ah2,'HS [m]')
    le1=legend(ah1,'Location','southoutside','Orientation','vertical');
    set(gcf,'Color','white')
    le1.Position=[0.0399 0.4315 0.9360 0.0964];
    le1.NumColumns=3;
    set(findall(gcf,'-property','FontSize'),'FontSize',11)
    le1.FontSize=8;

    save_name=fullfile(bf_out,'snow_station_comp',['HS_diff_' name_sta '.png']);
    save_name2=fullfile(bf_out,'snow_station_comp',['HS_diff' name_sta '.pdf']);
    exportgraphics(fh,save_name,'Resolution',700)
    exportgraphics(fh,save_name2)
    close(fh)

    %%
    timing_check=cell(size(years_u,1)-2,1);
    for siix=1:size(years_u,1)-2
        year_int=double(string(years_u(siix)));
        year_end=double(string(years_u(siix+1)));
        if leapyear(year_end)
            idss=time_overall>datetime(year_int,8,31) & time_overall<datetime(year_end,9,1) & time_overall~=datetime(year_end,2,29);
        else
            idss=time_overall>datetime(year_int,8,31) & time_overall<datetime(year_end,9,1);
        end
        timing_check{siix}=time_overall(idss);
        mat_all_oshd_newSurf(:,tix,siix)=diff_overall(idss);
        mat_all_oshd_origSurf(:,tix,siix)=diff_overall_orig(idss);
        mat_all_jim(:,tix,siix)=diff_jim(idss);
        mat_all_measured(:,tix,siix)=HS_meas_overall(idss);
    end

    years_crujra=categorical(year(daily_all_crujra.time_HS_meas));
    time_overall=daily_all_crujra.time_HS_meas;
    diff_overall=daily_all_crujra.diff_all_newSurf;
    diff_overall_origSurf=daily_all_crujra.diff_all_origSurf;
    diff_overall_nolapse=daily_all_crujra.diff_all_noLapse_newSurf;
    diff_overall_nolapse_origSurf=daily_all_crujra.diff_all_noLapse_origSurf;
    years_u_crujra=unique(years_crujra);

    timing_check_crujra=cell(size(years_u_crujra,1)-2,1);
    for siix=1:size(years_u_crujra,1)-2
        year_int=double(string(years_u_crujra(siix)));
        year_end=double(string(years_u_crujra(siix+1)));
        if leapyear(year_end)
            idss=time_overall>datetime(year_int,8,31) & time_overall<datetime(year_end,9,1) & time_overall~=datetime(year_end,2,29);
        else
            idss=time_overall>datetime(year_int,8,31) & time_overall<datetime(year_end,9,1);
        end
        timing_check_crujra{siix}=time_overall(idss);
        mat_all_crujra_newSurf(:,tix,siix)=diff_overall(idss);
        mat_all_crujra_origSurf(:,tix,siix)=diff_overall_origSurf(idss);
        mat_all_crujra_noLapse_newSurf(:,tix,siix)=diff_overall_nolapse(idss);
        mat_all_crujra_noLapse_origSurf(:,tix,siix)=diff_overall_nolapse_origSurf(idss);
    end
end

%%
save(fullfile(bf_out,'mat_for_wiggle.mat'),'timing_check','timing_check_crujra',...
    'elev_inter','time_int','comp_HS_all','mat_all_crujra_newSurf','mat_all_crujra_origSurf',...
    'mat_all_crujra_noLapse_newSurf','mat_all_crujra_noLapse_origSurf','mat_all_oshd_newSurf',...
    'mat_all_oshd_origSurf','mat_all_jim','daily_all','tt_all','tt_all_crujra','daily_all_crujra','locs');
%load(fullfile(bf_out,'mat_for_wiggle.mat'));
%% wiggle plot
[~,id_sort] = sort(elev_inter);
for jiix=2:4%size(mat_all_jim,3) %2014 not actually proper data
    %%
    timing_check1=timing_check{jiix};
    timing_check_crujra1=timing_check_crujra{jiix+4};

    mat_plot=mat_all_jim(:,id_sort,jiix);
    mat_plot_crujra_newSurf=mat_all_crujra_newSurf(:,id_sort,jiix+4);
    mat_plot_oshd_newSurf=mat_all_oshd_newSurf(:,id_sort,jiix);
    mat_plot_crujra_origSurf=mat_all_crujra_origSurf(:,id_sort,jiix+4);
    mat_plot_oshd_origSurf=mat_all_oshd_origSurf(:,id_sort,jiix);
    mat_plot_crujra_NoLapse_origSurf=mat_all_crujra_noLapse_origSurf(:,id_sort,jiix+4);
    mat_plot_crujra_NoLapse_newSurf=mat_all_crujra_noLapse_newSurf(:,id_sort,jiix+4);

    %wiggle plot does not work properly if there are any NANs -> only very
    %few, just set to 0 for now
    mat_plot(isnan(mat_plot))=0;
    mat_plot_crujra_newSurf(isnan(mat_plot_crujra_newSurf))=0;
    mat_plot_oshd_newSurf(isnan(mat_plot_oshd_newSurf))=0;
    mat_plot_crujra_origSurf(isnan(mat_plot_crujra_origSurf))=0;
    mat_plot_oshd_origSurf(isnan(mat_plot_oshd_origSurf))=0;
    mat_plot_crujra_NoLapse_newSurf(isnan(mat_plot_crujra_NoLapse_newSurf))=0;
    mat_plot_crujra_NoLapse_origSurf(isnan(mat_plot_crujra_NoLapse_origSurf))=0;

    labels=cell(1,size(mat_plot,2));
    elev_plot=elev_inter(id_sort);
    locs_plot=locs(id_sort);

    for piip=1:size(labels,2)
        if piip==size(labels,2)
            lab_in=[cell2mat(locs_plot(piip)) ' ' num2str(elev_plot(piip)) 'm'];
        else
            lab_in=[cell2mat(locs_plot(piip)) ' ' num2str(elev_plot(piip))];
        end
        labels{1,piip}=lab_in;
    end
    labels = cellfun(@(x) strrep(x,' ','\newline'), labels,'UniformOutput',false);

    id_1000=find(elev_plot<=1000);
    mat_plot_1000=[mat_plot(:,id_1000)];
    mat_plot_oshd_newSurf_1000=[mat_plot_oshd_newSurf(:,id_1000)];
    mat_plot_crujra_newSurf_1000=[mat_plot_crujra_newSurf(:,id_1000)];
    mat_plot_oshd_origSurf_1000=[mat_plot_oshd_origSurf(:,id_1000)];
    mat_plot_crujra_origSurf_1000=[mat_plot_crujra_origSurf(:,id_1000)];
    mat_plot_crujra_NoLapse_newSurf_1000=[mat_plot_crujra_NoLapse_newSurf(:,id_1000)];
    mat_plot_crujra_NoLapse_origSurf_1000=[mat_plot_crujra_NoLapse_origSurf(:,id_1000)];
    labels_1000=labels(id_1000);

    id_2000=find(elev_plot>1000 & elev_plot<=2000);
    mat_plot_2000=[mat_plot(:,id_2000)];
    mat_plot_oshd_newSurf_2000=[mat_plot_oshd_newSurf(:,id_2000)];
    mat_plot_crujra_newSurf_2000=[mat_plot_crujra_newSurf(:,id_2000)];
    mat_plot_oshd_origSurf_2000=[mat_plot_oshd_origSurf(:,id_2000)];
    mat_plot_crujra_origSurf_2000=[mat_plot_crujra_origSurf(:,id_2000)];
    mat_plot_crujra_NoLapse_newSurf_2000=[mat_plot_crujra_NoLapse_newSurf(:,id_2000)];
    mat_plot_crujra_NoLapse_origSurf_2000=[mat_plot_crujra_NoLapse_origSurf(:,id_2000)];
    labels_2000=labels(id_2000);

    id_3000=find(elev_plot>2000);
    mat_plot_3000=[mat_plot(:,id_3000)];
    mat_plot_oshd_newSurf_3000=[mat_plot_oshd_newSurf(:,id_3000)];
    mat_plot_crujra_newSurf_3000=[mat_plot_crujra_newSurf(:,id_3000)];
    mat_plot_oshd_origSurf_3000=[mat_plot_oshd_origSurf(:,id_3000)];
    mat_plot_crujra_origSurf_3000=[mat_plot_crujra_origSurf(:,id_3000)];
    mat_plot_crujra_NoLapse_newSurf_3000=[mat_plot_crujra_NoLapse_newSurf(:,id_3000)];
    mat_plot_crujra_NoLapse_origSurf_3000=[mat_plot_crujra_NoLapse_origSurf(:,id_3000)];
    labels_3000=labels(id_3000);

    %% wiggle 1000
    min_all_plot=min([min(mat_plot_1000,[],'all'),min(mat_plot_oshd_newSurf_1000,[],'all'),...
        min(mat_plot_crujra_newSurf_1000,[],'all'),min(mat_plot_crujra_NoLapse_newSurf_1000,[],'all')]);

    max_all_plot=max([max(mat_plot_1000,[],'all'),max(mat_plot_oshd_newSurf_1000,[],'all'),...
        max(mat_plot_crujra_newSurf_1000,[],'all'),max(mat_plot_crujra_NoLapse_newSurf_1000,[],'all')]);

    max_overall=max([abs(min_all_plot),abs(max_all_plot)]);
    fh=figure('Position', [491 13 481 667]);
    ah1=subplot(411);hold on;grid on
    mat_plot_1000(end,end)=max_overall;
    wiggle(round(mat_plot_1000,2),'BR2')
    xticks([]);

    id_first=find(day(timing_check1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels(datestr(timing_check1([id_first(1:2:end)' 365],1)))

    ah2=subplot(412);hold on;grid on
    mat_plot_oshd_newSurf_1000(end,end)=max_overall;
    wiggle(mat_plot_oshd_newSurf_1000,'BR2')
    xticks([]);
    yticks([id_first(1:2:end)' 365])
    yticklabels(datestr(timing_check1([id_first(1:2:end)' 365],1)))

    ah3=subplot(413);hold on;grid on
    mat_plot_crujra_newSurf_1000(end,end)=max_overall;
    wiggle(mat_plot_crujra_newSurf_1000,'BR2')
    xticks([]);

    id_first=find(day(timing_check_crujra1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels(datestr(timing_check_crujra1([id_first(1:2:end)' 365],1)))

    ah4=subplot(414);hold on;grid on
    an4=annotation('textbox',[0.6692 0.0938 0.2210 0.0545],'String',{['min \Delta = ' num2str(round(min_all_plot,2)) 'm'],...
        ['max \Delta = ' num2str(round(max_all_plot,2)) 'm']},'FitBoxToText','on','BackgroundColor',rgb('light red'));
    an5=annotation('textbox',[0.0241 0.1060 0.3758 0.0348],'String',{'blue=too much snow in model'},...
        'FitBoxToText','on','BackgroundColor',rgb('light blue'));
    mat_plot_crujra_NoLapse_newSurf_1000(end,end)=max_overall;
    wiggle(mat_plot_crujra_NoLapse_newSurf_1000,'BR2')
    xticks(1:size(mat_plot_crujra_NoLapse_newSurf_1000,2))
    a = gca;
    a.XTickLabel = labels_1000;
    id_first=find(day(timing_check_crujra1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels(datestr(timing_check_crujra1([id_first(1:2:end)' 365],1)))
    ylim([ah1 ah2 ah3 ah4],[1 id_first(end)]);
    xlim([ah1 ah2 ah3 ah4],[0 size(labels_1000,2)+0.9])
    set(gcf,'Color','white')
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    set([an4 an5],'FontSize',8)

    ti(1)=title(ah1,'Observed - FSM2');
    ti(2)=title(ah2,'Observed - CLM5  {\itClim_{OSHD 1km} LU_{HR 1km}}','interpreter','tex');
    ti(3)=title(ah3,'Observed - CLM5 {\itClim_{CRU* 1km} LU_{HR 1km}}','interpreter','tex');
    ti(4)=title(ah4,'Observed - CLM5 {\itClim_{CRU 1km} LU_{HR 1km}}','interpreter','tex');
    set(ti,'FontSize',8,'Color',rgb('red'),'FontWeight','Normal')

    ah1.Position=[0.1489 0.7673 0.7561 0.1577];
    ah2.Position=[0.1489 0.5800 0.7561 0.1488];
    ah3.Position=[0.1489 0.3882 0.7561 0.1488];
    ah4.Position=[0.1489 0.1964 0.7561 0.1488];

    fh.Position = [486 105 517 575];
    an5.Position = [0.0669 0.0924 0.3569 0.0400];
    an4.Position = [0.6794 0.0736 0.2099 0.0626];
    save_fig=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_all1000.png']);
    save_fig1=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_all1000.pdf']);
    exportgraphics(fh,save_fig,'Resolution',700)
    exportgraphics(fh,save_fig1)
    %% wiggle 2000
    min_all_plot=min([min(mat_plot_2000,[],'all'),min(mat_plot_oshd_newSurf_2000,[],'all'),...
        min(mat_plot_crujra_newSurf_2000,[],'all'),min(mat_plot_crujra_NoLapse_newSurf_2000,[],'all')]);

    max_all_plot=max([max(mat_plot_2000,[],'all'),max(mat_plot_oshd_newSurf_2000,[],'all'),...
        max(mat_plot_crujra_newSurf_2000,[],'all'),max(mat_plot_crujra_NoLapse_newSurf_2000,[],'all')]);

    max_overall=max([abs(min_all_plot),abs(max_all_plot)]);
    fh=figure('Position', [491 13 481 667]);
    ah1=subplot(411);hold on;grid on
    mat_plot_2000(end,end)=max_overall;
    wiggle(round(mat_plot_2000,2),'BR2')
    xticks([]);
    id_first=find(day(timing_check1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels([])

    ah2=subplot(412);hold on;grid on
    mat_plot_oshd_newSurf_2000(end,end)=max_overall;
    wiggle(mat_plot_oshd_newSurf_2000,'BR2')
    xticks([]);
    yticks([id_first(1:2:end)' 365]);
    yticklabels([]);

    ah3=subplot(413);hold on;grid on
    mat_plot_crujra_newSurf_2000(end,end)=max_overall;
    wiggle(mat_plot_crujra_newSurf_2000,'BR2')
    xticks([]);

    id_first=find(day(timing_check_crujra1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels([])

    ah4=subplot(414);hold on;grid on
    an4=annotation('textbox',[0.6692 0.0938 0.2210 0.0545],'String',{['min \Delta = ' num2str(round(min_all_plot,2)) 'm'],...
        ['max \Delta = ' num2str(round(max_all_plot,2)) 'm']},'FitBoxToText','on','BackgroundColor',rgb('cerulean'));
    mat_plot_crujra_NoLapse_newSurf_2000(end,end)=max_overall;
    wiggle(mat_plot_crujra_NoLapse_newSurf_2000,'BR2')
    xticks(1:size(mat_plot_crujra_NoLapse_newSurf_2000,2))
    a = gca;
    a.XTickLabel = labels_2000;
    id_first=find(day(timing_check_crujra1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels([]);
    ylim([ah1 ah2 ah3 ah4],[1 id_first(end)]);
    xlim([ah1 ah2 ah3 ah4],[0 size(labels_2000,2)+0.9])
    set(gcf,'Color','white')
    set(findall(gcf,'-property','FontSize'),'FontSize',7)
    set(an4,'FontSize',8)

    ti(1)=title(ah1,'Observed - FSM2');
    ti(2)=title(ah2,'Observed - CLM5  {\itClim_{OSHD 1km} LU_{HR 1km}}','interpreter','tex');
    ti(3)=title(ah3,'Observed - CLM5 {\itClim_{CRU* 1km} LU_{HR 1km}}','interpreter','tex');
    ti(4)=title(ah4,'Observed - CLM5 {\itClim_{CRU 1km} LU_{HR 1km}}','interpreter','tex');
    set(ti,'FontSize',8,'Color',rgb('blue'),'FontWeight','Normal');

    ah1.Position=[0.1489 0.7673 0.7561 0.1577];
    ah2.Position=[0.1489 0.5800 0.7561 0.1488];
    ah3.Position=[0.1489 0.3882 0.7561 0.1488];
    ah4.Position=[0.1489 0.1964 0.7561 0.1488];

    fh.Position = [486 105 517 575];
    an5.Position = [0.0669 0.0924 0.3569 0.0400];
    an4.Position = [0.6794 0.0736 0.2099 0.0626];
    save_fig=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_all2000.png']);
    save_fig1=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_all2000.pdf']);
    exportgraphics(fh,save_fig,'Resolution',700)
    exportgraphics(fh,save_fig1)
    %% wiggle 3000
    min_all_plot=min([min(mat_plot_3000,[],'all'),min(mat_plot_oshd_newSurf_3000,[],'all'),...
        min(mat_plot_crujra_newSurf_3000,[],'all'),min(mat_plot_crujra_NoLapse_newSurf_3000,[],'all')]);

    max_all_plot=max([max(mat_plot_3000,[],'all'),max(mat_plot_oshd_newSurf_3000,[],'all'),...
        max(mat_plot_crujra_newSurf_3000,[],'all'),max(mat_plot_crujra_NoLapse_newSurf_3000,[],'all')]);

    max_overall=max([abs(min_all_plot),abs(max_all_plot)]);
    %%
    fh=figure('Position', [491 13 481 667]);
    ah1=subplot(411);hold on;grid on
    mat_plot_3000(end,end)=max_overall;
    wiggle(round(mat_plot_3000,2),'BR2')
    xticks([]);
    id_first=find(day(timing_check1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels([])

    ah2=subplot(412);hold on;grid on
    mat_plot_oshd_newSurf_3000(end,end)=max_overall;
    wiggle(mat_plot_oshd_newSurf_3000,'BR2')
    xticks([]);
    yticks([id_first(1:2:end)' 365])
    yticklabels([])

    ah3=subplot(413);hold on;grid on
    mat_plot_crujra_newSurf_3000(end,end)=max_overall;
    wiggle(mat_plot_crujra_newSurf_3000,'BR2')
    xticks([]);

    id_first=find(day(timing_check_crujra1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels([])

    ah4=subplot(414);hold on;grid on
    an4=annotation('textbox',[0.6692 0.0938 0.2210 0.0545],'String',{['min \Delta = ' num2str(round(min_all_plot,2)) 'm'],...
        ['max \Delta = ' num2str(round(max_all_plot,2)) 'm']},'FitBoxToText','on','BackgroundColor',rgb('light green'));
    mat_plot_crujra_NoLapse_newSurf_3000(end,end)=max_overall;
    wiggle(mat_plot_crujra_NoLapse_newSurf_3000,'BR2')
    xticks(1:size(mat_plot_crujra_NoLapse_newSurf_3000,2))
    a = gca;
    a.XTickLabel = labels_3000;
    id_first=find(day(timing_check_crujra1)==1);
    yticks([id_first(1:2:end)' 365])
    yticklabels([])

    ylim([ah1 ah2 ah3 ah4],[1 id_first(end)]);
    xlim([ah1 ah2 ah3 ah4],[0 size(labels_3000,2)+0.9])
    set(gcf,'Color','white')
    set(findall(gcf,'-property','FontSize'),'FontSize',6)
    set(an4,'FontSize',8)

    ti(1)=title(ah1,'Observed - FSM2');
    ti(2)=title(ah2,'Observed - CLM5  {\itClim_{OSHD 1km} LU_{HR 1km}}','interpreter','tex');
    ti(3)=title(ah3,'Observed - CLM5 {\itClim_{CRU* 1km} LU_{HR 1km}}','interpreter','tex');
    ti(4)=title(ah4,'Observed - CLM5 {\itClim_{CRU 1km} LU_{HR 1km}}','interpreter','tex');
    set(ti,'FontSize',8,'Color',rgb('green'),'FontWeight','Normal');
    ah1.Position=[0.1489 0.7673 0.7561 0.1577];
    ah2.Position=[0.1489 0.5800 0.7561 0.1488];
    ah3.Position=[0.1489 0.3882 0.7561 0.1488];
    ah4.Position=[0.1489 0.1964 0.7561 0.1488];

    fh.Position = [486 105 517 575];
    an5.Position = [0.0669 0.0924 0.3569 0.0400];
    an4.Position = [0.6794 0.0736 0.2099 0.0626];
    save_fig=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_all3000.png']);
    save_fig1=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_all3000.pdf']);
    exportgraphics(fh,save_fig1)
    exportgraphics(fh,save_fig,'Resolution',700)

    %%
    fh=figure('Position',[290 424 1560 526]);
    annotation('textbox',[0.7944598062954,0.95,0.095036316960689,0.04892966245414],'String',{['min \Delta = ' num2str(round(min(min(mat_plot)),2)) 'm'],...
        ['max \Delta = ' num2str(round(max(max(mat_plot)),2)) 'm']},'FitBoxToText','on','BackgroundColor',rgb('light green'));
    annotation('textbox',[0.15 0.95 0.3758 0.0348],'String',{'blue==too much snow in model'},...
        'FitBoxToText','on','BackgroundColor',rgb('light blue'));
    %mat_plot(end,end)=max_abs;
    wiggle(round(mat_plot,2),'BR1')
    xticks(1:size(mat_plot,2))
    a = gca;
    a.XTickLabel = labels;
    id_first=find(day(timing_check1)==1);
    yticks([id_first' 365])
    yticklabels(datestr(timing_check1([id_first' 365],1)))
    title([num2str(double(string(years_u(jiix)))) '/' num2str(double(string(years_u(jiix)))+1) ': Observed - FSM2'])
    set(gcf,'Color','white')
    a.FontSize=12;
    a.XAxis.FontSize=8;
    ylim([1 id_first(end)]);
    xlim([0 size(labels,2)+0.9])
    save_fig=fullfile(bf_out,['wiggle_' num2str(double(string(years_u(jiix)))) '_fsm2.png']);
    exportgraphics(fh,save_fig,'Resolution',700)

end

end
