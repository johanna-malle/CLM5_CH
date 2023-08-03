function make_taylor_plots_new
bf_taylor ='/home/lud11/malle/CLM5_CH/taylor';
load(fullfile(bf_taylor,'taylor_stats_v3.mat'));

%% make split taylor diagram

fh=figure('Position',[565 461 1199 420]);hold on;%('Position',[1993 195 636 1115]);

colors_in = [[105/255,105/255,105/255];rgb('pinkish purple');rgb('leafy green');rgb('azure');rgb('mauve');rgb('pale green');rgb('very light blue');...
    [217/255,95/255,2/255];[117/255,112/255,179/255];[27/255,158/255,119/255];[217/255,95/255,2/255];[117/255,112/255,179/255];[27/255,158/255,119/255]];


tiledlayout(1,3,'TileSpacing','none','Padding','none');
ax = nexttile ;
ax.FontName = 'Serif';
statm_all=[HS_dec_ALL_stats;HS_feb_ALL_stats;HS_apr_ALL_stats];
max_std_all=max(statm_all(:,2));
[pp tt xx1 kk ww] = taylordiag(squeeze(HS_dec_ALL_stats(:,2)),squeeze(HS_dec_ALL_stats(:,3)),squeeze(HS_dec_ALL_stats(:,4)),...
    'tickRMS',[],'titleRMS',0,'tickRMSangle',87,'showlabelsRMS',1,'widthRMS',0.9,...
    'tickSTD',[0:0.1:max_std_all],'limSTD',round(max_std_all,2)+0.08,'npan',1,...
    'tickCOR',[.1:.2:.8 0.8 0.9 .95 .99 1],'showlabelsCOR',1,'titleCOR',1,'colCOR',rgb('silver'));
set(xx1,'Color',[105/255,105/255,105/255],'LineWidth',1.9);


for ii = 1 : length(tt)
    set(tt(ii),'fontsize',10,'fontweight','bold','FontName','SansSerif')
    set(pp(ii),'markersize',12)
    if ii == 1
        set(tt(ii),'String','FSM2 1 Dez','Color',[105/255,105/255,105/255],'fontsize',11,'verticalalignment','top','horizontalalignment','right','Rotation',10);
        set(pp(ii),'Color',[105/255,105/255,105/255],'MarkerSize',30)
    elseif ismember(ii,[8,9,10])
        set(tt(ii),'String',column_label_plots{ii},'Color',colors_in(ii,:),'horizontalalignment','left','Interpreter','latex');
        set(pp(ii),'Color',colors_in(ii,:),'marker','d','MarkerFaceColor',colors_in(ii,:),'MarkerEdgeColor',colors_in(1,:),'MarkerSize',7.5)
        set(kk(ii),'Color',colors_in(ii,:),'LineWidth',0.9);
        set(ww(ii),'Color',colors_in(ii,:),'FontSize',10,'verticalalignment','top','horizontalalignment','left');
    else
        set(tt(ii),'String',column_label_plots{ii},'Color',colors_in(ii,:),'horizontalalignment','left','Interpreter','latex');
        set(pp(ii),'Color',colors_in(ii,:),'marker','o','MarkerFaceColor',colors_in(ii,:),'MarkerEdgeColor',colors_in(1,:),'MarkerSize',7.5)
        set(kk(ii),'Color',colors_in(ii,:),'LineWidth',0.9);
        set(ww(ii),'Color',colors_in(ii,:),'FontSize',10,'verticalalignment','top','horizontalalignment','left');
   end
end
title('(a) Early accumulation period')
ax.TitleHorizontalAlignment = 'left';
set(gca,'Units','normalized')
titleHandle = get( gca ,'Title' );
pos  = get( titleHandle , 'position' );
pos1 = pos - [0 0.05 0] ;
set( titleHandle , 'position' , pos1 );

ax = nexttile ;
ax.FontName = 'Serif';

[pp2 tt2 xx2 kk ww] = taylordiag(squeeze(HS_feb_ALL_stats(:,2)),squeeze(HS_feb_ALL_stats(:,3)),squeeze(HS_feb_ALL_stats(:,4)),...
    'tickRMS',[],'titleRMS',0,'tickRMSangle',87,'showlabelsRMS',1,'widthRMS',0.9,...
    'tickSTD',[0:0.1:max_std_all],'limSTD',round(max_std_all,2)+0.2,'npan',1,...
    'tickCOR',[.1:.2:.8 0.8 0.9 .95 .99 1],'showlabelsCOR',1,'titleCOR',1,'colCOR',rgb('silver'));
set(xx2,'Color',[105/255,105/255,105/255],'LineWidth',1.9);


for ii = 1 : length(tt)
    set(tt2(ii),'fontsize',10,'fontweight','bold','FontName', 'SansSerif')
    set(pp2(ii),'markersize',12)
    if ii == 1
        set(tt2(ii),'String','FSM2 1 Feb','Color',[105/255,105/255,105/255],'fontsize',11,'verticalalignment','top','horizontalalignment','right','Rotation',10);
        set(pp2(ii),'Color',[105/255,105/255,105/255],'MarkerSize',30)
    elseif ismember(ii,[8,9,10])
        set(tt2(ii),'String',column_label_plots{ii},'Color',colors_in(ii,:),'horizontalalignment','left','Interpreter','latex');
        set(pp2(ii),'Color',colors_in(ii,:),'marker','d','MarkerFaceColor',colors_in(ii,:),'MarkerEdgeColor',colors_in(1,:),'MarkerSize',7.5)
        set(kk(ii),'Color',colors_in(ii,:),'LineWidth',0.9);
        set(ww(ii),'Color',colors_in(ii,:),'FontSize',10,'verticalalignment','top','horizontalalignment','left');
    else
        set(tt2(ii),'String',column_label_plots{ii},'Color',colors_in(ii,:),'horizontalalignment','left','Interpreter','latex');
        set(pp2(ii),'Color',colors_in(ii,:),'marker','o','MarkerFaceColor',colors_in(ii,:),'MarkerEdgeColor',colors_in(1,:),'MarkerSize',7.5)
        set(kk(ii),'Color',colors_in(ii,:),'LineWidth',0.9);
        set(ww(ii),'Color',colors_in(ii,:),'FontSize',10,'verticalalignment','top','horizontalalignment','left');
    end
end
title('(b) Mid-accumulation period')
ax.TitleHorizontalAlignment = 'left';
titleHandle = get( gca ,'Title' );
pos  = get( titleHandle , 'position' );
pos1 = pos - [0 0.05 0] ;
set( titleHandle , 'position' , pos1 );
ax = nexttile ;
ax.FontName = 'Serif';

[pp3 tt3 xx3 kk ww] = taylordiag(squeeze(HS_apr_ALL_stats(:,2)),squeeze(HS_apr_ALL_stats(:,3)),squeeze(HS_apr_ALL_stats(:,4)),...
    'tickRMS',[],'titleRMS',0,'tickRMSangle',87,'showlabelsRMS',1,'widthRMS',0.9,...
    'tickSTD',[0:0.1:max_std_all],'limSTD',round(max_std_all,2)+0.2,'npan',1,...
    'tickCOR',[.1:.2:.8 0.8 0.9 .95 .99 1],'showlabelsCOR',1,'titleCOR',1,'colCOR',rgb('silver'));
set(xx3,'Color',[105/255,105/255,105/255],'LineWidth',1.9);

for ii = 1 : length(tt3)
    set(tt3(ii),'fontsize',10,'fontweight','bold','FontName', 'SansSerif')
    set(pp3(ii),'markersize',12)
    if ii == 1
        set(tt3(ii),'String','FSM2 1 Apr','Color',[105/255,105/255,105/255],'fontsize',11,'verticalalignment','top','horizontalalignment','right','Rotation',10);
        set(pp3(ii),'Color',[105/255,105/255,105/255],'MarkerSize',30)
    elseif ismember(ii,[8,9,10])
        set(tt3(ii),'String',column_label_plots{ii},'Color',colors_in(ii,:),'horizontalalignment','left','Interpreter','latex');
        set(pp3(ii),'Color',colors_in(ii,:),'Marker','d','MarkerFaceColor',colors_in(ii,:),'MarkerEdgeColor',colors_in(1,:),'MarkerSize',7.5)
        set(kk(ii),'Color',colors_in(ii,:),'LineWidth',0.9);
        set(ww(ii),'Color',colors_in(ii,:),'FontSize',10,'verticalalignment','top','horizontalalignment','left');
    else
        set(tt3(ii),'String',column_label_plots{ii},'Color',colors_in(ii,:),'horizontalalignment','left','Interpreter','latex');
        set(pp3(ii),'Color',colors_in(ii,:),'marker','o','MarkerFaceColor',colors_in(ii,:),'MarkerEdgeColor',colors_in(1,:),'MarkerSize',7.5)
        set(kk(ii),'Color',colors_in(ii,:),'LineWidth',0.9);
        set(ww(ii),'Color',colors_in(ii,:),'FontSize',10,'verticalalignment','top','horizontalalignment','right');
    end
end
title('(c) Ablation period')
ax.TitleHorizontalAlignment = 'left';
titleHandle = get( gca ,'Title' );
pos  = get( titleHandle , 'position' );
pos1 = pos - [0 0.05 0] ;
set( titleHandle , 'position' , pos1 );
set(gcf,'Color','white');
set(gca, 'FontName', 'Serif')

exportgraphics(gcf,fullfile(bf_taylor,'HS_all_seasons_split_v2.pdf'),'ContentType','vector');
exportgraphics(gcf,fullfile(bf_taylor,'HS_all_seasons_split_v2.eps'),'ContentType','vector');