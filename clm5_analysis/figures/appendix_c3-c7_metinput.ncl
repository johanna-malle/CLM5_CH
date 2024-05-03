load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

setfileoption("nc","Format","NetCDF4Classic")
month_abbr = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
              "Oct","Nov","Dec"/)

bf_precip = "/home/malle/Documents/input_checks/precip/"
bf_temp = "/home/malle/Documents/input_checks/temp/"

pldir = "/home/malle/Documents/input_checks/plots_precip/"
pldir_temp = "/home/malle/Documents/input_checks/plots_temp/"

;need file to copy lat long from:
bf_OSHD1 = "/media/malle/LaCie1/CLM5_SP_TILES/OSHD_FILES/"   ; input directory
FSNO_OSHD = addfile(bf_OSHD1 + "FSNO_SP_OSHD.nc","r")
FSNO_OSHD_data=FSNO_OSHD->DATA
copy_latlon = FSNO_OSHD_data(lat|:,lon|:,time|:)  ; needed in order to copy lat lon later
copy_latlon@units="mm"
ismissing_oshd = (ismissing(FSNO_OSHD_data(lat|:,lon|:,time|0)))
;; for plotting later:
dir       = "shp_CH/" ;this needs to be in the same directory as this script for now
filename  = "gadm36_CHE_0.shp"
lnres                  = True
lnres@gsLineColor      = "orangered4"
lnres@gsLineThicknessF = 3
lbcolors_transparent = conform_dims(50,-1,-1)

color_plot12   = "precip_11lev";"WhiteBlueGreenYellowRed"     ; color map for 1st plot
color_plot3   = "nrl_sirkes"  ; color map for 2nd plot
color_plot4 = "precip2_17lev"     ; color map for 3rd and 4th plots

color_temp_diff = "temp_diff_18lev"
color_temp = "hotcold_18lev"

cmap = read_colormap_file("MPL_gist_earth")
cmap = cmap(::-1,:) ; reverse the color map

;-----------------------------------------------------------
;CRUJRA
;-----------------------------------------------------------
;start with monthly plots precip
fils_cru = systemfunc ("ls "+bf_precip+"month*CRUJRA*.nc") ; file paths
f_cru    = addfiles (fils_cru, "r")
X_cru    = addfiles_GetVar(f_cru,fils_cru,"precip")

;compute climatology over these 7 years
tClm_crujra2 = clmMonTLL(X_cru(time|:, lat|:, lon|:))
tClm_crujra = tClm_crujra2(lat|:, lon|:, month|:)
copy_VarCoords_2(copy_latlon,tClm_crujra)

wks = gsn_open_wks("png",pldir+"CRUJRA_month_avg_2015-2019")            ; open a png file
plot = new(13,graphic)
res          = True
res@gsnDraw  = False                          ; don't draw
res@gsnFrame = False                          ; don't advance frame
res@cnInfoLabelOn = False                     ; turn off cn info label

res@cnFillOn            = True
res@cnLinesOn           = False

res@mpLimitMode         = "Corners"
res@mpLeftCornerLatF    = 45.63320611929412
res@mpLeftCornerLonF    = 5.92;5.840585642273973
res@mpRightCornerLatF   = 48.0;48.05736089870589
res@mpRightCornerLonF   = 10.50507112172603
res@lbLabelBarOn        = False            ; turn off individual cb's

res@mpProjection        = "LambertConformal"
res@mpCenterLonF        = 8.1728
res@mpCenterLatF        = 46.8452
res@mpLambertParallel1F = 47
res@mpLambertParallel2F = 47
res@mpLambertMeridianF  = 8
res@cnLevelSelectionMode = "ManualLevels"
res@cnLineLabelsOn = False                  ; show line labels

res@pmTickMarkDisplayMode = "Always"
res@gsnRightString    = ""
res@gsnLeftString     = ""
res@gsnAddCyclic         = False              ; regional data
res@mpOutlineOn           = False         ; Use outlines from shapefile

res@gsnDraw              = False
res@gsnFrame             = False
res@cnLineLabelsOn       = False     ; Turn off contour line labels
res@cnFillOn             = True      ; Turn on contour fill
res@cnInfoLabelOn        = False     ; Turn off info label

sz_cmap=dimsizes(cmap)
res@cnFillPalette = cmap(5:sz_cmap(0)-1,:)

plot(0) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP                  = True                   ; modify the panel plot
resP@gsnFrame         = False                  ; don't advance the frame, so we can use gsn_text_ndc
resP@lbLabelFontHeightF = 0.012    ; set font height of Label Bar labels
resP@gsnPanelLabelBar = True                   ; add common colorbar
resP@lbTitleFontHeightF= .015                 ; make title smaller
resP@gsnPanelBottom   = 0.05      ; shrink panel plot by setting bottom edge of plot
resP@gsnMaximize    = True
resP@gsnFrame       = False
resP@gsnPanelFigureStrings= month_abbr ; add strings to panel
resP@gsnPanelMainString = "Monthly sum of precip for ClimCRU (2014-2019 mean)"
resP@lbBoxEndCapStyle  = "TriangleBothEnds"     ; default is "RectangleEnds", other option "TriangleBothEnds"
resP@lbTitleString = tClm_crujra@units
resP@lbTitlePosition            = "Bottom"            ; title position
resP@lbTitleDirection           = "Across"            ; title direction
resP@pmLabelBarOrthogonalPosF  = -0.02
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

;now loop through individual years:
years = (/2014,2015,2016,2017,2018,2019/)

do idx = 0,71,12
precip_int=X_cru(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

wks = gsn_open_wks("png",pldir+"CRUJRA_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly sum of P (CRUJRA): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

delete(precip_int)
end do

;-----------------------------------------------------------
;OSHD
;-----------------------------------------------------------
;start with monthly plots precip
fils_oshd = systemfunc ("ls "+bf_precip+"month*OSHD*.nc") ; file paths
f_oshd    = addfiles (fils_oshd, "r")
X_oshd    = addfiles_GetVar(f_oshd,fils_oshd,"precip")

;compute climatology over these 6 years
tClm_oshd1 = clmMonTLL(X_oshd(time|:, lat|:, lon|:))
tClm_oshd = tClm_oshd1(lat|:, lon|:, month|:)
copy_VarCoords_1(copy_latlon,tClm_oshd)

do check_nan = 0,11
   tClm_oshd(:,:,check_nan)=where(ismissing_oshd.eq.0,tClm_oshd(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir+"OSHD_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly sum of precip for ClimOSHD (2014-2019 mean)"
resP@lbTitleString = tClm_oshd@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

;now loop through individual years:
do idx = 0,71,12
precip_int=X_oshd(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

do check_nan = 0,11
   precip_int(:,:,check_nan)=where(ismissing_oshd.eq.0,precip_int(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir+"OSHD_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly sum of P (OSHD): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)
delete(precip_int)
end do

;this is just needed for diff plots later:
tClm_oshd_p=tClm_oshd
tClm_crujra_p=tClm_crujra
X_oshd_p=X_oshd
X_cru_p=X_cru

delete(tClm_oshd)
delete(tClm_crujra)
delete(X_oshd)
delete(X_cru)

;----------------------------------------------------------------
; NOW THE SAME FOR TEMP
;---------------------------------------------------------------
;-----------------------------------------------------------
;CRUJRA
;-----------------------------------------------------------
;start with monthly plots precip
fils_cru = systemfunc ("ls "+bf_temp+"month*CRUJRA*wL.nc") ; file paths
f_cru    = addfiles (fils_cru, "r")
X_cru    = addfiles_GetVar(f_cru,fils_cru,"temp")

;compute climatology over these 6 years
tClm_crujra1 = clmMonTLL(X_cru(time|:, lat|:, lon|:))
tClm_crujra = tClm_crujra1(lat|:, lon|:, month|:)
copy_VarCoords_1(copy_latlon,tClm_crujra)

wks = gsn_open_wks("png",pldir_temp+"CRUJRA_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

sz_cmap=dimsizes(cmap)
delete(res@cnFillPalette)
res@cnFillPalette = color_temp

plot(0) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,tClm_crujra(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelFigureStrings= month_abbr ; add strings to panel
resP@gsnPanelMainString = "Monthly mean temp for ClimCRU* (2014-2019 mean)"
resP@lbBoxEndCapStyle  = "TriangleBothEnds"     ; default is "RectangleEnds", other option "TriangleBothEnds"
resP@lbTitleString = tClm_crujra@units
resP@lbTitlePosition            = "Bottom"            ; title position
resP@lbTitleDirection           = "Across"            ; title direction
resP@pmLabelBarOrthogonalPosF  = -0.02
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

do idx = 0,71,12
precip_int=X_cru(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

wks = gsn_open_wks("png",pldir_temp+"CRUJRA_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly avg. T (CRUJRA, lapse): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

delete(precip_int)
end do
;-----------------------------------------------------------
;CRUJRA - no lapse
;-----------------------------------------------------------
;start with monthly plots precip
fils_cru_nL = systemfunc ("ls "+bf_temp+"month*CRUJRA*nL.nc") ; file paths
f_cru_nL    = addfiles (fils_cru_nL, "r")
X_cru_nL    = addfiles_GetVar(f_cru_nL,fils_cru_nL,"temp")

;compute climatology over these 6 years
tClm_crujra1 = clmMonTLL(X_cru_nL(time|:, lat|:, lon|:))
tClm_crujra_nL = tClm_crujra1(lat|:, lon|:, month|:)
copy_VarCoords_1(copy_latlon,tClm_crujra_nL)

wks = gsn_open_wks("png",pldir_temp+"CRUJRAnL_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,tClm_crujra_nL(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly mean temp for ClimCRU (2014-2019 mean)"
resP@lbTitleString = tClm_crujra@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

do idx = 0,71,12
precip_int=X_cru_nL(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

wks = gsn_open_wks("png",pldir_temp+"CRUJRAnL_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly avg. T (CRUJRA, noL): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

delete(precip_int)
end do
;-----------------------------------------------------------
;OSHD
;-----------------------------------------------------------
;start with monthly plots precip
fils_oshd = systemfunc ("ls "+bf_temp+"month*OSHD*.nc") ; file paths
f_oshd    = addfiles (fils_oshd, "r")
X_oshd    = addfiles_GetVar(f_oshd,fils_oshd,"temp")

;compute climatology over these 6 years
tClm_oshd1 = clmMonTLL(X_oshd(time|:, lat|:, lon|:))
tClm_oshd = tClm_oshd1(lat|:, lon|:, month|:)
copy_VarCoords_1(copy_latlon,tClm_oshd)

do check_nan = 0,11
   tClm_oshd(:,:,check_nan)=where(ismissing_oshd.eq.0,tClm_oshd(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir_temp+"OSHD_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,tClm_oshd(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly mean temp for ClimOSHD (2014-2019 mean)"
resP@lbTitleString = tClm_oshd@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

;now loop through individual years:
do idx = 0,71,12
precip_int=X_oshd(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

do check_nan = 0,11
   precip_int(:,:,check_nan)=where(ismissing_oshd.eq.0,precip_int(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir_temp+"OSHD_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly avg. P (OSHD): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)
delete(precip_int)
end do

;-----------------------------------------------------------
;DIFF - precip!
;-----------------------------------------------------------

diff_tClm=tClm_oshd_p
diff_tClm=tClm_oshd_p - tClm_crujra_p

diff_X=X_oshd_p
diff_X=X_oshd_p - X_cru_p

wks = gsn_open_wks("png",pldir+"Diff_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

delete(res@cnFillPalette)
res@cnFillPalette = color_plot3
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF     = ((max(abs(decimalPlaces(diff_tClm,0,True) )))*2)/20
res@cnMinLevelValF      = -(max(abs(decimalPlaces(diff_tClm,0,True))))
res@cnMaxLevelValF      = (max(abs(decimalPlaces(diff_tClm,0,True))))

plot(0) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Diff. in monthly precip sum (ClimOSHD-ClimCRU, 2014-2019 mean)"
resP@lbTitleString =  "~F33~D~F~"+diff_tClm@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

;now loop through individual years:
do idx = 0,71,12
precip_int=diff_X(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

do check_nan = 0,11
   precip_int(:,:,check_nan)=where(ismissing_oshd.eq.0,precip_int(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir+"Diff_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly sum of P (Diff, OSHD-CRU): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)
delete(precip_int)
end do

;-----------------------------------------------------------
;DIFF - temp
;-----------------------------------------------------------

diff_tClm=tClm_oshd
diff_tClm=tClm_oshd - tClm_crujra

diff_X=X_oshd
diff_X=X_oshd - X_cru

wks = gsn_open_wks("png",pldir_temp+"Diff_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

delete(res@cnFillPalette)
res@cnFillPalette = color_temp_diff
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF     = ((max(abs(decimalPlaces(diff_tClm,0,True) )))*2)/20
res@cnMinLevelValF      = -(max(abs(decimalPlaces(diff_tClm,0,True))))
res@cnMaxLevelValF      = (max(abs(decimalPlaces(diff_tClm,0,True))))

plot(0) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Diff. in monthly mean temp (ClimOSHD-ClimCRU*, 2014-2019 mean)"
resP@lbTitleString =  "~F33~D~F~"+diff_tClm@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

;now loop through individual years:
do idx = 0,71,12
precip_int=diff_X(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

do check_nan = 0,11
   precip_int(:,:,check_nan)=where(ismissing_oshd.eq.0,precip_int(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir_temp+"Diff_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly avg. T (Diff, OSHD-CRU): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)
delete(precip_int)
end do

;-----------------------------------------------------------
;DIFF - temperature, no lapse rate
;-----------------------------------------------------------

diff_tClm=tClm_oshd
diff_tClm=tClm_oshd - tClm_crujra_nL

diff_X=X_oshd
diff_X=X_oshd - X_cru_nL

wks = gsn_open_wks("png",pldir_temp+"Diff_nL_month_avg_2014-2019")            ; open a png file
plot = new(12,graphic)

delete(res@cnFillPalette)
res@cnFillPalette = color_temp_diff
res@cnLevelSelectionMode = "ManualLevels"
res@cnLevelSpacingF     = ((max(abs(decimalPlaces(diff_tClm,0,True) )))*2)/20
res@cnMinLevelValF      = -(max(abs(decimalPlaces(diff_tClm,0,True))))
res@cnMaxLevelValF      = (max(abs(decimalPlaces(diff_tClm,0,True))))

plot(0) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,diff_tClm(lat|:,lon|:,month|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Diff. in monthly mean temp (ClimOSHD-ClimCRU, 2014-2019 mean)"
resP@lbTitleString =  "~F33~D~F~"+diff_tClm@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)

;now loop through individual years:
do idx = 0,71,12
precip_int=diff_X(lat|:,lon|:,time|idx:idx+11)
if idx.eq.0 then
    year_int=years(idx)
else
    year_int=years((idx+1)/12)
end if
copy_VarCoords_1(copy_latlon,precip_int)

do check_nan = 0,11
   precip_int(:,:,check_nan)=where(ismissing_oshd.eq.0,precip_int(:,:,check_nan),copy_latlon@_FillValue)
end do

wks = gsn_open_wks("png",pldir_temp+"Diff_nL_month_"+year_int)            ; open a png file
plot = new(12,graphic)

plot(0) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|0),res)
poly0 = gsn_add_shapefile_polylines(wks,plot(0),dir+filename,lnres)
plot(1) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|1),res)
poly1 = gsn_add_shapefile_polylines(wks,plot(1),dir+filename,lnres)
plot(2) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|2),res)
poly2 = gsn_add_shapefile_polylines(wks,plot(2),dir+filename,lnres)
plot(3) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|3),res)
poly3 = gsn_add_shapefile_polylines(wks,plot(3),dir+filename,lnres)
plot(4) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|4),res)
poly4 = gsn_add_shapefile_polylines(wks,plot(4),dir+filename,lnres)
plot(5) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|5),res)
poly5 = gsn_add_shapefile_polylines(wks,plot(5),dir+filename,lnres)
plot(6) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|6),res)
poly6 = gsn_add_shapefile_polylines(wks,plot(6),dir+filename,lnres)
plot(7)  = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|7),res)
poly7 = gsn_add_shapefile_polylines(wks,plot(7),dir+filename,lnres)
plot(8) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|8),res)
poly8 = gsn_add_shapefile_polylines(wks,plot(8),dir+filename,lnres)
plot(9) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|9),res)
poly9 = gsn_add_shapefile_polylines(wks,plot(9),dir+filename,lnres)
plot(10) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|10),res)
poly10 = gsn_add_shapefile_polylines(wks,plot(10),dir+filename,lnres)
plot(11) = gsn_csm_contour_map(wks,precip_int(lat|:,lon|:,time|11),res)
poly11 = gsn_add_shapefile_polylines(wks,plot(11),dir+filename,lnres)

resP@gsnPanelMainString = "Monthly avg. T (Diff, OSHD-CRU nL): "+year_int
resP@lbTitleString = precip_int@units
gsn_panel(wks,(/plot/),(/4,3/),resP)             ; now draw as one plot
frame(wks)
delete(precip_int)
end do

delete(X_cru)
;-----------------------------------------------------------------------
; YEARLY SUM/YEARLY avg - precip
;----------------------------------------------------------------------

fils = systemfunc ("ls "+bf_precip+"yr*OSHD*.nc") ; file paths
f    = addfiles (fils, "r")
X    = addfiles_GetVar(f,fils,"precip")
OSHD_xAVG = dim_avg_Wrap( X(lat|:,lon|:,year|:) )
copy_VarAtts(copy_latlon, OSHD_xAVG)
copy_VarCoords_1(copy_latlon,OSHD_xAVG)
OSHD_xAVG@units="mm"
OSHD_xAVG=where(ismissing_oshd.eq.0,OSHD_xAVG,OSHD_xAVG@_FillValue)

fils_cru = systemfunc ("ls "+bf_precip+"yr*CRUJRA*.nc") ; file paths
f_crujra    = addfiles (fils_cru, "r")
X_cru    = addfiles_GetVar(f_crujra,fils_cru,"precip")
CRUJRA_xAVG = dim_avg_Wrap( X_cru(lat|:,lon|:,year|:) )
copy_VarAtts(copy_latlon, CRUJRA_xAVG)
copy_VarCoords_1(copy_latlon,CRUJRA_xAVG)
CRUJRA_xAVG@units="mm"

combi_yr_all=array_append_record(OSHD_xAVG,CRUJRA_xAVG,0)

diff_sum_all = CRUJRA_xAVG
diff_sum_all = OSHD_xAVG - CRUJRA_xAVG

chg_perc_all        =  OSHD_xAVG                        ;-- copy variable and retain metadata
chg_perc_all            =  abs(CRUJRA_xAVG-OSHD_xAVG)/where(CRUJRA_xAVG.ne.0,CRUJRA_xAVG,CRUJRA_xAVG@_FillValue) * 100       ;-- compute the difference

wks = gsn_open_wks("png",pldir+"comp_all_yrs")            ; open a png file
delete(plot)
delete(res)
delete(resP)

plot = new(4,graphic)
res          = True
res@gsnDraw  = False                          ; don't draw
res@gsnFrame = False                          ; don't advance frame
res@cnInfoLabelOn = False                     ; turn off cn info label

res@cnFillOn            = True
res@cnLinesOn           = False

res@mpLimitMode         = "Corners"
res@mpLeftCornerLatF    = 45.63320611929412
res@mpLeftCornerLonF    = 5.92;5.840585642273973
res@mpRightCornerLatF   = 48.0;48.05736089870589
res@mpRightCornerLonF   = 10.50507112172603
res@lbLabelBarOn        = True            ; turn off individual cb's

res@mpProjection        = "LambertConformal"
res@mpCenterLonF        = 8.1728
res@mpCenterLatF        = 46.8452
res@mpLambertParallel1F = 47
res@mpLambertParallel2F = 47
res@mpLambertMeridianF  = 8
res@cnLevelSelectionMode = "ManualLevels"
res@cnLineLabelsOn = False                  ; show line labels

res@pmTickMarkDisplayMode = "Always"
res@gsnRightString    = ""
res@gsnLeftString     = ""
res@gsnAddCyclic         = False              ; regional data
res@mpOutlineOn           = False         ; Use outlines from shapefile

res@gsnDraw              = False
res@gsnFrame             = False
res@cnLineLabelsOn       = False     ; Turn off contour line labels
res@cnFillOn             = True      ; Turn on contour fill
res@cnInfoLabelOn        = False     ; Turn off info label

res12=res
res12@cnFillPalette = color_plot12
res12@cnMinLevelValF      = (min(combi_yr_all))
res12@cnMaxLevelValF      = (max(combi_yr_all))
res12@cnLevelSpacingF     = ((max(combi_yr_all))-(min(combi_yr_all)))/12
res12@lbLabelFontColor    = "white"                ; this makes the labelbar labels invisible
res12@lbBoxLinesOn        = False
res12@lbFillColors        = lbcolors_transparent   ; the colors are fully transparent

res3=res
res3@cnFillPalette = color_plot3
res3@lbLabelBarOn        = True
res3@cnLevelSpacingF     = ((max(abs(diff_sum_all)))*2)/20
res3@cnMinLevelValF      = -(max(abs(diff_sum_all)))
res3@cnMaxLevelValF      = (max(abs(diff_sum_all)))
res3@lbTitleString       = "~F33~D~F~"+OSHD_xAVG@units
res3@lbTitlePosition     = "Bottom"            ; title position
res3@lbTitleDirection    = "Across"            ; title direction

res4=res
res4@cnFillPalette = color_plot4
res4@lbLabelBarOn        = True
res4@cnMinLevelValF      = 5
res4@cnMaxLevelValF      = 95
res4@cnLevelSpacingF     = 10
res4@lbTitleString       = "%"
res4@lbTitlePosition            = "Bottom"            ; title position
res4@lbTitleDirection           = "Across"            ; title direction

plot0 = gsn_csm_contour_map(wks,OSHD_xAVG,res12)
poly0 = gsn_add_shapefile_polylines(wks,plot0,dir+filename,lnres)
plot1 = gsn_csm_contour_map(wks,CRUJRA_xAVG,res12)
poly1 = gsn_add_shapefile_polylines(wks,plot1,dir+filename,lnres)
plot2 = gsn_csm_contour_map(wks,diff_sum_all,res3)
poly2 = gsn_add_shapefile_polylines(wks,plot2,dir+filename,lnres)
plot3 = gsn_csm_contour_map(wks,chg_perc_all,res4)
poly3 = gsn_add_shapefile_polylines(wks,plot3,dir+filename,lnres)

resP                  = True                   ; modify the panel plot
resP@gsnFrame         = False                  ; don't advance the frame, so we can use gsn_text_ndc
resP@lbLabelFontHeightF = 0.012    ; set font height of Label Bar labels
resP@lbTitleFontHeightF= .015                 ; make title smaller
resP@gsnPanelBottom   = 0.01      ; shrink panel plot by setting bottom edge of plot
resP@gsnMaximize    = True
resP@gsnFrame       = False
resP@gsnPanelFigureStrings= (/"OSHD","CRUJRA","Diff (OSDH-CRUJRA","%change crujra->OSHD"/) ; add strings to panel
resP@gsnPanelMainString = "Yearly sum of precip. - avg. 2014-2019"

gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),resP)             ; now draw as one plot

getvalues plot0
"vpXF"               : vpx3     ; this is for positioning the common labelbar
end getvalues
getvalues plot0@contour
"lbLabelStrings"     : labels
end getvalues
getvalues plot1
"vpXF"               : vpx4     ; this is for positioning
"vpWidthF"           : vpw4     ; the common lablebar
end getvalues
nlabels = dimsizes(labels)
txres               = True
txres@txFontHeightF = 0.015

;---Set some labelbar resources
lbcolors_fill            = span_color_rgba(res12@cnFillPalette,nlabels+1)
lbres                    = True
lbres@vpWidthF           = (vpx4 + vpw4) - vpx3
lbres@vpHeightF          = 0.08
lbres@lbAutoManage       = False
lbres@lbOrientation      = "horizontal"
lbres@lbJustification    = "CenterCenter"
lbres@lbMonoFillPattern  = True
lbres@lbFillColors       = lbcolors_fill
lbres@lbLabelFontHeightF = 0.012
lbres@lbPerimOn          = False
lbres@lbTitleString      = OSHD_xAVG@units           ; title string
lbres@lbTitleFontHeightF         = .015                 ; make title smaller
lbres@lbTitlePosition            = "Bottom"            ; title position
lbres@lbTitleDirection           = "Across"            ; title direction
lbres@lbLabelStride              = 2
lbres@pmLabelBarHeightF          = 0.1
lbres@lbBoxEndCapStyle           = "TriangleBothEnds"

gsn_labelbar_ndc(wks,nlabels,labels,vpx3,0.565,lbres)
frame(wks)

delete(f_crujra)

do yrs_int = 2014,2019
    print(yrs_int)
    fils_precip_OSHD = systemfunc ("ls "+bf_precip+"yr_"+"*OSHD*"+yrs_int+"*.nc") ; file paths
    f_OSHD = addfile(fils_precip_OSHD,"r")
    precip_sum_OSHD=f_OSHD->precip
    oshd_plot = precip_sum_OSHD(lat|:,lon|:,year|0)  ; needed in order to copy lat lon later
    copy_VarAtts(copy_latlon, oshd_plot)
    copy_VarCoords_1(copy_latlon,oshd_plot)
    oshd_plot=where(ismissing_oshd.eq.0,oshd_plot,oshd_plot@_FillValue)

    fils_precip_crujra = systemfunc ("ls "+bf_precip+"yr_"+"*CRUJRA*"+yrs_int+"*.nc") ; file paths
    print(fils_precip_crujra)
    f_crujra = addfile(fils_precip_crujra,"r")
    precip_sum_crujra=f_crujra->precip
    crujra_plot = precip_sum_crujra(lat|:,lon|:,year|0)  ; needed in order to copy lat lon later
    copy_VarAtts(copy_latlon, crujra_plot)
    copy_VarCoords_1(copy_latlon,crujra_plot)

    combi_yr=array_append_record(crujra_plot,oshd_plot,0)

    diff_sum = oshd_plot
    diff_sum = oshd_plot - crujra_plot

    chg_perc        =  oshd_plot                        ;-- copy variable and retain metadata
    chg_perc            =  abs(crujra_plot-oshd_plot)/where(crujra_plot.ne.0,crujra_plot,crujra_plot@_FillValue) * 100       ;-- compute the difference
    chg_perc@long_name  =  oshd_plot@long_name+" - %change" ;-- change long_name attribute

    wks = gsn_open_wks("png",pldir+"comp_yr_"+yrs_int)            ; open a png file
    plot = new(4,graphic)                          ; create a plot array

    res12@cnMinLevelValF      = (min(combi_yr))
    res12@cnMaxLevelValF      = (max(combi_yr))
    res12@cnLevelSpacingF     = ((max(combi_yr))-(min(combi_yr)))/12

    res3@cnMinLevelValF      = -(max(abs(diff_sum)))
    res3@cnMaxLevelValF      = (max(abs(diff_sum)))
    res3@cnLevelSpacingF     = ((max(abs(diff_sum)))*2)/20

    plot0 = gsn_csm_contour_map(wks,oshd_plot,res12)
    poly0 = gsn_add_shapefile_polylines(wks,plot0,dir+filename,lnres)
    plot1 = gsn_csm_contour_map(wks,crujra_plot,res12)
    poly1 = gsn_add_shapefile_polylines(wks,plot1,dir+filename,lnres)
    plot2 = gsn_csm_contour_map(wks,diff_sum,res3)
    poly2 = gsn_add_shapefile_polylines(wks,plot2,dir+filename,lnres)
    plot3 = gsn_csm_contour_map(wks,chg_perc,res4)
    poly3 = gsn_add_shapefile_polylines(wks,plot3,dir+filename,lnres)

    resP@gsnPanelMainString = "Yearly sum of precip. - "+yrs_int

    gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),resP)             ; now draw as one plot

    getvalues plot0
    "vpXF"               : vpx3     ; this is for positioning the common labelbar
    end getvalues
    getvalues plot0@contour
    "lbLabelStrings"     : labels
    end getvalues
    getvalues plot1
    "vpXF"               : vpx4     ; this is for positioning
    "vpWidthF"           : vpw4     ; the common lablebar
    end getvalues
    nlabels = dimsizes(labels)

    ;---Set some labelbar resources
    lbcolors_fill            = span_color_rgba(res12@cnFillPalette,nlabels+1)
    lbres@vpWidthF           = (vpx4 + vpw4) - vpx3
    gsn_labelbar_ndc(wks,nlabels,labels,vpx3,0.565,lbres)

    frame(wks)

    delete(oshd_plot)
    delete(crujra_plot)
    delete(combi_yr)
    delete(diff_sum)

end do

delete(OSHD_xAVG)
delete(CRUJRA_xAVG)
delete(f_crujra)
delete(X_cru)
delete(X)
delete(f)
;-----------------------------------------------------------------------
; YEARLY SUM/YEARLY avg - temp
;----------------------------------------------------------------------

fils = systemfunc ("ls "+bf_temp+"yr*OSHD*.nc") ; file paths
f    = addfiles (fils, "r")
X    = addfiles_GetVar(f,fils,"temp")
OSHD_xAVG = dim_avg_Wrap( X(lat|:,lon|:,year|:) )
copy_VarAtts(copy_latlon, OSHD_xAVG)
copy_VarCoords_1(copy_latlon,OSHD_xAVG)
;OSHD_xAVG@units="K"
OSHD_xAVG=where(ismissing_oshd.eq.0,OSHD_xAVG,OSHD_xAVG@_FillValue)

fils_cru = systemfunc ("ls "+bf_temp+"yr*CRUJRA*wL.nc") ; file paths
f_crujra    = addfiles (fils_cru, "r")
X_cru    = addfiles_GetVar(f_crujra,fils_cru,"temp")
CRUJRA_xAVG = dim_avg_Wrap( X_cru(lat|:,lon|:,year|:) )
copy_VarAtts(copy_latlon, CRUJRA_xAVG)
copy_VarCoords_1(copy_latlon,CRUJRA_xAVG)
CRUJRA_xAVG@units="K"

combi_yr_all=array_append_record(OSHD_xAVG,CRUJRA_xAVG,0)

diff_sum_all = CRUJRA_xAVG
diff_sum_all = OSHD_xAVG - CRUJRA_xAVG

chg_perc_all        =  OSHD_xAVG                        ;-- copy variable and retain metadata
chg_perc_all            =  abs(CRUJRA_xAVG-OSHD_xAVG)/where(CRUJRA_xAVG.ne.0,CRUJRA_xAVG,CRUJRA_xAVG@_FillValue) * 100       ;-- compute the difference

wks = gsn_open_wks("png",pldir_temp+"comp_all_yrs_wL")            ; open a png file
plot = new(4,graphic)

res12=res
res12@cnFillPalette = color_plot12
res12@cnMinLevelValF      = (min(combi_yr_all))
res12@cnMaxLevelValF      = (max(combi_yr_all))
res12@cnLevelSpacingF     = ((max(combi_yr_all))-(min(combi_yr_all)))/12

res12@lbLabelFontColor    = "white"                ; this makes the labelbar labels invisible
res12@lbBoxLinesOn        = False
res12@lbFillColors        = lbcolors_transparent   ; the colors are fully transparent

res3=res
res3@cnFillPalette = color_plot3
res3@lbLabelBarOn        = True
res3@cnLevelSpacingF     = ((max(abs(diff_sum_all)))*2)/20
res3@cnMinLevelValF      = -(max(abs(diff_sum_all)))
res3@cnMaxLevelValF      = (max(abs(diff_sum_all)))
res3@lbTitleString       = "~F33~D~F~"+OSHD_xAVG@units
res3@lbTitlePosition     = "Bottom"            ; title position
res3@lbTitleDirection    = "Across"            ; title direction

res4=res
res4@cnFillPalette = color_plot4
res4@lbLabelBarOn        = True
res4@cnMinLevelValF      = 5
res4@cnMaxLevelValF      = 95
res4@cnLevelSpacingF     = 10
res4@lbTitleString       = "%"
res4@lbTitlePosition            = "Bottom"            ; title position
res4@lbTitleDirection           = "Across"            ; title direction

plot0 = gsn_csm_contour_map(wks,OSHD_xAVG,res12)
poly0 = gsn_add_shapefile_polylines(wks,plot0,dir+filename,lnres)
plot1 = gsn_csm_contour_map(wks,CRUJRA_xAVG,res12)
poly1 = gsn_add_shapefile_polylines(wks,plot1,dir+filename,lnres)
plot2 = gsn_csm_contour_map(wks,diff_sum_all,res3)
poly2 = gsn_add_shapefile_polylines(wks,plot2,dir+filename,lnres)
plot3 = gsn_csm_contour_map(wks,chg_perc_all,res4)
poly3 = gsn_add_shapefile_polylines(wks,plot3,dir+filename,lnres)

resP                  = True                   ; modify the panel plot
resP@gsnFrame         = False                  ; don't advance the frame, so we can use gsn_text_ndc
resP@lbLabelFontHeightF = 0.012    ; set font height of Label Bar labels
resP@lbTitleFontHeightF= .015                 ; make title smaller
resP@gsnPanelBottom   = 0.01      ; shrink panel plot by setting bottom edge of plot
resP@gsnMaximize    = True
resP@gsnFrame       = False
resP@gsnPanelFigureStrings= (/"OSHD","CRUJRA","Diff (OSDH-CRUJRA","%change crujra->OSHD"/) ; add strings to panel
resP@gsnPanelMainString = "Yearly mean T - avg. 2014-2019"

gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),resP)             ; now draw as one plot

getvalues plot0
"vpXF"               : vpx3     ; this is for positioning the common labelbar
end getvalues
getvalues plot0@contour
"lbLabelStrings"     : labels
end getvalues
getvalues plot1
"vpXF"               : vpx4     ; this is for positioning
"vpWidthF"           : vpw4     ; the common lablebar
end getvalues
nlabels = dimsizes(labels)
txres               = True
txres@txFontHeightF = 0.015

;---Set some labelbar resources
lbres@vpWidthF           = (vpx4 + vpw4) - vpx3
lbres@lbTitleString      = OSHD_xAVG@units           ; title string


gsn_labelbar_ndc(wks,nlabels,labels,vpx3,0.565,lbres)
frame(wks)
delete(f_crujra)

do yrs_int = 2014,2019
    print(yrs_int)
    fils_precip_OSHD = systemfunc ("ls "+bf_temp+"yr_"+"*OSHD*"+yrs_int+"*.nc") ; file paths
    f_OSHD = addfile(fils_precip_OSHD,"r")
    precip_sum_OSHD=f_OSHD->temp
    oshd_plot = precip_sum_OSHD(lat|:,lon|:,year|0)  ; needed in order to copy lat lon later
    copy_VarAtts(copy_latlon, oshd_plot)
    copy_VarCoords_1(copy_latlon,oshd_plot)
    oshd_plot=where(ismissing_oshd.eq.0,oshd_plot,oshd_plot@_FillValue)

    fils_precip_crujra = systemfunc ("ls "+bf_temp+"yr_"+"*CRUJRA*"+yrs_int+"*wL.nc") ; file paths
    f_crujra = addfile(fils_precip_crujra,"r")
    precip_sum_crujra=f_crujra->temp
    crujra_plot = precip_sum_crujra(lat|:,lon|:,year|0)  ; needed in order to copy lat lon later
    copy_VarAtts(copy_latlon, crujra_plot)
    copy_VarCoords_1(copy_latlon,crujra_plot)

    combi_yr=array_append_record(crujra_plot,oshd_plot,0)

    diff_sum = oshd_plot
    diff_sum = oshd_plot - crujra_plot

    chg_perc        =  oshd_plot                        ;-- copy variable and retain metadata
    chg_perc            =  abs(crujra_plot-oshd_plot)/where(crujra_plot.ne.0,crujra_plot,crujra_plot@_FillValue) * 100       ;-- compute the difference
    chg_perc@long_name  =  oshd_plot@long_name+" - %change" ;-- change long_name attribute

    wks = gsn_open_wks("png",pldir_temp+"comp_yr_wL_"+yrs_int)            ; open a png file
    plot = new(4,graphic)                          ; create a plot array

    res12@cnMinLevelValF      = (min(combi_yr))
    res12@cnMaxLevelValF      = (max(combi_yr))
    res12@cnLevelSpacingF     = ((max(combi_yr))-(min(combi_yr)))/12

    res3@cnMinLevelValF      = -(max(abs(diff_sum)))
    res3@cnMaxLevelValF      = (max(abs(diff_sum)))
    res3@cnLevelSpacingF     = ((max(abs(diff_sum)))*2)/20

    plot0 = gsn_csm_contour_map(wks,oshd_plot,res12)
    poly0 = gsn_add_shapefile_polylines(wks,plot0,dir+filename,lnres)
    plot1 = gsn_csm_contour_map(wks,crujra_plot,res12)
    poly1 = gsn_add_shapefile_polylines(wks,plot1,dir+filename,lnres)
    plot2 = gsn_csm_contour_map(wks,diff_sum,res3)
    poly2 = gsn_add_shapefile_polylines(wks,plot2,dir+filename,lnres)
    plot3 = gsn_csm_contour_map(wks,chg_perc,res4)
    poly3 = gsn_add_shapefile_polylines(wks,plot3,dir+filename,lnres)

    resP@gsnPanelMainString = "Yearly avg. T - "+yrs_int

    gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),resP)             ; now draw as one plot

    getvalues plot0
    "vpXF"               : vpx3     ; this is for positioning the common labelbar
    end getvalues
    getvalues plot0@contour
    "lbLabelStrings"     : labels
    end getvalues
    getvalues plot1
    "vpXF"               : vpx4     ; this is for positioning
    "vpWidthF"           : vpw4     ; the common lablebar
    end getvalues
    nlabels = dimsizes(labels)

    ;---Set some labelbar resources
    lbcolors_fill            = span_color_rgba(res12@cnFillPalette,nlabels+1)
    lbres@vpWidthF           = (vpx4 + vpw4) - vpx3
    gsn_labelbar_ndc(wks,nlabels,labels,vpx3,0.565,lbres)

    frame(wks)

    delete(oshd_plot)
    delete(crujra_plot)
    delete(combi_yr)
    delete(diff_sum)

end do

delete(OSHD_xAVG)
delete(CRUJRA_xAVG)
delete(f_crujra)
delete(X_cru)
delete(X)
delete(f)
;-----------------------------------------------------------------------
; YEARLY SUM/YEARLY avg - temp - no lapse
;----------------------------------------------------------------------

fils = systemfunc ("ls "+bf_temp+"yr*OSHD*.nc") ; file paths
f    = addfiles (fils, "r")
X    = addfiles_GetVar(f,fils,"temp")
OSHD_xAVG = dim_avg_Wrap( X(lat|:,lon|:,year|:) )
copy_VarAtts(copy_latlon, OSHD_xAVG)
copy_VarCoords_1(copy_latlon,OSHD_xAVG)
;OSHD_xAVG@units="K"
OSHD_xAVG=where(ismissing_oshd.eq.0,OSHD_xAVG,OSHD_xAVG@_FillValue)

fils_cru = systemfunc ("ls "+bf_temp+"yr*CRUJRA*nL.nc") ; file paths
f_crujra    = addfiles (fils_cru, "r")
X_cru    = addfiles_GetVar(f_crujra,fils_cru,"temp")
CRUJRA_xAVG = dim_avg_Wrap( X_cru(lat|:,lon|:,year|:) )
copy_VarAtts(copy_latlon, CRUJRA_xAVG)
copy_VarCoords_1(copy_latlon,CRUJRA_xAVG)
CRUJRA_xAVG@units="K"

combi_yr_all=array_append_record(OSHD_xAVG,CRUJRA_xAVG,0)

diff_sum_all = CRUJRA_xAVG
diff_sum_all = OSHD_xAVG - CRUJRA_xAVG

chg_perc_all        =  OSHD_xAVG                        ;-- copy variable and retain metadata
chg_perc_all            =  abs(CRUJRA_xAVG-OSHD_xAVG)/where(CRUJRA_xAVG.ne.0,CRUJRA_xAVG,CRUJRA_xAVG@_FillValue) * 100       ;-- compute the difference

wks = gsn_open_wks("png",pldir_temp+"comp_all_yrs_nL")            ; open a png file
plot = new(4,graphic)

res12=res
res12@cnFillPalette = color_plot12
res12@cnMinLevelValF      = (min(combi_yr_all))
res12@cnMaxLevelValF      = (max(combi_yr_all))
res12@cnLevelSpacingF     = ((max(combi_yr_all))-(min(combi_yr_all)))/12

res12@lbLabelFontColor    = "white"                ; this makes the labelbar labels invisible
res12@lbBoxLinesOn        = False
res12@lbFillColors        = lbcolors_transparent   ; the colors are fully transparent

res3=res
res3@cnFillPalette = color_plot3
res3@lbLabelBarOn        = True
res3@cnLevelSpacingF     = ((max(abs(diff_sum_all)))*2)/20
res3@cnMinLevelValF      = -(max(abs(diff_sum_all)))
res3@cnMaxLevelValF      = (max(abs(diff_sum_all)))
res3@lbTitleString       = "~F33~D~F~"+OSHD_xAVG@units
res3@lbTitlePosition     = "Bottom"            ; title position
res3@lbTitleDirection    = "Across"            ; title direction

res4=res
res4@cnFillPalette = color_plot4
res4@lbLabelBarOn        = True
res4@cnMinLevelValF      = 5
res4@cnMaxLevelValF      = 95
res4@cnLevelSpacingF     = 10
res4@lbTitleString       = "%"
res4@lbTitlePosition            = "Bottom"            ; title position
res4@lbTitleDirection           = "Across"            ; title direction

plot0 = gsn_csm_contour_map(wks,OSHD_xAVG,res12)
poly0 = gsn_add_shapefile_polylines(wks,plot0,dir+filename,lnres)
plot1 = gsn_csm_contour_map(wks,CRUJRA_xAVG,res12)
poly1 = gsn_add_shapefile_polylines(wks,plot1,dir+filename,lnres)
plot2 = gsn_csm_contour_map(wks,diff_sum_all,res3)
poly2 = gsn_add_shapefile_polylines(wks,plot2,dir+filename,lnres)
plot3 = gsn_csm_contour_map(wks,chg_perc_all,res4)
poly3 = gsn_add_shapefile_polylines(wks,plot3,dir+filename,lnres)

resP                  = True                   ; modify the panel plot
resP@gsnFrame         = False                  ; don't advance the frame, so we can use gsn_text_ndc
resP@lbLabelFontHeightF = 0.012    ; set font height of Label Bar labels
resP@lbTitleFontHeightF= .015                 ; make title smaller
resP@gsnPanelBottom   = 0.01      ; shrink panel plot by setting bottom edge of plot
resP@gsnMaximize    = True
resP@gsnFrame       = False
resP@gsnPanelFigureStrings= (/"OSHD","CRUJRA","Diff (OSDH-CRUJRA","%change crujra->OSHD"/) ; add strings to panel
resP@gsnPanelMainString = "Yearly mean T - avg. 2014-2019 (nL)"

gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),resP)             ; now draw as one plot

getvalues plot0
"vpXF"               : vpx3     ; this is for positioning the common labelbar
end getvalues
getvalues plot0@contour
"lbLabelStrings"     : labels
end getvalues
getvalues plot1
"vpXF"               : vpx4     ; this is for positioning
"vpWidthF"           : vpw4     ; the common lablebar
end getvalues
nlabels = dimsizes(labels)
txres               = True
txres@txFontHeightF = 0.015

;---Set some labelbar resources
lbres@vpWidthF           = (vpx4 + vpw4) - vpx3
lbres@lbTitleString      = OSHD_xAVG@units           ; title string


gsn_labelbar_ndc(wks,nlabels,labels,vpx3,0.565,lbres)
frame(wks)
delete(f_crujra)

do yrs_int = 2014,2019
    print(yrs_int)
    fils_precip_OSHD = systemfunc ("ls "+bf_temp+"yr_"+"*OSHD*"+yrs_int+"*.nc") ; file paths
    f_OSHD = addfile(fils_precip_OSHD,"r")
    precip_sum_OSHD=f_OSHD->temp
    oshd_plot = precip_sum_OSHD(lat|:,lon|:,year|0)  ; needed in order to copy lat lon later
    copy_VarAtts(copy_latlon, oshd_plot)
    copy_VarCoords_1(copy_latlon,oshd_plot)
    oshd_plot=where(ismissing_oshd.eq.0,oshd_plot,oshd_plot@_FillValue)

    fils_precip_crujra = systemfunc ("ls "+bf_temp+"yr_"+"*CRUJRA*"+yrs_int+"*nL.nc") ; file paths
    f_crujra = addfile(fils_precip_crujra,"r")
    precip_sum_crujra=f_crujra->temp
    crujra_plot = precip_sum_crujra(lat|:,lon|:,year|0)  ; needed in order to copy lat lon later
    copy_VarAtts(copy_latlon, crujra_plot)
    copy_VarCoords_1(copy_latlon,crujra_plot)

    combi_yr=array_append_record(crujra_plot,oshd_plot,0)

    diff_sum = oshd_plot
    diff_sum = oshd_plot - crujra_plot

    chg_perc        =  oshd_plot                        ;-- copy variable and retain metadata
    chg_perc            =  abs(crujra_plot-oshd_plot)/where(crujra_plot.ne.0,crujra_plot,crujra_plot@_FillValue) * 100       ;-- compute the difference
    chg_perc@long_name  =  oshd_plot@long_name+" - %change" ;-- change long_name attribute

    wks = gsn_open_wks("png",pldir_temp+"comp_yr_nL_"+yrs_int)            ; open a png file
    plot = new(4,graphic)                          ; create a plot array

    res12@cnMinLevelValF      = (min(combi_yr))
    res12@cnMaxLevelValF      = (max(combi_yr))
    res12@cnLevelSpacingF     = ((max(combi_yr))-(min(combi_yr)))/12

    res3@cnMinLevelValF      = -(max(abs(diff_sum)))
    res3@cnMaxLevelValF      = (max(abs(diff_sum)))
    res3@cnLevelSpacingF     = ((max(abs(diff_sum)))*2)/20

    plot0 = gsn_csm_contour_map(wks,oshd_plot,res12)
    poly0 = gsn_add_shapefile_polylines(wks,plot0,dir+filename,lnres)
    plot1 = gsn_csm_contour_map(wks,crujra_plot,res12)
    poly1 = gsn_add_shapefile_polylines(wks,plot1,dir+filename,lnres)
    plot2 = gsn_csm_contour_map(wks,diff_sum,res3)
    poly2 = gsn_add_shapefile_polylines(wks,plot2,dir+filename,lnres)
    plot3 = gsn_csm_contour_map(wks,chg_perc,res4)
    poly3 = gsn_add_shapefile_polylines(wks,plot3,dir+filename,lnres)

    resP@gsnPanelMainString = "Yearly avg. T (nL) - "+yrs_int

    gsn_panel(wks,(/plot0,plot1,plot2,plot3/),(/2,2/),resP)             ; now draw as one plot

    getvalues plot0
    "vpXF"               : vpx3     ; this is for positioning the common labelbar
    end getvalues
    getvalues plot0@contour
    "lbLabelStrings"     : labels
    end getvalues
    getvalues plot1
    "vpXF"               : vpx4     ; this is for positioning
    "vpWidthF"           : vpw4     ; the common lablebar
    end getvalues
    nlabels = dimsizes(labels)

    ;---Set some labelbar resources
    lbcolors_fill            = span_color_rgba(res12@cnFillPalette,nlabels+1)
    lbres@vpWidthF           = (vpx4 + vpw4) - vpx3
    gsn_labelbar_ndc(wks,nlabels,labels,vpx3,0.565,lbres)

    frame(wks)

    delete(oshd_plot)
    delete(crujra_plot)
    delete(combi_yr)
    delete(diff_sum)

end do


end
