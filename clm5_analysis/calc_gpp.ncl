load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/cd_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/contrib/ut_string.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

;-----------------------------------------------------------------------
; make yearly and monthly aggregates of daily CLM5 GPP (FPSN)
;----------------------------------------------------------------------

setfileoption("nc","Format","NetCDF4Classic")

month_abbr = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
              "Oct","Nov","Dec"/)

bf_file_in = "/home/lud11/malle/CLM5_CH/"
g_carbon=12.0107;

reg_nm = (/"CRUJRA_FILES","CRUJRA_FILES_OLD","OSHD_FILES","OSHD_FILES_OLD","CRUJRA_FILES_noLapse","CRUJRA_FILES_noLapse_OLD"/)
var_nm = (/"CRUJRA","CRUJRA","OSHD","OSHD","CRUJRA_noLapse","CRUJRA_noLapse"/)

do i=0,5
    name_str=bf_file_in+reg_nm(i)+"/FPSN_SP_"+ var_nm(i)+".nc"
    print(name_str)
    file_in = addfile(name_str,"r")
    var_in_orig = file_in->DATA ;read in values
    copy_latlon = var_in_orig(lat|:,lon|:,time|1)
    printVarSummary(var_in_orig)

    var_in = var_in_orig
    var_in = var_in_orig*g_carbon*3600*24*0.000001 ;convert to gC per day from umol m-2 s-1

    MonthSum = calculate_monthly_values(var_in, "sum", 0, False)
    copy_VarAtts(copy_latlon, MonthSum)
    copy_VarCoords_1(copy_latlon,MonthSum)
    MonthSum@units="gC m-2 mo-1"
    printVarSummary(MonthSum)

    yearly_sum = month_to_annual(MonthSum(time|0:119,lat|:,lon|:),0)
    yearly_sum&year = ispan(2010,2019,1)
    copy_VarAtts(copy_latlon, yearly_sum)
    yearly_sum@units="gC m-2 y-1"
    printVarSummary(yearly_sum)

    outfile_yr= bf_file_in+reg_nm(i)+"/FPSN_SP_YR_SUM.nc"
    print ("Write output file year: " + outfile_yr)
    if (isfilepresent (outfile_yr)) then ; overwrite any previous file
        system ("rm " + outfile_yr)
    end if
    out = addfile (outfile_yr, "c")      ; create new output file
    copy_VarAtts (file_in, out)     ; write global attributes
    out->DATA = yearly_sum          ; write the main data array last,to improve speed
    delete(out)

    outfile_mo= bf_file_in+reg_nm(i)+"/FPSN_SP_MONTH_SUM.nc"
    print ("Write output file month: " + outfile_mo)
    if (isfilepresent (outfile_mo)) then ; overwrite any previous file
        system ("rm " + outfile_mo)
    end if
    out = addfile (outfile_mo, "c")      ; create new output file
    copy_VarAtts (file_in, out)     ; write global attributes
    out->DATA = MonthSum          ; write the main data array last,to improve speed
    delete(out)
    delete(file_in)
    delete(var_in)
end do

delete(var_in_orig)
delete(MonthSum)
delete(yearly_sum)
delete(reg_nm)
delete(copy_latlon)

; now coarse scale
reg_nm = (/"CRUJRA_FILES_025deg_cru_new","CRUJRA_FILES_05deg_cru_new","OSHD_FILES_025_new","OSHD_FILES_05_new"/)

do i=0,3
    name_str=bf_file_in+reg_nm(i)+"/FPSN.nc"
    print(name_str)
    file_in = addfile(name_str,"r")
    var_in_orig = file_in->DATA ;convert to mm/day which can then be summed to daily mm
    copy_latlon = var_in_orig(lat|:,lon|:,time|0)
    printVarSummary(var_in_orig)

    var_in = var_in_orig
    var_in = var_in_orig*g_carbon*3600*24*0.000001; ;convert to gC per day

    MonthSum = calculate_monthly_values(var_in, "sum", 0, False)
    copy_VarAtts(copy_latlon, MonthSum)
    copy_VarCoords_1(copy_latlon,MonthSum)
    MonthSum@units="gC m-2 mo-1"
    printVarSummary(MonthSum)

    yearly_sum = month_to_annual(MonthSum(time|0:119,lat|:,lon|:),0)
    yearly_sum&year = ispan(2010,2019,1)
    copy_VarAtts(copy_latlon, yearly_sum)
    ;copy_VarCoords_1(copy_latlon,yearly_sum)
    yearly_sum@units="gC m-2 y-1"
    printVarSummary(yearly_sum)

    outfile_yr= bf_file_in+reg_nm(i)+"/FPSN_SP_YR_SUM.nc"
    print ("Write output file: " + outfile_yr)
    if (isfilepresent (outfile_yr)) then ; overwrite any previous file
            system ("rm " + outfile_yr)
    end if
    out = addfile (outfile_yr, "c")      ; create new output file
    copy_VarAtts (file_in, out)     ; write global attributes
    out->DATA = yearly_sum          ; write the main data array last,to improve speed
    delete(out)

    outfile_mo= bf_file_in+reg_nm(i)+"/FPSN_SP_MONTH_SUM.nc"
    print ("Write output file: " + outfile_mo)
    if (isfilepresent (outfile_mo)) then ; overwrite any previous file
            system ("rm " + outfile_mo)
    end if
    out = addfile (outfile_mo, "c")      ; create new output file
    copy_VarAtts (file_in, out)     ; write global attributes
    out->DATA = MonthSum          ; write the main data array last,to improve speed
    delete(out)
    delete(file_in)
    delete(var_in)
    delete(var_in_orig)
    delete(copy_latlon)
    delete(MonthSum)
    delete(yearly_sum)
end do

end