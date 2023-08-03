;---------------------------------------------------------------------------
;
; NCL Program to calc evapotranspiration
;
; Usage:
;
;    ncl FILENAME.ncl
;
;---------------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
prog_id = "join.clm5_tiles.ncl"

begin
setfileoption ("nc", "Format", "64BitOffset")  ; use Netcdf large file format
header_space = 5000     ; pad size to optimize Netcdf-3 output
limit_2d_offset = 1.0e-08   ; offset limit for 2-D coordinates to trigger alignment warning

bf_OSHD="/media/malle/LaCie1/CLM5_SP_TILES/CRUJRA_FILES_05deg_cru_new/"
bf_crujra="/media/malle/LaCie1/CLM5_SP_TILES/CRUJRA_FILES_1deg_cru_new/"

infile_qvege_OSHD = bf_OSHD + "QVEGE.nc"
f_qvege_OSHD = addfile (infile_qvege_OSHD, "r")
QVEGE_OSHD=f_qvege_OSHD->DATA

infile_qvegt_OSHD = bf_OSHD + "QVEGT.nc"
f_qvegt_OSHD = addfile (infile_qvegt_OSHD, "r")
QVEGT_OSHD=f_qvegt_OSHD->DATA

infile_qsoil_OSHD = bf_OSHD + "QSOIL.nc"
f_qsoil_OSHD = addfile (infile_qsoil_OSHD, "r")
QSOIL_OSHD=f_qsoil_OSHD->DATA

;crujra
infile_qvege_CRUJRA = bf_crujra + "QVEGE.nc"
f_qvege_CRUJRA = addfile (infile_qvege_CRUJRA, "r")
QVEGE_CRUJRA=f_qvege_CRUJRA->DATA

infile_qvegt_CRUJRA = bf_crujra + "QVEGT.nc"
f_qvegt_CRUJRA = addfile (infile_qvegt_CRUJRA, "r")
QVEGT_CRUJRA=f_qvegt_CRUJRA->DATA

infile_qsoil_CRUJRA = bf_crujra + "QSOIL.nc"
f_qsoil_CRUJRA = addfile (infile_qsoil_CRUJRA, "r")
QSOIL_CRUJRA=f_qsoil_CRUJRA->DATA


;evapotranspiration=qflx_evap_soi + qflx_evap_can + qflx_tran_veg=QVEGT+QVEGE+QSOIL
;-------------------------------------------------------------------
;-- compute evapotranspiration and transpiration fraction
;-------------------------------------------------------------------
QFLX_EVAP_TOT_OSHD            =  QVEGE_OSHD                        ;-- copy variable and retain metadata
QFLX_EVAP_TOT_OSHD            =  QVEGE_OSHD+QVEGT_OSHD+QSOIL_OSHD                ;-- compute the difference
QFLX_EVAP_TOT_OSHD@long_name  = "Evapotranspiration (Qvege+qvegt+qsoil)" ;-- change long_name attribute

QFLX_EVAP_TOT_OSHD_no0            = QFLX_EVAP_TOT_OSHD
QFLX_EVAP_TOT_OSHD_no0            = where(QFLX_EVAP_TOT_OSHD.eq.0,QFLX_EVAP_TOT_OSHD@_FillValue,QFLX_EVAP_TOT_OSHD)

QFLX_EVAP_TOT_CRUJRA          =  QVEGE_CRUJRA                        ;-- copy variable and retain metadata
QFLX_EVAP_TOT_CRUJRA          =  QVEGE_CRUJRA+QVEGT_CRUJRA+QSOIL_CRUJRA    ;-- compute the difference
QFLX_EVAP_TOT_CRUJRA@long_name= "Evapotranspiration (Qvege+qvegt+qsoil)" ;-- change long_name attribute

QFLX_EVAP_TOT_CRUJRA_no0            = QFLX_EVAP_TOT_CRUJRA
QFLX_EVAP_TOT_CRUJRA_no0            = where(QFLX_EVAP_TOT_CRUJRA.eq.0,QFLX_EVAP_TOT_CRUJRA@_FillValue,QFLX_EVAP_TOT_CRUJRA)

TRANSP_FRAC_OSHD            =  QVEGT_OSHD                        ;-- copy variable and retain metadata
TRANSP_FRAC_OSHD            =  (QVEGT_OSHD/QFLX_EVAP_TOT_OSHD_no0)*100                ;-- compute the difference
TRANSP_FRAC_OSHD@long_name  = "Veg. Transpiration fraction" ;-- change long_name attribute
TRANSP_FRAC_OSHD@units      = "%" ;-- change long_name attribute

TRANSP_FRAC_CRUJRA          =  QVEGT_CRUJRA                        ;-- copy variable and retain metadata
TRANSP_FRAC_CRUJRA          =  (QVEGT_CRUJRA/QFLX_EVAP_TOT_CRUJRA_no0)*100    ;-- compute the difference
TRANSP_FRAC_CRUJRA@long_name= "Veg. Transpiration fraction" ;-- change long_name attribute
TRANSP_FRAC_CRUJRA@units      = "%" ;-- change long_name attribute

QSOIL_FRAC_OSHD            =  QSOIL_OSHD                        ;-- copy variable and retain metadata
QSOIL_FRAC_OSHD            =  (QSOIL_OSHD/QFLX_EVAP_TOT_OSHD_no0)*100          ;-- compute the difference
QSOIL_FRAC_OSHD@long_name  = "Soil evaporation fraction" ;-- change long_name attribute
QSOIL_FRAC_OSHD@units      = "%" ;-- change long_name attribute

QSOIL_FRAC_CRUJRA          =  QSOIL_CRUJRA                        ;-- copy variable and retain metadata
QSOIL_FRAC_CRUJRA          =  (QSOIL_CRUJRA/QFLX_EVAP_TOT_CRUJRA_no0)*100    ;-- compute the difference
QSOIL_FRAC_CRUJRA@long_name= "Soil evaporation fraction"   ;-- change long_name attribute
QSOIL_FRAC_CRUJRA@units      = "%" ;-- change long_name attribute

QVEGE_FRAC_OSHD            =  QVEGE_OSHD                        ;-- copy variable and retain metadata
QVEGE_FRAC_OSHD            =  (QVEGE_OSHD/QFLX_EVAP_TOT_OSHD_no0)*100               ;-- compute the difference
QVEGE_FRAC_OSHD@long_name  = "Veg. evaporation fraction" ;-- change long_name attribute
QVEGE_FRAC_OSHD@units      = "%" ;-- change long_name attribute


QVEGE_FRAC_CRUJRA          =  QVEGE_CRUJRA                        ;-- copy variable and retain metadata
QVEGE_FRAC_CRUJRA          =  (QVEGE_CRUJRA/QFLX_EVAP_TOT_CRUJRA_no0)*100    ;-- compute the difference
QVEGE_FRAC_CRUJRA@long_name= "Veg. evaporation fraction"   ;-- change long_name attribute
QSOIL_FRAC_CRUJRA@units      = "%" ;-- change long_name attribute
;-------------------------------------------------------------------
;-- define output files
;-------------------------------------------------------------------
outfile_evapo_OSHD=bf_OSHD + "QFLX_EVAP_TOT.nc"
outfile_evapo_CRUJRA=bf_crujra + "QFLX_EVAP_TOT.nc"
outfile_transFrac_OSHD=bf_OSHD + "TRANSP_FRAC.nc"
outfile_transFrac_CRUJRA=bf_crujra + "TRANSP_FRAC.nc"

outfile_QSoilFrac_OSHD=bf_OSHD + "QSOIL_FRAC.nc"
outfile_QSoilFrac_CRUJRA=bf_crujra + "QSOIL_FRAC.nc"
outfile_QVegeFrac_OSHD=bf_OSHD + "QVEGE_FRAC.nc"
outfile_QVegeFrac_CRUJRA=bf_crujra + "QVEGE_FRAC.nc"

if (isfilepresent (outfile_evapo_OSHD)) then   ; overwrite any previous file
  system ("rm " + outfile_evapo_OSHD)
end if
if (isfilepresent (outfile_evapo_CRUJRA)) then   ; overwrite any previous file
  system ("rm " + outfile_evapo_CRUJRA)
end if
if (isfilepresent (outfile_transFrac_OSHD)) then   ; overwrite any previous file
  system ("rm " + outfile_transFrac_OSHD)
end if
if (isfilepresent (outfile_transFrac_CRUJRA)) then   ; overwrite any previous file
  system ("rm " + outfile_transFrac_CRUJRA)
end if

if (isfilepresent (outfile_QSoilFrac_OSHD)) then   ; overwrite any previous file
  system ("rm " + outfile_QSoilFrac_OSHD)
end if
if (isfilepresent (outfile_QSoilFrac_CRUJRA)) then   ; overwrite any previous file
  system ("rm " + outfile_QSoilFrac_CRUJRA)
end if
if (isfilepresent (outfile_QVegeFrac_OSHD)) then   ; overwrite any previous file
  system ("rm " + outfile_QVegeFrac_OSHD)
end if
if (isfilepresent (outfile_QVegeFrac_CRUJRA)) then   ; overwrite any previous file
  system ("rm " + outfile_QVegeFrac_CRUJRA)
end if
;-------------------------------------------------------------------
;-- make output files
;-------------------------------------------------------------------
out = addfile (outfile_evapo_OSHD, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_OSHD, out)       ; write global attributes
out->DATA = QFLX_EVAP_TOT_OSHD            ; write the main data array last,to improve speed
;
delete(out)
out = addfile (outfile_evapo_CRUJRA, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_CRUJRA, out)       ; write global attributes
out->DATA = QFLX_EVAP_TOT_CRUJRA            ; write the main data array last,to improve speed

delete(out)
out = addfile (outfile_transFrac_OSHD, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_OSHD, out)       ; write global attributes
out->DATA = TRANSP_FRAC_OSHD            ; write the main data array last,to improve speed

delete(out)
out = addfile (outfile_transFrac_CRUJRA, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_CRUJRA, out)       ; write global attributes
out->DATA = TRANSP_FRAC_CRUJRA            ; write the main data array last,to improve speed

delete(out)
out = addfile (outfile_QSoilFrac_OSHD, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_OSHD, out)       ; write global attributes
out->DATA = QSOIL_FRAC_OSHD            ; write the main data array last,to improve speed
;
delete(out)
out = addfile (outfile_QSoilFrac_CRUJRA, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_CRUJRA, out)       ; write global attributes
out->DATA = QSOIL_FRAC_CRUJRA            ; write the main data array last,to improve speed

delete(out)
out = addfile (outfile_QVegeFrac_OSHD, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_OSHD, out)       ; write global attributes
out->DATA = QVEGE_FRAC_OSHD            ; write the main data array last,to improve speed

delete(out)
out = addfile (outfile_QVegeFrac_CRUJRA, "c")        ; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (infile_qvege_CRUJRA, out)       ; write global attributes
out->DATA = QVEGE_FRAC_CRUJRA            ; write the main data array last,to improve speed

end