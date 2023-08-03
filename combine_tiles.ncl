;---------------------------------------------------------------------------
;
; NCL Program to combine tiles over the extent of Switzerland
;
; Usage:
;
; 1.  Make a text file with tile numbers to be included.

; 2. Check which variables are in the respective file (h0, h1,...)
;    and adjust vars text file as well as add_yr parameter below accordingly
;
; 3.  Adjust the input and output parameters at the start of the
; program below, as needed.
;
; 4.  Run this script with no command arguments:
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

;OSHD or CRUJRA?
met_data= "CRUJRA"

;which h? - loop over h0-h3
h_num_all=asciiread ("h_num_file", -1, "string")
hfiles = dimsizes (h_num_all)
do hi = 0, hfiles-1
 h_num = h_num_all(hi);which h?

 ;this is needed for overall swiss grid:
 coord_file = "/media/malle/LaCie1/SPIN_UP_TILES/ch_1km_OSHD_BGC_test.clm2.h0.2010-02-01-00000.nc" ; supplemental coordinate file
 surf_file = "/media/malle/LaCie/CLM5_input/surfdata_1km_CH_v3_hist_16pfts_Irrig_CMIP6_simyr2000_std25.nc" ; supplemental coordinate file

 ;overall basefolder
 bf_all="/media/malle/LaCie1/CLM5_SP_TILES/all_simulations_oldSurf/"
 ;together with the tiles the simulation name will be determined, needing the following:
 var1    = "ch_1km_"
 run_ID    = "_SP_"+met_data;CRUJRA or OSHD?
 add_yr = ".clm2."+h_num+".2010-01-01-00000.nc"
 ;followed by:
 add_file="/"

 ;where to save out files:
 bf_out= "/media/malle/LaCie1/CLM5_SP_TILES/"+met_data+"_FILES_OLD/"
 ;-----------------------------------------------------------
 ; 1.  Read tile numbers and variables of interest from list file.
 ;-----------------------------------------------------------
 tiles=asciiread ("list_tiles", -1, "string")
 nfiles = dimsizes (tiles)
 variables_all=asciiread ("vars_"+h_num, -1, "string")
 nvars = dimsizes (variables_all)

 ;write first varible out, the script will just use this for overall dimensions
 var        = variables_all(0)
 ;-------------------------------------------------------------
 ; 2.  Read master coordinate arrays from supplemental file.
 ;-------------------------------------------------------------
 print ("")
 print ("Read coordinate file: " + coord_file)
 cf = addfile (coord_file, "r")
 sf = addfile (surf_file, "r")

 lat_main = sf->LATIXY; read master 2-D coordinates
 lon_main = sf->LONGXY

 xmain    = cf->lon    ; read master 1-D coordinates
 ymain    = cf->lat    ; read master 1-D coordinates
 nyc = dimsizes (ymain)
 nxc = dimsizes (xmain)
 print ("  Coordinate grid size lat, lon = " + nyc + " x " + nxc)
 ngrid = nxc * nyc
 ;-------------------------------------------------------------------------
 ; 3.  Scan through all input files, and check alignment of coordinates.
 ;-------------------------------------------------------------------------
 print ("")
 print ("Scan through input files to check coordinates.")

 isub = new ((/ nxc, 2 /), integer)	; arrays to save tile subscripts
 jsub = new ((/ nyc, 2 /), integer)

 lat_max_offset  = new (nfiles, double)
 lon_max_offset  = new (nfiles, double)
 lat_count_exact = new (nfiles, double)
 lon_count_exact = new (nfiles, double)

 do va = 0, nvars-1
 var_desc = variables_all(va)
 outfile    = bf_out+var_desc +run_ID+".nc"

 do fi = 0, nfiles-1
   tile_int = tiles(fi)
   bf_full = bf_all+var1+tile_int+run_ID+add_file

   infile = bf_full + var1 + tile_int + run_ID + add_yr
   f = addfile (infile, "r")
   x1 = f->$var$&lon    ; read 1-D coordinate vectors
   y1 = f->$var$&lat

   i1 = ind (x1(0) .eq. xmain)   ; find lower left grid point within
   j1 = ind (y1(0) .eq. ymain)   ; the master coordinate grid
   i2 = i1 + dimsizes (x1) - 1   ; find upper right subscripts
   j2 = j1 + dimsizes (y1) - 1

   isub(fi,:) = (/ i1, i2 /)   ; save subscripts of current tile
   jsub(fi,:) = (/ j1, j2 /)   ; within master grid

   ; Also check for matched time coordinates.
   time = f->$var$&time    ; read time coords for current file
   nt = dimsizes (time)
   time_dim_name = time!0    ; time dimension name, for later

   if (fi .eq. 0) then
    sample = f->$var$(0,0,0)    ; capture attributes from first file
    time1  = time     ; save first file time coordinates
    ntimes = nt
  end if

  delete (time)
  delete (y1)
  delete (x1)
end do

imin = min (isub(:,0))		; compute min and max tile subscripts
imax = max (isub(:,1))		; within master grid
jmin = min (jsub(:,0))
jmax = max (jsub(:,1))

lat_min = min (lat_main(jmin:jmax, imin:imax))   ; actual ranges for
lat_max = max (lat_main(jmin:jmax, imin:imax))   ; output 2-D coordinates
lon_min = min (lon_main(jmin:jmax, imin:imax))
lon_max = max (lon_main(jmin:jmax, imin:imax))
;-----------------------------------------------------------
; 5.  Create super array.
;-----------------------------------------------------------
print ("")
print ("Create super array.")

nxout = imax - imin + 1		; output dimensions
nyout = jmax - jmin + 1

out_dims   = (/ ntimes, nyout, nxout /)

type_size  = sizeof (sample)
array_size = type_size * product (todouble (out_dims))

print ("  Data type             = " + typeof (sample) + " (" + type_size \
       + " bytes)")
print ("  Requested dimensions  = (" + ntimes + ", " + nyout + ", " \
       + nxout + ")")
print ("  Requested memory size = " + sprintf ("%0.0f bytes", array_size) \
       + sprintf (" (%0.2f Gb)", array_size / 1e9))

; If there is to be a memory overflow error, it will probably happen here.
xout = new (out_dims, typeof (sample), sample@_FillValue)
;-----------------------------------------------------------
; 6.  Overlay input arrays.
;-----------------------------------------------------------
print ("")
print ("Overlay input arrays.")

isub = isub(:,:) - imin		; convert all tile subscripts
jsub = jsub(:,:) - jmin		; within master grid, to overlay subscripts within output grid

 do fi = 0, nfiles-1
  tile_int = tiles(fi)
  bf_full = bf_all+var1+tile_int+run_ID+add_file
  infile = bf_full + var1 + tile_int + run_ID + add_yr
  ;print(infile)

  f = addfile (infile, "r")
  i1 = isub(fi,0)			; recover i, j subscripts for
  i2 = isub(fi,1)			; current tile within output grid
  j1 = jsub(fi,0)
  j2 = jsub(fi,1)

  ; First tile:  Overlay main array; copy all metadata for data var.
  if (fi .eq. 0) then
    var_int=f->$var_desc$
    printMinMax(var_int, False)
    if (dimsizes(dimsizes(var_int)).eq.4) then ;necessary for Albedo values...
      var_int2=dim_avg_n_Wrap(var_int,1)
          printMinMax(var_int2, False)
    else
      var_int2=var_int
    end if

    xout(:,j1:j2,i1:i2) = var_int2(0:ntimes-1,:,:)
    delete (var_int)
    delete (var_int2)

    ; All other tiles:  Overlay main array; use mask to preserve existing data.
  else
    var_int=f->$var_desc$
    if (dimsizes(dimsizes(var_int)).eq.4) then
      var_int2=dim_avg_n_Wrap(var_int,1)
    else
      var_int2=var_int
    end if
    xout(:,j1:j2,i1:i2) =  where (ismissing (xout(:,j1:j2,i1:i2)), \
                                  var_int2(0:ntimes-1,:,:), xout(:,j1:j2,i1:i2))	; function suppresses metadata copy
    delete(var_int)
    delete(var_int2)
  end if

end do
;-----------------------------------------------------------
; 7.  Fix up metadata for output file.
;-----------------------------------------------------------
; Prepare global attributes.
f = addfile (infile, "r")		; re-open first input file to copy
; misc. vars and metadata
gatts = True				; copy most global attributes
copy_VarAtts (f, gatts)		; indirect for NCL 6.0.0 workaround

if (isatt (gatts, "tileid")) then	; remove tile ID attribute,
  delete (gatts@tileid)		;    no longer singular
end if

time_stamp  = systemfunc ("date")	; add history attribute
new_history = time_stamp + ": Tiles joined by " + prog_id

if (isatt (gatts, "history")) then	; append to previous history attribute
  gatts@history = new_history + gatts@history
else
  gatts@history = new_history
end if

; Clean up the X and Y coordinate variables.
xout&lon = xmain(imin:imax)		; apply only to super array, because
xout&lat = ymain(jmin:jmax)		;   this is where the gaps might be

delete (xout&lon@_FillValue)		; remove vestigial attributes
delete (xout&lat@_FillValue)

; Prepare supplemental variables.
tiles!0 = "ntiles"			; include list of tiles as data variable

; Compute actual_range attributes.
xout&lon@actual_range    = (/ min (xout&lon), max (xout&lon) /)   ; coordinate vars
xout&lat@actual_range    = (/ min (xout&lat), max (xout&lat) /)
xout&time@actual_range = (/ min (xout&time), max (xout&time) /)

xmin = min (xout)
xmax = max (xout)
xout@actual_range      = (/ xmin,xmax /)   ; main data array

; Here lat_main and lon_main are used only as carriers for these output attributes.
lat_main@actual_range  = (/ lat_min,      lat_max  /)   ; 2-D coordinates
lon_main@actual_range  = (/ lon_min,      lon_max  /)
;-----------------------------------------------------------
; 8.  Write output file.
;-----------------------------------------------------------
print ("")
print ("Write output file: " + outfile)

if (isfilepresent (outfile)) then	; overwrite any previous file
  system ("rm " + outfile)
end if

out = addfile (outfile, "c")		; create new output file
setfileoption (out, "HeaderReserveSpace", header_space) ; optimize Netcdf-3 output; see NCL docs
copy_VarAtts (gatts, out)		; write global attributes

dummy_size = -1			; make the time dimension unlimited
opt_unlimited = True
filedimdef (out, time_dim_name, dummy_size, opt_unlimited)

out->lat = ymain(jmin:jmax)   ; write 2-D output coordinates
out->lon = xmain(imin:imax)
out->DATA = xout			; write the main data array last,to improve speed
print ("Done.")

delete(var_desc)
delete(outfile)
;delete(var_int)
;delete(var_int2)
end do
delete(variables_all)
end do


end

exit