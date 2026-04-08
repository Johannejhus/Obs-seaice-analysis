load "ncl_funcs.ncl"

begin

;===============================
; Spatial domain
;===============================
latmin = -70.
latmax = -67.
lonmin = 0.
lonmax = 30.

; ice-sector longitude in 0–360
xmin = 0.
xmax = 30.

;===============================
; Time settings
;===============================
clst = 1982
clen = 2008

yrst = 1982
yren = 2025

;===============================
; SALINITY data (ORAS5)
;===============================
sNam     = "vosaline"
osrc     = "ORAS5"
diri_sal = "/g/data/gv90/wrh581/ORAS5/"
; CHANGE THIS to the correct file once you confirm it
fili_sal = "vosaline.nc"

;===============================
; Ice data
;===============================
diri_ice = "/g/data/gv90/wrh581/NSIDC/"
fili_ice = "NSIDC_SH_icearea_v06r00_197811-202512.nc"
iNam     = "SIA"

;===============================
; Output
;===============================
dirp  = "plots/"
filp  = "ORAS5_upper50m_salinity_regionalSIE_timeseries"
pType = "png"

;===============================
; READ SALINITY
;===============================
fo = addfile(diri_sal+fili_sal,"r")

;--------------------------------------------
; Read full salinity field
; Expecting dimensions like:
;   time, depth, lat, lon
; or
;   time_counter, deptht, y, x
;--------------------------------------------
salt_all = fo->$sNam$

print(" ")
print("Salinity variable summary:")
printVarSummary(salt_all)

;--------------------------------------------
; Try to identify coordinate names
; You may need to adjust these depending on ORAS5 file
;--------------------------------------------
time_name  = salt_all!0
depth_name = salt_all!1
lat_name   = salt_all!2
lon_name   = salt_all!3

print("Detected dimension names:")
print("time  = " + time_name)
print("depth = " + depth_name)
print("lat   = " + lat_name)
print("lon   = " + lon_name)

;--------------------------------------------
; Subset box first
;--------------------------------------------
salt_box = salt_all(:, :, {latmin:latmax}, {lonmin:lonmax})

;--------------------------------------------
; Identify depths <= 50 m
;--------------------------------------------
z = salt_box&$depth_name$
zind = ind(z.le.50.)

if (.not.any(ismissing(zind))) then
  salt_upper = dim_avg_n_Wrap(salt_box(:,zind,:,:), 1)
else
  print("ERROR: no depth levels <= 50 m found.")
  exit
end if

;--------------------------------------------
; Area average over lat/lon
;--------------------------------------------
wgty = cos(deg2rad(salt_upper&$salt_upper!1$))
salt_upper_reg = wgt_areaave(salt_upper, wgty, 1., 0)

salt_upper_reg!0 = "time"

;--------------------------------------------
; Attach time
;--------------------------------------------
if (isvar("time_counter")) then
  salt_upper_reg&time = cd_calendar(fo->time_counter, -1)
else if (isvar("time")) then
  salt_upper_reg&time = cd_calendar(fo->time, -1)
else
  print("ERROR: could not find time variable in salinity file.")
  exit
end if
end if

; keep monthly salinity series
salt_mon = salt_upper_reg({yrst*100+1:yren*100+12})

;===============================
; SALINITY anomalies
;===============================
salt_anoms = salt_upper_reg({yrst*100+1:yren*100+12})
salt_clim  = new(12, typeof(salt_upper_reg))

do i = 0, 11
  stdat = clst*100 + i+1
  endat = clen*100 + i+1
  salt_clim(i) = avg(salt_upper_reg({stdat:endat:12}))
  salt_anoms(i::12) = salt_anoms(i::12) - salt_clim(i)
end do

;===============================
; READ REGIONAL ICE
;===============================
fi  = addfile(diri_ice+fili_ice, "r")
ICE = fi->$iNam$({yrst*100+1:yren*100+12},:)
ICE = ICE * 1.e-6

; correct missing values
if (yrst.lt.1987 .and. yren.ge.1988) then
  ICE({198712:198801},:) = ICE@_FillValue
end if

; regional integral
if (xmin.lt.xmax) then
  ice_reg_raw = dim_sum_Wrap(ICE(:,{xmin:xmax}))
else
  ice_reg_raw = dim_sum_Wrap(ICE(:,{xmin:}))
  ice_reg_raw = ice_reg_raw + (/ dim_sum(ICE(:,{:xmax})) /)
end if

ice_reg = ice_reg_raw

;===============================
; ICE monthly anomalies
;===============================
do m = 1, 12
  clim = avg(ice_reg_raw({clst*100+m:clen*100+m:12}))
  ice_reg({yrst*100+m::12}) = ice_reg_raw({yrst*100+m::12}) - clim
end do

sie = ice_reg

;===============================
; PLOT
;===============================
pType@wkWidth  = 2500
pType@wkHeight = 2500

wks = gsn_open_wks(pType, dirp+filp)

resSAL                       = True
resSAL@gsnDraw               = False
resSAL@gsnFrame              = False
resSAL@gsnMaximize           = True
resSAL@vpHeightF             = 0.5
resSAL@vpWidthF              = 1.

resSAL@tiYAxisString         = "Upper-ocean salinity anomaly (0-50 m)"
resSAL@gsnYRefLine           = 0.
resSAL@xyLineThicknessF      = 2.5
resSAL@xyLineColor           = "black"
resSAL@gsnBelowYRefLineColor = "deepskyblue"
resSAL@gsnAboveYRefLineColor = "orangered"
resSAL@tiMainString          = "Regional upper-ocean salinity and sea-ice anomaly timeseries"

tdum = ispan(1,dimsizes(salt_anoms),1)

resSAL@tmXBMode              = "Explicit"
resSAL@tmXBValues            = tdum(::60)
resSAL@tmXBLabels            = ispan(yrst,yren,5)
resSAL@tmXBMinorValues       = tdum(::12)

resSAL@trYMaxF               = 0.3
resSAL@trYMinF               = -0.3

resSIE                  = True
resSIE@gsnDraw          = False
resSIE@gsnFrame         = False
resSIE@xyLineThicknessF = 7.
resSIE@xyLineColor      = "black"
resSIE@tiYAxisString    = "Regional SIA anomaly (10~S~6~N~ km~S~2~N~)"
resSIE@trYMaxF          = 3.
resSIE@trYMinF          = -3.
resSIE@trYReverse       = True

plot = gsn_csm_xy2(wks, tdum, salt_anoms, sie, resSAL, resSIE)

;========================================
; Add shaded periods
;========================================
pres = True
pres@gsEdgesOn   = False

nsh = 10
polys = new(nsh,graphic)

meas_years = (/2021, 2022, 2023, 2024, 2025/)

y0 = resSAL@trYMinF
y1 = resSAL@trYMaxF
ybox = (/y0, y0, y1, y1/)

k = 0
do i = 0, dimsizes(meas_years)-1
  yr = meas_years(i)
  prevyr = yr - 1

  ; Grey box: Jul-Aug in year before
  pres@gsFillColor = "gray85"
  xJul = (prevyr-yrst)*12 + 7
  xAug = (prevyr-yrst)*12 + 8
  xbox = (/xJul-0.5, xAug+0.5, xAug+0.5, xJul-0.5/)
  polys(k) = gsn_add_polygon(wks, plot, xbox, ybox, pres)
  k = k + 1

  ; Yellow box: Dec-Jan measurement season
  pres@gsFillColor = "khaki1"
  xDec = (prevyr-yrst)*12 + 12
  xJan = (yr-yrst)*12 + 1
  xbox = (/xDec-0.5, xJan+0.5, xJan+0.5, xDec-0.5/)
  polys(k) = gsn_add_polygon(wks, plot, xbox, ybox, pres)
  k = k + 1
end do

draw(plot)
frame(wks)

;=========================================================
; PRINT MEAN UPPER-OCEAN SALINITY FOR JANUARY AND WINTER
;=========================================================
print(" ")
print("======================================================")
print("MEAN UPPER-OCEAN SALINITY (0-50 m): JANUARY AND WINTER")
print("======================================================")

target_years_sal = (/2021,2022,2023,2024,2025/)

do i = 0, dimsizes(target_years_sal)-1
  yr = target_years_sal(i)
  prevyr = yr - 1

  sal_jan = salt_mon({yr*100+1})
  sal_winter = avg( salt_mon({prevyr*100+7:prevyr*100+8}) )

  print("Year " + yr)
  print("  January upper-ocean salinity         = " + sal_jan)
  print("  Previous winter upper-ocean salinity = " + sal_winter + "   (Jul-Aug " + prevyr + ")")
end do

;=========================================================
; PRINT SEA-ICE ANOMALY FOR JANUARY AND PREVIOUS WINTER
;=========================================================
print(" ")
print("======================================================")
print("SEA-ICE ANOMALY: JANUARY AND PREVIOUS WINTER (JUL-AUG)")
print("======================================================")

target_years_ice = (/2021,2022,2023,2024,2025/)

do i = 0, dimsizes(target_years_ice)-1
  yr = target_years_ice(i)
  prevyr = yr - 1

  ice_jan = sie({yr*100+1})
  ice_winter = avg( sie({prevyr*100+7:prevyr*100+8}) )

  print("Year " + yr)
  print("  January ice anomaly         = " + ice_jan)
  print("  Previous winter ice anomaly = " + ice_winter + "   (Jul-Aug " + prevyr + ")")
end do

end
