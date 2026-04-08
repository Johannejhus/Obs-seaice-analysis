load "ncl_funcs.ncl"

begin

;===============================
; Spatial domain
;===============================
latmin = -70.
latmax = -60.
lonmin = -10.
lonmax = 30.

; ice-sector longitude in 0–360
xmin = 350.
xmax = 30.

;===============================
; Time settings
;===============================
clst = 1982
clen = 2008

yrst = 1982
yren = 2026      ; ice only available to 2025 in current files

;===============================
; SST data
;===============================
oNam  = "sst"
osrc  = "OISSTv2"
diri_sst = "/g/data/gv90/wrh581/OISST/"
fili_sst = "OIv2_SST_1982-2026.nc"

;===============================
; Ice data
;===============================
diri_ice = "/g/data/gv90/wrh581/NSIDC/"
fili_ice = "NSIDC_SH_icearea_"+(/ \
             "v05r00_197811-202412.nc", \
             "v03r00_202501-202603.nc"/)

iNam = "SIA"

;===============================
; Output
;===============================
dirp  = "plots/"
filp  = "OISST_regionalSIE_timeseries"
pType = "png"

;===============================
; READ SST
;===============================
fo    = addfile(diri_sst+fili_sst,"r")
sst   = fo->$oNam$(:,{latmin:latmax},{lonmin:lonmax})

wgty  = cos(deg2rad(sst&$sst!1$))
sst  := wgt_areaave(sst,wgty,1.,0)

sst!0    = "time"
sst&time = cd_calendar(fo->time, -1)

; keep monthly SST series
sst_mon = sst({yrst*100+1:2026*100+12})

;===============================
; SST anomalies
;===============================
sst_anoms = sst({yrst*100+1:yren*100+12})
sst_clim  = new(12, typeof(sst))

do i = 0, 11
  stdat = clst*100 + i+1
  endat = clen*100 + i+1
  sst_clim(i) = avg(sst({stdat:endat:12}))
  sst_anoms(i::12) = sst_anoms(i::12) - sst_clim(i)
end do

;===============================
; READ REGIONAL ICE
;===============================
fi = addfiles(diri_ice+fili_ice,"r")

ICE = fi[:]->$iNam$({yrst*100+1:yren*100+12},:)
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

resSST                       = True
resSST@gsnDraw               = False
resSST@gsnFrame              = False
resSST@gsnMaximize           = True
resSST@vpHeightF             = 0.5
resSST@vpWidthF              = 1.

resSST@tiYAxisString         = "SST anomaly (°C)"
resSST@gsnYRefLine           = 0.
resSST@xyLineThicknessF      = 2.5
resSST@xyLineColor           = "black"
resSST@gsnBelowYRefLineColor = "deepskyblue"
resSST@gsnAboveYRefLineColor = "orangered"
resSST@tiMainString          = "Regional SST and sea-ice anomaly timeseries"

tdum = ispan(1,dimsizes(sst_anoms),1)

resSST@tmXBMode              = "Explicit"
resSST@tmXBValues            = tdum(::60)
resSST@tmXBLabels            = ispan(yrst,yren,5)
resSST@tmXBMinorValues       = tdum(::12)

resSST@trYMaxF               = 1.
resSST@trYMinF               = -1.

resSIE                  = True
resSIE@gsnDraw          = False
resSIE@gsnFrame         = False
resSIE@xyLineThicknessF = 7.
resSIE@xyLineColor      = "black"
resSIE@tiYAxisString    = "Regional SIA anomaly (10~S~6~N~ km~S~2~N~)"
resSIE@trYMaxF          = 3.
resSIE@trYMinF          = -3.
resSIE@trYReverse       = True

plot = gsn_csm_xy2(wks, tdum, sst_anoms, sie, resSST, resSIE)

;========================================
; Add shaded months
;========================================
pres = True
pres@gsFillColor = "gray85"
pres@gsEdgesOn   = False

nsh = 8
polys = new(nsh,graphic)

; x = (year-yrst)*12 + month

xJan21 = (2021-yrst)*12 + 1
xJan22 = (2022-yrst)*12 + 1
xJan23 = (2023-yrst)*12 + 1
xJan24 = (2024-yrst)*12 + 1
xJan25 = (2025-yrst)*12 + 1

xJul20 = (2020-yrst)*12 + 7
xAug20 = (2020-yrst)*12 + 8

xJul21 = (2021-yrst)*12 + 7
xAug21 = (2021-yrst)*12 + 8

xJul22 = (2022-yrst)*12 + 7
xAug22 = (2022-yrst)*12 + 8

y0 = resSST@trYMinF
y1 = resSST@trYMaxF

; winter before Jan 2021
xbox = (/xJul20-0.5, xAug20+0.5, xAug20+0.5, xJul20-0.5/)
ybox = (/y0, y0, y1, y1/)
polys(0) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

; winter before Jan 2022
xbox = (/xJul21-0.5, xAug21+0.5, xAug21+0.5, xJul21-0.5/)
polys(1) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

; winter before Jan 2023
xbox = (/xJul22-0.5, xAug22+0.5, xAug22+0.5, xJul22-0.5/)
polys(2) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

; January measurement periods
pres@gsFillColor = "khaki1"

xbox = (/xJan21-0.5, xJan21+0.5, xJan21+0.5, xJan21-0.5/)
polys(3) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xJan22-0.5, xJan22+0.5, xJan22+0.5, xJan22-0.5/)
polys(4) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xJan23-0.5, xJan23+0.5, xJan23+0.5, xJan23-0.5/)
polys(5) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xJan24-0.5, xJan24+0.5, xJan24+0.5, xJan24-0.5/)
polys(6) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xJan25-0.5, xJan25+0.5, xJan25+0.5, xJan25-0.5/)
polys(7) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

draw(plot)
frame(wks)

;=========================================================
; PRINT MEAN SST FOR JANUARY AND PREVIOUS WINTER (JUL-AUG)
;=========================================================
print(" ")
print("==============================================")
print("MEAN SST: JANUARY AND PREVIOUS WINTER (JUL-AUG)")
print("==============================================")

target_years_sst = (/2021,2022,2023,2024,2025,2026/)

do i = 0, dimsizes(target_years_sst)-1
  yr = target_years_sst(i)
  prevyr = yr - 1

  sst_jan = sst_mon({yr*100+1})
  sst_winter = avg( sst_mon({prevyr*100+7:prevyr*100+8}) )

  print("Year " + yr)
  print("  January SST            = " + sst_jan)
  print("  Previous winter SST    = " + sst_winter + "   (Jul-Aug " + prevyr + ")")
end do

;=========================================================
; PRINT SEA-ICE ANOMALY FOR JANUARY AND PREVIOUS WINTER
;=========================================================
print(" ")
print("======================================================")
print("SEA-ICE ANOMALY: JANUARY AND PREVIOUS WINTER (JUL-AUG)")
print("======================================================")

; only possible for years supported by ice file
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

print(" ")
print("NOTE:")
print("Sea-ice file only extends to May 2025 in this script.")
print("So January 2026 and winter 2025 sea-ice anomalies cannot be calculated")
print("unless you add a newer NSIDC file.")

end
