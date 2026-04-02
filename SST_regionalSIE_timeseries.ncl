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
yren = 2024

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
             "v03r00_202501-202505.nc"/)

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
  ice_reg = dim_sum_Wrap(ICE(:,{xmin:xmax}))
else
  ice_reg = dim_sum_Wrap(ICE(:,{xmin:}))
  ice_reg = ice_reg + (/ dim_sum(ICE(:,{:xmax})) /)
end if

; monthly anomalies
do m = 1, 12
  clim = avg(ice_reg({clst*100+m:clen*100+m:12}))
  ice_reg({yrst*100+m::12}) = ice_reg({yrst*100+m::12}) - clim
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

; ICE overlay
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

nsh = 6
polys = new(nsh,graphic)

; x index helper:
; x = (year-yrst)*12 + month
; since tdum starts at 1

xJan23 = (2023-yrst)*12 + 1
xSep22 = (2022-yrst)*12 + 9

xJan24 = (2024-yrst)*12 + 1
xSep23 = (2023-yrst)*12 + 9

xJan25 = (2025-yrst)*12 + 1
xSep24 = (2024-yrst)*12 + 9

; use SST axis range for shading polygons
y0 = resSST@trYMinF
y1 = resSST@trYMaxF

; September before measurements
xbox = (/xSep22-0.5, xSep22+0.5, xSep22+0.5, xSep22-0.5/)
ybox = (/y0, y0, y1, y1/)
polys(0) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xSep23-0.5, xSep23+0.5, xSep23+0.5, xSep23-0.5/)
polys(1) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xSep24-0.5, xSep24+0.5, xSep24+0.5, xSep24-0.5/)
polys(2) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

; January measurements
pres@gsFillColor = "khaki1"

xbox = (/xJan23-0.5, xJan23+0.5, xJan23+0.5, xJan23-0.5/)
polys(3) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xJan24-0.5, xJan24+0.5, xJan24+0.5, xJan24-0.5/)
polys(4) = gsn_add_polygon(wks, plot, xbox, ybox, pres)

xbox = (/xJan25-0.5, xJan25+0.5, xJan25+0.5, xJan25-0.5/)
polys(5) = gsn_add_polygon(wks, plot, xbox, ybox, pres)
draw(plot)
frame(wks)

do yr = 2023, 2025
  sst_jan = sst({yr*100+1})
  print("Mean SST, January "+yr+":")
  print(sst_jan)
end do


end
