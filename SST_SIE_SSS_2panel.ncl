load "ncl_funcs.ncl"

begin

;===============================
; Spatial domain
;===============================
latmin = -70.
latmax = -55.
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
fili_ice = "NSIDC_SH_icearea_"+(/"v05r00_197811-202412.nc","v03r00_202501-202505.nc"/)
iNam = "SIA"

;===============================
; SSS data (ORAS5)
;===============================
diri_sss = "/g/data/gv90/wrh581/ORAS5/vosaline/"
fili_sss = systemfunc("ls "+diri_sss+"ORAS5_vosaline_monthly_SOcean_{"+str_join(ispan(yrst,yren,1),",")+"}.nc 2>/dev/null")
vNam_sss = "vosaline"

;===============================
; Output
;===============================
dirp  = "plots/"
filp  = "OISST_regionalSIE_SSS_2panel"
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

; SST anomalies
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

tice = fi[:]->time
if (typeof(tice).eq."int64") then
  ttmp = tice
  tice := toint(ttmp)
  copy_VarAtts(ttmp,tice)
  delete(ttmp)
end if

YYYYMM = cd_calendar(tice,-1)

ICE = fi[:]->$iNam$
ICE&time := YYYYMM

ICE = ICE({yrst*100+1:yren*100+12},:)
ICE = ICE * 1.e-6



if (xmin.lt.xmax) then
  ice_reg = dim_sum_Wrap(ICE(:,{xmin:xmax}))
else
  ice_reg = dim_sum_Wrap(ICE(:,{xmin:}))
  ice_reg = ice_reg + (/dim_sum(ICE(:,{:xmax}))/)
end if

ice_reg!0    = "time"
ice_reg&time = ICE&time

do m = 1, 12
  clim = avg(ice_reg({clst*100+m:clen*100+m:12}))
  ice_reg({yrst*100+m::12}) = ice_reg({yrst*100+m::12}) - clim
end do

sie = ice_reg

;===============================
; READ SSS (surface level)
;===============================
fs = addfiles(fili_sss,"r")

; assumes ORAS5 tripolar grid, top level = first depth
sss_tmp = fs[:]->$vNam_sss$(:,0,:,:1439)

sss_tmp!1 = "lat"
sss_tmp&lat := fs[0]->nav_lat(:,0)
sss_tmp!2 = "lon"
sss_tmp&lon := fs[0]->nav_lon(0,:1439)

sss_tmp := sss_tmp(:,{latmin:latmax},:)

if (lonmin.lt.lonmax) then
  sss_tmp := sss_tmp(:,:,{lonmin:lonmax})
else
  sss_tmp(:,:,{lonmax:lonmin}) = sss_tmp@_FillValue
end if

; approximate area weighting
dlat = 0.12
dlon = 0.25
dy   = yresCalc(dlat)
dx   = xresCalc(dlon, sss_tmp&lat)
A    = conform(sss_tmp(0,:,:), dx*dy, 0)

sss = wgt_areaave2(sss_tmp, A, 0)

sss!0 = "time"
sss&time = cd_calendar(fs[:]->time_counter, -1)

sss_anoms = sss({yrst*100+1:yren*100+12})
sss_clim  = new(12, typeof(sss))

do i = 0, 11
  stdat = clst*100 + i+1
  endat = clen*100 + i+1
  sss_clim(i) = avg(sss({stdat:endat:12}))
  sss_anoms(i::12) = sss_anoms(i::12) - sss_clim(i)
end do

;===============================
; COMMON X
;===============================
nt = dimsizes(sst_anoms)
tdum = ispan(1, nt, 1)

;===============================
; OPEN WORKSTATION
;===============================
pType@wkWidth  = 2500
pType@wkHeight = 3200
wks = gsn_open_wks(pType, dirp+filp)

plots = new(2,graphic)

;===============================
; TOP PANEL: SST + SIE
;===============================
resSST = True
resSST@gsnDraw               = False
resSST@gsnFrame              = False
resSST@vpHeightF             = 0.32
resSST@vpWidthF              = 0.82
resSST@tiMainString          = "Regional SST, SIE and SSS anomalies"
resSST@tiYAxisString         = "SST anomaly (°C)"
resSST@gsnYRefLine           = 0.
resSST@xyLineThicknessF      = 2.5
resSST@xyLineColor           = "black"
resSST@gsnBelowYRefLineColor = "deepskyblue"
resSST@gsnAboveYRefLineColor = "orangered"
resSST@tmXBMode              = "Explicit"
resSST@tmXBValues            = tdum(::60)
resSST@tmXBLabels            = ispan(yrst,yren,5)
resSST@tmXBMinorValues       = tdum(::12)
resSST@trYMaxF               = 1.
resSST@trYMinF               = -1.

resSIE = True
resSIE@gsnDraw          = False
resSIE@gsnFrame         = False
resSIE@xyLineThicknessF = 6.
resSIE@xyLineColor      = "black"
resSIE@tiYAxisString    = "Regional SIE anomaly (10~S~6~N~ km~S~2~N~)"
resSIE@trYMaxF          = 3.
resSIE@trYMinF          = -3.
resSIE@trYReverse       = True

plots(0) = gsn_csm_xy2(wks, tdum, sst_anoms, sie, resSST, resSIE)

;===============================
; BOTTOM PANEL: SSS
;===============================
resSSS = True
resSSS@gsnDraw           = False
resSSS@gsnFrame          = False
resSSS@vpHeightF         = 0.22
resSSS@vpWidthF          = 0.82
resSSS@xyLineThicknessF  = 3.
resSSS@xyLineColor       = "navy"
resSSS@tiYAxisString     = "SSS anomaly"
resSSS@gsnYRefLine       = 0.
resSSS@tmXBMode          = "Explicit"
resSSS@tmXBValues        = tdum(::60)
resSSS@tmXBLabels        = ispan(yrst,yren,5)
resSSS@tmXBMinorValues   = tdum(::12)

plots(1) = gsn_csm_xy(wks, tdum, sss_anoms, resSSS)

;===============================
; PANEL
;===============================
pres = True
pres@gsnMaximize = True
pres@gsnPanelYWhiteSpacePercent = 5.

gsn_panel(wks, plots, (/2,1/), pres)

end
