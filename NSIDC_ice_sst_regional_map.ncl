;=========================================================
; NSIDC_ice_sst_regional_map.ncl
;
; Makes one Southern Ocean map of:
;   - SST anomaly (or absolute SST)
;   - monthly/mean sea ice concentration
;   - climatological ice edge
;   - selected period ice edge
;   - study region box
;
; Choose either:
;   1) January map for one year
;   2) Previous winter mean map (Jul-Aug)
;=========================================================

begin

;=========================================================
; USER SETTINGS
;=========================================================

  ;---------------------------------
  ; Choose map type
  ;---------------------------------
  map_type = "january"     ; options: "january" or "winter"

  ;---------------------------------
  ; Target year
  ; january: plots Jan of this year
  ; winter : plots Jul-Aug of previous year
  ;---------------------------------
  targyr = 2023

  ;---------------------------------
  ; Plot SST as anomaly or absolute
  ;---------------------------------
  isAnom = True

  ;---------------------------------
  ; Climatology period
  ;---------------------------------
  clst = 1982
  clen = 2011

  ;---------------------------------
  ; Max latitude shown on polar map
  ;---------------------------------
  latmax_plot = -40.

  ;---------------------------------
  ; Study box to mark on map
  ;---------------------------------
  box_lon = (/ -10., 30., 30., -10., -10. /)
  box_lat = (/ -70., -70., -60., -60., -70. /)

;=========================================================
; PATHS / FILENAMES
;=========================================================

  dirb = "/g/data/gv90/wrh581/"

  ;-------------------------
  ; Sea-ice data
  ;-------------------------
  isrc  = "NSIDC"
  diroi = dirb + isrc + "/monthly_sic/"
  filoi = diroi + "sic_pss25_197811-202512_v06r00.nc"
  oiNam = "cdr_seaice_conc_monthly"

  ;-------------------------
  ; SST data
  ;-------------------------
  osrc = "OISSTv2"
  diro = dirb + "/OISST/"
  filo = "OIv2_SST_1982-2026.nc"
  oNam = "sst"

  ;-------------------------
  ; Output
  ;-------------------------
  dirp  = "plots/"
  pType = "png"

;=========================================================
; DEFINE TIME RANGE
;=========================================================

  if (map_type.eq."january") then
    stdat = targyr*100 + 1
    endat = targyr*100 + 1
    tlabel = "January " + targyr
    filp = isrc + "_regional_" + osrc + "_january_map_" + targyr
  end if

  if (map_type.eq."winter") then
    stdat = (targyr-1)*100 + 7
    endat = (targyr-1)*100 + 8
    tlabel = "Jul-Aug " + (targyr-1) + " (winter before " + targyr + ")"
    filp = isrc + "_regional_" + osrc + "_winter_map_" + targyr
  end if

  print("Making map for: " + tlabel)
  print("Start date: " + stdat)
  print("End date  : " + endat)

;=========================================================
; READ SEA ICE
;=========================================================

  fi = addfile(filoi, "r")

  itime = fi->time
  if (typeof(itime).eq."int64") then
    ttmp  = itime
    itime := toint(ttmp)
    copy_VarAtts(ttmp, itime)
    delete(ttmp)
  end if
  itime := cd_calendar(itime, -1)

  sic = fi->$oiNam$    ; (time,y,x)

  if (isbyte(sic).or.isubyte(sic)) then
    sic := byte2flt(sic)
  end if

  sic&time = itime

  ; selected period SIC
  sic_sel = sic({stdat:endat}, :, :)

  ; mean over selected months
  if (dimsizes(sic_sel(:,0,0)).gt.1) then
    sic_plot = dim_avg_n_Wrap(sic_sel, 0)
  else
    sic_plot = sic_sel(0,:,:)
  end if

  ; monthly climatological SIC for same month(s)
  if (map_type.eq."january") then
    clim_ind = ind(itime%100.eq.1 .and. itime/100.ge.clst .and. itime/100.le.clen)
    sic_clim = dim_avg_n_Wrap(sic(clim_ind,:,:), 0)
  end if

  if (map_type.eq."winter") then
    clim_ind1 = ind(itime%100.eq.7 .and. itime/100.ge.clst .and. itime/100.le.clen)
    clim_ind2 = ind(itime%100.eq.8 .and. itime/100.ge.clst .and. itime/100.le.clen)

    sic_clim_jul = dim_avg_n_Wrap(sic(clim_ind1,:,:), 0)
    sic_clim_aug = dim_avg_n_Wrap(sic(clim_ind2,:,:), 0)
    sic_clim = (sic_clim_jul + sic_clim_aug)/2.
  end if

  ice_lat2d = fi->latitude
  ice_lon2d = fi->longitude

;=========================================================
; READ SST
;=========================================================

  fo = addfile(diro+filo, "r")

  sst_all = fo->$oNam$(:, {:latmax_plot}, :)
  sst_all&time := cd_calendar(fo->time, -1)

  sst_sel = sst_all({stdat:endat}, :, :)

  if (dimsizes(sst_sel(:,0,0)).gt.1) then
    sst_plot = dim_avg_n_Wrap(sst_sel, 0)
  else
    sst_plot = sst_sel(0,:,:)
  end if

  if (isAnom) then

    if (map_type.eq."january") then
      clim_ind = ind(sst_all&time%100.eq.1 .and. sst_all&time/100.ge.clst .and. sst_all&time/100.le.clen)
      sst_clim = dim_avg_n_Wrap(sst_all(clim_ind,:,:), 0)
    end if

    if (map_type.eq."winter") then
      clim_ind1 = ind(sst_all&time%100.eq.7 .and. sst_all&time/100.ge.clst .and. sst_all&time/100.le.clen)
      clim_ind2 = ind(sst_all&time%100.eq.8 .and. sst_all&time/100.ge.clst .and. sst_all&time/100.le.clen)

      sst_clim_jul = dim_avg_n_Wrap(sst_all(clim_ind1,:,:), 0)
      sst_clim_aug = dim_avg_n_Wrap(sst_all(clim_ind2,:,:), 0)
      sst_clim = (sst_clim_jul + sst_clim_aug)/2.
    end if

    sst_plot = sst_plot - sst_clim
  end if

;=========================================================
; ATTACH ICE GRID TO CLIM / SELECTED ICE
;=========================================================

  sic_plot@lat2d = ice_lat2d
  sic_plot@lon2d = ice_lon2d
  sic_clim@lat2d = ice_lat2d
  sic_clim@lon2d = ice_lon2d

;=========================================================
; OPEN WORKSTATION
;=========================================================

  wks = gsn_open_wks(pType, dirp+filp)

;=========================================================
; SST RESOURCES
;=========================================================

  res                         = True
  res@gsnDraw                 = False
  res@gsnFrame                = False
  res@gsnMaximize             = True
  res@gsnPolar                = "SH"
  res@mpMaxLatF               = latmax_plot

  res@cnLinesOn               = False
  res@cnFillOn                = True
  res@cnFillDrawOrder         = "PreDraw"
  res@cnFillMode              = "RasterFill"
  res@cnInfoLabelOn           = False
  res@cnLineLabelsOn          = False
  res@lbLabelBarOn            = True
  res@gsnLeftString           = ""
  res@tiMainString            = tlabel

  if (isAnom) then
    res@gsnRightString        = "(" + clst + "-" + clen + " climatology)"
    res@lbTitleString         = "SST anomaly (°C)"
    res@lbTitleFontHeightF    = 0.015
    res@lbTitlePosition       = "Bottom"
    res@cnLevelSelectionMode  = "ManualLevels"
    res@cnMaxLevelValF        = 3.
    res@cnMinLevelValF        = -3.
    res@cnLevelSpacingF       = 0.3
    res@cnFillPalette         = "cmocean_balance"
  else
    res@gsnRightString        = ""
    res@lbTitleString         = "SST (°C)"
    res@lbTitleFontHeightF    = 0.015
    res@lbTitlePosition       = "Bottom"
    res@cnLevelSelectionMode  = "ManualLevels"
    res@cnMaxLevelValF        = 1.
    res@cnMinLevelValF        = -1.5
    res@cnLevelSpacingF       = 0.1
    res@cnFillPalette         = "cmocean_thermal"
  end if

;=========================================================
; ICE EDGE RESOURCES
;=========================================================

  ires                        = True
  ires@gsnDraw                = False
  ires@gsnFrame               = False
  ires@trGridType             = "TriangularMesh"
  ires@cnLineLabelsOn         = False
  ires@cnInfoLabelOn          = False
  ires@cnLineThicknessF       = 4.
  ires@cnLevelSelectionMode   = "ExplicitLevels"
  ires@cnLevels               = 0.15
  ires@gsnLeftString          = ""
  ires@gsnRightString         = ""

;=========================================================
; SIC SHADE RESOURCES
;=========================================================

;=========================================================
; SIC SHADE RESOURCES
;=========================================================

  cres                        = True
  cres@gsnDraw                = False
  cres@gsnFrame               = False
  cres@trGridType             = "TriangularMesh"
  cres@cnLineLabelsOn         = False
  cres@cnInfoLabelOn          = False
  cres@gsnLeftString          = ""
  cres@gsnRightString         = ""

  cres@cnLinesOn              = False
  cres@cnFillOn               = True
  cres@cnLevelSelectionMode   = "ExplicitLevels"
  cres@cnLevels               = (/0.15,0.3,0.4,0.5,0.6,0.7,0.8,0.9/)

  nlev                        = dimsizes(cres@cnLevels)
  imap                        = new((/nlev,4/), float)
  imap(:,0)                   = 1.
  imap(:,1)                   = 1.
  imap(:,2)                   = 1.
  imap(:,3)                   = fspan(0.,1.,nlev)

  cres@cnFillPalette          = imap

;=========================================================
; DRAW MAP
;=========================================================

  plot = gsn_csm_contour_map_polar(wks, sst_plot, res)

  ; climatological ice edge
  ires@cnLineColor = "black"
  iclim = gsn_csm_contour(wks, sic_clim, ires)
  overlay(plot, iclim)

  ; selected-period ice edge
  ires@cnLineColor = "violet"
  isel = gsn_csm_contour(wks, sic_plot, ires)
  overlay(plot, isel)

  ; sea-ice concentration shading
  aplot = gsn_csm_contour(wks, sic_plot, cres)
  overlay(plot, aplot)

  ; regional box
  lres = True
  lres@gsLineColor      = "black"
  lres@gsLineThicknessF = 5.0
  dum = gsn_add_polyline(wks, plot, box_lon, box_lat, lres)

  draw(plot)
  frame(wks)

  print("Saved plot to: " + dirp + filp + "." + pType)

end
