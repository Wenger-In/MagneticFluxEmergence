; 读取一个文件夹下的Bt,Br,Bp文件，计算helicity flux, energy flux的时间序列

; written by Zheng Sun , 2022.9.4
compile_opt idl2

HARP_num = '1026'

files = '/Volumes/WD1/Helicity/SHARP_data/HARP' + HARP_num + '/'
Br_series = file_search(files, '*Br.fits')
Bt_series = file_search(files, '*Bt.fits')
Bp_series = file_search(files, '*Bp.fits')

read_sdo, Br_series, index_r, data_r
read_sdo, Bt_series, index_t, data_t
read_sdo, Bp_series, index_p, data_p

data_r = data_r[150 : 450, *, *]
data_t = data_t[150 : 450, *, *]
data_p = data_p[150 : 450, *, *]

nn = n_elements(index_r)
dHm_shear = fltarr(nn - 1)
dHm_emerge = fltarr(nn - 1)
dHm = fltarr(nn - 1)
Hm = fltarr(nn - 1)
dE = fltarr(nn - 1)
E = fltarr(nn - 1)
jz_total = fltarr(nn - 1)
Hc = fltarr(nn - 1)
time = fltarr(nn - 1)

DX = 1
DY = 1
; factor2=(0.5*725.*1000.*100.*0.5*725.*1000.*100.)^2  ; for DX=1 pixel (provided by QuanWang)
factor2 = 0.5 * 700.0 / 960. * 1e6 ; (provided by TingLi)

for i = 0, nn - 2 do begin
  tt1 = index_r[i].date_d$obs
  tt2 = index_r[i + 1].date_d$obs
  dd1 = strmid(tt1, 8, 2)
  dd2 = strmid(tt2, 8, 2)
  hh1 = strmid(tt1, 11, 2)
  hh2 = strmid(tt2, 11, 2)
  mm1 = strmid(tt1, 14, 2)
  mm2 = strmid(tt2, 14, 2)
  ss1 = strmid(tt1, 17, 2)
  ss2 = strmid(tt2, 17, 2)

  t_start = fix(dd1) * 24.0 * 3600.0 + fix(hh1) * 3600.0 + fix(mm1) * 60.0 + fix(ss1)
  t_stop = fix(dd2) * 24.0 * 3600.0 + fix(hh2) * 3600.0 + fix(mm2) * 60.0 + fix(ss2)

  day1 = float(dd1)
  day2 = float(dd2)
  day0 = day1 - day2
  if day0 gt 0. then t_stop = t_start + 720.

  ; ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; 用DAVE4VM计算

  time[i] = (t_start + t_stop) / 2.
  time[i] = time[i] / 3600.

  bx_start = data_p[*, *, i]
  by_start = -1. * data_t[*, *, i]
  bz_start = data_r[*, *, i]

  bx_stop = data_p[*, *, i + 1]
  by_stop = -1. * data_t[*, *, i + 1]
  bz_stop = data_r[*, *, i + 1]

  dobx = (bx_start + bx_stop) / 2.
  doby = (by_start + by_stop) / 2.
  dobz = (bz_start + bz_stop) / 2.

  dodave4vm, bx_start, by_start, bz_start, t_start, bx_stop, by_stop, bz_stop, t_stop, DX, DY, vel4vm, veldave, windowsize = 19

  dovx = vel4vm.u0
  dovy = vel4vm.v0
  dovz = vel4vm.w0

  dovx1 = dovx - (dovx * dobx + dovy * doby + dovz * dobz) * dobx / (dobx ^ 2 + doby ^ 2 + dobz ^ 2)
  dovy1 = dovy - (dovx * dobx + dovy * doby + dovz * dobz) * doby / (dobx ^ 2 + doby ^ 2 + dobz ^ 2)
  dovz1 = dovz - (dovx * dobx + dovy * doby + dovz * dobz) * dobz / (dobx ^ 2 + doby ^ 2 + dobz ^ 2)

  doAp = vectorpoten(dobz)

  vx = dovx1
  vy = dovy1
  vz = dovz1
  bx = dobx
  by = doby
  bz = dobz
  sz = size(vx)
  Ap_x = fltarr(sz[1], sz[2])
  Ap_y = fltarr(sz[1], sz[2])
  Ap_x[*, *] = doAp[0, *, *]
  Ap_y[*, *] = doAp[1, *, *]

  ; 计算磁螺度
  dHm_shear[i] = total((Ap_x[*, *] * vx[*, *] + Ap_y[*, *] * vy[*, *]) * bz[*, *], /nan)
  dHm_shear[i] = -2. * dHm_shear[i]
  dHm_emerge[i] = total((Ap_x[*, *] * bx[*, *] + Ap_y[*, *] * by[*, *]) * vz[*, *], /nan)
  dHm_emerge[i] = 2. * dHm_emerge[i]
  dHm[i] = dHm_shear[i] + dHm_emerge[i]

  dHm[i] = dHm[i] * factor2
  dHm_shear[i] = dHm_shear[i] * factor2
  dHm_emerge[i] = dHm_emerge[i] * factor2
  Hm[i] = total(dHm / 1e20, /nan)
  print, 'i=', i, '  dhdt=', dHm[i], dHm_shear[i], dHm_emerge[i]

  ; 计算能量
  dE[i] = (total((bx[*, *] * bx[*, *] + by[*, *] * by[*, *]) * vz[*, *], /nan) - total((bx[*, *] * vx[*, *] + by[*, *] * vy[*, *]) * bz[*, *], /nan)) / 4 / !pi * factor2
  E[i] = total(dE, /nan)

  ; 计算电流螺度
  jz = 1.0 / (4 * !pi * 1e-3) * ((by_stop - by_start) / (2. * factor2) - (bx_stop - bx_start) / (2. * factor2))
  jz1 = median(jz, 3)
  jz_total[i] = total(jz, /nan)
  Hc[i] = total(abs((4 * !pi * 0.001) * jz * bz), /nan)

  save, vx, vy, vz, Ap_x, Ap_y, filename = 'vm_' + HARP_num + int2str(i) + '.sav'
endfor

save, time, dHm_shear, dHm_emerge, dHm, Hm, dE, E, jz_total, Hc, filename = 'hf_' + HARP_num + '.sav'
end