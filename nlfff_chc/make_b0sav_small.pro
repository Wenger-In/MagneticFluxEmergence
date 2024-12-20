compile_opt idl2
cd, '/Users/zhengsun/Documents/code/Extrapolation/NLFFF_chc/nlfff_chc'

; .r '/Users/zhengsun/Documents/code/Extrapolation/NLFFF_chc/nlfff_chc/cea2map.pro'
ceaz = file_search('cea_raw/*.Br.fits')
ceax = file_search('cea_raw/*.Bp.fits')
ceay = file_search('cea_raw/*.Bt.fits', count = n)
ceaxerr = file_search('cea_raw/*.Bp_err.fits')
ceayerr = file_search('cea_raw/*.Bt_err.fits', count = n)
path = 'input/'

cea2map, ceaz[0], m
plot_image, m.data
cursor, x0, y0, /down, /data
plots, [x0, x0], [y0 - 2, y0 + 2], /data
plots, [x0 - 2, x0 + 2], [y0, y0], /data
cursor, x1, y1, /down, /data
plots, [x1, x1], [y1 - 2, y1 + 2], /data
plots, [x1 - 2, x1 + 2], [y1, y1], /data
plots, [x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], /data, col = 255

x0 = round(x0)
y0 = round(y0)
x1 = round(x1)
y1 = round(y1)
x00 = (x1 - x0) - [(x1 - x0) mod 8]
x1 = x00 + x0 - 1
y00 = (y1 - y0) - [(y1 - y0) mod 8]
y1 = y00 + y0 - 1
print, x0, y0, x1, y1
print, x1 - x0 + 1, y1 - y0 + 1

for i = 0, n - 1 do begin
  cea2map, ceaz[i], mapz
  cea2map, ceax[i], mapx
  cea2map, ceay[i], mapy
  cea2map, ceaxerr[i], mapxerr
  cea2map, ceayerr[i], mapyerr

  file1 = ceaz[i]
  if strpos(file1, '.20') lt 0 then message, 'note the filename of input and output' else $
    outname = strmid(file1, strpos(file1, '.20'), 17)

  xrange = [x0, x1]
  yrange = [y0, y1]
  ; 694        1213 ;chc2018-4-17
  ; 413          20

  ; 715        1234
  ; 5         420

  print, xrange, yrange
  ; gx=x1-x0+1
  ; gy=y1-y0+1

  gx = (xrange[1] - xrange[0] + 1) / 2
  gy = (yrange[1] - yrange[0] + 1) / 2
  mapz = rebin_map(get_sub_map(mapz, xrange = xrange, yrange = yrange, /pixel), gx, gy)
  mapx = rebin_map(get_sub_map(mapx, xrange = xrange, yrange = yrange, /pixel), gx, gy)
  mapy = rebin_map(get_sub_map(mapy, xrange = xrange, yrange = yrange, /pixel), gx, gy)
  xerr = rebin_map(get_sub_map(mapxerr, xrange = xrange, yrange = yrange, /pixel), gx, gy)
  yerr = rebin_map(get_sub_map(mapyerr, xrange = xrange, yrange = yrange, /pixel), gx, gy)

  xerr = xerr.data
  yerr = yerr.data
  err0 = (xerr ^ 2 + yerr ^ 2) ^ 0.5
  err = 1 / err0 / (max(1 / err0))

  bz = mapz.data
  bx = mapx.data
  by = -mapy.data
  xp = (get_map_xp(mapz))[*, 0]
  yp = rotate((get_map_yp(mapz))[0, *], 3)
  time = mapz.time
  save, bx, by, bz, err, xp, yp, time, filename = path + 'bfield' + outname + '.sav'
endfor
end