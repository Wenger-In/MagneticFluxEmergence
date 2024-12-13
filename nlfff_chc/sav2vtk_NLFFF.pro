; pro sav2vtk_chc,vecx=vecx,vecy=vecy,vecz=vecz,cx=cx,cy=cy,cz=cz,savname=savname
compile_opt idl2

fn = findfile('/Volumes/Elements/Bout_.20240503_000000_.sav', count = nn)
restore, fn[0], /v
vecx = TWBOX.bx
vecy = TWBOX.by
vecz = TWBOX.bz
; vecx=bx & vecy=by & vecz=bz
savname = strmid(fn, 0, 20)
print, savname
; stop

savname = 'lfff_during2'

if not keyword_set(cx) then cx = 1l
if not keyword_set(cy) then cy = 1l
if not keyword_set(cz) then cz = 1l
if not keyword_set(savname) then savname = '2021.1028.083600'
sz = size(vecz)
nx = sz[1] / cx
ny = sz[2] / cy
nz = sz[3] / cz
jvec = fltarr(nx, ny, nz, 3)
jvec[*, *, *, 0] = float(congrid(vecx, nx, ny, nz))
jvec[*, *, *, 1] = float(congrid(vecy, nx, ny, nz))
jvec[*, *, *, 2] = float(congrid(vecz, nx, ny, nz))
; bvec[*,*,*,0]=FLOAT(b3dx[irange[0]:irange[1],irange[2]:irange[3],0:nz-1])
; bvec[*,*,*,1]=FLOAT(b3dy[irange[0]:irange[1],irange[2]:irange[3],0:nz-1])
; bvec[*,*,*,2]=FLOAT(b3dz[irange[0]:irange[1],irange[2]:irange[3],0:nz-1])

x = findgen(nx)
y = findgen(ny)
z = findgen(nz)
dx = 1.
dy = dx
dz = dx

obj = create_struct( $
  ['nx', 'ny', 'nz', 'dx', 'dy', 'dz', $
  'x', 'y', 'z', $
  'jvec'], $
  nx, ny, nz, dx, dy, dz, $
  x, y, z, $
  jvec)

; Determine the coordinates and the precision.
; if N_ELEMENTS(file) eq 0 then file = 'work.vtk'
file = savname + '.vtk'
data_type = strlowcase(size(x, /tname))

; Find the dimensions of the data.
nx = n_elements(x)
ny = n_elements(y)
nz = n_elements(z)
nxny = long(nx * ny)
ntot = long(nx * ny * nz)
dimensions = fix([nx, ny, nz])
origin = fix([x[0], y[0], z[0]])
ndim = n_elements(where(dimensions gt 1))

sdim = strtrim(dimensions)
snx = strmid(sdim[0], strpos(sdim[0], ' ', /reverse_search) + 1, 5)
sny = strmid(sdim[1], strpos(sdim[1], ' ', /reverse_search) + 1, 5)
snz = strmid(sdim[2], strpos(sdim[2], ' ', /reverse_search) + 1, 5)

sori = strtrim(origin)
sox = strmid(sori[0], strpos(sori[0], ' ', /reverse_search) + 1, 5)
soy = strmid(sori[1], strpos(sori[1], ' ', /reverse_search) + 1, 5)
soz = strmid(sori[2], strpos(sori[2], ' ', /reverse_search) + 1, 5)

; Open the VTK file for write.
print, ' f2vtk: Writing .vtk file constituents ... '
openw, lun, file, /get_lun, /swap_if_big_endian

; Write the header information.
printf, lun, '# vtk DataFile Version 2.0'
printf, lun, 'NLFF volume data'
printf, lun, 'ASCII'
printf, lun, 'DATASET RECTILINEAR_GRID'
printf, lun, 'DIMENSIONS ', snx, ' ', sny, ' ', snz

printf, lun, 'X_COORDINATES ', snx, ' ', data_type
printf, lun, x
; PRINTF, lun, ' '
printf, lun, 'Y_COORDINATES ', sny, ' ', data_type
printf, lun, y
; PRINTF, lun, ' '
printf, lun, 'Z_COORDINATES ', snz, ' ', data_type
printf, lun, z
; PRINTF, lun, ' '
; PRINTF, lun, 'ORIGIN ', sox, ' ', soy, ' ', soz
printf, lun, 'POINT_DATA ', ntot

; Write out each data field.
tags = tag_names(obj)

for i = 0, n_tags(obj) - 1 do begin
  ; Skip non-data tags.
  if tags[i] eq 'X' or tags[i] eq 'Y' or tags[i] eq 'Z' or $
    tags[i] eq 'NX' or tags[i] eq 'NY' or tags[i] eq 'NZ' or $
    tags[i] eq 'DX' or tags[i] eq 'DY' or tags[i] eq 'DZ' then $
    continue

  print, ' ' + tags[i], ' : SIZE(obj.(i)) = ', (size(obj.(i)))
  data_type = strlowcase(size(obj.(i), /tname))

  data = reform(obj.(i))

  if (size(obj.(i)))[0] le ndim then begin
    if (size(obj.(i)))[0] lt ndim then $
      ntot2 = nxny $
    else $
      ntot2 = ntot
    printf, lun, ' SCALARS ', strlowcase(tags[i]), ' ', data_type
    printf, lun, ' LOOKUP_TABLE default'
    for j = 0l, ntot2 - 1 do $
      printf, lun, data[j]
    print, ' Scalar field written.'
  endif else begin
    data_x = reform(data[*, *, *, 0])
    data_y = reform(data[*, *, *, 1])
    data_z = reform(data[*, *, *, 2])
    printf, lun, ' VECTORS ', strlowcase(tags[i]), ' ', data_type
    for j = 0l, ntot - 1 do $
      printf, lun, data_x[j], data_y[j], data_z[j]
    print, ' Vector field written.'
  endelse
endfor

close, lun
free_lun, lun

print, ' f2vtk: done.'
print, ' File saved as: ', file

end