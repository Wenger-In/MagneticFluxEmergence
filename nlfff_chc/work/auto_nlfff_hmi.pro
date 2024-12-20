compile_opt idl2

dir_work = './' ; directory of nlfff C-programs
dir_data = '../input/' ; directory of boundary magnetograms
dir_sav = '../input/' ; directory of output 3d magnegograms

fn = findfile(dir_data + 'bfield.2011*_.sav') ; input boundary magnetograms
num = n_elements(fn)

; for i=160, 219 do begin
for i = 0, num - 1 do begin
  print, 'i=', i

  spawn, 'rm -f grid*.ini'
  spawn, 'rm -f allboundar*.dat'

  string_number = strmid(fn[i], strpos(fn[i], '_2'), 17)

  file1 = fn[i]
  if strpos(file1, '.20') lt 0 then message, 'note the filename of input and output' else $
    outname = strmid(file1, strpos(file1, '.20'), 17)
  def_nz = 200
  def_nd = 32
  def_level = 3
  folder = dir_work
  prep = 0
  ; ====== restore bx,by & bz ===============
  ; auto_load_mag,file1
  restore, file1

  ; bx=congrid(bx,256,256)
  ; by=congrid(by,256,256)
  ; bz=congrid(bz,256,256)
  s = size(bz)
  nx = s[1]
  ny = s[2]
  nz = min([nx, ny])
  nd = def_nd
  ; =========== multigrid ===================
  ; multigrid, def_level, 0
  level = def_level ; 3
  if n_elements(nz) eq 0 then nz = min([nx, ny]) ;
  if n_elements(nd) eq 0 then nd = min([nx, ny]) / 8 ; 32
  ; help, folder,nx,ny,nz,nd
  reducer = 2 ^ (level - 1) ; 4
  ;
  help, nx, ny, nz, nd
  ; test if nx,ny,nz,nd are multiples of reducer
  dummy = (nx mod reducer) + (ny mod reducer) + $
    (nz mod reducer) + (nd mod reducer) ; 4*
  if (dummy ne 0) then print, 'nx,ny,nz,nd must me multiples of', reducer

  if (dummy eq 0) then begin
    ; print,'Multigrid on', level, ' grids.'
    get_lun, u
    openw, u, folder + 'nd.ini'
    printf, u, 'nd'
    printf, u, nd
    printf, u, 'nxmax'
    printf, u, nx
    printf, u, 'nymax'
    printf, u, ny
    printf, u, 'nzmax'
    printf, u, nz
    close, u
    free_lun, u
    for ii = level, 1, -1 do begin
      print, 'Processing grid', nx, ny, nz, nd
      ; Preprocessing yes or no? (prep=1 means yes)
      if (prep eq 1) then begin
        bxorig = bx
        byorig = by
        bzorig = bz
        correct_magnetogram
      endif
      ; fullmag
      get_lun, u
      gridstring = 'grid' + strcompress(string(ii), /remove_all) + '.ini'
      openw, u, folder + gridstring
      printf, u, 'nx'
      printf, u, nx
      printf, u, 'ny'
      printf, u, ny
      printf, u, 'nz'
      printf, u, nz
      printf, u, 'mu'
      printf, u, 0.001
      printf, u, 'nd'
      printf, u, nd
      close, u
      allstring = 'allboundaries' + strcompress(string(ii), /remove_all) + '.dat'
      print, 'Saving: ', gridstring
      print, 'Saving: ', allstring
      openw, u, folder + allstring
      for iy = 0, ny - 1 do begin
        for ix = 0, nx - 1 do begin
          printf, u, bx[ix, iy]
          printf, u, by[ix, iy]
          printf, u, bz[ix, iy]
        endfor
      endfor
      close, u
      maskhmi, bz, u, 'mask' + strcompress(string(ii), /remove_all) + '.dat'
      close, u
      free_lun, u
      savstring = 'mag' + strcompress(string(ii), /remove_all) + '.sav'
      ; save,filename=folder+savstring,bx,by,bz,nd

      if (prep eq 1) then begin
        bx = bxorig
        by = byorig
        bz = bzorig
      endif
      ; Rebin grid
      if (ii gt 1) then begin
        nx = nx / 2
        ny = ny / 2
        nz = nz / 2
        nd = nd / 2
        bx = rebin(bx, nx, ny)
        by = rebin(by, nx, ny)
        bz = rebin(bz, nx, ny)
        Lx = 1.0 * (nx - 1)
        Ly = 1.0 * (ny - 1)
        Lz = 1.0 * (nz - 1)
      endif
    endfor
  endif

  ; ============preprocessing==========================
  ; file2='./multiprepro '+strcompress(string(def_level),/remove_all)
  ; spawn,file2

  ; =========nlfff & potential extrapolation===========

  file3 = './multigrid ' + strcompress(string(def_level), /remove_all)
  print, file3
  print, 'NLFFF extrapolating, i=', i
  spawn, file3
  print, 'Bout.bin -> Bout.sav, i=', i
  auto_convert, 1, file1, './', outname, filename = filename ; Bout.bin -> Bout.sav
  ; convert_bout2brtp, filename
  ; spawn,['mv',name,dir_sav],/noshell
  print, 'potential extrapolating, i=', i
  spawn, './relax2 23 0'
  spawn, './fftpot4'
  print, 'B0.bin -> B0.sav, i=', i
  auto_convert, 0, file1, './', outname, filename = filename ; B0.bin -> B0.sav
  ; convert_bout2brtp,filename
  ; spawn,['mv',name,dir_sav],/noshell
endfor

end