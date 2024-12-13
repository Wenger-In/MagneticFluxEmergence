;================================================================================
;auto_convert,1,file1,dir_sav,string_number  ; Bout.bin -> Bout.sav
PRO auto_convert, dummy, file1,dir_sav,string_number ,filename=filename


if n_elements(dummy) eq 0 then dummy=1
if n_elements(file1) eq 0 then file1=''
if n_elements(dir_sav) eq 0 then dir_sav=''
if n_elements(string_number) eq 0 then string_number=''
if dummy eq 0 then file4='B0' else file4='Bout'
Openr, 1,dir_sav+file4+'.bin'
fileInfo=FStat(1)
hilfe=fileinfo.size
dummyfeld=DblArr(hilfe/8)
readu,1, dummyfeld
close,1
param=strarr(10)
openr, 1,dir_sav+'grid.ini'
readf,1,param
close, 1
nx=long(param[1])
ny=long(param[3])
nz=long(param[5])
notneeded=float(param[7])
nd=long(param[9])
nynz=ny*nz
nxnynz=nx*ny*nz
;
b3dx=fltarr(nx,ny,nz)
b3dy=fltarr(nx,ny,nz)
b3dz=fltarr(nx,ny,nz)
;
for ix=0, nx-1 do begin
for iy=0, ny-1 do begin
for iz=0, nz-1 do begin
i=ix*nynz+iy*nz+iz
i2=i+nxnynz
i3=i2+nxnynz
b3dx[ix,iy,iz]=dummyfeld[i]
b3dy[ix,iy,iz]=dummyfeld[i2]
b3dz[ix,iy,iz]=dummyfeld[i3]
endfor
endfor
endfor
;;
b3dabs=sqrt(b3dx^2+b3dy^2+b3dz^2)
bx=reform(b3dx[*,*,0])
by=reform(b3dy[*,*,0])
bz=reform(b3dz[*,*,0])
text1a=strarr(10)
text1a[0]='Standard output'
text1a[1]='created by'
text1a[2]='auto_convert'
text1a[3]='3D magnetic fields'
text1a[4]=file1

restore,file1
twbox={bx:b3dx,by:b3dy,bz:b3dz,info:text1a,xp:xp,yp:yp,time:time}
filename=dir_sav+file4+'_'+string_number+'.sav'
save,filename=filename,twbox
print,'B-field in box saved'
print,filename
END

;================================================================================
;
