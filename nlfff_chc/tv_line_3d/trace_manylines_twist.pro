function trace_manylines_twist,bx,by,bz,sx0,sy0,sz0,nd=nd
s=size(bz)
curl_xyz, bx, by, bz, findgen(s[1]), findgen(s[2]), findgen(s[3]), jx, jy, jz
bb=(bx^2+by^2+bz^2)
jj=(jx^2+jy^2+jz^2)
jb=(jx*bx+jy*by+jz*bz)/bb


time=systime(1)
np=double(n_elements(sx0))
if not keyword_set(nd) then nd=double(1e6)

line={n:np,x:fltarr(np,2),y:fltarr(np,2),z:fltarr(np,3),t:fltarr(np),tw:fltarr(np),tw1:fltarr(np),tw2:fltarr(np),fl:fltarr(np)}

for1=findgen(long(np/nd)+1)*nd
for2=(findgen(long(np/nd)+1)*nd+nd-1)<(np-1) 

for i=0,n_elements(for1)-1 do begin
sx1=sx0[for1[i]:for2[i]]
sy1=sy0[for1[i]:for2[i]]
sz1=sz0[for1[i]:for2[i]]
l1=supper_trace_flines_xyz_bot_endpoint(bx,by,bz,sx1,sy1,sz1)

vw=where(ptr_valid(l1.x) eq 1,vn)
for j=0L,vn-1L do begin
k=vw[j]

;line.x[for1[i]:for2[i],*]=l1.x
;line.y[for1[i]:for2[i],*]=l1.y
;line.z[for1[i]:for2[i],*]=l1.z
line.x[k+for1[i],*]=[  (*(l1.x[k]))[0],(*(l1.x[k]))[l1.l[k]-1] ]
line.y[k+for1[i],*]=[  (*(l1.y[k]))[0],(*(l1.y[k]))[l1.l[k]-1] ]
line.z[k+for1[i],*]=[  (*(l1.z[k]))[0],(*(l1.z[k]))[l1.l[k]-1] ,max(*(l1.z[k]))]
 
xl=*l1.x[k] & yl=*l1.y[k] & zl=*l1.z[k]  &nl=l1.l[k]
if nl gt 20 then line.t[k+for1[i]]=c_twists(xl,yl,zl)
if zl[nl-1] lt 1 then begin
dl=((xl-shift(xl,1))^2+(yl-shift(yl,1))^2+(zl-shift(zl,1))^2)^0.5
dl[0]=0 
      ljb=interpolate(jb,xl,yl,zl)
      lb=interpolate(bb,xl,yl,zl)
line.tw[k+for1[i]]=total(ljb/(1.*n_elements(lb)));/(4*!pi)
line.tw1[k+for1[i]]=total(ljb*lb/total(lb));*total(dl)/(4*!pi)
line.tw2[k+for1[i]]=total(ljb*lb^2/total(lb^2));*total(dl)/(4*!pi)
line.fl[k+for1[i]]=total(dl)
endif
endfor
endfor

print,'trace_manylines_endpoint.pro used time ',systime(1)-time
return, line
end
