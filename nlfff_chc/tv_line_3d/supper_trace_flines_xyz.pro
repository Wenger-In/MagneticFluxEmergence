function supper_trace_flines_xyz,bx,by,bz,sx0,sy0,sz0,def_sign=def_sign,endpoint=endpoint,noreverse=noreverse,keepvalid=keepvalid
st=systime(1)
maxiter=3e4
s=size(bx) & nx=s[1] & ny=s[2] & nz=s[3]
x=findgen(nx)
y=findgen(ny)
z=findgen(nz)
xmin = min(x, max = xmax)
ymin = min(y, max = ymax)
zmin = min(z, max = zmax)

np=n_elements(sx0)
line={n:np,x:ptrarr(np),y:ptrarr(np),z:ptrarr(np)}
if not keyword_set(def_sign) then sign=1 else sign=def_sign

sx=[[findgen(np)],[reform(sx0)]]
sy=[[findgen(np)],[reform(sy0)]]
sz=[[findgen(np)],[reform(sz0)]]
endp=1;(size(sx))[2]-
ok=where((sx[*,endp] ge 0) And (sy[*,endp] ge 0) And (sz[*,endp] ge 0) $ 
        And (sx[*,endp] Le nx-1) And (sy[*,endp] Le ny-1) And (sz[*,endp] Le nz -1),nok)
 count = 0l    & cok=-1
      while  nok gt 0 $
      and  count lt maxiter do begin
                                ; RHS of ODEs determined by linear
                                ; interpolation
if cok[0] ne -1 then begin
sx=sx[ok,*]
sy=sy[ok,*]
sz=sz[ok,*]
endif 
       count = count+1l
        bxinterp=interpolate(bx,sx[*,endp],sy[*,endp],sz[*,endp])
byinterp=interpolate(by,sx[*,endp],sy[*,endp],sz[*,endp])
bzinterp=interpolate(bz,sx[*,endp],sy[*,endp],sz[*,endp])
binterp=sqrt(bxinterp^2+byinterp^2+bzinterp^2);>1.e-9
    dydx=sign*[[bxinterp/binterp],[byinterp/binterp],bzinterp/binterp]
sx=[[sx],[sx[*,endp]+sign*bxinterp/binterp]]
sy=[[sy],[sy[*,endp]+sign*byinterp/binterp]]
sz=[[sz],[sz[*,endp]+sign*bzinterp/binterp]]
endp=endp+1l

ok=where((sx[*,endp] ge 0) And (sy[*,endp] ge 0) And (sz[*,endp] ge 0) $ 
        And (sx[*,endp] Le nx-1) And (sy[*,endp] Le ny-1) And (sz[*,endp] Le nz -1),nok,complement=cok)
if cok[0] ne -1 then begin
for j=0L,n_elements(cok)-1 do begin

if endp gt 2 then begin
if keyword_set(endpoint) then begin
line.x[sx[cok[j],0] ]=ptr_new(reform(sx[cok[j],[1,endp-1]]))
line.y[sx[cok[j],0] ]=ptr_new(reform(sy[cok[j],[1,endp-1]]))
line.z[sx[cok[j],0] ]=ptr_new([reform(sz[cok[j],[1,endp-1]]),max(sz[cok[j],1:endp-1])])
endif else begin
line.x[sx[cok[j],0] ]=ptr_new(reform(sx[cok[j],1:endp-1]))
line.y[sx[cok[j],0] ]=ptr_new(reform(sy[cok[j],1:endp-1]))
line.z[sx[cok[j],0] ]=ptr_new(reform(sz[cok[j],1:endp-1]))
endelse
endif
endfor

endif

        ;result = rk4(coords, dydx, s, h, 'differential2')
         ;   dydx2= differential2(s, coords+dydx*h*0.5)
          ; dydx3= differential2(s, coords+dydx2*h*0.5)
           ; dydx4= differential2(s, coords+dydx3*h)
        ;coords =coords+h*( 2*dydx+dydx2+dydx3+dydx4)/5.
         ;xp=[xp,coords[0]] & yp=[yp,coords[1]] & zp=[zp,coords[2]]
      endwhile
;=======================sign =1===========================
if (not keyword_set(def_sign)) and (not  keyword_set(endpoint)) then begin
pv=ptr_valid(line.x)
if keyword_set(noreverse) then rnum=0 else rnum=2
sign=-1
sx=[[findgen(np)],[reform(sx0)]]
sy=[[findgen(np)],[reform(sy0)]]
sz=[[findgen(np)],[reform(sz0)]]
endp=1;(size(sx))[2]-
ok=where((sx[*,endp] ge 0) And (sy[*,endp] ge 0) And (sz[*,endp] ge 0) $ 
        And (sx[*,endp] Le nx-1) And (sy[*,endp] Le ny-1) And (sz[*,endp] Le nz -1),nok)
 count = 0l    & cok=-1
      while  nok gt 0 $
      and  count lt maxiter do begin
                                ; RHS of ODEs determined by linear
                                ; interpolation
if cok[0] ne -1 then begin
sx=sx[ok,*]
sy=sy[ok,*]
sz=sz[ok,*]
endif 
       count = count+1l
        bxinterp=interpolate(bx,sx[*,endp],sy[*,endp],sz[*,endp])
byinterp=interpolate(by,sx[*,endp],sy[*,endp],sz[*,endp])
bzinterp=interpolate(bz,sx[*,endp],sy[*,endp],sz[*,endp])
binterp=sqrt(bxinterp^2+byinterp^2+bzinterp^2);>1.e-9
    dydx=sign*[[bxinterp/binterp],[byinterp/binterp],bzinterp/binterp]
sx=[[sx],[sx[*,endp]+sign*bxinterp/binterp]]
sy=[[sy],[sy[*,endp]+sign*byinterp/binterp]]
sz=[[sz],[sz[*,endp]+sign*bzinterp/binterp]]
endp=endp+1l

ok=where((sx[*,endp] ge 0) And (sy[*,endp] ge 0) And (sz[*,endp] ge 0) $ 
        And (sx[*,endp] Le nx-1) And (sy[*,endp] Le ny-1) And (sz[*,endp] Le nz -1),nok,complement=cok)
if cok[0] ne -1 then begin
for j=0L,n_elements(cok)-1 do begin

if endp gt 2 then begin
if pv[sx[cok[j],0]] eq 0 then begin
line.x[sx[cok[j],0] ]=ptr_new(rotate(reform(sx[cok[j],1:endp-1]),rnum))
line.y[sx[cok[j],0] ]=ptr_new(rotate(reform(sy[cok[j],1:endp-1]),rnum))
line.z[sx[cok[j],0] ]=ptr_new(rotate(reform(sz[cok[j],1:endp-1]),rnum))
endif 
if pv[sx[cok[j],0]] eq 1 and rnum eq 2 then begin
line.x[sx[cok[j],0] ]=ptr_new([rotate(reform(sx[cok[j],1:endp-1]),rnum),*line.x[sx[cok[j],0] ]])
line.y[sx[cok[j],0] ]=ptr_new([rotate(reform(sy[cok[j],1:endp-1]),rnum),*line.y[sx[cok[j],0] ]])
line.z[sx[cok[j],0] ]=ptr_new([rotate(reform(sz[cok[j],1:endp-1]),rnum),*line.z[sx[cok[j],0] ]])
endif
endif
endfor

endif

        ;result = rk4(coords, dydx, s, h, 'differential2')
         ;   dydx2= differential2(s, coords+dydx*h*0.5)
          ; dydx3= differential2(s, coords+dydx2*h*0.5)
           ; dydx4= differential2(s, coords+dydx3*h)
        ;coords =coords+h*( 2*dydx+dydx2+dydx3+dydx4)/5.
         ;xp=[xp,coords[0]] & yp=[yp,coords[1]] & zp=[zp,coords[2]]
      endwhile

print,count
endif
if not keyword_set(keepvalid) then begin
valid=where(ptr_valid(line.x) eq 1,nvalid)
print,nvalid,'/',np
if nvalid lt np then line={n:nvalid,x:line.x[valid],y:line.y[valid],z:line.z[valid]}
endif
print,'trace lines: time used',systime(1)-st
return,line
end
