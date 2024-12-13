function trace_flines_xyz,bx_in,by_in,bz_in,sx,sy,sz,def_sign=def_sign
common shared, bx, by, bz, x, y, z, sign
bx=bx_in & by=by_in & bz=bz_in
maxiter=3e4
s=size(bx) & nx=s[1] & ny=s[2] & nz=s[3]
x=findgen(nx)
y=findgen(ny)
z=findgen(nz)
xmin = min(x, max = xmax)
ymin = min(y, max = ymax)
zmin = min(z, max = zmax)

np=n_elements(sx)
line={n:np,x:ptrarr(np),y:ptrarr(np),z:ptrarr(np)}

; determine whether field is up or down
; h is stepsize, by default 1/8 of minimum gridsize
for i=0L,np-1L do begin
x0=sx[i] & y0=sy[i] & z0=sz[i]

    If(x0 ge 0) And (y0 ge 0) And (z0 ge 0) $ 
        And (x0 Le nx-1) And (y0 Le ny-1) And (z0 Le nz -1) Then Begin
      h = 1/4.0 & s = 0.0 & coords = [x0, y0, z0] & xp=x0 & yp=y0 & zp=z0
      sign = float(0)
      bznearby = lin3interp(bz, x, y, z, x0, y0, z0)
                                ; always integrate into the box
         
      if (bznearby gt 0.0) then sign = 1.0 else sign = -1.0 
          if keyword_set(def_sign) then sign =def_sign
                                ; integrate until field line leaves box
      ;If(keyword_set(linecolor)) Then flcol = linecolor $
      ;Else flcol = bzb[ix, iy]
      ;print, bz[ix, iy, 0], flcol 
      count = 0l
      while (coords[0] ge xmin) and (coords[0] le xmax) and $
        (coords[1] ge ymin) and (coords[1] le ymax) and $
        (coords[2] ge zmin) and (coords[2] le zmax) and $
        count lt maxiter do begin
                                ; RHS of ODEs determined by linear
                                ; interpolation
        count = count+1l
        dydx = differential2(s, coords)
        result = rk4(coords, dydx, s, h, 'differential2')
                    
        coords = result
         xp=[xp,result[0]] & yp=[yp,result[1]] & zp=[zp,result[2]]
      endwhile
 
           if z0 gt 0 and (keyword_set(def_sign) eq 0)  then begin 
            ;  sign=-sign
               ;x0=xp[n_elements(xp)-2]
               ;y0=yp[n_elements(xp)-2]
               ;z0=zp[n_elements(xp)-2]
             h = 1/4.0 &  s = 0.0 & coords = [x0, y0, z0] 
         ;bznearby = lin3interp(bz, x, y, z, x0,y0,z0)
                                ; always integrate into the box
            ;  xp=x0 & yp=y0 & zp=z0
     if (bznearby gt 0.0) then sign = 1.0 else sign = -1.0 
               sign=-sign  
    
       count = 0l
      while (coords[0] ge xmin) and (coords[0] le xmax) and $
        (coords[1] ge ymin) and (coords[1] le ymax) and $
        (coords[2] ge zmin) and (coords[2] le zmax) and $
        count lt maxiter do begin
                        ; RHS of ODEs determined by linear
                                ; interpolation
        count = count+1l
        dydx = differential2(s, coords)
        result = rk4(coords, dydx, s, h, 'differential2')
                                ; plot increment in field line
        ;plots, [coords[0], result[0]], [coords[1], result[1]], $
         ; [coords[2], result[2]], /t3d, /data, thick = thick_x, color = flcol
        coords = result
         xp=[result[0],xp] & yp=[result[1],yp] & zp=[result[2],zp]
      endwhile

             endif
    Endif
line.x[i]=ptr_new(xp)
line.y[i]=ptr_new(yp)
line.z[i]=ptr_new(zp)
  Endfor
return,line
end
