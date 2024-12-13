function c_twists,x0,y0,z0,tm=tm

 ;y=x0
 ;x=y0
 ;z=z0

;if n_elements(x) gt 10 then begin
x=smooth(x0,5)
y=smooth(y0,5)
z=smooth(z0,5)
;x=smooth(x,100)
;y=smooth(y,100,/e)
;z=smooth(z,100,/e)
;endif

d1=((x-shift(x,1))^2+(y-shift(y,1))^2+(z-shift(z,1))^2)^0.5
d1[0]=0 & d=d1*0
for i=0,n_elements(d)-1 do d[i]=total(d1[0:i])
tx=deriv(d, x) & ty=deriv(d, y) & tz=deriv(d, z)
vx=deriv(d,tx) & vy=deriv(d,ty) & vz=deriv(d,tz)
tl=(tx^2+ty^2+tz^2)^0.5
vl=(vx^2+vy^2+vz^2)^0.5
tx=tx/tl & ty=ty/tl & tz=tz/tl
vx=vx/vl & vy=vy/vl & vz=vz/vl
dx=deriv(d,vx) & dy=deriv(d,vy) & dz=deriv(d,vz)

   cx = vy*dz - vz*dy
   cy = vz*dx - vx*dz
   cz = vx*dy - vy*dx
tt1=(d1*((tx*cx+ty*cy+tz*cz)/2./!pi))>(-0.02)<0.02
n=n_elements(tt1)

;tt=total(tt1[n/8.:7*n/8.])*8/6
tt=total(tt1[2:n-3])
return,tt

end

