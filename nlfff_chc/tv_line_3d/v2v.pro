function vv
end
pro rotz,x,y,z,ang,xr,yr,zr
ang=ang*!dtor
xr=x*cos(ang)-y*sin(ang)
yr=x*sin(ang)+y*cos(ang)
zr=z
end

pro rotx,x,y,z,ang,xr,yr,zr
ang=ang*!dtor
yr=y*cos(ang)-z*sin(ang)
zr=y*sin(ang)+z*cos(ang)
xr=x
end

pro roty,x,y,z,ang,xr,yr,zr
ang=ang*!dtor
zr=z*cos(ang)-x*sin(ang)
xr=z*sin(ang)+x*cos(ang)
yr=y
end

pro make_cir,qqrin,cirx,ciry
cirx=fltarr(360) & ciry=cirx
for ic=0,359 do begin
cirx[ic]=qqrin[0]+qqrin[2]*cos(ic*!dtor)
ciry[ic]=qqrin[1]+qqrin[2]*sin(ic*!dtor)
endfor
end

pro plot_cline,x,y,col,nop=nop,_extra=extra
tvlct,r0,g0,b0,/get
r=r0 & b=b0 & g=g0
r[253:255]=[0  ,255,0]
g[253:255]=[0  ,255,255]
b[253:255]=[255,0,0]
tvlct,r,g,b
if col eq 'y' then col0=254
if col eq 'g' then col0=255
if col eq 'b' then col0=253
if not exist(col0) then col0=0
 plots,x,y,col=col0,/data,_extra=extra
tvlct,r0,g0,b0
end

pro l_v2v,x,y,map,v1,refa,v2,xr,yr,zr
Carr =GET_STEREO_CARR_ROT(refa.time,v2)-GET_STEREO_CARR_ROT(map.time,v1)
carr=carr*360.
b0=map.b0
z=(map.rsun^2-x^2-y^2)^0.5
rotx,x,y,z,-b0,x1,y1,z1
roty,x1,y1,z1,carr,xt,yt,zt

rotx,xt,yt,zt,refa.b0,xr,yr,zr
xr=xr/map.rsun*refa.rsun
yr=yr/map.rsun*refa.rsun
end


pro v2v,map,mapv,goal,goalv,pbr=pbr
if keyword_set(pbr) then begin
stereo=mapv 
if not (stereo eq 'A' or stereo eq 'B') then stereo=''
pb0r=pb0r(map.time,/arcsec,stere=stereo)
map.b0=pb0r[1]
map.rsun=pb0r[2]
stereo=goalv 
if not (stereo eq 'A' or stereo eq 'B') then stereo=''
pb0r=pb0r(goal.time,/arcsec,stereo=stereo)
goal.b0=pb0r[1]
goal.rsun=pb0r[2]
endif

win9,0
plot_map,map,/log,/limb
cursor,x,y,/data
make_cir,[x,y,50.],cirx,ciry
                ;map.rsun=(pb0r(map.time,/arcsec))[2]
;gg=(findgen(101)/50.-1)*!pi/2.
;cirx=(map.rsun-1)*cos(gg)
;ciry=(map.rsun-1)*sin(gg)
;x=cirx[50]  &  y=ciry[50]

plot_cline,cirx,ciry,'g' ; plot_cline,cxa,cya,cola,noclip=0
helsdo=arcmin2hel(x/60.,y/60.,b0=map.b0,rsun=map.rsun,p=map.roll_angle) 

 l_v2v,x,y,map,mapv,goal,goalv,xa,ya,za
l_v2v,cirx,ciry,map,mapv,goal,goalv,cxa,cya,cza
 if za lt 0 then begin
  hela=arcmin2hel(xa/60.,ya/60.,b0=goal.b0,rsun=goal.rsun,p=goal.roll_angle,/backside) 
  if hela[1] lt 0 then hela[1]=-180-hela[1]
  if hela[1] ge 0 then hela[1]=180-hela[1]
  hela=[hela[0],hela[1]]
  linesta=1 & cola='y'
 endif else begin
  hela=arcmin2hel(xa/60.,ya/60.,b0=goal.b0,rsun=goal.rsun,p=goal.roll_angle)
  linesta=0 & cola='g'
 endelse
 win9,1
 plot_map,goal,/log,/limb
 plot_cline,cxa,cya,cola,noclip=0



print,'the location in map: '
print,' arcsec: ',strtrim(x,2)+','+strtrim(y,2)
print,' lat-lon: ',helsdo
print,'the location in goal: '
print,' arcsec: ',strtrim(xa,2)+','+strtrim(ya,2)
print,' lat-lon: ',hela
end
