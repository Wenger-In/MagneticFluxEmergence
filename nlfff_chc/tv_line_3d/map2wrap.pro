function l2r,m
;jsoc2map,'../fd/m720/hmi_m_2011.0606.031025_720s.fits',m,/fast
;m=rot_map(m,-m.roll_angle)
s=size(m.data)
xc=0.5*s[1]-m.xc/m.dx
yc=0.5*s[2]-m.yc/m.dy
rs=m.rsun/m.dx
y=long(dindgen(s[1],s[2])/s[1])*1.

x=dindgen(s[1],s[2]) mod s[1]

co=((x-xc)^2.+(y-yc)^2.)/rs^2
out=where(co gt sin(80*!dtor))
if n_elements(out) gt 1 then co[out]=0
co=1/(1-co)^0.5

;m.data=co
;m.rsun=0.99*m.rsun
;plot_map,m,grid=10
return,co
end


pro map2wrap,fn,savepath=savepath,skip=skip,log=log,verbose=verbose,prep=prep,lon=lon,lat=lat,los2r=los2r,step=step,threshold=threshold,smooth=smooth,jsoc=jsoc


fn=findfile(fn)
n=n_elements(fn)
;fits2map,fn[0],map
;xr=get_map_xrange(map)
;yr=get_map_yrange(map)
;lon1=((arcmin2hel(xr[0]/60.,mean(yr)/60.,b0=map.b0,rsun=map.rsun,p=map.roll_angle))[0]) 
;lon2=((arcmin2hel(xr[1]/60.,mean(yr)/60.,b0=map.b0,rsun=map.rsun,p=map.roll_angle))[0]) 
;lat1=((arcmin2hel(mean(xr)/60.,yr[0]/60.,b0=map.b0,rsun=map.rsun,p=map.roll_angle))[1]) 
;lat2=((arcmin2hel(mean(xr)/60.,yr[1]/60.,b0=map.b0,rsun=map.rsun,p=map.roll_angle))[1]) 
if not keyword_set(lon) then begin
lon1=-90. & lon2=90.
endif else begin
lon1=lon[0] & lon2=lon[1]
endelse

if not keyword_set(lat) then begin
lat1=-90. & lat2=90.
endif else begin
lat1=lat[0] & lat2=lat[1]
endelse

;lat1=-90. & lat2=90.

;lon1=-44  &  lon2=-26
;lat1=9   &  lat2=22

lat0=lat1 & lon0=lon1 & dlat=lat2-lat1 & dlon=lon2-lon1
xl=long(dlon/0.036)*1.
yl=long(dlat/0.036)*1.
lon=(dindgen(xl*yl) mod xl)/(xl-1)

lon=dlon*lon+lon0

;lat=asin(2.*long(dindgen(xl*yl)/xl)/(yl-1)-1.)/!dtor
lat=long(dindgen(xl*yl)/xl)/(yl-1.)

lat=dlat*lat+lat0
print,'Longitude: ',[min(lon),max(lon)]
print,'Latitude:  ',[min(lat),max(lat)]



ll=fltarr(2,xl*yl)


rtime=anytim2tai('2011-07-09 00:00:00')
if not keyword_set(step) then l=1 else l=step
forbegin=0
forend=n-1
for k=forbegin,forend,l do begin
nameb=rstrpos(fn[k],'/')
namel=strlen(fn[k])
if not keyword_set(savepath) then savepath='./'  else $
if not file_exist(savepath,/dir) then file_mkdir,savepath

wn=savepath+'wrap_'+strmid(fn[k],nameb+1,namel-nameb)

if findfile(wn) eq '' or not keyword_set(skip) then begin

if keyword_set(jsoc) then $
jsoc2map,fn[k],m   ,prep=prep  $ ; <============== raw fits
else fits2map,fn[k],m  
m=rot_map(m,-m.roll_angle)
;m.rsun=(pb0r(m.time,/arcsec))[2]
;help,m,/st
if keyword_set(los2r) then m.data=m.data*l2r(m)

dt=0;get_map_time(m,/tai)-rtime
;print,'dt=',dt/3600.

ll[0,*]=rotate(lon+diff_rot(dt/86400., lat, /synodic,_extra=extra),1)
ll[1,*]=rotate(lat,1)
xy=lonlat2xy(ll,radius=m.rsun,b0=m.b0*!dtor,behind=behind)
s=size(m.data)
;xy=0.5*s[1]-m.xc/m.dx+xy/m.dx
xy[0,*]=xy[0,*]/m.dx+0.5*s[1]-m.xc/m.dx
xy[1,*]=xy[1,*]/m.dy+0.5*s[2]-m.yc/m.dy

bxy=INTERPOLATE(m.data,xy[0,*],xy[1,*], MISSING = -1) 
bxy=bxy-(bxy+3e+4)*behind
w=fltarr(xl,yl)
for i=0,(yl-1) do begin
w[*,i]=bxy[i*xl:i*xl+(xl-1)]
endfor




xc=(lon2+lon1)/2.  & yc=(lat2+lat1)/2.
dx=(lon2-lon1)/xl
dy=(lat2-lat1)/yl
if keyword_Set(smooth) then w=smooth(w,smooth)
if keyword_set(threshold) then begin
wlt=where(abs(w) lt threshold,nwlt)
if nwlt gt 0 then w[wlt]=0
endif


map=make_map(w,xc=xc,yc=yc,dx=dx,dy=dy,time=m.time,l0=0,b0=0,rsun=960)
;help,xc,yc,dx,dy
stop
map2fits,map,wn
if keyword_set(verbose) or forbegin eq forend then begin
window,1,xs=800,ys=800
if keyword_set(log) then begin
plot_map,map,/log
endif else begin
plot_map,map,drange=[-300,300]
endelse

window,0,xs=800,ys=800
plot_map,m,grid=10,drange=[-300,300]



endif

endif
print,'writefits  ',wn,' ',strtrim(long(k),2),'/',strtrim(long(n-1),2)

endfor
end
