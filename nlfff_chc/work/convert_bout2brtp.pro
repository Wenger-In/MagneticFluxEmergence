pro convert_bout2brtp,f,name
f=file_search(f,count=n)
for i=0,n-1 do begin 
restore,f[i]
br=twbox.bz
bth=-twbox.by
bph=twbox.bx
s=size(br)
nr=s[3]      
nlon=s[1]  
nlat=s[2]  

lon=twbox.xp ;  lontitude degree
lat=twbox.yp ;
theta=(90-lat)*!dtor
phi=(lon mod 360)*!dtor
rix=findgen(nr)*(max(lon)-min(lon))/nlon*2.*!pi/360.+1  ;  unit:rsun   e.g. [1,1.01,1.02.........]
name=strmid(f[i],0,strpos(f[i],'.sav'))+'.rtp.sav'
save,br,bth,bph,nr,nlat,nlon,rix,lat,lon,theta,phi,name=name
endfor
end
