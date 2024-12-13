ceaz=file_search('./*.Br.fits')
ceax=file_search('./*.Bp.fits')
ceay=file_search('./*.Bt.fits',count=n)
ceaxerr=file_search('./*.Bp_err.fits')
ceayerr=file_search('./*.Bt_err.fits',count=n)
path='/data1/hd/hd1/20150203/3d_line/input'

for i=0,n-1 do begin
;fits2map,ceaz[i],mapz
;fits2map,ceax[i],mapx
;fits2map,ceay[i],mapy
;fits2map,ceaxerr[i],mapxerr
;fits2map,ceayerr[i],mapyerr

read_sdo,ceaz[i],in,data
wcs=fitshead2wcs(struct2fitshead(in)) & wcs.simple=1
wcs2map,data,wcs,mapz

read_sdo,ceax[i],in,data
wcs=fitshead2wcs(struct2fitshead(in)) & wcs.simple=1
wcs2map,data,wcs,mapx

read_sdo,ceay[i],in,data
wcs=fitshead2wcs(struct2fitshead(in)) & wcs.simple=1
wcs2map,data,wcs,mapy

read_sdo,ceaxerr[i],in,data
wcs=fitshead2wcs(struct2fitshead(in)) & wcs.simple=1
wcs2map,data,wcs,mapxerr

read_sdo,ceayerr[i],in,data
wcs=fitshead2wcs(struct2fitshead(in)) & wcs.simple=1
wcs2map,data,wcs,mapyerr

mapz=get_sub_map(mapz,xrange=[350,989],yrange=[100,499],/pixel)
mapx=get_sub_map(mapx,xrange=[350,989],yrange=[100,499],/pixel)
mapy=get_sub_map(mapy,xrange=[350,989],yrange=[100,499],/pixel)
xerr=get_sub_map(mapxerr,xrange=[350,989],yrange=[100,499],/pixel)
yerr=get_sub_map(mapyerr,xrange=[350,989],yrange=[100,499],/pixel)
xerr=xerr.data
yerr=yerr.data
err0=(xerr^2+yerr^2)^0.5
err=1/err0/(max(1/err0))


bz=mapz.data
bx=mapx.data
by=-mapy.data
xp=(get_map_xp(mapz))[*,0]
yp=rotate((get_map_yp(mapz))[0,*],3)
time=mapz.time

time=sxpar(headfits(ceaz[i],/ext),'T_rec')
if strpos(time,'_TAI') gt 0 then time=strmid(time,0,strpos(time,'_TAI') )
save,bx,by,bz,err,xp,yp,time,filename=path+'bfield_'+anytime2name(time)+'.sav'
endfor
end
