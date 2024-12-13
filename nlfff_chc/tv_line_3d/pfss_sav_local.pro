pro pfss_sav_local,hmir,rebin=rebin
;pfss_sav_local,'../vm/cea/hmi.sharp_cea_720s.3321.20131104_*4800_TAI.Br.fits'
@pfss_data_block
;hmir='/data1/biyi/vm/11429/rtp/hmi.B_720s_e15w1332_cutout.11429.20120307_000000_TAI.br.fits'
;hmir='../vm/cea/hmi.sharp_cea_720s.3321.20131104_*4800_TAI.Br.fits'
hmir=file_search(hmir,count=n)

if not keyword_set(rebin) then rebin=4
smooth=0

for j=0,n-1 do begin
;addrtp2synotic.pro <=======================================
;synoticn='down_pfss_sav/Bfield_20120306_000400.sav' & restore,synoticn


fits2map,hmir[j],mapr
;mapr=get_sub_map(mapr,xrange=[18,39],yrange=[-28,-7])

;synoticn='down_pfss_sav/'+path2name(pfss_time2file(mapr.time,/ssw_cat,/url))
synoticn='Bfield_20120909_120400.h5'
print,mapr.time,synoticn
 pfss_restore,synoticn,/refresh
slat=get_map_yrange(mapr)
slon=get_map_xrange(mapr)
imr=mapr.data

synotic=br[*,*,0]
s=size(synotic)
;rebin=3
nlat2=s[2]*rebin
nlon2=nlat2*2
 gaussquad_legendre,nlat2,cth,weights
    theta2=acos(cth)  ;  radians
    lat2=90-theta2*180/!dpi
lon2=(linrange(nlon2+1,0,360))(0:nlon2-1)
phi2=lon2*!dtor

nlat1=s[2]
nlon1=nlat1*2
 gaussquad_legendre,nlat1,cth,weights
    theta1=acos(cth)  ;  radians
    lat1=90-theta1*180/!dpi
lon1=(linrange(nlon1+1,0,360))(0:nlon1-1)
phi1=lon1*!dtor

;dlatix=asin(linrange(180,89.5,-89.5)/90.0)*180/!dpi  ;  180 slat bins
;    dlonix=linrange(360,1,360)  ;  360 longitude bins

    dlatinterp=get_interpolation_index(lat1,lat2)
    dloninterp=get_interpolation_index(lon1,lon2)
    synotic=interpolate(synotic,dloninterp,dlatinterp,/grid)


s=size(synotic)
;carr=index.CRLN_OBS 
Carr=GET_STEREO_CARR_ROT(mapr.time) &carr=(long(carr)+1-carr)*360.; carr=ir.CRLN_OBS;

dcarr=(180-carr)/360.*s[1]
synotic=shift(synotic,dcarr)
;lon2=shift(lon2,dcarr)
;slat=[ir.minlat,ir.maxlat]
;slon=[ir.minlon,ir.maxlon]-carr
sm=size(mapr.data)
xpr=get_map_xrange(mapr) & ypr=get_map_yrange(mapr)
xp=linrange(sm[1],xpr[0],xpr[1]);-carr
yp=linrange(sm[2],ypr[0],ypr[1])
;slon=slon-carr

;map2wrap,map, wrap, lat=slat,lon=slon
;map2wrap,mapx,wrapx,lat=slat,lon=slon
;map2wrap,mapy,wrapy,lat=slat,lon=slon
slatn=get_interpolation_index(lat2,slat)
slonn=get_interpolation_index(lon2,slon+180)
slatn[0]=long(slatn[0]+1-1e-5)  &  slatn[1]=long(slatn[1])
slonn[0]=long(slonn[0]+1-1e-5)  &  slonn[1]=long(slonn[1]); sub of synotic (pixel)
dslonn=slonn[1]-slonn[0] & dslatn=slatn[1]-slatn[0]
x1=slonn[0] & x2=slonn[1] & y1=slatn[0] & y2=slatn[1]
;if dslonn gt dslatn then begin
;x1=x1+(dslonn-dslatn)/2.
;x2=x1+dslatn
;endif
;if dslonn lt dslatn then begin
;y1=y1-(dslonn-dslatn)/2.
;y2=y1+dslonn
;endif
slats=linrange(y2-y1+1,lat2[y1],lat2[y2])
slons=linrange(x2-x1+1,lon2[x1]-180,lon2[x2]-180);   lat-lon of sub of synotic  
dlatinterp=get_interpolation_index(yp,slats)  
dloninterp=get_interpolation_index(xp,slons)
sub =interpolate(imr, dloninterp,dlatinterp,/grid)

;stop
thindex=(90-slats)*!dpi/180.
phindex=(carr+slons+360 mod 360)*!dtor


;fff_temp_pre_s, subt,subp,subr, thindex, phindex, suby, subx, sub,te=50,pe=50,re=20,/noplot
;fff_temp_pre1, subp,subt,subr, phindex,thindex, subx, suby, sub

stop1=0
if stop1 then begin
ys=(y2-y1)/((x2-x1)*3)*900
window,1,xs=900,ys=ys
tvscl,congrid(subt,300,ys)>(-500)<500,0
tvscl,congrid(subp,300,ys)>(-500)<500,1
tvscl,congrid(subr,300,ys)>(-500)<500,2
window,2,xs=900,ys=ys
tvscl,congrid(suby,300,ys)>(-500)<500,0
tvscl,congrid(subx,300,ys)>(-500)<500,1
tvscl,congrid(sub ,300,ys)>(-500)<500,2
stop
endif

print,'mean',mean(synotic)
print,mean(abs(sub)),mean(abs(synotic[x1:x2,y1:y2]))
factor=(mean(abs(sub))/mean(abs(synotic[x1:x2,y1:y2])))
print,factor
synotic=synotic*factor
print,mean(synotic)
help,sub
verb=1
if verb then begin
win9,1,xs=500,ys=500,title='sub'
plot_image,sub>(-500)<500
win9,2,xs=500,ys=500,title='synotic[x1:x2,y1:y2]'
plot_image,synotic[x1:x2,y1:y2]>(-500)<500
endif
synotic0=synotic
synotic[x1:x2,y1:y2]=sub
synotic=shift(synotic,-dcarr)
synotic0=shift(synotic0,-dcarr)
if verb then begin
win9,3,title='synotic0'
plot_image,synotic0>(-400)<400
plots,[x1,x2,x2,x1,x1]-dcarr,[y1,y1,y2,y2,y1],/data

win9,4,title='synotic'
plot_image,synotic>(-400)<400
plots,[x1,x2,x2,x1,x1]-dcarr,[y1,y1,y2,y2,y1],/data
stop
endif
theta=theta2
phi=phi2
;;lat2=90-theta2*180/!dpi
;;lon2=(linrange(nlon2+1,0,360))(0:nlon2-1)
;;phi2=lon2*!dtor

;=======define rindex=======================
rgrid=1
dr0=(!dpi/nlat2) ;  r grid spacing at r=1, make it half avg lat grid spacing
rra=[1d0,2.5d0]  ;  range of r
case rgrid of
  2: begin  ;  radial gridpoint separation is proportional to r^2
    rix=[rra(0)]
    lastr=rra(0)
    repeat begin 
      nextr=lastr+dr0*(lastr/rra(0))^2
      rix=[rix,nextr]
      lastr=nextr
    endrep until nextr ge rra(1)
    rix2=rix/((max(rix)-rra(0))/(rra(1)-rra(0)))
    rix=rix2+(rra(0)-rix2(0))
    nr=n_elements(rix)
   ; rix=[linrange(5,rix[0],rix[1]),rix[2:nr-1]]
    nr=n_elements(rix)
    end
else: begin  ;  radial gridpoints uniformly spaced
    nr=round((rra(1)-rra(0))/dr0)
    rix=linrange(nr,rra(0),rra(1))
  ; rix=[linrange(5,rix[0],rix[1]),rix[2:nr-1]]
 nr=n_elements(rix)
    end
endcase

rindex=rix
bbox=[long(min(phindex)/!dtor+1-1e-5),long(slat[0]+1-1e-5),long(max(phindex)/!dtor),long(slat[1])]
;==============================
;filename=path2name(hmir)+'_'+strtrim(long(rebin),2)+'.sav'
;print,'write sav file: '+filename
;save,synotic,theta,phi,rindex,thindex,phindex,bbox,dcarr,x1,x2,y1,y2,sub,subx,suby,filename=filename
;end
sav=path2name(hmir[j])+'.sav';'_'+strtrim(long(rebin),2)+'.sav'









;pfss_sav_local.pro <=============================================
;sav='hmi.B_720s_e15w1332_cutout.11429.20120306_000000_TAI.br.fits_2.sav'
;restore,sav,/verb
theta0=theta
phi0=phi

mag=synotic
;  next get PFSS coefficients
rss=2.5  ;  source surface radius
pfss_get_potl_coeffs,mag,rtop=rss

;  next reconstruct the coronal field in a spherical shell between 1 and rss

pfss_potl_field,2.5,3,rindex=rindex,thindex=thindex,phindex=phindex;,/trunc;,/quiet

;  record some diagnostic info from this extrapolation
cth=cos(theta2)
;monopole=(mean_dtheta(total(mag,1),cth)/n_elements(lon2))(0)
cth=cos(theta)
monopole=(mean_dtheta(total((br(*,*,0)),1),cth)/nlon)(0)
surfflux=mean_dtheta(total((abs(br(*,*,0))),1),cth)*nlat*1e18*rix(0)^2
openflux=mean_dtheta(total((abs(br(*,*,nr-1))),1),cth)*nlat*1e18*rix(nr-1)^2
monfrac=monopole*nlon*nlat*1e18/[surfflux,openflux]
r2inv=rebin(reform(1./rix^2,1,1,nr),nlon,nlat,nr)
;br=br-float(monopole*r2inv)  ;  removes monopole from coronal field

;  print some of the diagnostic info
print,'monopole = ',monopole
print,'unsigned flux = ',surfflux
print,'open flux =     ',openflux
print,'monopole fractions = ',monfrac

;  at this point field lines can be drawn, etc., as in pfss_sample1.pro
;=============================
srindex=rindex & stheta=theta & sphi=phi &
 brp=br & bthp=bth & bphp=bph 
filename='pfss_'+sav
print,filename
save,brp,bthp,bphp,rindex,theta,phi,carr,mapr,mapt,mapp,filename='pfss_'+strtrim(rebin)+'.sav'
;pfss2nlff2,sav=filename,rsize=250,tsize=200,vsize=200;=============================================================
print,j,'/',n-1
endfor
end
;nlat=n_elements(theta0)
;nlon=n_elements(phi0)
;nr=n_elements(rindex)
;br0=fltarr(nlon,nlat,nr)
;bth0=br0 & bph0=br0

;for i=0,nr-1 do begin
 ;br0[x1:x2,y1:y2,i]=br[*,*,i]  & br0[*,*,i]=shift(br0[*,*,i],-dcarr)
;bth0[x1:x2,y1:y2,i]=bth[*,*,i] & bth0[*,*,i]=shift(bth0[*,*,i],-dcarr)
;bph0[x1:x2,y1:y2,i]=bph[*,*,i] & bph0[*,*,i]=shift(bph0[*,*,i],-dcarr)

;endfor

;br=br0 & bth=bth0 & bph=bph0 & theta=theta0 & phi=phi0
;lat=90-theta*180./!dpi  
;lon=phi/!dtor
;filename='pfss_'+path2name(sav)
;print,'write sav file: '+filename
;save,br,bth,bph,nr,nlat,nlon,rix,lat,lon,theta,phi,phiat,phibt,l0,b0,bbox,$
;theta0,phi0,subx,suby,srindex,stheta,sphi,brp,bthp,bphp,x1,x2,y1,y2,dcarr,carr,filename=filename
;endfor

;end
