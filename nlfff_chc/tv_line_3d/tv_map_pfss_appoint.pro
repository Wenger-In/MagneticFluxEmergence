@pfss_data_block
;===============================================[
!p.multi=[0,2,2] & !p.position=[0,0,1,1]
;pfss_restore,'Bfield_20110427_000400.sav'
;pfss_restore,'Bfield_20120909_180400.h5'
restore,'/data2/biyi/null_0909/pfss_wrap_sub_hmi_m_2010.1115.143028.fits.sav',/ver
;br =brp
;bth=bthp
;bph=bphp
s=size(br)
nr=s[3]
nlon=s[1]
nlat=s[2]
lat=(!pi/2-theta)/!dtor
lon=phi/!dtor
;rix=rindex
;theta=theta+158./!dtor
en='../dr/211/dr_aia_2011.0427.052548_0211.lev1.fits'
fits2map,'/data2/biyi/null_0909/sub_hmi_m_2010.1115.143028.fits',map    
;ref_time=map.time     
ref_time='2010-11-15 1430:00'                               
                             
;map=get_sub_map(map,xrange=[-1100,-500],yrange=[0,800])

ratio=0.7    ;tv,congrid(map.data,ratio*xs[1],ratio*xs[2]) 
log='n'        
vmin=-100  & vmax=100

;about pfss lines
width=1;  image out to 2.5 R_sun   \  together these keywords produce
mag  =2 ;  magnification factor     /  a 720x720 image (below, in outim)
imsc =100  ;  data values at which image of background magnetogram saturates
radstart=1.5
corl=1.1

fieldtype=2 & invdens =7;  factor inverse to line density, i.e. lower values = more lines
fieldtype=6 & invdens =400

;  usage:  pfss_field_start_coord,fieldtype,spacing,radstart=radstart,/top,
;            bbox=bbox,/add
;          where fieldtype = 1=starting points fall along the equator
;                            2=uniform grid, with a random offset
;                            3=points are distributed randomly in both
;                              latitude and longitude
;                            4=read in from a file
;                            5=uniform grid (default)
;                            6=points are weighted by flux
;                            7=rectangular grid
;===============================================]

;-------------------------
ang=pb0r(ref_time,/arcsec)
   map.roll_angle=ang[0] & map.b0=ang[1] & map.rsun=ang[2]
sout=map.data
Carr = GET_STEREO_CARR_ROT(ref_time) &
 lcent=(long(carr)+1-carr)*360;.+20
;lcent=h.CRLN_OBS
bcent=map.b0 
;----------------------------------------------------------------------
xr=get_map_xrange(map)/60.  & yr=get_map_yrange(map)/60.
hel1=arcmin2hel(xr[0],mean(yr),rsun=map.rsun,b0=0)
hel2=arcmin2hel(mean(xr),yr[0],rsun=map.rsun,b0=0)
hel3=arcmin2hel(xr[1],mean(yr),rsun=map.rsun,b0=0)
hel4=arcmin2hel(mean(xr),yr[1],rsun=map.rsun,b0=0)
bbox=[lcent+hel1[1],hel2[0],lcent+hel3[1],hel4[0]]
;bbox=[105,0,128+10,45]

;---------------------------------
down_color=[0b,255b,0b]
PIL_color=[255b,0b,0b]
up_color=[0b,0b,255b]
low_color=[255b,255b,0b]
high_color=[255b,255b,255b]
;----------------------



sub=[0.,0,0,0];[left,right,down,up]
rsun=1.0*map.rsun 
xr=1.0*get_map_xrange(map)/rsun
yr=1.0*get_map_yrange(map)/rsun
sub[0]=0.5*(width+xr[0])/width
sub[1]=0.5*(width+xr[1])/width
sub[2]=0.5*(width+yr[0])/width
sub[3]=0.5*(width+yr[1])/width


radstart=[1.05,2]
fieldtype=[2,2] & invdens =[4,4]
space=[[7.68,118.42,110.07] ,[10.79,119.76,119.35]]
space=[8.35,241.84,242.55];[7,41,12]+[0,200,230]+1
for itf=0,0 do begin

;[7.68,118.42,110.07]   ^_^
;pfss_field_start_coord1,2,10,radstart=1.2,bbox=bbox
;restore,'nullpos4.sav'
;space=nullpos
;space[0,*]=nullpos[2,*]
;space[1,*]=nullpos[1,*]
;space[2,*]=nullpos[0,*]
;space=space[*,6]
;space=[8.35,241.84,242.55]
pfss_field_start_coord1,4,space,radstart=1.2,bbox=bbox
restore,'line/null.sav'
pfss_trace_field
if n_elements(radstart) gt 1 then begin
if itf eq 0 then begin
ptr0=ptr  &  ptth0=ptth  & ptph0=ptph  & nstep0=nstep
endif else begin
s0=size(ptr0) & s1=size(ptr)
ptr1=fltarr(max([s0[1],s1[1]]),s0[2]+s1[2])
nstep1=[nstep0,nstep]
ptth1=ptr1 & ptph1=ptr1 
for jtf=0,s0[2]-1 do begin
 ptr1[0:nstep1[jtf]-1,jtf]= ptr0[0:nstep1[jtf]-1,jtf]
ptth1[0:nstep1[jtf]-1,jtf]=ptth0[0:nstep1[jtf]-1,jtf]
ptph1[0:nstep1[jtf]-1,jtf]=ptph0[0:nstep1[jtf]-1,jtf]
endfor
for jtf=s0[2],s0[2]+s1[2]-1 do begin
 ptr1[0:nstep1[jtf]-1,jtf]= ptr[0:nstep[jtf-s0[2]]-1,jtf-s0[2]]
ptth1[0:nstep1[jtf]-1,jtf]=ptth[0:nstep[jtf-s0[2]]-1,jtf-s0[2]]
ptph1[0:nstep1[jtf]-1,jtf]=ptph[0:nstep[jtf-s0[2]]-1,jtf-s0[2]]
endfor
ptr=ptr1  &  ptth=ptth1  & ptph=ptph1  & nstep=nstep1
endelse
endif
endfor


                                            
pfss_draw_field_new,bcent=bcent,lcent=lcent,width=width,imsc=imsc,colr=corl,mag=mag,lxyz=lxyz,$
                                                 /inv,height_inv=2.5,drawclose=1,drawopen=1


loadct,0  ;  loadct,3 also looks nice too
tvlct,re,gr,bl,/get
re(251:255)=[down_color[0],PIL_color[0],up_color[0],low_color[0],high_color[0]]
gr(251:255)=[down_color[1],PIL_color[1],up_color[1],low_color[1],high_color[1]]
bl(251:255)=[down_color[2],PIL_color[2],up_color[2],low_color[2],high_color[2]]
tvlct,re,gr,bl
nax=size(outim1,/dim)
window,0,xsiz=nax(0),ysiz=nax(1)
tv,outim1
plots,[sub[0],sub[1],sub[1],sub[0],sub[0]]*nax(0),[sub[2],sub[2],sub[3],sub[3],sub[2]]*nax(1),/device


xpd=lxyz.xpd  &ypd=lxyz.ypd &zpd=lxyz.zpd &ppns=lxyz.ppns &nppns=lxyz.nppns &ncol=lxyz.ncol 
 width=lxyz.width &wopen=lxyz.wopen &wp=lxyz.wp &xpil=lxyz.xpil & zpil=lxyz.zpil 
                                                  xpil0=lxyz.xpil0 & zpil0=lxyz.zpil0 
xs=size(sout)
sout=congrid(sout,ratio*xs[1],ratio*xs[2]) 
xs=size(sout)
mx=xs[1]
my=xs[2]

xpd=mx*(xpd-2.*width*(sub[0]-0.5))/(2.*width*(sub[1]-sub[0]))
zpd=my*(zpd-2.*width*(sub[2]-0.5))/(2.*width*(sub[3]-sub[2]))


xpil[1:(size(xpil))[1]-1,*]= $
mx*(xpil[1:(size(xpil))[1]-1,*]-2.*width*(sub[0]-0.5))/(2.*width*(sub[1]-sub[0]))
len=40
zpil=my*(zpil-2.*width*(sub[2]-0.5))/(2.*width*(sub[3]-sub[2]))
window,1,xs=xs[1],ys=xs[2]
tvlct,re,gr,bl,/get



if log eq 'y' then begin
sout=alog10(1.*sout>10);--------------------------
sout=bytscl(sout,min=min(sout),max=max(sout),top=248-0.)+0.
endif else $
sout=bytscl(sout,min=vmin,max=vmax,top=248-0.)+0.
;tv,sout     
map=rep_tag_value(map,sout,'data') 
plot_map,map        
;t3d,/reset


pixel2data,xpd,zpd,map,xpd,zpd
for i=0L,nppns do begin
xpp=xpd(0:ppns(i)-1,i)
ypp=ypd(0:ppns(i)-1,i)
zpp=zpd(0:ppns(i)-1,i)
col=ncol(i)
;pixel2data,xpp,zpp,map,xpp,zpp
plots,xpp,zpp,col=col,thick=thick,/data
endfor
sizexpil=size(xpil)
if sizexpil[0] eq 1 then sizexpil[2]=1
for j=0L,sizexpil[2]-1 do begin
n1=strtrim(long(j),2)
for i=1L,xpil[0,j],10 do begin
n2=strmid(strtrim(i/10000.,2),2,4)
xyouts,xpil[i,j],zpil[i,j],n1+'.'+n2,col=252b,/device
plots,xpil[1:xpil[0,j],j],zpil[1:xpil[0,j],j],col=252b,psym=1,/device
endfor
endfor


;----------------------------
 ;re=plinep(map,'/data1/biyi/2011.0606/sav/run_scc/0606_p.dat',0)
; re=plinep(map,'/data1/biyi/2011.0606/sav/run_scc/all_p.dat',0)
;-----------------------------
;window,2,xs=xs[1],ys=xs[2]
;outim2=outim1(nax[0]*sub[0]:nax[0]*sub[1],nax[1]*sub[2]:nax[1]*sub[3]);
;tv,congrid(outim2,xs[1],xs[2])
;for i=0,nppns do begin
;xpp=xpd(0:ppns(i)-1,i)
;ypp=ypd(0:ppns(i)-1,i)
;zpp=zpd(0:ppns(i)-1,i)

;col=ncol(i)
;plots,xpp,zpp,col=255b,thick=0.1,/device


;endfor
print,'lcent:',lcent
print,'bcent:',bcent
;pfss_to_spherical,pfss_data
;spherical_trackball_widget,pfss_data
;twidget,pfss_data
end
