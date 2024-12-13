;;@pfss_data_block
;===============================================[


;restore,'Bout__2015.0203.104800.rtp.sav',/verb
restore,'Bout__2015.0221.180000.rtp.sav',/verb
;restore,'Bout__2015.0221.140000.rtp.sav',/v


window,2,xs=700,ys=600
;sav=['null_fan.sav','null_spine.sav'];,'nullv2.sav']
;sav=['start_rtp_ht_1048.sav']
sav=['start_rtp_ht_140000.sav']
im=bytscl(br[*,*,0],min=-500,max=500)
xrange=[1,800]
yrange=[1,384]
gx=200
gy=96



radstart=1.00 ;unit :rsun
fieldtype=2 & invdens =1.5  ;number
fieldtype=6 & invdens =100 ;number
;fieldtype=5 & invdens =1.5  ;density
bbox=[min(lon),min(lat),max(lon),max(lat)];===================

map=make_map(im,xc=mean(lon),yc=mean(lat),dx=lon[1]-lon[0],dy=abs(lat[1]-lat[0]))


bbox=[-38,5,-23,10]

nlat=n_elements(lat) & nlon=n_elements(lon) & nr=n_elements(rix)
ax=25 & az=35;======================================

;ax=40 & az=160
;ax=80 & az=80
ax=90 & az=0
;ax=90 & az=48
scale3,ax=ax,az=az,xrange=lon[[0,nlon-1]],yrange=lat[[0,nlat-1]],zrange=rix[[0,nr*0.5-1]]

c1=convert_coord(min(lon),min(lat),min(rix),/data,/to_device,/t3d)
c2=convert_coord(max(lon),min(lat),min(rix),/data,/to_device,/t3d)
c3=convert_coord(max(lon),max(lat),min(rix),/data,/to_device,/t3d)
c4=convert_coord(min(lon),max(lat),min(rix),/data,/to_device,/t3d)
x0=([min(lon),max(lon),max(lon),min(lon)]-min(lon))/(max(lon)*1.-min(lon))*nlon 
y0=([min(lat),min(lat),max(lat),max(lat)]-min(lat))/(max(lat)*1.-min(lat))*nlat
x1=[c1[0],c2[0],c3[0],c4[0]] & x1=(x1-min(x1))/(max(x1)*1.-min(x1))*nlon
y1=[c1[1],c2[1],c3[1],c4[1]] & y1=(y1-min(y1))/(max(y1)*1.-min(y1))*nlat
 polywarp,x0,y0,x1,y1,1,p,q
 b=poly_2d(im,p,q,mis=255)
loadct,0

;plots,lon[[0,nlon-1,nlon-1,0,0]],lat[[0,0,nlat-1,nlat-1,0]],[1,1,1,1,1],/t3d,col=255,thick=5
if az ge   0 and az lt  90 then size_xy=[c2[0]-c4[0],c3[1]-c1[1],c4[0],c1[1]]
if az gt  90 and az lt 180 then size_xy=[c1[0]-c3[0],c2[1]-c4[1],c3[0],c4[1]]
if az gt 180 and az lt 270 then size_xy=[c4[0]-c2[0],c1[1]-c3[1],c2[0],c3[1]]
if az gt 270 and az lt 360 then size_xy=[c3[0]-c1[0],c4[1]-c2[1],c1[0],c2[1]]
imb=congrid(b,size_xy[0],size_xy[1])
print,size_xy

tv,imb,size_xy[2],size_xy[3],/dev
;----------------------------------------------------------------------------


nt=0
colors=[110,55,0]
thicks=[1,2,3]
linestyles=[0,1,1]
for jf=0,n_elements(sav)-1 do begin

;if jf eq 0 then begin
fieldtype=4 & invdens =sav[jf]
;endif 
col=100+jf*100
loadct,5;18-jf*13


restore,invdens
;===========
sb=size(br)
fl=10
fl=fl*1.
gx=linrange(10,10,sb[1]-10);findgen(sb[1]/fl)*fl
gy=linrange(5,10,sb[2]-10);findgen(sb[2]/fl)*fl

;sx=reform(     gx#(fltarr(n_elements(gy))+1) ,n_elements(gx#(fltarr(n_elements(gy))+1)))
;sy=reform(    (fltarr(n_elements(gx))+1 )#gy,n_elements(sx))
;sz=sx*0



;=============
sxp=interpolate(lon,sx)
syp=interpolate(lat,sy)
szp=interpolate(rix,sz)



plots,sxp,syp,szp,psym=1,/t3d,col=110



;line=trace_flines_xyz(bph,-bth,br,sx,sy,sz)
line=supper_trace_flines_xyz(bph,-bth,br,sx,sy,sz)
if jf eq 1 then line=trace_flines_xyz(bph,-bth,br,sx,sy,sz,def_sign=-1);(1 or -1);======================note!!!!
newsx=sx*0



if jf eq 0 then flag=0     
for i=0L,line.n-1 do begin
lx=*line.x[i]
ly=*line.y[i]
lz=*line.z[i]
ln=n_elements(lx)
x=interpolate(lon,*line.x[i])
y=interpolate(lat,*line.y[i])
z=interpolate(rix,*line.z[i])
np=n_elements(x)-1 ;&z=1+(z-1)*5

bzp1=interpolate(br,lx[0],ly[0],lz[0])
bzp2=interpolate(br,lx[ln-1],ly[ln-1],lz[ln-1])



maxr=max(z)
;plot1
if max(z) gt 1.05 then begin 
if  bzp2 lt -500 then begin
plots,x,y,z,/t3d,color='00ff00'xl
endif
endif
;plot2
if max(z) lt 1.01 then begin 
if  bzp2 gt -200 && bzp2 le -150 then begin
plots,x,y,z,/t3d
endif
endif


endfor
endfor

end
