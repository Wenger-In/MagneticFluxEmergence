;goto,s2
;file=file_search('Bout__2015.0203.104800.rtp.sav',count=fn)
;file=file_search('Bout__2015.0221.180000.rtp.sav',count=fn)
file=file_search('Bout__2015.0221.140000.rtp.sav',count=fn)

for fi=0,fn-1 do begin

restore,file[fi],/ver & bx=bph & by=-bth & bz=br;==note1
s=size(bz)
;bz=bz[*,*,3:s[3]-1]
;bx=bx[*,*,3:s[3]-1]
;by=by[*,*,3:s[3]-1]
s=size(bz)
xp=[0,s[1]-1]
yp=[0,s[2]-1]
;fl=10.
fl=1.
xp=linrange(fl*(xp[1]-xp[0])+1,xp[0],xp[1]) &nx=n_elements(xp)
yp=linrange(fl*(yp[1]-yp[0])+1,yp[0],yp[1]) &ny=n_elements(yp)
dx=1./fl & dy=1./fl
sx0=xp#replicate(1,ny)
sy0=replicate(1,nx)#yp
sz0=sx0*0
sx=reform(sx0,n_elements(sx0))
sy=reform(sy0,n_elements(sx0))
sz=reform(sz0,n_elements(sx0))
line=trace_manylines_twist(bx,by,bz,sx,sy,sz,nd=1e7)
;line=supper_trace_flines_xyz_bot(bph,-bth,br,sx,sy,sz,/endpoint)
sx=sx0*0 & sy=sx & sz=sx
ex=sx0*0 & ey=ex & ez=ex &zmax=ex-1
tw=sx & tw1=sx & tw2=sx & t=sx & fl=sx
for i=0L,line.n-1 do begin
 sx[i]=line.x[i,0] & sy[i]=line.y[i,0] & sz[i]=line.z[i,0]
 ex[i]=line.x[i,1] & ey[i]=line.y[i,1] & ez[i]=line.z[i,1] & zmax[i]=line.z[i,2]
tw[i]=line.tw[i] & tw1[i]=line.tw1[i] & tw2[i]=line.tw2[i] & t[i]=line.t[i] & fl[i]=line.fl[i]
endfor

s2:
yx=shift(ey,-1,0)
yy=shift(ey,0,-1)
xx=shift(ex,-1,0)
xy=shift(ex,0,-1)

yx1=shift(ey,1,0)
yy1=shift(ey,0,1)
xx1=shift(ex,1,0)

xy1=shift(ex,0,1)

bs=interpolate(bz[*,*,0],sx0,sy0)
be=interpolate(bz[*,*,0],ex,ey)

   ok0=where(abs(bs) eq 0 or abs(be) eq 0 or   $
            ez gt 2 or $
shift(ez,-1,0) gt 2 or $
shift(ez,0,-1) gt 2 or $
shift(ez,1,0)  gt 2 or $
shift(ez,0,1)  gt 2 or $
            zmax le 3 or $
shift(zmax,-1,0) le 3  or $
shift(zmax,0,-1) le 3  or $
shift(zmax,1,0)  le 3  or $
shift(zmax,0,1) le 3  or $
abs(bs) lt 10 or shift(abs(bs),1,0) lt 10 or shift(abs(bs),0,1) lt 10 or  $
                 shift(abs(bs),-1,0) lt 10 or shift(abs(bs),0,-1) lt 10 or $ 
abs(be) lt 10 or shift(abs(be),1,0) lt 10 or shift(abs(be),0,1) lt 10 or $
                  shift(abs(be),-1,0) lt 10 or shift(abs(be),0,-1) lt 10  );maybe should consider the rmax if not keyword_set(endpoint), it means discard the lines that two low! 
fac=bs*0+1 & fac[[0,nx-1],*]=0 & fac[*,[0,ny-1]]=0
fac[ok0]=0 & okf=where(fac ne 0 )
bb=fac*0 & bb[okf]=abs(be[okf]/bs[okf])
n=(((xx -ex)/dx)^2+((xy -ex)/dx)^2+$
   ;((xx1-ex)/dx)^2+((xy1-ex)/dx)^2+$
   ((yx -ey)/dy)^2+((yy -ey)/dy)^2);+$
   ;((yx1-ey)/dy)^2+((yy1-ey)/dy)^2)
n=n*fac
imq=n*bb

wn='imq/40/imq_'+path2name(file[fi]) 
wn='./imq_'+path2name(file[fi]) 
save,xp,yp,sx,sy,sz,ex,ey,ez,zmax,bs,be,imq,t,tw,tw1,tw2,fl,filename=wn
print,'save,    ',wn
endfor
end

