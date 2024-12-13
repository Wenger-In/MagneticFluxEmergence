;restore,'imq_Bout__2015.0203.104800.rtp.sav',/v
;restore,'imq_Bout__2015.0221.180000.rtp.sav',/v
restore,'imq_Bout__2015.0221.140000.rtp.sav',/v

im=imq
im=tw2

rt=1.
wxs=(size(im))[1]/rt & wys=(size(im))[2]/rt
window,xs=3*wxs,ys=3*wys
;tvscl,alog10(congrid(im,wxs,wys)>2)>2
tvscl,congrid(im,wxs*3,wys*3)>(-1.5)<1.5

roi=0
if roi eq 1 then begin
tmpre=defroi(wxs,wys,tmpx,tmpy,/nofill)
tmpim=im*0
tmpim[congrid(tmpre,(size(im))[1],(size(im))[2])]=1
endif else tmpim=im*0+1

nx=n_elements(xp)
ny=n_elements(yp)
sx0=xp#replicate(1,ny)
sy0=replicate(1,nx)#yp

sx=sx0[where(im gt -10.0 and tmpim eq 1)];
sy=sy0[where(im gt -10.0 and tmpim eq 1)]
sz=replicate(1,n_elements(sx))


save,sx,sy,sz,filename='start_rtp_ht_140000.sav'

end
