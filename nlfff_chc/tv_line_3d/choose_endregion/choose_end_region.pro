;pro choose_end_region
restore,'../Bout__2015.0221.180000.rtp.sav',/verb
im=BR[*,*,0]

plot_image,im,min=-100,max=100
s=size(im) & roi=fltarr(s[1],s[2])

sub_region=defroi_chc(s[1],s[2],x,y,/nofill,/RESTORE)
roi[sub_region]=1


plot_image,im*roi,min=-100,max=100

end
