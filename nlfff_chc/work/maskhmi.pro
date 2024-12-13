PRO maskhmi,br,u,filename

br2=reform(br)
help,br2
u0=where(br2 eq 0.,nu0)
u1=where(abs(br2) ge 150.)
u2=where((abs(br2) lt 150.) and (abs(br2) ne 0.))
print,'number of br2=',n_elements(br2)
print,'number of u0=',n_elements(u0)
print,'number of u1=',n_elements(u1)
print,'number of u2=',n_elements(u2)
print,'number of u0+u1+u2=',(n_elements(u0)+n_elements(u1)+n_elements(u2))
print,'max of br2=',max(br2)
help,u0
s=size(br2)
nt=s[1]
np=s[2]
w=fltarr(nt,np)
;w[u0]=0.0
w[u1]=1.0
w[u2]=(abs(br2[u2])/150.)^2
;;
;;;;;;;; File mask.dat;;;;;;;;;
openw,u,filename

FOR iy=0,np-1 DO BEGIN
      FOR ix=0,nt-1 DO BEGIN
        printf,u, w[ix,iy]
    
      ENDFOR
    ENDFOR
close,u
print,'File mask.dat generated!!'
;;;;;;;;;;;;;;;;;;;;;;;
;;
help,w
print,'========================='
print,'w_min=',min(w[u2])
print,'max_br2u2=',max(br2[u2])
print,'max_br2u2=',min(br2[u2])
x=findgen(abs(max(br2)))
;plot,br2,w

END
