function ymy,m
bilinear4,m.br00,m.br03,m.br13,m.br10,r1a,r1b,r1c,r1d
bilinear4,m.br10,m.br11,m.br12,m.br13,r2a,r2b,r2c,r2d
bilinear4,m.br01,m.br11,m.br12,m.br02,r3a,r3b,r3c,r3d
bilinear4,m.br00,m.br01,m.br02,m.br03,r4a,r4b,r4c,r4d
bilinear4,m.br00,m.br01,m.br11,m.br10,r5a,r5b,r5c,r5d
bilinear4,m.br03,m.br02,m.br12,m.br13,r6a,r6b,r6c,r6d

bilinear4,m.bt00,m.bt03,m.bt13,m.bt10,t1a,t1b,t1c,t1d
bilinear4,m.bt10,m.bt11,m.bt12,m.bt13,t2a,t2b,t2c,t2d
bilinear4,m.bt01,m.bt11,m.bt12,m.bt02,t3a,t3b,t3c,t3d
bilinear4,m.bt00,m.bt01,m.bt02,m.bt03,t4a,t4b,t4c,t4d
bilinear4,m.bt00,m.bt01,m.bt11,m.bt10,t5a,t5b,t5c,t5d
bilinear4,m.bt03,m.bt02,m.bt12,m.bt13,t6a,t6b,t6c,t6d

bilinear4,m.bp00,m.bp03,m.bp13,m.bp10,p1a,p1b,p1c,p1d
bilinear4,m.bp10,m.bp11,m.bp12,m.bp13,p2a,p2b,p2c,p2d
bilinear4,m.bp01,m.bp11,m.bp12,m.bp02,p3a,p3b,p3c,p3d
bilinear4,m.bp00,m.bp01,m.bp02,m.bp03,p4a,p4b,p4c,p4d
bilinear4,m.bp00,m.bp01,m.bp11,m.bp10,p5a,p5b,p5c,p5d
bilinear4,m.bp03,m.bp02,m.bp12,m.bp13,p6a,p6b,p6c,p6d

ra=[r1a,r2a,r3a,r4a,r5a,r6a]
rb=[r1b,r2b,r3b,r4b,r5b,r6b]
rc=[r1c,r2c,r3c,r4c,r5c,r6c]
rd=[r1d,r2d,r3d,r4d,r5d,r6d]

ta=[t1a,t2a,t3a,t4a,t5a,t6a]
tb=[t1b,t2b,t3b,t4b,t5b,t6b]
tc=[t1c,t2c,t3c,t4c,t5c,t6c]
td=[t1d,t2d,t3d,t4d,t5d,t6d]

pa=[p1a,p2a,p3a,p4a,p5a,p6a]
pb=[p1b,p2b,p3b,p4b,p5b,p6b]
pc=[p1c,p2c,p3c,p4c,p5c,p6c]
pd=[p1d,p2d,p3d,p4d,p5d,p6d]

res=fltarr(3,12)
for i=0,5 do begin
for j=0,2 do begin 
if j eq 0 then begin
a1=ra[i] & b1=rb[i] & c1=rc[i] & d1=rd[i]
a2=ta[i] & b2=tb[i] & c2=tc[i] & d2=td[i]
a3=pa[i] & b3=pb[i] & c3=pc[i] & d3=pd[i]
endif
if j eq 1 then begin
a1=ra[i] & b1=rb[i] & c1=rc[i] & d1=rd[i]
a2=pa[i] & b2=pb[i] & c2=pc[i] & d2=pd[i]
a3=ta[i] & b3=tb[i] & c3=tc[i] & d3=td[i]
endif
if j eq 2 then begin
a1=pa[i] & b1=pb[i] & c1=pc[i] & d1=pd[i]
a2=ta[i] & b2=tb[i] & c2=tc[i] & d2=td[i]
a3=ra[i] & b3=rb[i] & c3=rc[i] & d3=rd[i]
endif

if d1*c2-d2*c1 ne 0 then begin;
m=(d2*b1-d1*b2)/(d1*c2-d2*c1)
n=(d2*a1-d1*a2)/(d1*c2-d2*c1)
a=d1*m
if a ne 0 then begin;
b=b1+c1*m+d1*n
c=a1+c1*n
bac=b*b-4*a*c
if bac gt 0 then begin;
x1=(-b+bac^0.5)/a/2.
x2=(-b-bac^0.5)/a/2.
if (c1+d1*x1) eq 0 and (a1+b1*x1) ne 0 then y1=0
if (c1+d1*x2) eq 0 and (a1+b1*x2) ne 0 then y2=0
if (c1+d1*x1) ne 0 then y1=-(a1+b1*x1)/(c1+d1*x1)
if (c1+d1*x2) ne 0 then y2=-(a1+b1*x2)/(c1+d1*x2)
if x1 ge 0 and x1 le 1 and y1 ge 0 and y1 le 1 then begin
;res[j*3+0,i*2]=x1
;res[j*3+1,i*2]=y1
res[j,i*2]=a3+b3*x1+c3*y1+d3*x1*y1
endif
if x2 ge 0 and x2 le 1 and y2 ge 0 and y2 le 1 then begin
;res[j*3+0,i*2+1]=x2
;res[j*3+1,i*2+1]=y2
res[j,i*2+1]=a3+b3*x2+c3*y2+d3*x2*y2
endif

endif
endif
endif
endfor
endfor
fr=0 & ft=0 & fp=0
r=res[0,*]
wr=where(r ne 0,nwr)
if nwr eq 2 then if r[wr[0]]*r[wr[1]] lt 0 then fr=1

t=res[1,*]
wt=where(t ne 0,nwt)
if nwt eq 2 then if t[wt[0]]*t[wt[1]] lt 0 then ft=1

p=res[2,*]
wp=where(p ne 0,nwp)
if nwp eq 2 then if p[wp[0]]*p[wp[1]] lt 0 then fp=1

yy=0
if fr eq 1 and ft eq 1 and fp eq 1 then yy=1
return, yy
end
