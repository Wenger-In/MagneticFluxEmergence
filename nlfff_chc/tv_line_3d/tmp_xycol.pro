pro tmp_xycol,x0,y0,col0,x,y,col,n,xyc,colc

if n eq 0 then begin
x=x0
y=y0
col=col0
n=1
xyc=[-1,n_elements(x)-1]
colc=[-1,n_elements(col)-1]
endif

if n gt 0 then begin
x=[reform(x),reform(x0)]
y=[reform(y),reform(y0)]
col=[reform(col),reform(col0)]
n=n+1
xyc=[xyc,n_elements(x)-1]
colc=[colc,n_elements(colc)-1]
endif

end
;x[(xyc[i]+1):xyc[i+1]]
;col[(colc[i]+1):colc[i+1]]



;ps3d={n:nok,np:fltarr(nok),x:ptrarr(nok),y:ptrarr(nok),col:ptrarr(nok),im:img,imx:imx,imy:imy}
