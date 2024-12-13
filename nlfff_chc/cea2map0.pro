pro cea2map0,f1,map
read_sdo,f1,index,data,/UNCOMP_DELETE
xc=index.xcen  -index.crln_obs              ;*(index.londtmax+index.londtmin)/2.
yc=index.ycen                         ;*(index.latdtmax+index.latdtmin)/2.
dx=index.cdelt1                     ;*(index.londtmax-index.londtmin)/index.naxis1
dy=index.cdelt2                   ;*(index.latdtmax-index.latdtmin)/index.naxis2

map=make_map(data,xc=xc,yc=yc,dx=dx,dy=dy,time=index.T_REC,b0=index.crlt_obs,l0=index.crln_obs,rsun=960)
end
