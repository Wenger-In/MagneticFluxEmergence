pro cea2map,f1,map
read_sdo,f1,index,data,/UNCOMP_DELETE
wcs=fitshead2wcs(struct2fitshead(index)) & wcs.simple=1
wcs2map,data,wcs,map
end
