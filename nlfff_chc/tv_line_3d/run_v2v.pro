

jsoc2map,'/data1/biyi/20110427/fd/304/aia_2011.0427.020008_0304.lev1.fits',map1
fits2map,'/data1/biyi/20110427/stereo/euvi/sb_2011.0427.002530_195_s2048.fits',map2

map1v='earth'  ;  'earth', 'B' or 'A'
map2v='B'       ;  'earth', 'B' or 'A'
;4 inputs
v2v,map1,map1v,map2,map2v,/pbr;  view to view
end
