function path2name,name,ntime=ntime
re=name
ntime=re
for i=0,n_elements(name)-1 do begin
nameb=rstrpos(name[i],'/')
namel=strlen(name[i])
re[i]=strmid(name[i],nameb+1,namel-nameb)
if keyword_set(ntime) then ntime[i]=strmid(re[i],strpos(re[i],'_20'),17)
endfor
return,re
end
