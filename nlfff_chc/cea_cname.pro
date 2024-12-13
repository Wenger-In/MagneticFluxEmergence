function anytime2name,at
if strpos(at,'_TAI') gt 0 then at=strmid(at,0,strpos(at,'_TAI') )
tai=anytim2tai(at)
utc=tai2utc(tai,/ext)
return,strtrim(utc.year,2)+'.'+$
    strmid(strtrim(utc.month/100.,2),2,2)+$
    strmid(strtrim(utc.day/100.,2),2,2)+'.'+$
    strmid(strtrim(utc.hour/100.,2),2,2)+$
    strmid(strtrim(utc.minute/100.,2),2,2)+$
    strmid(strtrim(utc.second/100.,2),2,2)
end

pro cea_cname,fn,mv=mv,step=step
if not keyword_set(step) then step=1
fn=file_search(fn,count=n)
for i=0L,n-1,step do begin
;h=headfits(fn[i],/ext)
;tt=sxpar(h,'T_REC')
read_sdo,fn[i],index
yymm=anytime2name(index.T_REC)
f=fn[i]
name_end=strmid(f,strpos(strmid(f,0,strpos(f,'.',/reverse_s)),'.',/reverse_s),50)
name=strmid(f,0,rstrpos(f,'/')+1)
name=name+'cea_'+yymm+name_end;---------------------<<<  note
;===============================================================

if keyword_set(mv) then spawn,['mv',fn[i],name],/noshell
print,name
cc:
print,'now',i,'   /',n-1
endfor
end
