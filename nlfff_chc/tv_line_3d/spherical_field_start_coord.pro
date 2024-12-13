;+
;  NAME: spherical_field_start_coord
; 
;  PURPOSE:
;    This procedure chooses points through which to trace fieldlines
;    of the vector field.
;  
;  CALLING SEQUENCE:
;    spherical_field_start_coord,sph_data,fieldtype,spacing,bbox=bbox,
;      radstart=radstart,/add
; 
;  INPUTS: 
;    sph_data = a structure of type spherical_field_data (see
;      spherical_field_data__define.pro) with the following fields defined on
;      input: br,bth,bph,nlat,nlon,nr,rix,thix,phix,lat,lon,str,stth,stph.
;      Basically, one needs the vector field (br,bth,bph), its dimension
;      (nr,nlat,nlon), and its indexing (rix,thix,phix,lat,lon).  If a set of
;      starting points is already defined in the structure (str,stth,stph),
;      the new fieldlines will either overwrite or be added on to the existing
;      set, depending on whether the /add flag is set.
;    fieldtype = 1 = starting points fall along the equator
;                2 = uniform grid, with a random offset
;                3 = points are distributed randomly in latitude and longitude
;                4 = read in from a file (not implemented)
;                5 = uniform grid (default)
;                6 = points are weighted by radial flux at the start radius
;                7 = rectangular grid
;    spacing = controls density of points, does different things depending on 
;              fieldtype
;    bbox = either [lon1,lon2] (for fieldtype1) or [lon1,lat1,lon2,lat2] (for
;           fieldtypes 2,3,5,6,7) defining bounding box (in degrees) outside
;           of which no fieldline starting points lie
;    radstart = a scalar equal to the radius at which all fieldlines should
;               start (default=minimum radius in the domain)
;    add = set this flag if the new starting points are to be added to the
;          existing set already defined in sph_data (default is to overwrite)
; 
;  OUTPUTS:
;    sph_data = a structure of type spherical_field_data (see
;      spherical_field_data__define.pro) with the following fields set on
;      output: str,stth,stph.  These are the r-, theta-, and phi-components of
;      the fieldline starting points.
;
;  NOTES:
;    1.  fieldtype=4 option is not implemented.
; 
;  MODIFICATION HISTORY:
;    M.DeRosa - 13 Dec 2005 - created from guts of pfss_field_start_coord.pro
;               19 Dec 2005 - now handles "wraparound" bounding boxes,
;                             i.e. bounding boxes that span the 0 meridian
;               24 Jan 2006 - now deals with bounded datasets
;               14 Dec 2009 - added fieldtype 7
;               23 Dec 2009 - added latitudinal boundary check
; 
;-

pro spherical_field_start_coord,sph_data,fieldtype,spacing,bbox=bbox, $
  radstart=radstart,add=add
help,*sph_data.lat
;  print usage message
if n_params() eq 0 then begin
  print,'  ERROR in spherical_field_start_coord: no input data provided'
  print,'  Calling sequence: spherical_field_start_coord,sph_data,'+ $
    'fieldtype,spacing,bbox=bbox,radstart=radstart,/add'
  return
endif

;  get rmin,rmax    
rmin=min(*sph_data.rix,max=rmax)

;  set radstart if not set
if n_elements(radstart) eq 0 then radstart=rmin else radstart=radstart(0)

case fieldtype of

  1: begin  ;  points along the equator

    npt=round(spacing(0))  ;  spacing = number of points
    str=replicate(radstart,npt)
    stth=replicate(!dpi/2,npt)
    if n_elements(bbox) eq 2 then begin
      phmin=bbox(0)
      phmax=bbox(1)
      if phmax lt phmin then phmax=phmax+360
      stph=linrange(npt,phmin,phmax)
      stph=((stph+360) mod 360) * !dpi/180
    endif else begin
      stph=(linrange(npt+1,0,360))(0:npt-1)
      stph=stph*!dpi/180
    endelse
     
    end

  2: begin  ;  uniform grid, with a random offset

    ;  set up binning
    nlatbin=round(double(sph_data.nlat)/spacing(0))
    dlatbin=!dpi/nlatbin
    latbin=linrange(nlatbin,dlatbin/2,!dpi-dlatbin/2)
    nlonbin=round(nlatbin*2*sin(latbin))
    dlonbin=2*!dpi/nlonbin
    nloncum=lonarr(nlatbin)
    for i=0l,nlatbin-1 do nloncum(i)=total(nlonbin(0:i))
    npt=round(total(nlonbin))  ;  total number of points

    ;  add random offsets
    stth=randomu(seed,npt)-0.5
    stph=randomu(seed,npt)-0.5
    for i=0l,npt-1 do begin
      lonbinix=(where(i lt nloncum))(0)
      stth(i)=latbin(lonbinix)+stth(i)*dlatbin/2
      stph(i)=(i-(nloncum(lonbinix)-nlonbin(lonbinix)))*dlonbin(lonbinix) + $
        stph(i)*dlonbin(lonbinix)/2
    endfor
    stph=(stph+2*!dpi) mod (2*!dpi)
    
    if n_elements(bbox) eq 4 then begin  ;  remove all points outside box

      ;  filter in latitude
      whlat=where((stth ge (90-bbox(3))*!dpi/180) $
        and (stth le (90-bbox(1))*!dpi/180),nwh)
      if nwh gt 0 then begin
        stth=stth(whlat)
        stph=stph(whlat)
      endif

      ;  filter in longitude
      lon1=((bbox(0)+360) mod 360)*!dpi/180
      lon2=((bbox(2)+360) mod 360)*!dpi/180
      if (lon2 gt lon1) then begin
        whlon=where((stph ge lon1) and (stph le lon2),nwh)
      endif else begin
        whlon=where((stph le lon2) or (stph ge lon1),nwh)
      endelse
      if nwh gt 0 then begin
        stth=stth(whlon)
        stph=stph(whlon)
      endif

    endif

    ;  set radial starting points
    str=replicate(radstart,n_elements(stth))

    end

  3: begin  ;  choose locations at random

    npt=round(spacing(0))  ;  number of points
    str=replicate(radstart,npt)
    if n_elements(bbox) eq 4 then begin
      thmin=(90-bbox(3))*!dpi/180
      thmax=(90-bbox(1))*!dpi/180
      stth=randomu(seed,npt)*(thmax-thmin)+thmin
      phmin=bbox(0)
      phmax=bbox(2)
      if phmax lt phmin then phmax=phmax+360
      stph=((randomu(seed,npt)*(phmax-phmin)+phmin) mod 360) * !dpi/180
    endif else begin
      thmin=min(*sph_data.theta,max=thmax)
      stth=randomu(seed,npt)*(thmax-thmin)+thmin
      stph=randomu(seed,npt)*2*!dpi
    endelse

    end

  4: begin  ;  read in locations from a file

     ; print,'  WARNING in spherical_field_start_coord: fieldtype=4 not '+ $
      ;  'implemented.'
ads=[0,3,3,3]
sp=[-1,1]*2.
adz=linrange(ads[3],sp[0],sp[1])
adz3=rebin(adz,ads[3],ads[1],ads[2]) &adz3=transpose(adz3,[1,2,0])

ady=linrange(ads[2],sp[0],sp[1])
ady3=rebin(ady,ads[2],ads[1],ads[3]) &ady3=transpose(ady3,[1,0,2])

adx=linrange(ads[1],sp[0],sp[1])
adx3=rebin(adx,ads[1],ads[2],ads[3])

adx3=rebin(adx3,ads[1]*ads[2]*ads[3])
adz3=rebin(adz3,ads[1]*ads[2]*ads[3])
ady3=rebin(ady3,ads[1]*ads[2]*ads[3])
help,spacing
if (size(spacing))[0] eq 1 then begin
str0 =spacing[0]+adx3
stth0=spacing[1]+ady3
stph0=spacing[2]+adz3
endif
if (size(spacing))[0] eq 2 then begin
str0 =spacing[0,0]+adx3
stth0=spacing[1,0]+ady3
stph0=spacing[2,0]+adz3
for j=1,(size(spacing))[2]-1 do begin
str0 =[str0,  spacing[0,j]+adx3]
stth0=[stth0, spacing[1,j]+ady3]
stph0=[stph0, spacing[2,j]+adz3]
endfor
endif
help,str0
str=interpolate(*sph_data.rix,str0)
stth=(90-interpolate(*sph_data.lat,stth0))/!radeg
stph=interpolate(*sph_data.lon,stph0)/!radeg
    end

  5: begin  ;  uniform grid

    ;  set up binning
    nlatbin=round(double(sph_data.nlat)/spacing(0))
    dlatbin=!dpi/nlatbin
    latbin=linrange(nlatbin,dlatbin/2,!dpi-dlatbin/2)
    nlonbin=round(nlatbin*2*sin(latbin))
    dlonbin=2*!dpi/nlonbin
    nloncum=lonarr(nlatbin)
    for i=0l,nlatbin-1 do nloncum(i)=total(nlonbin(0:i))
    npt=round(total(nlonbin))  ;  total number of points

    ;  set stth,stph
    stth=dblarr(npt,/noz)
    stph=dblarr(npt,/noz)
    for i=0l,npt-1 do begin
      lonbinix=(where(i lt nloncum))(0)
      stth(i)=latbin(lonbinix)
      stph(i)=(i-(nloncum(lonbinix)-nlonbin(lonbinix)))*dlonbin(lonbinix)
    endfor

    if n_elements(bbox) eq 4 then begin  ;  remove all points outside box

      ;  filter in latitude
      whlat=where((stth ge (90-bbox(3))*!dpi/180) $
        and (stth le (90-bbox(1))*!dpi/180),nwh)
      if nwh gt 0 then begin
        stth=stth(whlat)
        stph=stph(whlat)
      endif

      ;  filter in longitude
      lon1=((bbox(0)+360) mod 360)*!dpi/180
      lon2=((bbox(2)+360) mod 360)*!dpi/180
      if (lon2 gt lon1) then begin
        whlon=where((stph ge lon1) and (stph le lon2),nwh)
      endif else begin
        whlon=where((stph le lon2) or (stph ge lon1),nwh)
      endelse
      if nwh gt 0 then begin
        stth=stth(whlon)
        stph=stph(whlon)
      endif

    endif

    ;  set radial starting points
    str=replicate(radstart,n_elements(stth))

    end

  6:  begin  ;  flux-based

    ;  get random starting points, oversample by a factor of 10
    npt=round(spacing(0))
    oversampling=10l  ;  oversampling factor
    if n_elements(bbox) eq 4 then begin
      thmin=(90-bbox(3))*!dpi/180
      thmax=(90-bbox(1))*!dpi/180
      stth=randomu(seed,npt*oversampling)*(thmax-thmin)+thmin
      phmin=bbox(0)*!dpi/180
      phmax=bbox(2)*!dpi/180
      stph=randomu(seed,npt*oversampling)*(phmax-phmin)+phmin
    endif else begin
      thmin=min(*sph_data.theta,max=thmax)
      stth=randomu(seed,npt*oversampling)*(thmax-thmin)+thmin
      stph=randomu(seed,npt*oversampling)*2*!dpi
    endelse
    str=replicate(radstart,npt*oversampling)

    ;  get br at the starting radius of the random points
    ir=get_interpolation_index(*sph_data.rix,str)
    ith=get_interpolation_index(*sph_data.lat,90-stth*180/!dpi)
    iph=get_interpolation_index(*sph_data.lon,stph*180/!dpi)
    br0=interpolate(*sph_data.br,iph,ith,ir)

    ;  reverse sort b values (should we take into account the unequal areas?)
    magix=reverse(sort(abs(br0)))
    str=str(magix(0:npt-1))
    stth=stth(magix(0:npt-1))
    stph=stph(magix(0:npt-1))

    end

  7:  begin  ;  rectangular grid 

    ;  set up binning
    nlatbin=round(double(sph_data.nlat)/spacing(0))
    dlatbin=!dpi/nlatbin
    latbin=linrange(nlatbin,dlatbin/2,!dpi-dlatbin/2)
    nlonbin=nlatbin*2
    lonbin=(linrange(nlonbin+1,0,2*!dpi))(0:nlonbin-1)
    npt=nlonbin*nlatbin  ;  total number of points

    ;  set stth,stph
    stph=reform(lonbin#replicate(1,nlatbin),npt,/overwrite)
    stth=reform(replicate(1,nlonbin)#latbin,npt,/overwrite)

    if n_elements(bbox) eq 4 then begin  ;  remove all points outside box

      ;  filter in latitude
      whlat=where((stth ge (90-bbox(3))*!dpi/180) $
        and (stth le (90-bbox(1))*!dpi/180),nwh)
      if nwh gt 0 then begin
        stth=stth(whlat)
        stph=stph(whlat)
      endif

      ;  filter in longitude
      lon1=((bbox(0)+360) mod 360)*!dpi/180
      lon2=((bbox(2)+360) mod 360)*!dpi/180
      if (lon2 gt lon1) then begin
        whlon=where((stph ge lon1) and (stph le lon2),nwh)
      endif else begin
        whlon=where((stph le lon2) or (stph ge lon1),nwh)
      endelse
      if nwh gt 0 then begin
        stth=stth(whlon)
        stph=stph(whlon)
      endif

    endif

    ;  set radial starting points
    str=replicate(radstart,n_elements(stth))

    end

  else: begin
    print,'  ERROR in spherical_field_start_coord:  invalid fieldtype'
    return
    end

endcase

;  if bounded, select only those points that lie inside longitude boundaries
if n_elements(sph_data.lonbounds) gt 0 then begin
  if sph_data.lonbounds(0) ge 0 then begin
    lonb2=sph_data.lonbounds*!dpi/180
    if lonb2(1) gt lonb2(0) then begin
      wh=where((stph ge lonb2(0)) and (stph le lonb2(1)),nwh)
    endif else begin
      wh=where((stph ge lonb2(0)) or (stph le lonb2(1)),nwh)
    endelse
    if nwh eq 0 then begin
      print,'ERROR in spherical_field_start_coord: no points in bounds'
      return
    endif else begin
      str=str(wh)
      stth=stth(wh)
      stph=stph(wh)
    endelse
  endif
endif

;  select only those points that lie inside latitudinal boundaries of domain
thmin=min(*sph_data.theta,max=thmax)
whlat=where((stth ge thmin) and (stth le thmax),nwh)
if nwh gt 0 then begin
  stph=stph(whlat)
  stth=stth(whlat)
  str=str(whlat)
endif

;  assign pointers, nullify if already valid
if keyword_set(add) then begin
  addvalid=1
  if ptr_valid(sph_data.str) then oldstr=*sph_data.str else addvalid=0
  if ptr_valid(sph_data.stth) then oldstth=*sph_data.stth else addvalid=0
  if ptr_valid(sph_data.stph) then oldstph=*sph_data.stph else addvalid=0
endif else addvalid=0
if ptr_valid(sph_data.str) then ptr_free,sph_data.str
if ptr_valid(sph_data.stth) then ptr_free,sph_data.stth
if ptr_valid(sph_data.stph) then ptr_free,sph_data.stph
if keyword_set(add) and addvalid then begin
  sph_data.str=ptr_new([oldstr,str])
  sph_data.stth=ptr_new([oldstth,stth])
  sph_data.stph=ptr_new([oldstph,stph])    
endif else begin
  sph_data.str=ptr_new(str)
  sph_data.stth=ptr_new(stth)
  sph_data.stph=ptr_new(stph)    
endelse

end
