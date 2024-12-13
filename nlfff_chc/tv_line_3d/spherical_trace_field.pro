;+
;  NAME: spherical_trace_field
; 
;  PURPOSE: 
;    Given vector field data (such as the magnetic field) on a spherical grid,
;    and a set of starting points, this procedure will trace the fieldlines
;    that pass through each of the starting points.
; 
;  CALLING SEQUENCE: 
;    spherical_trace_field,sph_data,stepmax=stepmax,safety=safety,
;      outfield=outfield,linelengths=linelengths,linekind=linekind,
;      /endpoints,/oneway,/noreverse,/trim,/quiet
; 
;  INPUTS: 
;    sph_data = a structure of type spherical_field_data (see
;      spherical_field_data__define.pro) with the following fields defined on
;      input: br,bth,bph,nlat,nlon,nr,rix,thix,phix,lat,lon,str,stth,stph.
;      Basically, one needs the vector field (br,bth,bph), its dimension
;      (nr,nlat,nlon), its indexing (rix,thix,phix,lat,lon), and the starting
;      points (str,stth,stph).  
;    stepmax = max number of steps per field line (default=3000)
;    safety = maximum ds along each field line, in units of minimum grid
;             spacing (default = 0.5)
;    endpoints = set if endpoints only are desired
;    oneway = set if loops are to be traced in one direction only, either
;             upward or downward depending on whether starting point is closer
;             to the bottom or to the top of the domain, and then if the
;             starting point is no more than 1% away from the upper or lower
;             boundary
;    noreverse = set if field lines are not to be reversed (by default, closed
;                field lines are oriented so that the end with negative
;                polarity comes first and open field lines are oriented so
;                that the photospheric end comes first
;    trim = trim fieldlines so that the trajectories do not extend beyond the
;           boundaries of the domain
;    quiet = set this keyword to inhibit screen output
; 
;  OUTPUTS:
;    sph_data = a structure of type spherical_field_data (see
;      spherical_field_data__define.pro) with the following fields set on
;      output: ptr,ptth,ptph,nstep.  These are the (r,theta,phi)-trajectories
;      of each fieldlines (ptr,ptth,ptph) and the number of points in each
;      line (nstep).
;    outfield = an array of dimension [nlines,3,3], containing the three
;      components of the vector field at both endpoints and the starting point
;      for each line.  The first dimension indexes the fieldlines, the second
;      dimension indexes the points of interest (endpoint #1, starting point,
;      endpoint #2), and the third dimension indexes the coordinate
;      (r,theta,phi).
;    linelengths = an array of dimension [nlines] containing the lengths of
;                  the fieldlines that were traced (in whatever units the
;                  radius variable is in)
;    linekind = an nlines-element array containing an indeger coded to
;      represent the type of fieldline, as follows:
;        -1=line starting point out of bounds
;         0=error of some sort?
;         1=maximum step limit reached
;         2=line intersects inner and outer radial boundaries, 
;         3=both endpoints of line lie on inner boundary, 
;         4=both endpoints of line lie on outer boundary,
;         5=line intersects inner and side boundaries,
;         6=line intersects outer and side boundaries,
;         7=both endpoints of line lie on side boundary/ies
;         8=one of the endpoints is a null
;
;  NOTES:  -The latitude array (sph_data.lat) MUST be in a monotonically
;           ascending order for this procedure to work properly.
;
;  MODIFICATION HISTORY:
;    M.DeRosa - 13 Dec 2005 - copied the guts of pfss_trace_field
;               24 Jan 2006 - now accommodates data that is bounded in phi
;               27 Jan 2006 - now accommodates starting points that are out of
;                             bounds, as well as fieldlines that hit nulls
;               20 Sep 2006 - added quiet keyword
;               18 Sep 2008 - added a ptr_free command prior to, and added the
;                             /no_copy flag during, the new pointer assignment
;                             statements to prevent memory leaks
;               16 Dec 2009 - added endpoints,oneway,noreverse keywords
;               18 Dec 2009 - added trim and linelengths keywords
;               27 Dec 2009 - fixed two issues with linelength calculation
;                             (put inside loop, error in cartesian conversion)
;               20 Jan 2010 - fixed concatenation error occurring when the
;                             starting point arrays are not one-dimensional
;               15 Mar 2010 - fixed issue with logic when /trim keyword is set
;               15 Nov 2010 - fixed minor issue for when nlines=1
;-

pro spherical_trace_field,sph_data,stepmax=stepmax,safety=safety, $
  outfield=outfield,linelengths=linelengths,linekind=linekind, $
  endpoints=endpoints,oneway=oneway,noreverse=noreverse,trim=trim,quiet=quiet

;  usage message
if n_elements(sph_data) eq 0 then begin
  print,'  ERROR in spherical_trace_field: no input data provided'
  print,'  Calling sequence: spherical_trace_field,sph_data,'+ $
    'stepmax=stepmax,safety=safety,outfield=outfield,linelengths=linelengths,'+$
    'linekind=linekind,/oneway,/endpoints,/noreverse,/trim,/quiet'
  return
endif

;  parameters
if n_elements(stepmax) eq 0 then stepmax=15000l else stepmax=long(stepmax(0))
if n_elements(safety) eq 0 then safety=0.1 else safety=float(safety(0))
nlines=n_elements(*sph_data.str)  ;  nlines = number of field lines to be drawn
rmin=min(*sph_data.rix,max=rmax)
thmin=min(*sph_data.theta,max=thmax)
if n_elements(sph_data.lonbounds) gt 0 then begin
  if sph_data.lonbounds(0) ge 0 then begin
    ph1=sph_data.lonbounds(0)*!dpi/180
    ph2=sph_data.lonbounds(1)*!dpi/180
    bounded=1b
  endif else bounded=0b
endif else bounded=0b

;  do some boundary checking
if bounded then begin
  if ph2 gt ph1 then begin
    inbounds=(*sph_data.stph ge ph1) and (*sph_data.stph le ph2)
  endif else begin
    inbounds=(*sph_data.stph gt ph2) or (*sph_data.stph lt ph1)
  endelse
endif else inbounds=replicate(1b,nlines)  ;  all 1's 
inbounds=inbounds and (*sph_data.str ge rmin) and (*sph_data.str le rmax)
inbounds=inbounds and (*sph_data.stth ge thmin) and (*sph_data.stth le thmax)

;  get deltas
deltar=(*sph_data.rix)(1)-(*sph_data.rix)(0)
deltath=(*sph_data.theta)(0)-(*sph_data.theta)(1)
deltaph=((*sph_data.lon)(1)-(*sph_data.lon)(0))*!dtor

;  initialize arrays
nstep=lonarr(nlines)
linekind=intarr(nlines)
linelengths=fltarr(nlines)
stbph=(stbth=(stbr=fltarr(nlines)))
if keyword_set(endpoints) then begin
  ptr=fltarr(2,nlines)
  ptth=fltarr(2,nlines)
  ptph=fltarr(2,nlines)
endif else begin
  ptr=fltarr(stepmax,nlines)
  ptth=fltarr(stepmax,nlines)
  ptph=fltarr(stepmax,nlines)
endelse
ptbr=fltarr(nlines,2)
ptbth=fltarr(nlines,2)
ptbph=fltarr(nlines,2)
ptr(0,*)=*sph_data.str
ptth(0,*)=*sph_data.stth
ptph(0,*)=*sph_data.stph
ir=fltarr(stepmax)
ith=fltarr(stepmax)
iph=fltarr(stepmax)

;  loop through each line
if not keyword_set(quiet) then print,$
  '  spherical_trace_field: tracing '+strcompress(nlines,/r)+' field lines'
for i=0l,nlines-1 do begin

  ;  print time left
  if not keyword_set(quiet) then $
    pfss_print_time,'  spherical_trace_field: ',i+1,nlines,tst,slen1 $
    else tst=long(systime(1))

  ;  initialize
  ir(0)=ptr(0,i)
  ith(0)=ptth(0,i)
  iph(0)=ptph(0,i)
  step=1l
;irc=get_interpolation_index(*sph_data.rix,ir(0))
 ;     ithc=get_interpolation_index(*sph_data.lat,90-ith(0)*!radeg)
  ;    iphc=get_interpolation_index(*sph_data.lon,(iph(0)*!radeg+360) mod 360)
;print,[irc,ithc,iphc]
  ;  trace out line
  if inbounds(i) then begin
    repeat begin

      ;  current point
      ptc=[ir(step-1),ith(step-1),iph(step-1)]  ;  (r,th,ph) of current point
      spth=sin(ptc(1))

      ;  calculate value of Br,Bth,Bph at current point
      irc=get_interpolation_index(*sph_data.rix,ptc(0))
      ithc=get_interpolation_index(*sph_data.lat,90-ptc(1)*!radeg)
      iphc=get_interpolation_index(*sph_data.lon,(ptc(2)*!radeg+360) mod 360)
      brc=interpolate(*sph_data.br,iphc,ithc,irc)
      bthc=interpolate(*sph_data.bth,iphc,ithc,irc)/ptc(0)
      bphc=interpolate(*sph_data.bph,iphc,ithc,irc)/(ptc(0)*spth)
      ds=reform([brc,bthc,bphc])

      ;  initialization
      if step eq 1 then begin

        ;  set initial fields
        ptbr(i,0)=brc
        ptbth(i,0)=bthc
        ptbph(i,0)=bphc

        ;  determine of loop is to be traced both ways or one way only
        if keyword_set(oneway) then begin

          ;  determine if starting point is at top, bottom or middle
          pct=(ptc(0)-rmin)/(rmax-rmin)
          case 1 of
            pct gt 0.99: begin
              top=1
              linekind(i)=4
              end
            ((pct lt 0.01) and (float((*sph_data.rix)(1)) ge 1.01)): begin
              top=-1
              linekind(i)=3
              end
            ((ptc(0) lt float((*sph_data.rix)(1))) or (pct lt 0.01)): begin
              top=-1
              linekind(i)=3
              end
            else: top=0
          endcase

        endif else top=0

        ;  set initial stepsize and direction for first step
        stbr(i)=brc  &  stbth(i)=bthc  &  stbph(i)=bphc
        bsign=sign_mld(brc)
        steplen=safety*deltar/(abs(brc)>1.0)
        if top ge 0 then bsign=-bsign

      endif
      if bsign lt 0 then ds=-ds

      ;  take a step, repeat if too big
      stepflag=0
      repeat begin

        ;  step forward by an amount steplen along direction ds
        result=forw_euler(step-1,ptc,steplen,ds)

        ;  evaluate if current step was too big or too small
        diff=result-ptc
        err=max(abs(diff/[deltar,deltath,deltaph]))/safety
        case 1 of
          err gt 1: begin  ;  stepsize too big, reduce and do over
            steplen=0.5*steplen
            stepflag=0
            end
          err lt 0.1: begin  ;  stepsize too small, increase for next time
            steplen=2*steplen
            stepflag=1
            end
          else: stepflag=1  ;  don't change stepsize, just set flag to exit
        endcase
        
      endrep until stepflag

      ;  store in ir,ith,iph arrays
      ir(step)=result(0)
      ith(step)=result(1)
      iph(step)=result(2)

      ;  set flags
      hitnull=total(diff^2) eq 0.0
      hitsides=(ith(step) lt thmin) or (ith(step) gt thmax)
      if bounded then begin
        if ph2 gt ph1 then begin
          hitsides=hitsides or ((iph(step) gt ph2) or (iph(step) lt ph1))
        endif else begin
          hitsides=hitsides or ((iph(step) gt ph1) and (iph(step) lt ph2))
        endelse
      endif
      hittop=(ir(step) gt rmax)
      hitbot=(ir(step) lt rmin)
      hitstepmax=((step+1) eq stepmax)

      ;  deal with flags
      if (hitsides or hittop or hitbot or hitstepmax or hitnull) then begin
        if top eq 0 then begin
          if hitstepmax then begin  ;  line may close on itself?
            linekind(i)=1
            flag=1
          endif else begin  ;  hit boundary

            ;  reset flag
            flag=0

            ;  set case to indicate which boundary
            case 1 of
              hittop: begin
                linekind(i)=4
                top=1
                end
              hitbot: begin
                linekind(i)=3
                top=-1
                end
              hitsides: if linekind(i) eq 7 then flag=1 else linekind(i)=7
              hitnull: if linekind(i) eq 8 then flag=1 else linekind(i)=8
              else: begin
                print,'  shouldn''t be able to get here'
                stop
                end
            endcase

            ;  reverse order of data and continue in other direction
            ir(0:step)=reverse(ir(0:step))
            ith(0:step)=reverse(ith(0:step))
            iph(0:step)=reverse(iph(0:step))
          
            ;  reset initial fields
            ptbr(i,0)=brc
            ptbth(i,0)=bthc
            ptbph(i,0)=bphc

            ;  reset sign, step length
            bsign=-bsign

          endelse
        endif else begin
          flag=1
          kindc=linekind(i)
          case 1 of
            hittop: begin
              case kindc of
                3: linekind(i)=2  ;  I'm not sure this can happen
                4: linekind(i)=4
                7: linekind(i)=6
                8: linekind(i)=8
                else: begin
                  print,'  shouldn''t be able to get here'
                  stop
                  end
              endcase
              end
            hitbot: begin
              case kindc of
                3: linekind(i)=3  ;  I'm not sure this can happen either 
                4: begin
                  linekind(i)=2
                  if not keyword_set(noreverse) then begin
                    ir(0:step)=reverse(ir(0:step))
                    ith(0:step)=reverse(ith(0:step))
                    iph(0:step)=reverse(iph(0:step))
                  endif
                  end
                7: linekind(i)=5
                8: linekind(i)=8
                else: begin
                  print,'  shouldn''t be able to get here'
                  stop
                  end
              endcase
              end
            hitsides: begin
              case kindc of
                3: linekind(i)=5
                4: linekind(i)=6
                8: linekind(i)=8
                else: begin
                  print,'  shouldn''t be able to get here'
                  stop
                  end
              endcase
              end
            hitnull: linekind(i)=8
            hitstepmax: linekind(i)=1
          endcase
        endelse
      endif else flag=0

      ;  increment counter
      step=step+1

    endrep until flag

    ;  set final fields
    ptbr(i,1)=brc
    ptbth(i,1)=bthc
    ptbph(i,1)=bphc

    ;  if closed, make field line go from negative to positive
    if not keyword_set(noreverse) then begin
      if (linekind(i) eq 3) and (ptbr(i,0) gt 0) then begin
        ir(0:step-1)=reverse(ir(0:step-1))
        ith(0:step-1)=reverse(ith(0:step-1))
        iph(0:step-1)=reverse(iph(0:step-1))
        ptbr(i,*)=reverse(ptbr(i,*))
        ptbth(i,*)=reverse(ptbth(i,*))
        ptbph(i,*)=reverse(ptbph(i,*))
      endif
    endif

    ;  trim lines if desired
    if keyword_set(trim) then begin

      ;  loop through all of the different types of boundaries
      for k=0,5 do begin

        ;  set field component variable and boundary values
        case k of
          0: begin  ;  fieldlines that cross inner radial boundary
            fcomp=ir(0:step-1)
            bval=rmin
            end
          1: begin  ;  fieldlines that cross outer radial boundary
            fcomp=ir(0:step-1)
            bval=rmax
            end
          2: begin  ;  fieldlines that cross lower latitudinal boundary
            fcomp=ith(0:step-1)
            bval=thmin
            end
          3: begin  ;  fieldlines that cross upper latitudinal boundary
            fcomp=ith(0:step-1)
            bval=thmax
            end
          4: begin  ;  fieldlines that cross lower longitudinal boundary
            fcomp=iph(0:step-1)
            if bounded then begin
              if ph2 gt ph1 then bval=ph1 else bval=ph2
            endif else bval=min(fcomp)-1  ;  rigged so that ncr=0 below
            end
          5: begin  ;  fieldlines that cross upper longitudinal boundary
            fcomp=ir(0:step-1)
            if bounded then begin
              if ph2 gt ph1 then bval=ph2 else bval=ph1
            endif else bval=max(fcomp)+1  ;  rigged so that ncr=0 below
            end
        endcase

        ;  see how many crossings of the boundary value there are
        crosstest=(fcomp(0:step-2)-bval)*(fcomp(1:step-1)-bval)
        whcross=where(crosstest lt 0,ncr)

        ;  interpolate between the two points on either side of the crossing
        if ncr gt 0 then begin
          for j=0l,ncr-1 do begin
            case whcross(j) of
              0: begin  ;  first point
                coeff=(bval-fcomp(0))/(fcomp(1)-fcomp(0))
                ir(0)=ir(0)+coeff*(ir(1)-ir(0))
                ith(0)=ith(0)+coeff*(ith(1)-ith(0))
                iph(0)=iph(0)+coeff*(iph(1)-iph(0))
                end
              step-2: begin  ;  last point
                coeff=(bval-fcomp(step-2))/(fcomp(step-1)-fcomp(step-2))
                ir(step-1)=ir(step-2)+coeff*(ir(step-1)-ir(step-2))
                ith(step-1)=ith(step-2)+coeff*(ith(step-1)-ith(step-2))
                iph(step-1)=iph(step-2)+coeff*(iph(step-1)-iph(step-2))
                end
              step-3: begin  ;  second to last point (sometimes occurs)
                step=step-1
                coeff=(bval-fcomp(step-2))/(fcomp(step-1)-fcomp(step-2))
                ir(step-1)=ir(step-2)+coeff*(ir(step-1)-ir(step-2))
                ith(step-1)=ith(step-2)+coeff*(ith(step-1)-ith(step-2))
                iph(step-1)=iph(step-2)+coeff*(iph(step-1)-iph(step-2))
                end
              else: begin  ;  likely bug if stopped here
                print,'  ERROR in spherical_trace_field: fieldline crosses '+$
                  'boundary value more than twice (for some reason)'
                stop,'  (type .c to return)'
                return
                end
            endcase
          endfor
        endif

        ;  check to see if there are endpoints right on the upper/outer bound
        if ((crosstest(step-3) eq 0) and (crosstest(step-2) eq 0)) then begin
          step=step-1  ;  (last point is eliminated)
        endif

        ;  check to see if there are endpoints right on the lower/inner bound
        if ((crosstest(0) eq 0) and (crosstest(1) eq 0)) then begin
          ir=shift(ir,-1)
          ith=shift(ith,-1)
          iph=shift(iph,-1)
          step=step-1  ;  (first point is eliminated)
        endif

      endfor

    endif

    ;  determine line length
    xpt=ir*sin(ith)*sin(iph)
    ypt=ir*sin(ith)*cos(iph)
    zpt=ir*cos(ith)
    ptdist=sqrt((shift(xpt,-1)-xpt)^2+(shift(ypt,-1)-ypt)^2+ $
      (shift(zpt,-1)-zpt)^2)
    linelengths(i)=total(ptdist(0:step-2))

    ;  fill pt arrays
    nstep(i)=step
    if keyword_set(endpoints) then begin
      ptr(*,i)=ir([0,step-1])
      ptth(*,i)=ith([0,step-1])
      ptph(*,i)=iph([0,step-1])
    endif else begin
      ptr(*,i)=ir
      ptth(*,i)=ith
      ptph(*,i)=iph
    endelse

  endif else begin  ;  line starting point is out of bounds 
    nstep(i)=-1
    linekind(i)=-1
  endelse

endfor

;  exit if all starting points out of bounds
if round(total(inbounds)) eq 0 then begin
  print,'  ERROR in spherical_trace_field: starting points out of bounds'
  return
endif

;  print total time
ttime=round((long(systime(1))-tst)/60.)
if not keyword_set(quiet) then $
  print,'  spherical_trace_field: total time = '+strcompress(ttime,/r)+' min'

;  save a little memory hopefully
maxnstep=max(nstep)
if (not keyword_set(endpoints)) and (max(nstep) lt stepmax) then begin
  ptr=ptr(0:maxnstep-1,*)
  ptth=ptth(0:maxnstep-1,*)
  ptph=ptph(0:maxnstep-1,*)
endif

;  get interpolates for outfield
naxptr=size(ptr,/dim)
nwidth=naxptr(0)
;nlines=naxptr(1)  ;  already defined, but causes problems for when nlines=1
outfield=fltarr(nlines,3,3,/noz)
endix=(nstep-1)+l64indgen(nlines)*nwidth
points=[[reform(ptph(0,*)),reform(*sph_data.stph,nlines),ptph(endix)], $
        [reform(ptth(0,*)),reform(*sph_data.stth,nlines),ptth(endix)], $
        [reform(ptr(0,*)),reform(*sph_data.str,nlines),ptr(endix)]]

;  now interpolate to get the field
iphc=get_interpolation_index(*sph_data.lon,points(*,0)*180./!dpi)
ithc=get_interpolation_index(*sph_data.lat,90-points(*,1)*180./!dpi)
irc=get_interpolation_index(*sph_data.rix,points(*,2))
outfield(*,*,0)=interpolate(*sph_data.br,iphc,ithc,irc)
outfield(*,*,1)=interpolate(*sph_data.bth,iphc,ithc,irc)
outfield(*,*,2)=interpolate(*sph_data.bph,iphc,ithc,irc)

;  free some memory (if needed) and assign pointers
if ptr_valid(sph_data.nstep) then ptr_free,sph_data.nstep
if ptr_valid(sph_data.ptr) then ptr_free,sph_data.ptr
if ptr_valid(sph_data.ptth) then ptr_free,sph_data.ptth
if ptr_valid(sph_data.ptph) then ptr_free,sph_data.ptph
sph_data.nstep=ptr_new(nstep,/no_copy)
sph_data.ptr=ptr_new(ptr,/no_copy)
sph_data.ptth=ptr_new(ptth,/no_copy)
sph_data.ptph=ptr_new(ptph,/no_copy)

end

