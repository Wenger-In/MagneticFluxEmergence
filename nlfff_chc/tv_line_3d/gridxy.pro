;+
; NAME:
;	gridxy
; PURPOSE:
;	Create a 2-D grid of X & Y coordinates, suitable for computations,
;	in the form of 2 images, stacked as 3-D array,
;	or 2 vectors (collapsed images, 2-D array result).
; CALLING:
;	gxy = gridxy( Ngx, Ngy, XRAN=[-1,1], YRAN=[-1,1] )
; INPUTS:
;	Ngx, Ngy = number of X & Y grid points, resp.
; KEYWORDS:
;	XRANGE, YRANGE = 2 element arrays specifying [min,max] of X & Y range.
;		Note that the ranges are saved in common MapLim, so that if
;		next call does not specify them, the common values are used.
;		Default values are [0,Ngx] and [0,Ngy].
;	DIMENSIONS = 2 elements array, alternative specification of [Ngx,Ngy].
;	/INCLUDE_MAXRAN causes the last coordinate in output array to be
;			the maximum of X & Y range (default is max - grid_size).
;	SHIFT = fraction of grid size to shift points so as to stay within
;			specified range, default is start at minimum of range.
;
;	/VECTOR : collapse 2-D to 1-D and return 2-D array : ( Ngx*Ngy, 2 ).
;	/TRANSPOSE applies only when /VECTOR, result is fltarr( 2, Ngx*Ngy )
;	/FLIP transposes the x and y coordinates.
;	POWER = scalar, optional exponent used to raise elements to a power
;		before creating matrix, thus faster than afterwards.
;		Power is not applied to zero values, they are kept zero.
; OUTPUTS:
;	Result is an array of X & Y coordinate points:
;	3-D array = fltarr( Ngx, Ngy, 2 )
;	or 2-D array = fltarr( Ngx*Ngy, 2 )  if /VECTOR set.
; COMMON BLOCKS:
;	common MapLim ,xminc,xmaxc ,yminc,ymaxc		;can specify limits.
; PROCEDURE:
;	Use the matrix outer product (#).
; HISTORY:
;	Written: Frank Varosi U.MD. 1987
;	F.V. 1990, updated for IDL version-2, added keywords.
;-

function gridxy, Ngx, Ngy, DIMENSIONS=dimen, XRANGE=xran, YRANGE=yran,  $
			INCLUDE_MAXRAN=include_maxran, SHIFT_FRACTION=shift, $
			FLIP=flip, TRANSPOSE=transp, VECTOR=vector, POWER=power
							

  common MapLim ,xminc,xmaxc ,yminc,ymaxc

	if N_elements( xran ) EQ 2 then begin
		xminc = xran(0)
		xmaxc = xran(1)
	  endif else if (N_elements( xminc ) EQ 1) AND $
			(N_elements( xmaxc ) EQ 1) then xran = [xminc,xmaxc]

	if N_elements( yran ) EQ 2 then begin
		yminc = yran(0)
		ymaxc = yran(1)
	  endif else if (N_elements( yminc ) EQ 1) AND $
			(N_elements( ymaxc ) EQ 1) then yran = [yminc,ymaxc]

	if N_elements( dimen ) EQ 2 then begin
		Ngx = dimen(0)
		Ngy = dimen(1)
	   endif

	if N_elements( xran ) NE 2 then xran = [0,Ngx]
	if N_elements( yran ) NE 2 then yran = [0,Ngy]
	xgs = float( xran(1) - xran(0) )/Ngx
	ygs = float( yran(1) - yran(0) )/Ngy

	if keyword_set( include_maxran ) then begin
		xgs = xgs * float( Ngx )/( Ngx - 1 )
		ygs = ygs * float( Ngy )/( Ngy - 1 )
	   endif

	Gx = xgs * findgen( Ngx ) + xran(0)
	Gy = ygs * findgen( Ngy ) + yran(0)

	if keyword_set( shift ) then begin
		Gx = Gx + shift * xgs
		Gy = Gy + shift * ygs
	   endif

	if keyword_set( power ) then begin
		w = where( Gx NE 0, nw )
		if (nw GT 0) then Gx(w) = Gx(w)^power
		w = where( Gy NE 0, nw )
		if (nw GT 0) then Gy(w) = Gy(w)^power
	   endif

	if keyword_set( flip ) then begin

		Gxm = replicate( 1,Ngy ) # Gx 	;generate X matrix by columns
		Gym = Gy # replicate( 1,Ngx ) 	;generate matrix of Y by rows

	  endif else begin

		Gxm = Gx # replicate( 1,Ngy )	;generate matrix of X by rows
		Gym = replicate( 1,Ngx ) # Gy	;generate Y matrix by columns
	   endelse

	if keyword_set( vector ) then begin	;collapse matrices into vectors
  						;concatenate to form 2-D array
		Ngxy = Ngx * Long( Ngy )
		grid = [ [reform( Gxm, Ngxy )] , [reform( Gym, Ngxy )] ]

		if keyword_set( transp ) then return, transpose( grid ) $
					 else return, grid

	  endif else   return, [ [[Gxm]] , [[Gym]] ]	;default is 3-D array
end
