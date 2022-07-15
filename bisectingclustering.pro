;----------------------------------------------------------------------
; This is an implementation of a Divisive Bisecting Clustering
; algorithm.
;
; This is the IDL implementation by Dr. Christiaan Boersma
; (Christiaan.Boersma@nasa.gov).
;
; Please refer to the following paper
;
; Boersma, C., Bregman, J., Allamandola, L.J., "The Charge State of
; Polycyclic Aromatic Hydrocarbons Across Reflection Nebulae: PAH
; Charge Balance and Calibration", 2016, ApJ, 832, 51
; https://doi.org/10.3847/0004-637X/832/1/51
;
;----------------------------------------------------------------------
; Permission to use, copy, or modify this software and its documentation
; for educational and research purposes only and without fee is hereby
; granted, provided that this copyright notice and the original authors'
; names appear on all copies and supporting documentation. This program
; shall not be used, rewritten, or adapted as the basis of a commercial
; software or hardware product without first obtaining permission of the
; authors. The authors make no representations about the suitability of
; this software for any purpose. It is provided "as is" without express
; or implied warranty.
;----------------------------------------------------------------------
;
;Input : (1) data: the data to cluster (samples x features).
;        (2) bisection_max: maximum number of bisections. The expected
;                           number of bins is 2 x bisection_max.
;        (3) min_samples: minimum number of samples needed to be
;            considered a bin.
;
;Output: (1) bin_labels: the label for each sample it belongs to
;                        (samples).
;
;Basic Usage:
;   Given 2 test data formatted as samples x features
;
;   bin_labels = BISECTINGCLUSTERING(data, 2, 10)
;
;Advanced Usage:
;   A spectral cube (n, m, w). For example
;
;   bin_labels = REFORM( $
;                 BISECTINGCLUSTERING( $
;                  REFORM(cube, n * m, w), 2, 10), n, m)
;========================================================================

PRO NODES,matrix,depth,max_depth,min_samples,bin_labels,parent

  COMPILE_OPT IDL2

  ON_ERROR, 2

  MESSAGE,STRING(FORMAT='("DEPTH:",X,I0)', depth),/INFORMATIONAL

  IF depth GT 0 THEN d = depth - 1 $
  ELSE d = depth

  k = WHERE(bin_labels[*,d] EQ parent, nbin)

  IF nbin LT min_samples THEN BEGIN

     MESSAGE,STRING(FORMAT='(A0,X,"(NBIN=",I0,")")', "SINGULAR DECOMPOSITION OR SPARSE BIN: DONE", nbin),/INFORMATIONAL

     bin_labels[k,depth:*] = parent

     RETURN

  ENDIF

  dummy = MIN(matrix, imin, /NAN, SUBSCRIPT_MAX=imax)

  bin1 = parent * 2 + 1 & bin2 = bin1 + 1

  ij = ARRAY_INDICES(matrix, imax)

  srt = SORT(ij)

  ij = ij[srt]

  MESSAGE,STRING(FORMAT='("bin_labels:",X,I0,X,"&",X,I0)',bin1,bin2),/INFORMATIONAL

  bin_labels[ij[0],depth] = bin1 & bin_labels[ij[1],depth] = bin2

  m1 = matrix & m2 = matrix

  m1[ij[1],*] = !VALUES.D_NAN & m1[*,ij[1]] = !VALUES.D_NAN

  m2[ij[0],*] = !VALUES.D_NAN & m2[*,ij[0]] = !VALUES.D_NAN

  FOR j = 0, nbin - 1 DO BEGIN

     IF k[j] EQ ij[0] OR k[j] EQ ij[1] THEN CONTINUE

     IF NOT FINITE(matrix[ij[0], k[j]]) OR NOT FINITE(matrix[ij[1], k[j]]) THEN CONTINUE

    b = matrix[ij[0], k[j]] ; bin1

    a = matrix[ij[1], k[j]] ; bin2

    ; if b<a then add k[j] to bin1 - i.e., a>b

     IF a GT b THEN BEGIN

        bin_labels[k[j],depth] = bin1

        m2[k[j],*] = !VALUES.D_NAN

        m2[*,k[j]] = !VALUES.D_NAN

     ENDIF ELSE BEGIN

        bin_labels[k[j],depth] = bin2

        m1[k[j],*] = !VALUES.D_NAN

        m1[*,k[j]] = !VALUES.D_NAN

     ENDELSE

  ENDFOR

  depth++

  IF depth EQ max_depth THEN BEGIN

     MESSAGE,'HIT MAX DEPTH: DONE',/INFORMATIONAL

     RETURN

  ENDIF

  d = depth

  NODES,m1,d,max_depth,min_samples,bin_labels,bin1

  NODES,m2,depth,max_depth,min_samples,bin_labels,bin2

END

FUNCTION BISECTINGCLUSTERING,data,max_depth,min_samples

  COMPILE_OPT IDL2

  ON_ERROR, 2

  IF N_PARAMS() NE 3 THEN $
     MESSAGE,"EXPECTING THREE PARAMETERS"

  IF SIZE(data, /N_DIMENSIONS) NE 2 THEN $
     MESSAGE,"EXPECTING TWO DIMENSIONS"

  ndata = N_ELEMENTS(data[*,0])

  matrix = DBLARR(ndata, ndata)

  FOR i = 0, ndata - 1 DO BEGIN

    FOR j = 0, i DO BEGIN

      IF i EQ j THEN BEGIN

         matrix[i,j] = !VALUES.D_NAN

         CONTINUE

      ENDIF

      matrix[i,j] = NORM(data[i,*]-data[j,*], LNORM=2)

      matrix[j,i] = matrix[i,j]

    ENDFOR

  ENDFOR

  m = matrix

  depth = 0 & parent = 0

  bin_labels = INTARR(ndata, max_depth)

  NODES,m,depth,max_depth,min_samples,bin_labels,parent

  FOR i = 0, max_depth - 1 DO BEGIN

     slice = bin_labels[*,i]

     u = slice[UNIQ(slice, SORT(slice))]

     nu = N_ELEMENTS(u)

     FOR j = 0, nu - 1 DO BEGIN

        sel = WHERE(bin_labels[*,i] EQ u[j])

        slice[sel] = u[j]

     ENDFOR

  ENDFOR

  u = UNIQ(slice, SORT(slice))

  nu = N_ELEMENTS(u)

  bin_labels = slice * 0

  srt = DBLARR(nu)

  offset = 1

  IF MIN(slice[u]) EQ 0 THEN offset = 0

  FOR i = 0, nu - 1 DO BEGIN

     sel = WHERE(slice EQ slice[u[i]])

     srt[i] = MIN(sel)

     IF slice[u[i]] EQ 0 THEN srt[i] = -1D

     bin_labels[sel] = i + offset

  ENDFOR

  srt = SORT(srt)

  srt = SORT(srt) ; yes, a second time!

  bin_labels = srt[bin_labels - offset] + offset

  RETURN, bin_labels

END
