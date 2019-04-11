;this is plr.pro, version 19 Oct 1998 by RHM. This program plots 
;the mags V and I as a function of log P for the cepheids analyzed
;with lafkin.pro, and then helps to determine the relative distance
;modulus (relative to LMC) using the cepheid data in file 'plr.dat'.
;the user is asked to estimate the distance modulus as many times as
;needed until he is satisfied with the fit.

logP=fltarr(50) &    logPLMC=logP
magv=logP       &    magvLMC=logP
magi=logP       &    coviLMC=logP
tit1=''
tit2=''
fname=''
on_ioerror, goahead
read, 'enter filename, e.g. n1365.dat, do not use quotes: ', fname
openr, 1, fname
readf, 1, tit1
readf, 1, tit2
jj=0
  while not eof(1) do begin
  readf, 1, cant0, cant1, cant2 
  logP(jj) = cant0         ; periods in days
  magv(jj) = cant1         ; mags V
  magi(jj) = cant2         ; mags I
  jj=jj+1 
  endwhile
goahead:close, 1

logP=alog10(logP(0:jj-1))  ;here we calculate the logs of periods
magv=magv(0:jj-1)
magi=magi(0:jj-1)

!x.style=1
!y.style=1
!x.thick=2
!y.thick=2
!p.charsize=2
!p.thick=2
!x.range=[0.6,1.9]
!x.title='!17log of period (days)'
window, 1

openr, 1, 'plr.dat'
readf, 1, tit1
readf, 1, tit2
jj=0
  while not eof(1) do begin
  readf, 1, cant0, cant1, cant2 
  logPLMC(jj) = cant0         ; logs of periods in days
  magvLMC(jj) = cant1         ; mags V
  coviLMC(jj) = cant2         ; colors V-I
  jj=jj+1 
  endwhile
close, 1

logPLMC=logPLMC(0:jj-1)
magvLMC=magvLMC(0:jj-1)
coviLMC=coviLMC(0:jj-1)

magiLMC=magvLMC-coviLMC

;now determine (interactively) the distance modulus differences V, I
read, 'enter first estimate of delta(m-M) in V and I: ', dmagV, dmagI

while (dmagV ne 0) do begin
magvL=magvLMC+dmagV
magiL=magiLMC+dmagI
wset,1
!y.range=[29.,24.]
!y.title='!17V magnitudes'
plot, logP, magv, psym=2
oplot, logPLMC, magvL, psym=6
set_plot, 'ps'
device, /landscape, filename='finv.ps'
plot, logP, magv, psym=2
oplot, logPLMC, magvL, psym=6
endpl
wset,0
!y.range=[28.,23.]
!y.title='!17I magnitudes'
plot, logP, magi, psym=2
oplot, logPLMC, magiL, psym=6
set_plot, 'ps'
device, /landscape, filename='fini.ps'
plot, logP, magi, psym=2
oplot, logPLMC, magiL, psym=6
endpl
read, 'enter delta(m-M) in V and I, zero to finish: ', dmagV, dmagI
endwhile

end
