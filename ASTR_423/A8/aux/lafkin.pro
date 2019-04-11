; this is lafkin.pro, version 16 Oct 1998 by RHM. Period search
; using the method by Lafler & Kinman 1965 ApJS 11, 216.

tit1=''
tit2=''
tit3=''
fnam=''
t=fltarr(30)  ;maximum number of data (verify before using file)
v=t
; select and read data
read, 'enter filename, e.g. cepv3.dat, do not use quotes: ', fnam
openr, 1, fnam
readf, 1, tit1
readf, 1, tit2
readf, 1, tit3
jj=0
  while not eof(1) do begin
  readf, 1, cant0, cant1 
  t(jj) = cant0         ; time of measurement
  v(jj) = cant1         ; measurement
  jj=jj+1 
  endwhile
close, 1
;print, jj
t=t(0:jj-1)
v=v(0:jj-1)
f=v   &   pha=v   &   fap=v   &   vv=v

IP=3.5               ; the initial period in days
FP=57.               ; the final period in days
DP=0.02                ; the increment

NP=(FP-IP)/DP      ; number of periods to be tested
th=fltarr(NP+1)
P=th        &         the=th
;now calculate sigma denominator
avv=total(v)/jj
sqs=(v-avv)^2
sig=total(sqs)
;now we need sigma numerator for each trial period
  for k=0, NP do begin
  P(k) = IP + k*DP
  pinv=1./P(k)
  f = (t-t(0))*pinv          ;the phases
  f = f - fix(f)             ;are calculated
  vv = v(sort(f))            ;and sorted
     for i=0, jj-2 do begin
     th(k) = th(k) + (vv(i) - vv(i+1))^2
     end
  th(k) = th(k) + (vv(0) - vv(jj-1))^2  ;numerator calculated
  end
the = th/sig    ;numerator/denominator for all tested periods

!x.style=1
!y.style=1
!x.thick=2
!y.thick=2
!p.charsize=2
!p.thick=2
!x.range=[1.,FP]
!y.range=[0.,2.5]
!x.title='!17period (days)'
!y.title='!17function theta'

wset, 0
plot, P, the, title=fnam
set_plot, 'ps'
device, /landscape, filename='pese.ps'
plot, P, the, title=fnam
endpl
indi=min(the)
pp=P(where(the EQ indi))
print, pp, '   is the best period.'

;print, where(the EQ indi)
;now we want to check the shape of the light (or velocity) curve,
;for each of several different periods.
window, 1      ;make another window
wset, 1        ;and activate it
psel=pp(0)  ;pp could be a vector in case of several best periods

;repeat begin   ;start loop testing different periods
while (psel ne 0) do begin

piv = 1./psel
    for kk=0,jj-1 do begin
    fap(kk) = (t(kk)-t(0)) * piv
    pha(kk) = fap(kk) - fix(fap(kk))      ;phases are calculated
    end
;print,pha

!x.range=[0.,2.]
!y.range=[29.,24.]
!x.title='!17phase'
!y.title='!17V'

pha2=1.+pha
plot, pha, v, psym=4, title=string(psel)
oplot,pha2,v, psym=4
set_plot, 'ps'
device, /landscape, filename='phas.ps'
plot, pha, v, psym=4, title=string(psel)
oplot,pha2,v, psym=4
endpl

read, 'enter new period, 0 to finish: ', psel
endwhile
;endrep until (psel eq 0)  ;end of loop testing different periods

end
