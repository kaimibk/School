;this is zeropo.pro, version 23 Oct 1998 by RHM. Here we enter our
;abs visual mag of S Normae and compare the other Galactic cepheids
;with the period-mag relation for the LMC.

Pgal=fltarr(19)   &   Mvga=Pgal
tit1=''
tit2=''
openr, 1, 'cepgal.dat'   ;read data for Galactic cepheids
readf, 1, tit1
readf, 1, tit2
for jj=0,18 do begin
  readf, 1, cant0, cant1
  Pgal(jj) = cant0       ; Periods (days) of Gal cephs
  Mvga(jj) = cant1       ; Absolute visual mags of Gal cephs
endfor
close,1
logPLMC=fltarr(26)  &  magvLMC=logPLMC   &  coviLMC=logPLMC
openr, 1, 'plr.dat'      ;now read data of LMC cepheids
readf, 1, tit1
readf, 1, tit2
for jj=0,25 do begin
  readf, 1, cant0, cant1, cant2
  logPLMC(jj) = cant0    ; logs of Periods (days) of LMC cephs
  magvLMC(jj) = cant1    ; Apparent visual mags of LMC cephs
  coviLMC(jj) = cant2    ; (V-I) colors of LMC cephs
endfor
close,1
logPLM=fltarr(18)  &  magvLM=logPLM  
openr, 1, 'spc.dat'      ;now read data of short-period LMC cepheids
readf, 1, tit1
readf, 1, tit2
for jj=0,17 do begin
  readf, 1, cant0, cant1
  logPLM(jj) = cant0    ; logs of Periods (days) of LMC cephs
  magvLM(jj) = cant1    ; Apparent visual mags of LMC cephs
endfor
close,1
!x.style=1
!y.style=1
!x.thick=2
!y.thick=2
!p.charsize=2
!p.thick=2
!x.range=[0.3,1.9]
!y.range=[-1.8,-7.2]
!x.title='!17log of period (days)'
!y.title='!17absolute V mags'
Pgal=alog10(Pgal)        ;here we calculate logs of periods
plot, Pgal,Mvga,psym=1
MvS=fltarr(1)
read, 'enter abs V mag (average) of cepheid S Normae:', cant3
logPsnor=[0.989]
MvS(0)=cant3
oplot, logPsnor, MvS, psym=4
read, 'enter E(B-V) for LMC cepheids: ', EBmV
AvLMC=3.1*EBmV
magvLMC=magvLMC-AvLMC    ;the ap mags corrected for extinction
magvLM =magvLM-AvLMC      ;the ap mags corrected for extinction
;now determine (interactively) the distance modulus of the LMC
read, 'enter first estimate of dist mod of LMC: ', mmM
while (mmM ge 0) do begin
MvLMC=magvLMC-mmM
MvLM =magvLM-mmM
plot, Pgal,Mvga,psym=1
oplot, logPsnor, MvS, psym=4
oplot, logPLMC,MvLMC,psym=6
oplot, logPLM,MvLM,psym=6
set_plot, 'ps'
device, /landscape, filename='lmc.ps'
plot, Pgal,Mvga,psym=1
oplot, logPsnor, MvS, psym=4
oplot, logPLMC,MvLMC,psym=6
oplot, logPLM,MvLM,psym=6
endpl
read, 'enter new (m-M) (negative to finish): ', mmM
endwhile

end