;this is clusters.pro, version 22 Oct 1998, by RHM. Construction
;of color-color and color-mag diagrams for extinction and m-M
;determinations.

V1=fltarr(50)  &  BmV1=V1   &   UmB1=V1
V2=fltarr(29)  &  BmV2=V2   &   UmB2=V2
tit1=''
tit2=''
tit3=''
fname=''
;read, 'enter filename, eg cluple.dat, do not use quotes: ', fname
openr, 1, 'cluple.dat'
readf, 1, tit1
readf, 1, tit2
readf, 1, tit3
for jj=0,49 do begin
  readf, 1, cant0, cant1, cant2 
  V1(jj) = cant0         ; mags V of Pleiades
  BmV1(jj) = cant1       ; colors B-V of Pleiades
  UmB1(jj) = cant2       ; colors U-B of Pleiades
endfor
close,1
openr, 1, 'clu6087.dat'
readf, 1, tit1
readf, 1, tit2
readf, 1, tit3
for jj=0,28 do begin
  readf, 1, cant0, cant1, cant2 
  V2(jj) = cant0         ; mags V of NGC 6087
  BmV2(jj) = cant1       ; colors B-V of NGC 6087
  UmB2(jj) = cant2       ; colors U-B of NGC 6087
endfor
close,1

!x.style=1
!y.style=1
!x.thick=2
!y.thick=2
!p.charsize=2
!p.thick=2
!x.range=[-0.3,0.5]
!x.title='!17B-V'
wset,0
!y.range=[13.,1.]
!y.title='!17V magnitudes'
plot, BmV1,V1,psym=1
oplot, BmV2,V2,psym=6
window,1
wset,1
!y.range=[0.8,-0.8]
!y.title='!17U-B'
plot, BmV1,UmB1,psym=1
oplot, BmV2,UmB2,psym=6

;now determine (interactively) the reddening correction
read, 'enter first estimate of deltaE(B-V): ', dEBmV
while (dEBmV ge 0) do begin
dEUmB=0.84*dEBmV
BmV2o=BmV2-dEBmV
UmB2o=UmB2-dEUmB
plot, BmV1,UmB1,psym=1
oplot, BmV2o,UmB2o,psym=6
set_plot,'ps'
device, /landscape, filename='clu1.ps'
plot, BmV1,UmB1,psym=1
oplot, BmV2o,UmB2o,psym=6
endpl
dAv=3.1*dEBmV
read, 'enter deltaE(B-V) (negative to finish): ', dEBmV
endwhile

V2o=V2-dAv
wset,0
!y.range=[13.,1.]
!y.title='!17V magnitudes'
plot, BmV1,V1,psym=1
oplot, BmV2o,V2o,psym=6
;now determine (interactively) the distance modulus diff
read, 'enter first estimate of delta(m-M): ', dmM
while (dmM ge 0) do begin
V2c=V2o-dmM
plot, BmV1,V1,psym=1
oplot, BmV2o,V2c,psym=6
set_plot,'ps'
device, /landscape, filename='clu2.ps'
plot, BmV1,V1,psym=1
oplot, BmV2o,V2c,psym=6
endpl
read, 'enter delta(m-M) (negative to finish): ', dmM
endwhile

end
