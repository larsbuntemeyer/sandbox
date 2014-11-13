;--------------------------------------------------------------------
;Copyright (c) 2008 Sebastian Knop
;All rights reserved.
;
;Redistribution and use in source and binary forms, with or without
;modification, are permitted provided that the following conditions
;are met:
;1. Redistributions of source code must retain the above copyright
;   notice, this list of conditions and the following disclaimer.
;2. Redistributions in binary form must reproduce the above copyright
;   notice, this list of conditions and the following disclaimer in the
;   documentation and/or other materials provided with the distribution.
;3. The name of the author may not be used to endorse or promote products
;   derived from this software without specific prior written permission.
;
;   THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
;   IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
;   OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
;   IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
;   INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
;   NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
;   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
;   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
;   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
;   THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
pro spacetimewindow
; this procedure redraws the plotwindow with the actual parameters in the
; common block
COMMON mass_data2, mass_m, mass_gramm, radius_rs, radius_cm, rschwarz, fac
window,0,xsize=700,ysize=700,title='Plotting window',xpos=520, ypos=40
;this array is just for displaying purposes.
;we want the sphere of tangency. so fill an array with
;2 pi/1000 angles with constant radius.....
layerrmu=dblarr(2,1000)
x=dblarr(1001)
y=dblarr(1001)
for i=0,999 do begin
layerrmu[0,i]=rschwarz
layerrmu[1,i]=2.0*3.1415927/1000*i
x[i]=cos(layerrmu[1,i])*rschwarz
y[i]=sin(layerrmu[1,i])*rschwarz
endfor
bndry=1;*rschwarz
; plot the circle of the innermost layer and define the plotsize
plot, layerrmu[0,*] , layerrmu[1,*] ,/polar,xrange=[-fac*bndry,fac*bndry], yrange=[-fac*bndry,fac*bndry],xstyle=9, ystyle=9,xtitle='!6X [cm]', ytitle='!6Y [cm]',xmargin=[12,6.25],ymargin=[4,4],charsize=1.5,thick=2 ,/isotropic,color=0,background=255
axis,xaxis=1,xrange=[!X.CRANGE[0]/rschwarz,!X.CRANGE[1]/rschwarz],xtitle='!6 X [1/R!I s!N]',xstyle=1,charsize=1.5,color=0
axis,yaxis=1,yrange=[!y.CRANGE[0]/rschwarz,!y.CRANGE[1]/rschwarz],ytitle='!6 Y [1/R!I s!N]',ystyle=1,charsize=1.5,color=0
x[1000]=x[0]
y[1000]=y[0]
polyfill,x,y,color=0
massstring=string(mass_gramm)
massenangabe=strcompress(strjoin(['!6Masse: ',massstring,' !6g'],/single))
xyouts,0.4,0.95 ,massenangabe , /normal, color=0,charsize=1.5
massstring=string(mass_m)
massenangabe=strcompress(strjoin(['!6Masse: ',massstring,' !61/M_solar'],/single))
xyouts,0.4,0.92 ,massenangabe , /normal, color=0,charsize=1.5
end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
function raydgl, x, uphi_rd
; differential equation for the photon orbits.
COMMON mass_data2, mass_m, mass_gramm, radius_rs, radius_cm, rschwarz, fac
return ,[uphi_rd[1],1.5d0*rschwarz*uphi_rd[0]^2-uphi_rd[0]]
end
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
pro plotray
;
; makes a simple Runge-Kutta integration of the photon orbits.
; uses the variables in the common blocks and overplots in the
; only active window.
;
COMMON mass_data2, mass_m, mass_gramm, radius_rs, radius_cm, rschwarz, fac
COMMON coords2, r, phi, drdphi
; define the local arrays
uv_plus=dblarr(2)
uv_minus=dblarr(2)
duvdphi_plus=dblarr(2)
duvdphi=dblarr(2)
rphi=dblarr(4,1000000)
; fill the stuff
uv_plus[0]=1.d0/r
uv_plus[1]=-1.d0/r^2*drdphi
uv_minus[0]=1.d0/r
uv_minus[1]=-1.d0/r^2*drdphi

; check if we have hit the unstable orbit!
if (r eq 1.5d0*rschwarz) then begin
for i=1,999 do begin
rphi[0,i]=r
rphi[1,i]=2.0*3.1415927/1000*i
endfor
indx=where(rphi[0,*] ne 0.d0)
length=n_elements(indx)
rphi[0:1,0:length-1]=rphi[0:1,indx]
oplot,rphi[0,0:length-1],rphi[1,0:length-1],color=0,/polar,thick=2
goto, end1
endif
;
path_plus=phi
path_minus=phi
;
; if the derivative is too steep the stepsize will be too big as there
; will only be one step. In order to avoid this we dynamically decrease
; the stepsize.
step_plus=1.d-4
if (abs(drdphi) gt 1.d5) then step_plus=step_plus/alog(1.d2*abs(drdphi))
step_minus=-1.d-4
if (abs(drdphi) gt 1.d5) then step_minus=step_minus/alog(1.d2*abs(drdphi))
;
i=0L
; copy the starting points
rphi[0,i]=1.d0/uv_plus[0]
rphi[1,i]=path_plus
rphi[2,i]=1.d0/uv_minus[0]
rphi[3,i]=path_minus
; integrate the positive branch
while(1.d0/uv_plus[0] gt rschwarz and 1.d0/uv_plus[0] lt 2.d0*fac) do begin
duvdphi_plus=raydgl(path_plus,uv_plus)
uv_plus=rk4(uv_plus,duvdphi_plus,path_plus,step_plus,'raydgl',/double)
;print, 1.d0/uv_plus[0],path_plus
i=i+1L
rphi[0,i]=1.d0/uv_plus[0]
rphi[1,i]=path_plus
path_plus=path_plus+step_plus
endwhile
; integrate the negative branch
i=0L
while(1.d0/uv_minus[0] gt rschwarz and 1.d0/uv_minus[0] lt 2.d0*fac) do begin
duvdphi_minus=raydgl(path_minus,uv_minus)
uv_minus=rk4(uv_minus,duvdphi_minus,path_plus,step_minus,'raydgl',/double)
i=i+1L
rphi[2,i]=1.d0/uv_minus[0]
rphi[3,i]=path_minus
path_minus=path_minus+step_minus
endwhile
; post processing, just use the nonzero entries for plotting and
; color code the rays.
;
; positive branch
indx=where(rphi[0,*] ne 0.d0)
length=n_elements(indx)
rphi[0:1,0:length-1]=rphi[0:1,indx]
if (rphi[0,1] - rphi[0,2] gt 0.d0 )  then begin
oplot,rphi[0,0:length-1],rphi[1,0:length-1],color=80,/polar
endif else begin
oplot,rphi[0,0:length-1],rphi[1,0:length-1],color=250,/polar
endelse
; negative branch
indx=where(rphi[2,*] ne 0.d0)
length=n_elements(indx)
rphi[2:3,0:length-1]=rphi[2:3,indx]
if (rphi[2,1] - rphi[2,2] gt 0.d0 )  then begin
oplot,rphi[2,0:length-1],rphi[3,0:length-1],color=80,/polar
endif else begin
oplot,rphi[2,0:length-1],rphi[3,0:length-1],color=250,/polar
endelse
end1:
end
;
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
pro spacetimeplot
print, '- - - - S T A R T I N G   U P   G E O D E S I C    P L O T T E R - - - -'
COMMON physconst2, M_solar,gravconst,clight,pi
COMMON mass_data2, mass_m, mass_gramm, radius_rs, radius_cm, rschwarz, fac
COMMON coords2, r, phi, drdphi
loadct, 39;, ncol=16
device,retain=2
device,true_color=24
device,direct_color=24
device,pseudo_color=24
device,decomposed=0
; strings needed in advance
ch=''
errorstring=''
; print the help
   		spawn, 'clear'
		spawn, 'echo "********************************"'
		spawn, 'echo Hilfe zum Photonenbahnen-Plotter'
		spawn, 'echo "********************************"'
		spawn, 'echo'
		spawn, 'echo "Dieses Programm stellt die Lichtwege (auch Photonenbahnen genant) um ein kompaktes  "'
		spawn, 'echo "Objekt grafisch dar. Die Lichtstrahlen kann man sich als einen Laserstrahl aus      "'
		spawn, 'echo "einem Laserpointer vorstellen. Im alltaeglichen Leben erscheint der Lichtweg eines  "'
		spawn, 'echo "solchen Lichtstrahls immer als gradlinig.                                           "'
		spawn, 'echo "Die allgemeine Relativitaeetstheorie sagt jedoch eine Ablenkung von Licht durch     "'
		spawn, 'echo "Massen vorraus. Dieser Effekt wird aber nur fuer sehr grosse Massen sichtbar, wes-  "'
		spawn, 'echo "halb er in der Astrophysik eine wichtige Rolle spielt.                              "'
		spawn, 'echo "Dieses Programm berechnet nun die Lichtwege um eine solche Masse, wobei die Para-   "'
		spawn, 'echo "meter des Systems nach Belieben veraendert werden koennen.                          "'
		spawn, 'echo "Das bedeutet, dass die vorhandene Masse wie auch der Startpunkt des Lichtstrahls    "'
		spawn, 'echo "veraendert werden koennen und die entprechenden Bahnen dargestellt werden.          "'
		spawn, 'echo "Die dafuer notwendigen Befehle werden im Folgenden kurz erlaeutert.                 "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Alle Befehle muessen mit der Tastatur eingegeben werden und mit RETURN bestaetigt  "'
		spawn, 'echo " werden!                                                                            "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Koordinatenwahl mit der Maus:                                                      "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Startpunkt im Fenster eintragen         :  s                              "'
		spawn, 'echo "          Durch einen Klick mit der linken Maustaste an der gewuenschten Position   "'
		spawn, 'echo "          im Fenster, wird der entsprechende Photonenorbit dargestellt.             "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Koordinatenwahl mit der Tastatur:                                                  "'
		spawn, 'echo "          Winkel Koordinate                       : p                               "'
		spawn, 'echo "          Radiale Koordinate                      : r                               "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Anstatt die Koordinaten des Startpunktes mit der Maus zu bestimmen        "'
		spawn, 'echo "          koennen die Koordinaten auch in Polarkoordinaten direkt ueber die         "'
		spawn, 'echo "          Tastatur bestimmt werden.                                                 "'
		spawn, 'echo "          Der Winkel wird dabei in Radian angegeben, waehrend die radiale           "'
		spawn, 'echo "          Koordinate in Schwarzschild-Radien gemessen wird (Zur Erinnerung:         "'
		spawn, 'echo "          Fuer eine Sonnenmasse ist ein Schwarzschild-Radius ungefaehr 3 km).       "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Weitere Einstellungen:                                                             "'
		spawn, 'echo "          lokale Ableitung dr/dphi                :  d                              "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          (Fortgeschrittener Einstellung!) Legt die initiale Ausbreitungsrichtung   "'
		spawn, 'echo "          des Lichtstrahles fest. == 0 bedeutet eine tangentiale Bahn.              "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Zentral-Masse                           :  m                              "'
		spawn, 'echo "          Die Masse des Zentralobjektes. Sie wird in Einheiten von der Masse unserer"'
		spawn, 'echo "          Sonne angegeben.                                                          "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Darstellungsgrenzen                     :  g                              "'
		spawn, 'echo "          Legt fest wie gross der Ausschnitt im Plotfenster sein soll. Wird in      "'
		spawn, 'echo "          cm angegeben. Werte koennen in wissenschaftlicher Notation: z.B: 1.d7     "'
		spawn, 'echo "          angegeben werden. Wird dieser Wert geaendert werden alle alten Bahnen     "'
		spawn, 'echo "          geloescht!                                                                "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Koordinaten zuruecksetzen               :  n                              "'
		spawn, 'echo "          Sollte kein vernuenftiges Bild mehr zustande kommen, da die Koordinaten   "'
		spawn, 'echo "          ausserhalb des vernuenftigen Rahmens liegen, kann hiermit die Ausgangs-   "'
		spawn, 'echo "          konfiguration hergestellt werden.                                         "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Alte Bahnen loeschen                             : l                               "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Hilfe                                            : hilfe                           "'
		spawn, 'echo "          Zeigt diese Hilfe an.                                                     "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Beenden                                          : e oder q                        "'
		spawn, 'echo "          Beendet das Programm                                                      "'
		spawn, 'echo "                                                                                    "'
		read,'hit enter',ch
;jump point
restart:
;filling the physconst common
M_solar=1.989d33
gravconst=6.6732d-8
clight=2.997924562d+10
pi=acos(-1.d0)
;print, pi
; here comes the mass_data common
mass_m = 1.d0
mass_gramm = mass_m*M_solar
radius_rs=2.d0
radius_cm=radius_rs*(2*gravconst*mass_gramm/clight^2)
rschwarz=2*gravconst*mass_gramm/clight^2
fac=4.d0*rschwarz
; taking care of the coords common
r=radius_cm
phi=0.5d0*pi
drdphi=0.d0
; create the plotting window
spacetimewindow
; the main loop 
while (ch ne 'ende' and ch ne 'e' and ch ne 'stop' and ch ne 'q' and ch ne 'quit') do begin

 spawn,'clear'
 print,'******************************************************************************'
 print,''
 print,'               ================='
 print,'               Photonen Orbits: '
 print,'               ================='
 print,''
 print,'Koordinatenwahl mit der Maus:'
 print,''
 print,'               Startpunkt im Fenster eintragen         :  s'
 print,''
 print,'Koordinatenwahl mit der Tastatur:'
 print,''
 print,'               Winkel Koordinate        (ändern mit p) : ',phi,'  rad,  '
 print,''
 print,'               Radiale Koordinate       (ändern mit r) : ',radius_rs,'  r/R_s,  ';, radius_cm ,'  cm'
 print,'                                                         ',radius_cm ,'  cm'
 print,'Weitere Einstellungen:'
 print,''
 print,'               lokale Ableitung dr/dphi (ändern mit d) : ',drdphi
 print,''
 print,'               Zentral-Masse            (ändern mit m) : ',mass_m,'  m/M_solar,  ';, mass_gramm,'  g'
 print,'                                                         ',mass_gramm,'  g'
 print,'               Darstellungsgrenzen      (ändern mit g) : ',fac,'  cm,  ';, radius_cm ,'  cm'
 print,''
; print,'  letzte Graphik drucken : d'
 print,'               Koordinaten zuruecksetzen               :  n'
 print,'               Alte Bahnen loeschen                    :  l'
 print,''
 print,''
 print,'  Hilfe                                                : hilfe'
 print,'  Beenden                                              : e oder q'
 print,''
 print,'********************************************************************************'
 if (errorstring ne '') then begin
  print,'[31;01mFEHLER[m: ',errorstring
  print,''
 endif
 errorstring=''

 read,'Geodesic-plotter> ',ch

 switch ch of
   'l': begin
   	  spacetimewindow
          break
        end

   'm': begin
          read,' neue Masse in solaren Einheiten :',ch
	  if (ch eq '') then break
            tmp=double(ch)
            if (tmp lt 0.d0 or tmp gt 1.d100) then begin
              errorstring=' Massenangabe ungültig!'
            endif else begin
              mass_m=tmp
	      mass_gramm = mass_m*M_solar
	      rschwarz=2*gravconst*mass_gramm/clight^2
            endelse
            spacetimewindow
          break
        end

   'r': begin
          read,' neue radiale Koordinate in Schwarzschild-Radien: ',ch
	  if (ch eq '') then break
            tmp=double(ch)
            if (tmp lt 1.d0 or tmp gt 2.d1) then begin
              errorstring=' Radiusangabe ungültig!'
            endif else begin
              radius_rs=tmp
	      radius_cm = radius_rs*(2*gravconst*mass_gramm/clight^2)
            endelse
              r=radius_cm
	      plotray
          break
        end

   'p': begin
          read,' neue Winkel Koordinate in Radian: ',ch
	  if (ch eq '') then break
              tmp=double(ch)
              phi=tmp
	      plotray
          break
        end

   's': begin
          cursor,x,y,/data
;	  print, x,y
;	  xyouts,x,y,'x',color=0
              radius_cm=sqrt(x^2+y^2)
;	      print, radius_cm
	      if (radius_cm lt rschwarz) then begin
              errorstring=' Radiusangabe ungültig!'
	      endif else begin
;	      radius_rs = radius_cm/(2*gravconst*mass_gramm/clight^2)
	      radius_rs = radius_cm/rschwarz
;	      print, radius_rs
	      endelse
              r=radius_cm
	      phi=acos(x/r)
	      if(y lt 0.d0) then phi=-phi
	      plotray
          break
        end


   'g': begin
          read,' neue Darstellungsgrenze in cm: ',ch
	  if (ch eq '') then break
	  tmp=double(ch)
	      if (tmp lt 1.1d0*rschwarz) then begin
              errorstring=' Angabe ungültig!'
	      endif else begin
	      fac=tmp
	      endelse
	      spacetimewindow
          break
        end

   'd': begin
          read,' neue lokale Ableitung: ',ch
	  if (ch eq '') then break
	  tmp=double(ch)
	      if (abs(tmp) gt 5.d9) then begin
	      tmp=5.d9*tmp/abs(tmp)
	      print,'Ableitung zu steil! Begrenzt auf: ',tmp
	      endif
	      drdphi=tmp
 	      plotray
          break
        end
	
   'n': begin
   	  goto, restart
          break
        end

	
   'hilfe': begin
   		spawn, 'clear'
		spawn, 'echo "********************************"'
		spawn, 'echo Hilfe zum Photonenbahnen-Plotter'
		spawn, 'echo "********************************"'
		spawn, 'echo'
		spawn, 'echo "Dieses Programm stellt die Lichtwege (auch Photonenbahnen genant) um ein kompaktes  "'
		spawn, 'echo "Objekt grafisch dar. Die Lichtstrahlen kann man sich als einen Laserstrahl aus      "'
		spawn, 'echo "einem Laserpointer vorstellen. Im alltaeglichen Leben erscheint der Lichtweg eines  "'
		spawn, 'echo "solchen Lichtstrahls immer als gradlinig.                                           "'
		spawn, 'echo "Die allgemeine Relativitaeetstheorie sagt jedoch eine Ablenkung von Licht durch     "'
		spawn, 'echo "Massen vorraus. Dieser Effekt wird aber nur fuer sehr grosse Massen sichtbar, wes-  "'
		spawn, 'echo "halb er in der Astrophysik eine wichtige Rolle spielt.                              "'
		spawn, 'echo "Dieses Programm berechnet nun die Lichtwege um eine solche Masse, wobei die Para-   "'
		spawn, 'echo "meter des Systems nach Belieben veraendert werden koennen.                          "'
		spawn, 'echo "Das bedeutet, dass die vorhandene Masse wie auch der Startpunkt des Lichtstrahls    "'
		spawn, 'echo "veraendert werden koennen und die entprechenden Bahnen dargestellt werden.          "'
		spawn, 'echo "Die dafuer notwendigen Befehle werden im Folgenden kurz erlaeutert.                 "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Alle Befehle muessen mit der Tastatur eingegeben werden und mit RETURN bestaetigt  "'
		spawn, 'echo " werden!                                                                            "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Koordinatenwahl mit der Maus:                                                      "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Startpunkt im Fenster eintragen         :  s                              "'
		spawn, 'echo "          Durch einen Klick mit der linken Maustaste an der gewuenschten Position   "'
		spawn, 'echo "          im Fenster, wird der entsprechende Photonenorbit dargestellt.             "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Koordinatenwahl mit der Tastatur:                                                  "'
		spawn, 'echo "          Winkel Koordinate                       : p                               "'
		spawn, 'echo "          Radiale Koordinate                      : r                               "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Anstatt die Koordinaten des Startpunktes mit der Maus zu bestimmen        "'
		spawn, 'echo "          koennen die Koordinaten auch in Polarkoordinaten direkt ueber die         "'
		spawn, 'echo "          Tastatur bestimmt werden.                                                 "'
		spawn, 'echo "          Der Winkel wird dabei in Radian angegeben, waehrend die radiale           "'
		spawn, 'echo "          Koordinate in Schwarzschild-Radien gemessen wird (Zur Erinnerung:         "'
		spawn, 'echo "          Fuer eine Sonnenmasse ist ein Schwarzschild-Radius ungefaehr 3 km).       "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Weitere Einstellungen:                                                             "'
		spawn, 'echo "          lokale Ableitung dr/dphi                :  d                              "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          (Fortgeschrittener Einstellung!) Legt die initiale Ausbreitungsrichtung   "'
		spawn, 'echo "          des Lichtstrahles fest. == 0 bedeutet eine tangentiale Bahn.              "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Zentral-Masse                           :  m                              "'
		spawn, 'echo "          Die Masse des Zentralobjektes. Sie wird in Einheiten von der Masse unserer"'
		spawn, 'echo "          Sonne angegeben.                                                          "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Darstellungsgrenzen                     :  g                              "'
		spawn, 'echo "          Legt fest wie gross der Ausschnitt im Plotfenster sein soll. Wird in      "'
		spawn, 'echo "          cm angegeben. Werte koennen in wissenschaftlicher Notation: z.B: 1.d7     "'
		spawn, 'echo "          angegeben werden. Wird dieser Wert geaendert werden alle alten Bahnen     "'
		spawn, 'echo "          geloescht!                                                                "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo "          Koordinaten zuruecksetzen               :  n                              "'
		spawn, 'echo "          Sollte kein vernuenftiges Bild mehr zustande kommen, da die Koordinaten   "'
		spawn, 'echo "          ausserhalb des vernuenftigen Rahmens liegen, kann hiermit die Ausgangs-   "'
		spawn, 'echo "          konfiguration hergestellt werden.                                         "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Alte Bahnen loeschen                             : l                               "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Hilfe                                            : hilfe                           "'
		spawn, 'echo "          Zeigt diese Hilfe an.                                                     "'
		spawn, 'echo "                                                                                    "'
		spawn, 'echo " Beenden                                          : e oder q                        "'
		spawn, 'echo "          Beendet das Programm                                                      "'
		spawn, 'echo "                                                                                    "'
		read,'hit enter',ch
		break
            end

 endswitch

endwhile
wdelete
retall
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
