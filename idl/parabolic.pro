;main program

dx = 0.1
t0 = 0.0
timesteps = 300
x0 = 0.0
x1 = 100.0
n  = 100
ni = 10

;
dx = (x1-x0)/n
s1 = dindgen(n)
i1 = dindgen(n)
ds1 = dindgen(n-1)
ds2 = dindgen(n-1)
si = dindgen((n-1)*ni)
tau = dindgen((n-1)*ni)
x1 = dindgen((n-1)*ni)
x = dindgen(n)
dtau = 1.e-3
tau(0) = 0.0

;prepare x-grid
x(0) = x0
for i=1,n-1 do begin
   x(i) = x(i-1) + dx
endfor

;intialize the source function at the grid points
for i=0,n-1 do begin
   if(x(i) le 40.0 or x(i) gt 60.0) then s1(i) = 0.1 else s1(i) = 10.0
   ;s1(i) = 10.0
endfor

;s1 = 5+5.0*randomn(2.0,n) + 5.0
;s1 = 5.0+2.0*sin(0.25*x)
i1(0) = s1(0)

;intialize the source function gradients
ds2(0) = 0.0
ds1(0) = (s1(1)-s1(0)) / dtau
for i=1,n-2 do begin
   ds2(i) = ((s1(i+1)-s1(i))*dtau + (s1(i-1)-s1(i))*dtau) / (dtau*dtau*dtau + dtau*dtau*dtau) 
   ds1(i) = (s1(i+1)-s1(i))/dtau - ((s1(i+1)-s1(i))*dtau + (s1(i-1)-s1(i))*dtau) / (dtau*dtau + dtau*dtau)
;   ds2(i) = 1.0 / (2.0*dtau) * ((s1(i+1)-s1(i)) / dtau - (s1(i)-s1(i-1)) / dtau) 
;   ds1(i) = (s1(i-1)-s1(i+1))/(2.0*dtau) 
endfor

;interpolate the source function between the grid points
for i=0,n-2 do begin
   for j=0,ni-1 do begin
      k = i*ni+j
      x1(k) = x(0) + dx/ni * k
      tau(k) = tau(0)+dtau/ni * k
      si(k) = s1(i) + ((dtau/ni)*j)*ds1(i) + ((dtau/ni)*j)^2*ds2(i)   
   endfor
endfor

;compute source integral
for i=1,n-2 do begin
   xp = exp(-dtau)
   e0 = 1.0 - xp
   e1 = dtau - e0
   e2 = dtau*dtau-2.0*e1
   ;third order
   u  = e0+(e2-(2.0*dtau+dtau)*e1)/(dtau*(dtau+dtau))
   v  = ((dtau+dtau)*e1-e2)/(dtau*dtau)
   w  = (e2-dtau*e1)/(dtau*(dtau+dtau))
   print, 'u',u,'   v',v,'   w',w
   print, 'sum 3:', u+v+w
   q3  = u*s1(i-1) + v*s1(i) + w*s1(i+1)
   ;second order
   u = (1.0-(1.0+dtau)*xp)/dtau
   v = (dtau-1.0+xp)/dtau
   w = 0.0
   print, 'u',u,'   v',v,'   w',w
   print, 'sum 2:', u+v+w
   q2 = u*s1(i-1)+v*s1(i) 
   ; q_max
   q_max1 = 0.5 *  (s1(i)+s1(i-1)) * dtau
;
   print, 'x(i)',x(i)
   print, 'q3:',q3
   print, 'q2:',q2
   print, 'q_max1:',q_max1
   print, 's1(i-1)',s1(i-1),'   s1(i)',s1(i),'   s1(i+1)',s1(i+1)
;   if(i gt 0) then begin
;     print, 'i', i
;     print, 'ds1(i-1)',ds1(i-1),'   ds1(i)',ds1(i),'   ds1(i+1)',ds1(i+1)
;     print, 'ds2(i-1)',ds2(i-1),'   ds2(i)',ds2(i),'   ds2(i+1)',ds2(i+1)
;   endif
   print, '-----------------------------------------------------------------------'
;
;   if(q3 gt q2) then begin
;      print, 'WARNING q3>q2'
;      print, 'dtau',dtau
;      print, 's1(i-1)',s1(i-1),'   s1(i)',s1(i),'   s1(i+1)',s1(i+1)
;   endif
   if(q3 lt 0.0) then begin
      print, 'WARNING q3<0.0'
      print, 'dtau',dtau
      print, 's1(i-1)',s1(i-1),'   s1(i)',s1(i),'   s1(i+1)',s1(i+1)
      print, 'ds1(i-1)',ds1(i-1),'   ds1(i)',ds1(i),'   ds1(i+1)',ds1(i+1)
      print, 'ds2(i-1)',ds2(i-1),'   ds2(i)',ds2(i),'   ds2(i+1)',ds2(i+1)
      q3=q2
   endif
;   if(q3 gt q2) then begin
;      print, 'WARNING q3>q2 '
;      print, 'q2',q2
;      print, 'dtau',dtau
;      print, 'u',u,'   v',v,'   w',w
;      print, 's1(i-1)',s1(i-1),'   s1(i)',s1(i),'   s1(i+1)',s1(i+1)
;      print, '-----------------------------------------------------------------------'
;     ; q3=q2
;   endif
;
   q = q3
   if(abs(s1(i)-s1(i-1)) lt 1.e-8) then begin
     q = q2
     print, '=========================='
     print, '==using 2nd order========='
     print, '=========================='
   endif
   dSdTauP = (s1(i)-s1(i-1))/dtau
   dSdTauN = (s1(i+1)-s1(i))/dtau
   d2Sd2Tau = (dSdTauN-dSdTauP)/(2.0*dtau)
;   q_max1 = 0.5 *  (s1(i)+s1(i-1)) * dtau
;   q_max2 = 0.5 * (1.0-xp)* (s1(i)+s1(i+1))
   print, 'x(i)',x(i),'  |  q',q,'  |  dSdTauP',dSdTauP,'  |  dSdTauN',dSdTauN,'  |  d2Sd2Tau',d2Sd2Tau
   print, '-----------------------------------------------------------------------'
;   if (abs(dSdTauP) lt 0.01) then begin
;      print, 'WARNING! linear interpolation'
;      i1(i) = i1(i-1)*xp + q1 
;   endif else begin
;      print, 'using parabolic interpolation'
;      i1(i) = i1(i-1)*xp + q 
;   endelse
;  q = max(min(q3,q_max1),0.0)
;  q = min(q3,q_max1)
  i1(i) = i1(i-1)*xp + q 
endfor

xp = exp(-dtau)
u = (1.0-(1.0+dtau)*xp)/dtau
v = (dtau-1.0+xp)/dtau
w = 0.0
i1(n-1) = i1(n-2)*xp + u*s1(n-2) + v*s1(n-1)

;plot
plot, x, s1, psym = 1, color=0, thick=4, yrange=[-1.0,12.0]
oplot, x1, si, psym = 0, color=120, thick=1
oplot, x, i1, psym = 0, color=200, thick=2
;oplot, x, ds1, psym = -1, color=0, thick=1
;oplot, x, ds2, psym = -2, color=120, thick=2
END

