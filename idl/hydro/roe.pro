

;simple gamma eos
pro eos,dens,pres,ener,eint,enth,velx,vely,velz,gamma

   eint = ener - 0.5*(velx^2+vely^2+velz^2)
   pres = dens*eint*(gamma-1.0)
   ener = eint + 0.5*(velx^2+vely^2+velz^2)
   enth = ener+pres/dens

end

;average state
pro average,q,qa,w,n

   for i=1,n-1 do begin
      qa(i) = (q(i-1)*w(i-1)+q(i)*w(i))/(w(i-1)+w(i)) 
   endfor
   qa(0) = qa(1)
   qa(n) = qa(n-1)

end

;main program
dt = 1.0
dt_max = 10.0
t0 = 0.0
tmax = 5000.0
timesteps = 10000
x0 = 0.0
x1 = 100.0
xc = 0.5*(x1-x0)
n  = 200
k  = 1.0
u0 = 0.0
gamma = 7.0/5.0
dens_left = 1.0e5
dens_right = 1.25e4
pres_left = 1.0
pres_right = 0.1

;dx
dx = (x1-x0)/n

;hydro variables at cell centers
dens = dindgen(n)
pres = dindgen(n)
velx = dindgen(n)
vely = dindgen(n)
velz = dindgen(n)
ener = dindgen(n)
eint = dindgen(n)
enth = dindgen(n)

;averaged hydro variables at cell faces
dens_a = dindgen(n+1)
pres_a = dindgen(n+1)
velx_a = dindgen(n+1)
vely_a = dindgen(n+1)
velz_a = dindgen(n+1)
ener_a = dindgen(n+1)
eint_a = dindgen(n+1)
enth_a = dindgen(n+1)

;conserved quantities at cell centers
q1 = dindgen(n)
q2 = dindgen(n)
q3 = dindgen(n)
q4 = dindgen(n)
q5 = dindgen(n)

;jumps in conserved quantities at cell faces
dq1 = dindgen(n+1)
dq2 = dindgen(n+1)
dq3 = dindgen(n+1)
dq4 = dindgen(n+1)
dq5 = dindgen(n+1)

;averaged states at cell interfaces
q1a = dindgen(n+1)
q2a = dindgen(n+1)
q3a = dindgen(n+1)
q4a = dindgen(n+1)
q5a = dindgen(n+1)

;flux at cell interfaces
f1 = dindgen(n+1)
f2 = dindgen(n+1)
f3 = dindgen(n+1)
f4 = dindgen(n+1)
f5 = dindgen(n+1)
f1l = dindgen(n+1)
f2l = dindgen(n+1)
f3l = dindgen(n+1)
f4l = dindgen(n+1)
f5l = dindgen(n+1)
f1r = dindgen(n+1)
f2r = dindgen(n+1)
f3r = dindgen(n+1)
f4r = dindgen(n+1)
f5r = dindgen(n+1)

;eigenvelocities
l1 = dindgen(n+1)
l2 = dindgen(n+1)
l3 = dindgen(n+1)
l4 = dindgen(n+1)
l5 = dindgen(n+1)

;eigenvelocities
a1 = dindgen(n+1)
a2 = dindgen(n+1)
a3 = dindgen(n+1)
a4 = dindgen(n+1)
a5 = dindgen(n+1)

;eigenvectors
K1 = dindgen(n+1,5)
K2 = dindgen(n+1,5)
K3 = dindgen(n+1,5)
K4 = dindgen(n+1,5)
K5 = dindgen(n+1,5)

;velocity at cell interaces
uf = dindgen(n+1)
;x-grid
x = dindgen(n)
;average weighting
w = dindgen(n)

;prepare x-grid
x(0) = x0
for i=1,n-1 do begin
   x(i) = x(i-1) + dx
endfor

;initial conditions
for i=0,n-1 do begin
   if(x(i) lt xc) then begin
      dens(i) = dens_left
      pres(i) = pres_left
   endif else begin 
      dens(i) = dens_right
      pres(i) = pres_right
   endelse
   velx(i) = 0.0
   vely(i) = 0.0
   velz(i) = 0.0
   eint(i) = pres(i)/(dens(i)*(gamma-1.0))
   ener(i) = eint(i) + 0.5*(velx(i)^2+vely(i)^2+velx(i)^2) 
endfor

;eos
eos,dens,pres,ener,eint,enth,velx,vely,velz,gamma

time = 0.0

for t=1,timesteps do begin

   print, 'timestep',t,'       time',time,'      timestep',dt

   ;conserved quantities
   for i=0,n-1 do begin
      q1(i) = dens(i) 
      q2(i) = dens(i)*velx(i) 
      q3(i) = dens(i)*vely(i) 
      q4(i) = dens(i)*velz(i) 
      q5(i) = ener(i)*dens(i)
   endfor
   ;jumps in conserved quantities
   dq1(0) = 0.0 
   dq2(0) = 0.0 
   dq3(0) = 0.0 
   dq4(0) = 0.0 
   dq5(0) = 0.0 
   dq1(n) = 0.0 
   dq2(n) = 0.0 
   dq3(n) = 0.0 
   dq4(n) = 0.0 
   dq5(n) = 0.0 
   for i=1,n-1 do begin
      dq1(i) = q1(i)-q1(i-1) 
      dq2(i) = q2(i)-q2(i-1) 
      dq3(i) = q3(i)-q3(i-1) 
      dq4(i) = q4(i)-q3(i-1) 
      dq5(i) = q5(i)-q5(i-1) 
   endfor

   ;average 
   average,velx,velx_a,sqrt(dens),n
   average,vely,vely_a,sqrt(dens),n
   average,velz,velz_a,sqrt(dens),n
   average,enth,enth_a,sqrt(dens),n
   vel_a = sqrt(velx_a^2+vely_a^2+velz_a^2)
   cs_a = sqrt((gamma-1.0)*((enth_a-0.5*vel_a^2)))

   ;average eigenvalues
   l1 = velx_a - cs_a
   l2 = velx_a
   l3 = velx_a
   l4 = velx_a
   l5 = velx_a + cs_a

   ;average eigenvectors
   K1(0) = 1.0  
   K1(1) = velx_a - cs_a  
   K1(2) = vely_a
   K1(3) = velz_a
   K1(4) = enth_a - velx_a*cs_a

   K2(0) = 1.0  
   K2(1) = velx_a
   K2(2) = vely_a
   K2(3) = velz_a
   K2(4) = 0.5*vel_a^2 

   K3(0) = 0.0  
   K3(1) = 0.0 
   K3(2) = 1.0 
   K3(3) = 0.0 
   K3(4) = vely_a 

   K4(0) = 0.0  
   K4(1) = 0.0 
   K4(2) = 0.0 
   K4(3) = 1.0 
   K4(4) = velz_a 

   K5(0) = 1.0  
   K5(1) = velx_a + cs_a 
   K5(2) = vely_a 
   K5(3) = velz_a
   K5(4) = enth_a + velx_a*cs_a 

   ;wavestrengths
   a3 = dq3-vely_a*dq1
   a4 = dq4-velz_a*dq1
   dq5a = dq5-(dq3-vely_a*dq1)*vely_a-(dq4-velz_a*dq1)*velz_a
   a2 = (gamma-1.0)/(cs_a^2) * (dq1*(enth_a-velx_a^2)+velx_a*dq2-dq5a)
   a1 = 1.0/(2.0*cs_a) * (dq1*(velx_a+cs_a)-dq2-cs_a*a2)
   a5 = dq1-(a1+a2)

   ;
   for i=1,n-1 do begin
      f1l(i) = q1(i-1)*velx(i-1)
      f2l(i) = q2(i-1)*velx(i-1)+pres(i-1) 
      f3l(i) = q3(i-1)*velx(i-1)
      f4l(i) = q4(i-1)*velx(i-1)
      f5l(i) = (q5(i-1)+pres(i-1))*velx(i-1)
      f1r(i) = q1(i)*velx(i)
      f2r(i) = q2(i)*velx(i)+pres(i) 
      f3r(i) = q3(i)*velx(i)
      f4r(i) = q4(i)*velx(i)
      f5r(i) = (q5(i)+pres(i))*velx(i)
   endfor
   f1l(0) = q1(0)*velx(0)
   f2l(0) = q2(0)*velx(0)+pres(0) 
   f3l(0) = q3(0)*velx(0)
   f4l(0) = q4(0)*velx(0)
   f5l(0) = (q5(0)+pres(0))*velx(0)
   f1r(0) = f1l(0)
   f2r(0) = f2l(0)
   f3r(0) = f3l(0)
   f4r(0) = f4l(0)
   f5r(0) = f5l(0)
   f1l(n) = f1l(n-1)
   f2l(n) = f2l(n-1)
   f3l(n) = f3l(n-1)
   f4l(n) = f4l(n-1)
   f5l(n) = f5l(n-1)
   f1r(n) = f1l(n-1) 
   f2r(n) = f2l(n-1)
   f3r(n) = f3l(n-1)
   f4r(n) = f4l(n-1)
   f5r(n) = f5l(n-1)


   ;compute fluxes, FINALLY!
   f1(0) = 0.0
   f2(0) = 0.0 
   f3(0) = 0.0 
   f4(0) = 0.0 
   f5(0) = 0.0 
   f1(n) = 0.0
   f2(n) = 0.0 
   f3(n) = 0.0 
   f4(n) = 0.0 
   f5(n) = 0.0 

   f1 = 0.5*(f1l+f1r) - 0.5 * (a1*abs(l1)*K1(0) + a2*abs(l2)*K2(0) + $
              a3*abs(l3)*K3(0) + a4*abs(l4)*K4(0) + a5*abs(l5)*K5(0))
   f2 = 0.5*(f2l+f2r) - 0.5 * (a1*abs(l1)*K1(1) + a2*abs(l2)*K2(1) + $
              a3*abs(l3)*K3(1) + a4*abs(l4)*K4(1) + a5*abs(l5)*K5(1))
   f3 = 0.5*(f3l+f3r) - 0.5 * (a1*abs(l1)*K1(2) + a2*abs(l2)*K2(2) + $
              a3*abs(l3)*K3(2) + a4*abs(l4)*K4(2) + a5*abs(l5)*K5(2))
   f4 = 0.5*(f4l+f4r) - 0.5 * (a1*abs(l1)*K1(3) + a2*abs(l2)*K2(3) + $
              a3*abs(l3)*K3(3) + a4*abs(l4)*K4(3) + a5*abs(l5)*K5(3))
   f5 = 0.5*(f5l+f5r) - 0.5 * (a1*abs(l1)*K1(4) + a2*abs(l2)*K2(4) + $
              a3*abs(l3)*K3(4) + a4*abs(l4)*K4(4) + a5*abs(l5)*K5(4))

   ;save old timestep quantities
   q1_old = q1
   q2_old = q2
   q3_old = q3
   q4_old = q4
   q5_old = q5

   ;
   ;advect with upstrem differencing
   ;
   for i=0,n-1 do begin
      
       q1(i) = q1_old(i) + dt/dx * (f1(i)-f1(i+1))
       q2(i) = q2_old(i) + dt/dx * (f2(i)-f2(i+1))
       q3(i) = q3_old(i) + dt/dx * (f3(i)-f3(i+1))
       q4(i) = q4_old(i) + dt/dx * (f4(i)-f4(i+1))
       q5(i) = q5_old(i) + dt/dx * (f5(i)-f5(i+1))
      
   endfor
    
   ;hydro quantities
   dens = q1
   velx = q2/q1
   vely = q3/q1
   velz = q4/q1
   ener = q5/q1

   ;eos
   eos,dens,pres,ener,eint,enth,velx,vely,velz,gamma

   dt = 0.1 * dx/max(abs(velx)) 
   if(dt gt dt_max) then dt = dt_max 
   plot,  x, dens, color=0, psym=-1
   time = time+dt
   if(time gt tmax) then break
  
endfor



END

