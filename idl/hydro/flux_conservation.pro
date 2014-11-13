;main program

dt = 0.01
t0 = 0.0
timesteps = 1000
x0 = 0.0
x1 = 100.0
n  = 3000
k  = 1.0
u0 = 0.0
wait_dt = 0.000


q11 = 0.0
;dx
dx = (x1-x0)/n
;quantity to be advected
q1 = dindgen(n)
;flux at cell interfaces
f1 = dindgen(n+1)
;velocity at cell interaces
u1 = dindgen(n+1)
;x-grid
x = dindgen(n)
;time 
t = dindgen(n)
t(0) = t0

;prepare x-grid
x(0) = x0
for i=1,n-1 do begin
   x(i) = x(i-1) + dx
endfor

;initial conditions and velocity field
for i=0,n-1 do begin
   if(x(i) gt 20.0 and x(i) lt 40) then q1(i) = 1.0 else q1(i) = 0.0
   if(x(i) lt 40.0) then u1(i) = 1.0
   if(x(i) ge 40.0) then u1(i) = 1.0 * ((cos(((x(i)-40.0)/20.0) *1.0*!pi)+1.0) * 0.25) + 0.5
   if(x(i) gt 60.0) then u1(i) = 0.5
endfor

;boundaries
u1(n) = u1(n-1)

;q1 will be advected with flux conservation and q2 without
q2 = q1

;initate flux with donor cell scheme
f1(0) = u1(0) * q1(0)
f1(n) = u1(n) * q1(n-1)
for i=1,n-1 do begin
   if (u1(i) gt 0.0) then q11 = q1(i-1) else q11 = q1(i) 
   f1(i) = u1(i) * q11
endfor

;advect with upstrem differencing
for t=1,timesteps do begin
   ;save old timestep quantities
   q_old  = q1
   q_old1 = q2
   ;update fluxes
   f1(0) = u1(0) * q1(0)
   f1(n) = u1(n) * q1(n-1)
   for i=1,n-1 do begin
      if (u1(i) gt 0.0) then q11 = q1(i-1) else q11 = q1(i) 
      f1(i) = u1(i) * q11 
   endfor
   ;update q with flux conservation
   for i=0,n-1 do begin
      q1(i) = q_old(i) + dt/dx * (f1(i)-f1(i+1))
   endfor
   ;update q without flux conservation
   for i=1,n-2 do begin
      q2(i) = q_old1(i) - dt/dx * u1(i) * (q_old1(i)-q_old1(i-1))
   endfor
   plot,  x, q1, yrange=[0.0,2.0], color=0, /nodata
   oplot, x, u1, psym=0, linestyle=2, color=0
   oplot, x, q1, psym=0, thick=2, color=30
   oplot, x, q2, psym=0, thick=2, color=120
   wait, wait_dt
endfor

print, "dx:",dx
print, "dt:",dt
print, "max(u1) :",max(u1)
print, "max CFL:", (max(u1)*dt)/dx

END

