;main program

dt = 0.1
t0 = 0.0
timesteps = 300
x0 = 0.0
x1 = 100.0
n  = 100
u  = 1.0



;
dx = (x1-x0)/n
q1 = dindgen(n)
x = dindgen(n)
t = dindgen(n)
t(0) = t0

;prepare x-grid
x(0) = x0
for i=1,n-1 do begin
   x(i) = x(i-1) + dx
endfor

;initial conditions
for i=0,n-1 do begin
   if(x(i) le 30.0) then q1(i) = 1.0 else q1(i) = 0.0
endfor
q2 = q1

;advect with centered differencing
for t=1,timesteps do begin
   ;save old timestep
   q_old = q1
   for i=1,n-2 do begin
      q1(i) = q_old(i) - dt/(2*dx) * u * (q_old(i+1)-q_old(i-1))
   endfor
endfor

;advect with upstrem differencing
for t=1,timesteps do begin
   ;save old timestep
   q_old = q2
   for i=1,n-2 do begin
      q2(i) = q_old(i) - dt/dx * u * (q_old(i)-q_old(i-1))
   endfor
   plot, x, q2, psym = -1, yrange=[0.0,1.4]
   wait, 0.001
endfor

;plot
;plot, x, q1, psym = 0, yrange=[0.0,1.4]
;oplot, x, q2, psym = -6

print, "dx:",dx
print, "dt:",dt
print, "u :",u
print, "CFL:", (u*dt)/dx

END

