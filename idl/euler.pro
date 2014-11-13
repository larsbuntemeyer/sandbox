;main program

n = 100
t0 = 0.0
tn = 10.0
q0 = 1.0

;
dt = (tn-t0)/n
t = dindgen(n)
q1 = dindgen(n)
q2 = dindgen(n)
x = dindgen(n)
print, t
t(0) = t0
q1(0) = q0
q2(0) = q0
advect(0) = q0

;prepare t
for i=1,n-1 do begin
   t(i) = t(i-1) + dt
   print,t(i)
endfor

print, 'initial: ',t(0),q(0)

;forward euler
for i=1,n-1 do begin
   q1(i) = q1(i-1) - dt * q1(i-1)^2
   print, 'numerical: ',t(i),q1(i)
endfor

;backward euler
for i=1,n-1 do begin
   q2(i) = q2(i-1)/(1+dt)
   print, 'numerical: ',t(i),q2(i)
endfor
print, 'dt:',dt

;plot
plot, t, q1, psym = 1
;oplot, t, q2, psym = 2
oplot, t, 1.0/(t+1)

END

