function diff_analytic,x,x0,q0,a,t,t0,dim

   exp_fact = exp(-1.d0*(x-x0)^2/(4.d0*a*(t)))
   diff = q0*(4.d0*!pi*a*(t))^(-0.5d0*dim)*exp_fact

   return, diff

end

c = 2.997d10
D = 1.d4
x0 = 0.d0
x1 = 1.d0
n = 100
dx = (x1-x0)/n
xc = 0.5d0*(x1-x0)
q0 = 1.d5
xrange=[x0,x1]
yrange=[0.1*q0,3.0*q0]

t0 =1.d-7
t = 1.d-6
x = dindgen(n)
q = dindgen(n)
for i=0,n-1 do begin
  x(i) = x0+0.5*dx+i*dx
endfor

for i=0,n-1 do begin
  q(i) = diff_analytic(x(i),xc,q0,D,t,t0,1)
endfor

plot,x,q,xrange=xrange,yrange=yrange

END
