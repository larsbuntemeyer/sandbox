;
;olson & kunasz quadrature
;---------------------------------------------------------
function qdr,s1,s2,s3,dtau1,dtau2,order
   xp = exp(-dtau1)

   case order of

   ;third order
   3: begin
      e0 = 1.0 - xp
      e1 = dtau1 - e0
      e2 = dtau1*dtau1-2.0*e1
      u  = e0+(e2-(2.0*dtau1+dtau2)*e1)/(dtau1*(dtau1+dtau2))
      v  = ((dtau1+dtau2)*e1-e2)/(dtau1*dtau2)
      w  = (e2-dtau1*e1)/(dtau2*(dtau1+dtau2))
      end

   ;second order
   2: begin
      u = (1.0-(1.0+dtau1)*xp)/dtau1
      v = (dtau1-1.0+xp)/dtau1
      w = 0.d0
      end

   ;first order
   1: begin
      u  = 0.5d0*(1.d0-xp)
      v  = 0.5d0*(1.d0-xp)
      w  = 0.d0
      end

   else: begin
      print, 'qdr order invalid!'
      stop
      end
   endcase

   dI = u*s1+v*s2+w*s3 

   return, dI

end

;main program
;---------------------------------------------------------
ps = 0
wf = 0
dt_init = 1.d-20
dt_min = 1.d-20
dt_fix = 1.e-14
dt = dt_init
t0 = 0.0
tmax = 1.e10
nplot = 100
e0 = 1.d4
timesteps = 200000
orderOfQdr = 3 
cfl = 0.8
x0 = -0.5
x1 = 0.5
n  = 64
xc = 0.0;0.5*(x1-x0)
kappa = 1000.0
c  = 2.99792458d10
D = c / (1.0*kappa)
dIdt = 1.0
dx = (x1-x0)/n

xrange = [x0,x1]
yrange = [1.e-10,1.2]

;
q = dindgen(n)
q0 = dindgen(n)
j = dindgen(n)
j0 = dindgen(n)
s = dindgen(n)
s_new = dindgen(n)
i_up       = dindgen(n)
i_down     = dindgen(n)
i_up_old   = dindgen(n)
i_down_old = dindgen(n)
x    = dindgen(n)
t = dindgen(n)
error = dindgen(n)
t(0) = t0


;prepare x-grid
x(0) = x0+0.5*dx
for i=1,n-1 do begin
   x(i) = x(i-1) + dx
endfor

;initial conditions
for i=0,n-1 do begin
   r_sqr = (x(i)-xc)^2
   q0(i) = e0 / ((4.d0*!pi*D*t0)^(0.5d0)) * exp(-r_sqr/(4.d00*D*t0))
   q(i) = q0(i) 
   s(i) = q0(i)
   j(i) = 0.0
   i_up_old(i)   = 0.0 
   i_down_old(i) = 0.0 
endfor
stop
irradiation_left = 0.0
 
;init timer
time = t0
plot_counter = 0

   ; compute time-independent solution

   dtau = kappa/(sqrt(1.0)) * dx
   ; integrate upwards
   s1 = s(0)
   s2 = s(0)
   s3 = s(1)
   i_up(0) = irradiation_left ;qdr(s1,s2,s3,dtau,dtau,1)
   for i=1,n-2 do begin
      s1 = s(i-1)
      s2 = s(i)
      s3 = s(i+1)
      dI = qdr(s1,s2,s3,dtau,dtau,orderOfQdr)
      if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,3)
      i_up(i) = i_up(i-1)*exp(-dtau) + dI
   endfor
   s1 = s(n-2)
   s2 = s(n-1)
   s3 = 0.0
   dI = qdr(s1,s2,s3,dtau,dtau,2)
   if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,1)
   i_up(n-1) = i_up(n-2)*exp(-dtau) + dI 

   ; integrate downwards
   s1 = s(n-1)
   s2 = s(n-1)
   s3 = s(n-2)
   i_down(n-1) = 0.0;s2;qdr(s1,s2,s3,dtau,dtau,1)
   for i=n-2,1,-1 do begin
      s1 = s(i+1)
      s2 = s(i)
      s3 = s(i-1)
      dI = qdr(s1,s2,s3,dtau,dtau,orderOfQdr)
      if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,3)
      i_down(i) = i_down(i+1)*exp(-dtau) + dI
   endfor
   s1 = s(1)
   s2 = s(0)
   s3 = 0.0
   dI = qdr(s1,s2,s3,dtau,dtau,2)
   if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,1)
   i_down(0) = i_down(1)*exp(-dtau) + dI 

   for i=0,n-1 do begin
       j0(i) = 0.5*(i_up(i)+i_down(i))
   endfor
;evolve
for t=0,timesteps do begin

   print, 'time:',t, time, dt
   dIdt_fact = dIdt * 1.0 / (c*dt)
   kappa_hat = (kappa + dIdt_fact)
   dtau = kappa_hat/(sqrt(1.0)) * dx

   ; integrate upwards
   s1 = kappa/kappa_hat*(s(0)+1.0/(kappa)*dIdt_fact*i_up_old(0))
   s2 = kappa/kappa_hat*(s(0)+1.0/(kappa)*dIdt_fact*i_up_old(0))
   s3 = kappa/kappa_hat*(s(1)+1.0/(kappa)*dIdt_fact*i_up_old(1))
   i_up(0) = irradiation_left ;qdr(s1,s2,s3,dtau,dtau,1)
   for i=1,n-2 do begin
      s1 = kappa/kappa_hat*(s(i-1)+1.0/(kappa)*dIdt_fact*i_up_old(i-1))
      s2 = kappa/kappa_hat*(s(i)+1.0/(kappa)*dIdt_fact*i_up_old(i))
      s3 = kappa/kappa_hat*(s(i+1)+1.0/(kappa)*dIdt_fact*i_up_old(i+1))
      dI = qdr(s1,s2,s3,dtau,dtau,orderOfQdr)
      if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,orderOfQdr-1)
      i_up(i) = i_up(i-1)*exp(-dtau) + dI
   endfor
   s1 = kappa/kappa_hat*(s(n-2)+1.0/(kappa)*dIdt_fact*i_up_old(n-2))
   s2 = kappa/kappa_hat*(s(n-1)+1.0/(kappa)*dIdt_fact*i_up_old(n-1))
   s3 = 0.0
   dI = qdr(s1,s2,s3,dtau,dtau,2)
   if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,1)
   i_up(n-1) = i_up(n-2)*exp(-dtau) + dI 

   ; integrate downwards
   s1 = kappa/kappa_hat*(s(n-1)+1.0/(kappa)*dIdt_fact*i_down_old(n-1))
   s2 = kappa/kappa_hat*(s(n-1)+1.0/(kappa)*dIdt_fact*i_down_old(n-1))
   s3 = kappa/kappa_hat*(s(n-2)+1.0/(kappa)*dIdt_fact*i_down_old(n-2))
   i_down(n-1) = 0.0;s2;qdr(s1,s2,s3,dtau,dtau,1)
   for i=n-2,1,-1 do begin
      s1 = kappa/kappa_hat*(s(i+1)+1.0/(kappa)*dIdt_fact*i_down_old(i+1))
      s2 = kappa/kappa_hat*(s(i)+1.0/(kappa)*dIdt_fact*i_down_old(i))
      s3 = kappa/kappa_hat*(s(i-1)+1.0/(kappa)*dIdt_fact*i_down_old(i-1))
      dI = qdr(s1,s2,s3,dtau,dtau,orderOfQdr)
      if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,orderOfQdr-1)
      i_down(i) = i_down(i+1)*exp(-dtau) + dI
   endfor
   s1 = kappa/kappa_hat*(s(1)+1.0/(kappa)*dIdt_fact*i_down_old(1))
   s2 = kappa/kappa_hat*(s(0)+1.0/(kappa)*dIdt_fact*i_down_old(0))
   s3 = 0.0
   dI = qdr(s1,s2,s3,dtau,dtau,2)
   if(dI lt 0.0) then dI = qdr(s1,s2,s3,dtau,dtau,1)
   i_down(0) = i_down(1)*exp(-dtau) + dI 

   ; compute mean intensity and save old intensities
   ds_max = 0.d0
   for i=0,n-1 do begin
       i_up_old(i) = i_up(i) 
       i_down_old(i) = i_down(i)
       j(i) = 0.5*(i_up(i)+i_down(i))
       s_new(i) = s(i); + kappa*c*dt*(j(i)-s(i))
       ;if(dIdt gt 0.0) then s_new(i) = j(i)
       ds = abs(s_new(i)-s(i))/s(i)
       if (s_new(i) lt 1.d-40) then begin
           s_new(i)=1.d-40
           ds = 0.0
       end
       if (ds gt ds_max) then ds_max = ds
       s(i) = s_new(i)
       error(i) = abs(s(i)-q(i))/q(i)
   endfor
   if (time gt tmax) then begin
      print, "tmax exceeded, stopping..."
      break
   end
   if(plot_counter gt nplot) then begin
      plot,  x, j0, xrange=xrange,yrange=yrange, linestyle=2,/ylog
      oplot,  x, j, psym=-1 
      plot_counter = 0.0
   end
   ;new timestep
   dt_old = dt
   ;if(ds_max gt 0.0) then dt = cfl*dt/ds_max
   dt = cfl*dx/c
   if(dt lt dt_min) then dt = dt_min
   if(dt gt 2.0*dt_old) then dt = 2.0*dt_old
   if(dt_fix gt 0.0) then dt = dt_fix
   ;if(dIdt gt 0.0) then dt = dt_init
   time = time + dt
   plot_counter = plot_counter + 1
endfor

;plot
print, 'plotting...'
plot, x, q, xrange=xrange,yrange=yrange, linestyle=2
oplot,  x, j, psym=1 

sn = string(n,FORMAT='(i4.4)')
if(ps) then begin
   items = ['numerically', 'analytically']
   psyms = [-15, -16]
   colors = ['red7', 'blu7']
   cgPS_Open, 'diffusion_1D_' + sn + '.ps'
   xmin = min(x)
   xmax = max(x)
   xshift = 0.1*(xmax-xmin)
   ymin = min(q0)
   ymax = max(q0)
   yshift = 0.1*(ymax-ymin)
   xrange=[xmin-xshift,xmax+xshift]
   yrange=[1.e-4,1.e6]
   xtitle='x'
   ytitle='S'
   cgplot,  x, q0, xrange=xrange,yrange=yrange, linestyle=2,  $
            xtitle=xtitle, ytitle=ytitle, /ylog
   cgplot, x, q, /overplot
   cgplot, x, s, psym = 1, /overplot
   cgPS_Close
end

if(wf) then begin
   outfile = 'diffusion_1D_' + sn + '.dat'
   openw,1, outfile
   printf,1, 'starting time', t0
   printf,1, 'total time', time
   printf,1, 'total timesteps', t
   printf,1, 'cfl', cfl
   printf,1, 'kappa', kappa
   printf,1, 'x,q0,q,s'
   for i=0,n-1 do begin
      printf,1,x(i),q0(i),q(i),s(i)
   endfor
   close,1
   print, "written file:",outfile
end

print, "initial time t0:",t0, " s"
print, "total time:",time
print, "evolved time:",time-t0
print, "evolved time/initial time:",(time-t0)/t0
print, "Diffusion Coefficient, D:",D," s^-1"
print, "D*evolved time:",D*(time-t0)
print, "D^-1:",1.0/D," s"
print, "e0:",e0," erg"
print, 'dIdt_fact:', dIdt_fact
print, 'dt:', dt 
print, 'dx/c', dx/c
print, 'dx/(c*dt)',dx/(c*dt)
print, 'dtau/(c*dt)',dtau/(c*dt)
print, 'kappa*c*dt',kappa*c*dt
;
end

