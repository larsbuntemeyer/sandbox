

c        = 2.997d10
trad     = 100.d0
it       = 3000
tgas     = 100.d0
const    = 1.8049444d-05
gasConst = 8.3145119843000d+07  
gamma    = 5.d0/3.d0
abar     = 1.d0
mean     = const * trad^4

urad     = 1.d12
ugas     = 1.d2
kappa    = 4.d-8
rho      = 1.d-7
trange   = 1.d-7
dtini    = 1.d-20
dtmax    = 1.d0
cfl      = 0.1
trange   = 1.e-3

tgas     = ugas*abar*(gamma-1.d0)/gasConst
opac     = kappa*rho
ei       = tgas*gasConst/(abar*(gamma-1.d0))
time     = dindgen(it+1)
eint     = dindgen(it+1)
erad     = dindgen(it+1)
temp     = dindgen(it+1)
cool     = dindgen(it+1)
eint[0]  = ei
temp[0]  = tgas
erad[0]  = urad
planck   = const*tgas^4
qrad     = dtini*(c*kappa*urad-4.d0*!pi*kappa*planck)
cool[0]  = qrad
dt       = cfl*dtini*(ugas+qrad)/abs(qrad)
time[0]  = dt
dt = trange/it

for n=1,1 do begin

for i=1,it do begin


    time[i] = time[i-1]+dt

    print, 'time:',time[i],'               dt:',dt

    planck = const*tgas^4

    qrad = dt*(c*kappa*urad-4.d0*!pi*kappa*planck)

    ugas_old = ugas 
    ugas = ugas + qrad 
    urad = urad - qrad
 
    tgas = ugas*abar*(gamma-1.d0)/gasConst

    temp[i] = tgas
    eint[i] = ugas 
    erad[i] = urad 
    cool[i] = qrad

;    dt = cfl * dt * ugas/abs(qrad)
;    if (dt gt dtmax) then dt = dtmax

 ;   if(abs(qrad/ugas) lt 1.e-20) then begin
 ;      dt = cfl * dt * ugas/abs(qrad)
 ;      if (dt gt dtmax) then dt = dtmax
 ;   endif else begin
 ;      dt = dtmax ;cfl * dt * ugas/abs(qrad)
 ;   endelse
    

endfor

xrange = [min(time),max(time)]
yrange = [0.9*min(temp),max(temp)*1.1]
plot, time,temp

endfor


end
