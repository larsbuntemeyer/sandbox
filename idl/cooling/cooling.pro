

c        = 2.997d10
trad     = 10.d0
it       = 1
tgas     = 12.d0
const    = 1.8049444d-05
gasConst = 8.3145119843000d+07  
gamma    = 5.d0/3.d0
abar     = 1.d0
mean     = const * trad^4

kappa    = 1.0
rho      = 1.d-10
trange   = 1.d8
dt       = trange/it


urad     = const*trad^4
ugas     = tgas*gasConst/(abar*(gamma-1.0))
uegas    = trad*gasConst/(abar*(gamma-1.0))
opac     = kappa*rho
time     = dindgen(it+1)
eint     = dindgen(it+1)
erad     = dindgen(it+1)
temp     = dindgen(it+1)
cool     = dindgen(it+1)
eint[0]  = ugas
temp[0]  = tgas
source   = const*tgas^4
qrad     = 4.0*!pi*opac*(mean-source)
emit     = 4.0*!pi*opac*source
abso     = 4.0*!pi*opac*mean
cool[0]  = qrad
time[0]  = 0
erad[0] = 4.d0*!pi*mean/c/rho 
t_resid  = 1.d99

if(emit gt abso) then begin
   t_cool = ugas*rho/(emit-abso)
   print, 'cell is cooling:',t_cool/dt
endif
if(emit lt abso) then begin
   t_heat = (uegas-ugas)*rho/(abso-emit)
   print, 'cell is heating:',t_heat/dt
endif

for n=1,1 do begin

for i=1,it do begin

    time[i] = time[i-1]+dt

    print,i, ' time:',time[i],'               dt:',dt

    source   = const*tgas^4

    qrad = 4.0*!pi*opac*(mean-source)

    ugas = ugas+dt*qrad/rho
;    mean = mean-dt*qrad

    tgas = ugas*abar*(gamma-1.d0)/gasConst

    if(t_resid lt 1.e-4) then break
 
    temp[i] = tgas
    eint[i] = ugas 
    erad[i] = 4.d0*!pi*mean/c/rho 
    cool[i] = qrad

endfor

xrange = [min(time),max(time)]
yrange = [0.5*min(temp),max(temp)*1.5]
;plot, time,temp,/xlog,/ylog
;oplot, time,temp

endfor


end
