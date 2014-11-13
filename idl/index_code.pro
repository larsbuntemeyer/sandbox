function index,ii,jj,kk,face

nxb = 8
nyb = 8
nzb = 8
k3d = 1
k2d = 1

case face of
  "x": faceIndex = jj+(kk-1)*(nyb+1)+(ii/(nxb+1))*(nyb+1)*(nzb+1*k3d)
  "y": faceIndex = 2*(nyb+1)*(nzb+1*k3d)+ii+(kk-1)*(nxb+1)+(jj/(nyb+1))*(nxb+1)*(nzb+1*k3d) 
  "z": faceIndex = 2*(nyb+1)*(nzb+1)+2*(nxb+1)*(nzb+1)+ii+(jj-1)*(nxb+1)+(kk/(nzb+1))*(nxb+1)*(nyb+1) 
  else: print, 'invalid face'
endcase

return, faceIndex

end

function boundaryIndex,ii,jj,kk,face

nxb = 8
nyb = 8
nzb = 8
k3d = 1
k2d = 1

case face of
  "x": faceIndex = jj+(kk-1)*(nyb+1)+(ii/(nxb+1))*(nyb+1)*(nzb+1*k3d)
  "y": faceIndex = 2*(nyb+1)*(nzb+1*k3d)+ii+(kk-1)*(nxb+1)+(jj/(nyb+1))*(nxb+1)*(nzb+1*k3d) 
  "z": faceIndex = 2*(nyb+1)*(nzb+1)+2*(nxb+1)*(nzb+1)+ii+(jj-1)*(nxb+1)+(kk/(nzb+1))*(nxb+1)*(nyb+1) 
  else: print, 'invalid face'
endcase

return, faceIndex

end

k3d = 1
k2d = 1
nxb = 8
nyb = 8
nzb = 8
nguard = 4

ib = nguard
ie = nguard+nxb
jb = nguard
je = nguard+nyb
kb = nguard
ke = nguard+nzb

nrOfFaceValues  = 2*((nxb+1)*(nzb+1*k3d)+(nyb+1)*(nzb+1*k3d)+(nxb+1)*(nyb+1)*k3d)

print, '-------'
print, 'x-faces'
print, '-------'
for i=ib,ie,nxb do begin
    ii=i+1-nguard
    for k=kb,ke do begin
        kk=k+(1-nguard)*k3d
        for j=jb,je do begin
            jj=j+(1-nguard)*k2d
            print, 'i,j,k',i,j,k,'   ii,jj,kk',ii,jj,kk,'         Index',index(ii,jj,kk,'x')
        endfor
    endfor
endfor
print, '-------'
print, 'y-faces'
print, '-------'
for j=jb,je,nyb do begin
    jj=j+1-nguard
    for k=kb,ke do begin
        kk=k+(1-nguard)*k3d
        for i=ib,ie do begin
            ii=i+1-nguard
            print, 'i,j,k',i,j,k,'   ii,jj,kk',ii,jj,kk,'         Index',index(ii,jj,kk,'y')
        endfor
    endfor
endfor
print, '-------'
print, 'z-faces'
print, '-------'
for k=kb,ke,nzb do begin
    kk=k+1-nguard
    for j=jb,je do begin
        jj=j+1-nguard
        for i=ib,ie do begin
            ii=i+1-nguard
            print, 'i,j,k',i,j,k,'   ii,jj,kk',ii,jj,kk,'         Index',index(ii,jj,kk,'z')
        endfor
    endfor
endfor
print, 'nrOfFaceValues:',nrOfFaceValues

print, '-------'
print, 'x-faces'
print, '-------'
for i=ib,ie+1,nxb+1 do begin
    ii=i+1-nguard
    for k=kb,ke do begin
        kk=k+(1-nguard)*k3d
        for j=jb,je do begin
            jj=j+(1-nguard)*k2d
            print, 'i,j,k',i,j,k,'   ii,jj,kk',ii,jj,kk,'         Index',boundaryIndex(ii,jj,kk,'x')
        endfor
    endfor
endfor
print, '-------'
print, 'y-faces'
print, '-------'
for j=jb,je+1,nyb+1 do begin
    jj=j+1-nguard
    for k=kb,ke do begin
        kk=k+(1-nguard)*k3d
        for i=ib,ie do begin
            ii=i+1-nguard
            print, 'i,j,k',i,j,k,'   ii,jj,kk',ii,jj,kk,'         Index',boundaryIndex(ii,jj,kk,'y')
        endfor
    endfor
endfor
print, '-------'
print, 'z-faces'
print, '-------'
for k=kb,ke+1,nzb+1 do begin
    kk=k+1-nguard
    for j=jb,je do begin
        jj=j+1-nguard
        for i=ib,ie do begin
            ii=i+1-nguard
            print, 'i,j,k',i,j,k,'   ii,jj,kk',ii,jj,kk,'         Index',boundaryIndex(ii,jj,kk,'z')
        endfor
    endfor
endfor
print, 'nrOfFaceValues:',nrOfFaceValues

end
