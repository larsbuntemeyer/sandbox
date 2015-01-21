!cccccccccccccccccccccccc     Program 1.1     cccccccccccccccccccccccccc

! Code converted using TO_F90 by Alan Miller
! Date: 2015-01-07  Time: 14:54:29

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                      c
! Please Note:                                                         c
!                                                                      c
! (1) This computer program is part of the book, "An Introduction to   c
!     Computational Physics," written by Tao Pang and published and    c
!     copyrighted by Cambridge University Press in 1997.               c
!                                                                      c
! (2) No warranties, express or implied, are made for this program.    c
!                                                                      c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

PROGRAM one_d_motion

! Program for the motion of a particle subject to an external
! force f(x) = -x.   We have divided the total time 2*pi into
! 10000 intervals with an equal time step.   The position and
! velocity of the particle are written out at every 500 steps.
! This program is written by Tao Pang in 1997.

PARAMETER (n=10001,in=500)
!     REAL T(N),V(N),X(N)

! Assign constants, initial position, and initial velocity

pi   = 4.0*ATAN(1.0)
dt   = 2.0*pi/FLOAT(n-1)
x(1) = 0.0
t(1) = 0.0
v(1) = 1.0

! Recursion for position and velocity at later time

DO        i = 1, n-1
  t(i+1) = dt*i
  x(i+1) = x(i)+v(i)*dt
  v(i+1) = v(i)-x(i)*dt
END DO

! Write the position and velocity every 500 steps

WRITE (6,999) (t(i),x(i),v(i),i=1,n,in)
STOP
999 FORMAT (3F16.8)
END PROGRAM one_d_motion
