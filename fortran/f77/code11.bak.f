ccccccccccccccccccccccccc     Program 1.1     cccccccccccccccccccccccccc
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c Please Note:                                                         c
c                                                                      c
c (1) This computer program is part of the book, "An Introduction to   c
c     Computational Physics," written by Tao Pang and published and    c
c     copyrighted by Cambridge University Press in 1997.               c
c                                                                      c
c (2) No warranties, express or implied, are made for this program.    c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      PROGRAM ONE_D_MOTION
C
C Program for the motion of a particle subject to an external
C force f(x) = -x.   We have divided the total time 2*pi into
C 10000 intervals with an equal time step.   The position and
C velocity of the particle are written out at every 500 steps.
C This program is written by Tao Pang in 1997.
C
      PARAMETER (N=10001,IN=500)
      REAL T(N),V(N),X(N)
C
C Assign constants, initial position, and initial velocity
C
      PI   = 4.0*ATAN(1.0)
      DT   = 2.0*PI/FLOAT(N-1)
      X(1) = 0.0
      T(1) = 0.0
      V(1) = 1.0
C
C Recursion for position and velocity at later time
C
      DO      100  I = 1, N-1
        T(I+1) = DT*I
        X(I+1) = X(I)+V(I)*DT
        V(I+1) = V(I)-X(I)*DT
  100 CONTINUE
C
C Write the position and velocity every 500 steps
C
      WRITE (6,999) (T(I),X(I),V(I),I=1,N,IN)
      STOP
  999 FORMAT (3F16.8)
      END
