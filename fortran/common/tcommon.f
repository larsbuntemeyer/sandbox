      program tcommon
C
C      implicit none
C
C _X(Here we use implicit declarations ... all REAL variables)
C
      include "block.h"
C
      logical logic
C
      call suba()
      call subb()
      call subc()
C
      logic = .true.
C
      write(6,'(1x,"a,b,c=",3f7.2,1L)') a,b,c,logic
      stop
      end
c=======================================================================
      subroutine suba()
      include "block.h"
C _X(must declare 'abc' exactly the same in each routine)
      integer a
      a = 1.0
      return
      end
c=======================================================================
      subroutine subb()
      include "block.h"
C _X(all variables are visible to each routine that accesses 'abc')
      b = 2.0
      return
      end
c=======================================================================
      subroutine subc()
      include "block.h"
C _X(Common source of problems, suppose implicitly)
C _X(declare 'C' variables complex ... yields)
C _X(different size common block)
      c = 3.0
      return
      end
