! -------------------------------------------------------------------
!
!     Solve the Kirchhoff-Master equations for arbitrary 
!     functions of the local voltage u. Newton solutions for the p, 
!     Monte Carlo solution for the u. Marcus production version.
!
!     metallic lead model
!
!     Numerical evaluation of the Jacobian matrix, the
!     KM equations have to be encoded only once (subroutine fxn) 
!
!     notation throughout the program:
!
!            A     C
!           ->    ->
!      n-1      n      n+1
!           <-    <-
!            B     D
!
!     A - B = C - D
!     A - B - C + D = 0
!
!     forward only motion implies that B = D = 0
!
!     T. Koslowski, M. Castellano, University of Freiburg 2023
!
!     Makes use of the following linpack subroutines and functions:
!     dgefa, dgesl, daxpy, dscal, ddot, idamax
!
!     requires a double precision random number generator (ranf)
!
! --------------------------------------------------------------------

      program newtmc

      parameter (nmax=100)

      implicit real*8 (a-h,o-z)

      common uin,uout,gin,gout,dgin,dgout,alin,alout,tk,beta

      real*8 p(nmax),q(nmax)      !  populations, holes
      real*8 pold(nmax)           !  old populations
      real*8 g(nmax,nmax)         !  conductivity matrix
      real*8 deltag(nmax,nmax)    !  driving forces
      real*8 alambda(nmax,nmax)   !  reorganization energies
      real*8 u(nmax)              !  voltages
      real*8 f(nmax)              !  functions (see SR fxn)

!     workspace for solving linear systems of equations
      real*8 amat(nmax,nmax),b(nmax)
      integer iwork(nmax)

!     voltage and population backups from old Monte Carlo cycles
      real*8 usave(nmax),psave(nmax)

!     control switches - use saved voltage and population
!     backups from old calculations
      logical use_old_u,use_old_p

! ------------------------------------------------------------------

      use_old_u=.false. 
      use_old_p=.false. 

!     set some matrices to zero
      do i=1,nmax
        do j=1,nmax
          g(i,j)=0.0d0
          deltag(i,j)=0.0d0
          alambda(i,j)=0.0d0
        enddo
      enddo

!     input of Marcus theory data
      read (5,*) n 
      do i=1,n-1
        read (5,*) g(i,i+1),deltag(i,i+1),alambda(i,i+1)
        g(i+1,i)=g(i,i+1)
        deltag(i+1,i)=-deltag(i,i+1)
        alambda(i+1,i)=alambda(i,i+1)
      enddo
      read (5,*) gin,gout
      read (5,*) dgin,dgout
      read (5,*) alin, alout

      tk=0.025d0       ! thermal energy
      mcstep=10000*n
      beta=1.0d0/tk
      iseed=91881

! -----------------------------------------------------

      open (unit=8,file='newtp.deb')
      open (unit=10,file='current.out')
      open (unit=14,file='pop.out')
      open (unit=16,file='mu.out')

      nc=0

!     voltage / chemical potential loop

      do iiuu=-600,600,10

      if (iiuu.eq.0) goto 200

      nc=nc+1

      write (6,*) iiuu, ' ------------------------ '

      uin = 0.01d0*iiuu
      uout= 0.0d0*iiuu

!     initial p
      if (nc.eq.1.or.(.not.use_old_p)) then
        do i=1,n
          p(i)=0.5d0
        enddo
      else
        do i=1,n
          p(i)=psave(i)
        enddo
      endif

!     initial u 
      du=(uin-uout)/(n+1)
      ustep=0.1d0*dabs(du)
      if (nc.eq.1.or.(.not.use_old_u)) then
        do i=1,n
          u(i)=uin-i*du
        enddo
      else
        do i=1,n
          u(i)=usave(i)
        enddo
      endif

      call newtsolv (n,nmax,amat,b,u,f,p,q,g,deltag,alambda,iwork,info)
      call fxn (n,nmax,u,p,q,g,deltag,alambda,f,2)
      cold=f(1)

!     Monte Carlo loop for u

      do istep=1,mcstep

 100    continue
        ii=int(ranf(iseed)*n)+1
        unew=uold+(0.5d0-ranf(iseed))*ustep

        uold=u(ii)
        do i=1,n
          pold(i)=p(i)
        enddo
        u(ii)=unew
        itest=1

        call newtsolv (n,nmax,amat,b,u,f,p,q, &
                       g,deltag,alambda,iwork,info)
        call fxn (n,nmax,u,p,q,g,deltag,alambda,f,2)
        cnew=f(1)
        if (info.ne.0) itest=0
        do i=1,n
          if (p(i).lt.0.0d0) itest=0
          if (p(i).gt.1.0d0) itest=0
        enddo

        dcold=dabs(cold)
        dcnew=dabs(cnew)
        if (dcnew.gt.dcold.and.itest.eq.1) then
          cold=cnew
          isucc=1
        else 
          u(ii)=uold
          do i=1,n
            p(i)=pold(i)
          enddo
          isucc=0
        endif

        if (isucc.eq.1) write (6,*) iiuu,istep,cold,p(1)

      enddo ! MC loop

      call newtsolv (n,nmax,amat,b,u,f,p,q, &
                     g,deltag,alambda,iwork,info)

!     make a final test on the populations
      do i=1,n
        if (p(i).lt.0.0d0.or.p(i).gt.1.0d0) itest=0
      enddo

      if (itest.eq.1) then

        izero=0
        write (6,*)
        write (6,*) 'i, p, U '
        write (8,*)
        write (8,*) iiuu,p(1),uin
        do i=1,n
          write (6,*) i,p(i),u(i)
          write (8,*) i,p(i),u(i)
        enddo
        write (8,*) n+1,p(n),uout

        write (6,*)
        write (6,*) ' current check: '
        call fxn (n,nmax,u,p,q,g,deltag,alambda,f,0)
        do i=1,n
          write (6,*) i,f(i)
        enddo

        write (10,*) uin-uout,f(1)
        write (6,*) uin-uout,f(1)
        write (14,99) uin-uout,(p(i),i=1,n)
        write (16,99) uin-uout,(u(i),i=1,n)

        do i=1,n
          usave(i)=u(i)
          psave(i)=p(i)
        enddo

      endif

!     jump here at zero voltage
 200  continue

      enddo ! voltage loop

      close (8)
      close (10)
      close (14)
      close (16)

 99   format (5f12.7)

      stop
      end

! -------------------------------------------------------------

      subroutine newtsolv (n,nmax,amat,b,u,f,p,q,g, &
                           deltag,alambda,iwork,info)

      implicit real*8 (a-h,o-z)

      common uin,uout,gin,gout,dgin,dgout,alin,alout,tk,beta

      real*8 p(*),q(*)
      real*8 g(nmax,*),deltag(nmax,*),alambda(nmax,*)
      real*8 u(*),f(*)
      real*8 amat(nmax,*),b(*)
      integer iwork(*)
      logical done

      epsit=1.0d-10
      maxit=100

!     Newton search for p
      do iter=1,maxit
        done=.true.
        call fxn (n,nmax,u,p,q,g,deltag,alambda,f,1)
        do i=1,n
          b(i)=-f(i)
        enddo
        call jmatpn (n,nmax,amat,u,p,q,g,deltag,alambda,f)
        call dgefa (amat,nmax,n,iwork,info)
        if (info.ne.0) then
          write (6,*) info,' = info, should be zero '
          exit
        endif
        call dgesl (amat,nmax,n,iwork,b,0)
        errmax=0.0d0
        do i=1,n
          pold=p(i)
          p(i)=p(i)+b(i)
          if (p(i).ne.0.0d0) err=dabs((p(i)-pold)/p(i))
          if (err.gt.epsit) done=.false.
          errmax=dmax1(err,errmax)
        enddo
        if (done) exit
      enddo

      return
      end

! -----------------------------------------------------------
!
!     Newton's solution for a Kirchhoff-Master equation
!     compute the function f(x_n) = rhs of each iteration step.
!     Also used to compute the currents.
!
!     a - b - c + d = 0
!
!     the iflag means:
!
!     0 - compute all currents
!     1 - compute the functions 
!     2 - compute the first current
!
! -----------------------------------------------------------
 
      subroutine fxn (n,nmax,u,p,q,g,deltag,alambda,f,iflag)

      implicit real*8 (a-h,o-z)

      common uin,uout,gin,gout,dgin,dgout,alin,alout,tk,beta

      real*8 p(*),q(*),u(*),f(*)
      real*8 g(nmax,*),deltag(nmax,*),alambda(nmax,*)

!     compute holes from charges
      do i=1,n
        q(i)=1.0d0-p(i)
      enddo

!     left lead
      aright=alambda(1,2)
      ce=0.25*beta*(aright-u(1)+u(2)+deltag(1,2))**2 / aright
      de=0.25*beta*(aright+u(1)-u(2)+deltag(2,1))**2 / aright

      dw=2.0d0*dsqrt(alin*tk)
      a=q(1)*gin*dlead((alin+u(1)-uin+dgin)/dw)
      b=p(1)*gin*dlead((alin-u(1)+uin-dgin)/dw)
      c=p(1)*q(2)*g(1,2)*dexp(-ce)
      d=p(2)*q(1)*g(1,2)*dexp(-de)

      f(1)=a-b
      if (iflag.eq.2) return
      if (iflag.eq.1) f(1)=f(1)-c+d

      do i=2,n-1

!       exponents
        aleft=alambda(i-1,i)
        aright=alambda(i,i+1)
        ae=0.25*beta*(aleft+u(i)-u(i-1)+deltag(i-1,i))**2 / aleft
        be=0.25*beta*(aleft-u(i)+u(i-1)+deltag(i,i-1))**2 / aleft
        ce=0.25*beta*(aright-u(i)+u(i+1)+deltag(i,i+1))**2 / aright
        de=0.25*beta*(aright+u(i)-u(i+1)+deltag(i+1,i))**2 / aright

!       contributions
        a=p(i-1)*q(i)*g(i,i-1)*dexp(-ae)
        b=q(i-1)*p(i)*g(i,i-1)*dexp(-be)
        c=q(i+1)*p(i)*g(i,i+1)*dexp(-ce)
        d=p(i+1)*q(i)*g(i,i+1)*dexp(-de)

        f(i)=a-b
        if (iflag.eq.1) f(i)=f(i)-c+d

      enddo

      aleft=alambda(n-1,n)
      ae=0.25*beta*(aleft+u(n)-u(n-1)+deltag(n-1,n))**2 / aleft
      be=0.25*beta*(aleft-u(n)+u(n-1)+deltag(n,n-1))**2 / aleft

      a=p(n-1)*q(n)*g(n,n-1)*dexp(-ae)
      b=q(n-1)*p(n)*g(n,n-1)*dexp(-be)
      dw=2.0d0*dsqrt(alout*tk)
      c=p(n)*gout*dlead((alout+uout-u(n)+dgout)/dw)
      d=q(n)*gout*dlead((alout-uout+u(n)-dgout)/dw)

      f(n)=a-b
      if (iflag.eq.1) f(n)=f(n)-c+d

      return
      end

! -----------------------------------------------------------
!
!      Newton's solution for a Kirchhoff-Master equation.
!      Compute the Jacobi matrix w.r.t. p numerically
!
!      this is not really efficient yet, but it works
!
!      f = a - b - c + d 
!
! -----------------------------------------------------------

      subroutine jmatpn (n,nmax,amat,u,p,q,g,deltag,alambda,f)

      implicit real*8 (a-h,o-z)

      common uin,uout,gin,gout,dgin,dgout,alin,alout,tk,beta

      real*8 uin,uout,gin,gout,dgin,dgout,alin,alout,tk,beta
      real*8 a,b,c,d,ae,be,ce,de

      real*8 p(*),q(*),u(*),f(*)
      real*8 g(nmax,*),deltag(nmax,*),alambda(nmax,*)
      real*8 amat(nmax,*)

      integer i,j,n,nmax

      deltap=1.0d-5

      do i=1,n
        do j=1,n
          call fxn (n,nmax,u,p,q,g,deltag,alambda,f,1)
          fi=f(i)
          pold=p(j)
          p(j)=pold+deltap
          call fxn (n,nmax,u,p,q,g,deltag,alambda,f,1)
          amat(i,j)=(f(i)-fi)/deltap
          p(j)=pold
        enddo
      enddo

      return
      end

! --------------------------------------------------------------

!     ideal metallic leads (constant DOS), input and output 
!     conductivity, Fermi function = step function
     
      double precision function dlead(x)

      implicit real*8 (a-h,o-z)

!     exact
      dlead=0.5d0*erfc(x)

!     approximation by a single Gaussian
!     dlead=exp(-0.518*(x+1.18)**2)
!     if (x.lt.-1.18) dlead=1.0d0

      return
      end

! -----------------------------------------------------------

      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
!
!     dgefa factors a double precision matrix by gaussian elimination.
!
!     dgefa is usually called by dgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
!
!     on entry
!
!        a       double precision(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that dgesl or dgedi will divide by zero
!                     if called.  use  rcond  in dgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,dscal,idamax
!
!     internal variables
!
      double precision t
      integer idamax,j,k,kp1,l,nm1
!
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (a(l,k) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!
!           compute multipliers
!
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
!
!     dgesl solves the double precision system
!     a * x = b  or  trans(a) * x = b
!     using the factors computed by dgeco or dgefa.
!
!     on entry
!
!        a       double precision(lda, n)
!                the output from dgeco or dgefa.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!        ipvt    integer(n)
!                the pivot vector from dgeco or dgefa.
!
!        b       double precision(n)
!                the right hand side vector.
!
!        job     integer
!                = 0         to solve  a*x = b ,
!                = nonzero   to solve  trans(a)*x = b  where
!                            trans(a)  is the transpose.
!
!     on return
!
!        b       the solution vector  x .
!
!     error condition
!
!        a division by zero will occur if the input factor contains a
!        zero on the diagonal.  technically this indicates singularity
!        but it is often caused by improper arguments or improper
!        setting of lda .  it will not occur if the subroutines are
!        called correctly and if dgeco has set rcond .gt. 0.0
!        or dgefa has set info .eq. 0 .
!
!     to compute  inverse(a) * c  where  c  is a matrix
!     with  p  columns
!           call dgeco(a,lda,n,ipvt,rcond,z)
!           if (rcond is too small) go to ...
!           do 10 j = 1, p
!              call dgesl(a,lda,n,ipvt,c(1,j),0)
!        10 continue
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas daxpy,ddot
!
!     internal variables
!
      double precision ddot,t
      integer k,kb,l,nm1
!
      nm1 = n - 1
      if (job .ne. 0) go to 50
!
!        job = 0 , solve  a * x = b
!        first solve  l*y = b
!
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
!
!        now solve  u*x = y
!
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
!
!        job = nonzero, solve  trans(a) * x = b
!        first solve  trans(u)*y = b
!
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
!
!        now solve trans(l)*x = y
!
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end

      double precision function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + &
         dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

      subroutine  dscal(n,da,dx,incx)
!
!     scales a vector by a constant.
!     uses unrolled loops for increment equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision da,dx(*)
      integer i,incx,m,mp1,n,nincx
!
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
!
!        code for increment equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

      integer function idamax(n,dx,incx)
!
!     finds the index of element having max. absolute value.
!     jack dongarra, linpack, 3/11/78.
!     modified 3/93 to return if incx .le. 0.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dmax
      integer i,incx,ix,n
!
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
!
!        code for increment not equal to 1
!
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
!
!        code for increment equal to 1
!
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end
!
!     #############################################################
!     ##                                                         ##
!     ##  function ranf  --  portable random number generator    ##
!     ##                                                         ##
!     #############################################################
!
!     "random" generates a random number on [0,1] via a long
!     period generator due to L'Ecuyer with Bays-Durham shuffle
!
!     literature references:
!
!     P. L'Ecuyer, Communications of the ACM, 31, 742-774 (1988)
!
!     W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P.
!     Flannery, Numerical Recipes (Fortran), 2nd Ed., Cambridge
!     University Press, 1992, Section 7-1
!
      function ranf (seed)

      implicit none

      integer im1,ia1,iq1,ir1
      integer im2,ia2,iq2,ir2
      integer big,ntable
      integer imm1,ndiv

      real*8 factor

      parameter (im1=2147483563)
      parameter (ia1=40014)
      parameter (iq1=53668)
      parameter (ir1=12211)
      parameter (im2=2147483399)
      parameter (ia2=40692)
      parameter (iq2=52774)
      parameter (ir2=3791)
      parameter (big=141803398)
      parameter (ntable=32)
      parameter (imm1=im1-1)
      parameter (ndiv=1+imm1/ntable)
      parameter (factor=1.0d0/im1)
      integer i,k,next,seed,seed2
      integer iy,itable(ntable)
      real*8 ranf
      logical initial
      save initial,seed2,iy,itable
      data initial /.true./
!
!     warm up and then load the shuffling table
!
      if (initial) then
         initial=.false.
         seed2 = seed
         do i = ntable+8, 1, -1
            k = seed / iq1
            seed = ia1 * (seed-k*iq1) - k*ir1
            if (seed .lt. 0)  seed = seed + im1
            if (i .le. ntable)  itable(i) = seed
         end do
         iy = itable(1)
      end if
!
!     get a new random number value each call
!
      k = seed / iq1
      seed = ia1*(seed-k*iq1) - k*ir1
      if (seed .lt. 0)  seed = seed + im1
      k = seed2 / iq2
      seed2 = ia2*(seed2-k*iq2) - k*ir2
      if (seed2 .lt. 0)  seed2 = seed2 + im2
      i = 1 + iy/ndiv
      iy = itable(i) - seed2
      itable(i) = seed
      if (iy .lt. 1)  iy = iy + imm1
      ranf = factor * iy

      return
      end
