! ******************************************************************
! ******************************************************************

program algencanma

  use modamerica

  implicit none
  
  ! PARAMETERS
  integer, parameter :: ndiff = 3, ntrials = 2000, nvb = 86
  integer, parameter :: diffnballs(ndiff) = (/ 50, 75, 100 /)

  ! LOCAL SCALARS
  logical :: checkder,inA
  integer :: allocerr,answer,hnnzmax,i,inform,itrial,jcnnzmax,k,m,n,nvparam
  real(kind=8) :: bestrad,cnorm,efacc,efstain,eoacc,eostain,epsfeas,epsopt,f,nlpsupn,seed,snorm
  
  ! LOCAL ARRAYS
  character(len=80) :: specfnm,outputfnm,vparam(10),solfilename
  logical :: coded(11)
  logical, pointer :: equatn(:),linear(:)
  integer :: vb(nvb),cuba(4),tf(4)
  real(kind=8), pointer :: l(:),lambda(:),u(:),x(:)
  real(kind=8) :: p(2)

  ! COMMON SCALARS
  integer :: nballs
  real(kind=8) :: h

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdata/ dbl,dtr,h,nballs

  ! FUNCTIONS
  real(kind=8) :: drand

  ! EXTERNAL SUBROUTINES
  external :: myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
              myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp

  h = 1.0d-03

  write(*,*) 'xxxxxxx'
  write(*,*) 'h = ',h
  write(*,*) 'xxxxxxx'

  ! Constraints

  m = 1

  allocate(equatn(m),linear(m),lambda(m),stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Allocation error in main program'
     stop
  end if
     
  equatn(1:m) = .true.
  linear(1:m) = .false.
     
  ! Coded subroutines

  coded(1:11) = .false.
  coded(1:2) = .true. ! fsub, gsub
  coded(4:5) = .true. ! csub, jacsub
     
  ! Checking derivatives?

  checkder = .false.

  ! Parameters setting
        
  epsfeas   = 0.1d0 * h
  epsopt    = 1.0d-01
        
  efstain   = 1.0d-08
  eostain   = 1.0d-08
        
  efacc     = 1.0d-08
  eoacc     = 1.0d-08
  
  outputfnm = ''
  specfnm   = ''
        
  nvparam = 2
  vparam(1) = 'OUTER-ITERATIONS-LIMIT 20'
  vparam(2) = 'LARGEST-PENALTY-PARAMETER-ALLOWED 1.0d+08'
  
  ! Problem data
  
  put(1:2, 1:np) = put(1:2, 1:np) / 20.0d0

  dbl(1) = minval(put(1,1:np))
  dbl(2) = minval(put(2,1:np))
  dtr(1) = maxval(put(1,1:np))
  dtr(2) = maxval(put(2,1:np))
  
  h = 1.0d-02
  
  ! Indices of the vertices (outer boundary)
  
  vb(1:nvb) = (/ 43, 80, 44, 78, 124, 45, 79, 46, 83, 35, 120, 56, 121, 33, 118, 104, 116,  &
                119, 105, 106, 117, 107, 86, 113, 114, 31, 92, 93, 30, 28, 27, 128, 59,     &
                130, 57, 58, 2, 69, 67, 66, 4, 3, 5, 68, 76, 48, 129, 47, 125, 17, 19, 122, &
                123, 24, 29, 91, 90, 32, 87, 115, 88, 103, 89, 85, 34, 94, 62, 98, 61, 60,  &
                96, 97, 36, 99, 100, 37, 38, 39, 101, 102, 40, 41, 42, 127, 81, 82 /)
                
  cuba(1:4) = (/ 49, 50, 51, 52 /)
  
  tf(1:4) = (/ 73, 77, 71, 72 /)

  do k = 1,ndiff

     nballs = diffnballs(k)
     
     ! Number of variables (x = (c_1,...,c_nballs,r))

     n = 2 * nballs + 1

     ! Set lower bounds, upper bounds, and initial guess

     allocate(x(n),l(n),u(n),stat=allocerr)

     if ( allocerr .ne. 0 ) then
        write(*,*) 'Allocation error in main program'
        stop
     end if

     do i = 1,nballs
        l(2*i-1) = dbl(1)
        l(2*i)   = dbl(2)
        u(2*i-1) = dtr(1)
        u(2*i)   = dtr(2)
     end do
     
     l(2*nballs+1) = 0.0d0
     u(2*nballs+1) = 1.0d+20

     ! Upper bounds on the number of sparse-matrices non-null elements
     
     jcnnzmax = n
     hnnzmax  = 0
        
     bestrad = 1.0d+20
     
     do itrial = 1,ntrials
     
        seed = 654321.0d0 * itrial

        do i = 1,nballs
10         continue
           p(1) = dbl(1) + ( dtr(1) - dbl(1) ) * drand(seed)
           p(2) = dbl(2) + ( dtr(2) - dbl(2) ) * drand(seed)
           
           call polygon_contains_point_2d_3 ( nvb, put(1:2, vb(1:nvb)), p, inA )
           
           if ( .not. inA ) then
              call polygon_contains_point_2d_convex( 4, put(1:2, cuba(1:4)), p, inA )
           end if
           
           if ( .not. inA ) then
              call polygon_contains_point_2d_convex( 4, put(1:2, tf(1:4)), p, inA )
           end if
           
           if ( .not. inA ) then
              go to 10
           end if
           x(2*i-1) = p(1)
           x(2*i)   = p(2)
        end do

        x(2*nballs+1) = ( 0.5d0 + 1.0d0 * drand(seed) ) / nballs

        lambda(1:m) = 0.0d0
     
        ! Optimize
        
        call algencan(myevalf,myevalg,myevalh,myevalc,myevaljac,myevalhc, &
             myevalfc,myevalgjac,myevalgjacp,myevalhl,myevalhlp,jcnnzmax, &
             hnnzmax,epsfeas,epsopt,efstain,eostain,efacc,eoacc,outputfnm, &
             specfnm,nvparam,vparam,n,x,l,u,m,lambda,equatn,linear,coded, &
             checkder,f,cnorm,snorm,nlpsupn,inform)

        write(*,*) 'inform = ',inform,' feasible = ',(cnorm .le. epsfeas),' cnorm = ',cnorm,' radius = ',x(2*nballs+1)

        if ( inform .eq. 0 .and. cnorm .le. epsfeas .and. x(2*nballs+1) .lt. bestrad ) then
           write(*,*) 'A new solution was found at trial ',itrial,' with r = ',x(2*nballs+1)
           bestrad = x(2*nballs+1)
           call drawsol(n,x)
        else
           write(*,*) 'The output of trial ',itrial,' did not improve the best known solution.'
        end if
     end do
  
     deallocate(x,l,u,stat=allocerr)

     if ( allocerr .ne. 0 ) then
        write(*,*) 'Deallocation error in main program'
        stop
     end if

  end do

  deallocate(lambda,equatn,linear,stat=allocerr)

  if ( allocerr .ne. 0 ) then
     write(*,*) 'Deallocation error in main program'
     stop
  end if

  stop

end program algencanma

! ******************************************************************
! ******************************************************************

subroutine myevalf(n,x,f,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  integer :: nballs
  real(kind=8) :: h

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdata/ dbl,dtr,h,nballs

  flag = 0

  f = x(2*nballs+1)

end subroutine myevalf

! ******************************************************************
! ******************************************************************

subroutine myevalg(n,x,g,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n
  integer, intent(out) :: flag

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n)

  ! COMMON SCALARS
  integer :: nballs
  real(kind=8) :: h

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdata/ dbl,dtr,h,nballs

  flag = 0

  g(1:2*nballs) = 0.0d0
  g(2*nballs+1) = 1.0d0

end subroutine myevalg

! ******************************************************************
! ******************************************************************

subroutine myevalh(n,x,hrow,hcol,hval,hnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,n
  integer, intent(out) :: flag,hnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hcol(lim),hrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hval(lim)

  flag = - 1

end subroutine myevalh

! ******************************************************************
! ******************************************************************

subroutine myevalc(n,x,ind,c,flag)

  use modamerica

  implicit none
  
  ! PARAMETERS
  integer, parameter :: nvb = 86

  ! SCALAR ARGUMENTS
  integer, intent(in)  :: ind,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: c

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in)  :: x(n)

  ! COMMON SCALARS
  integer :: nballs
  real(kind=8) :: h

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdata/ dbl,dtr,h,nballs

  ! LOCAL SCALARS
  logical :: inA,insomeball
  integer :: count,countinA,i,j,k,nx,ny
  real(kind=8) :: hx,hy

  ! LOCAL ARRAYS
  integer :: vb(nvb),cuba(4),tf(4)
  real(kind=8) :: p(2)
  
  ! Indices of the vertices (outer boundary)
  
  vb(1:nvb) = (/ 43, 80, 44, 78, 124, 45, 79, 46, 83, 35, 120, 56, 121, 33, 118, 104, 116,  &
                119, 105, 106, 117, 107, 86, 113, 114, 31, 92, 93, 30, 28, 27, 128, 59,     &
                130, 57, 58, 2, 69, 67, 66, 4, 3, 5, 68, 76, 48, 129, 47, 125, 17, 19, 122, &
                123, 24, 29, 91, 90, 32, 87, 115, 88, 103, 89, 85, 34, 94, 62, 98, 61, 60,  &
                96, 97, 36, 99, 100, 37, 38, 39, 101, 102, 40, 41, 42, 127, 81, 82 /)
                
  cuba(1:4) = (/ 49, 50, 51, 52 /)
  
  tf(1:4) = (/ 73, 77, 71, 72 /)
  
  flag = 0

  if ( ind .ne. 1 ) then
     write(*,*) 'Wrong ind in evalc.'
     flag = - 1
     return
  end if
  
  nx = ceiling( ( dtr(1) - dbl(1) ) / h )
  ny = ceiling( ( dtr(2) - dbl(2) ) / h )
    
  hx = ( dtr(1) - dbl(1) ) / nx
  hy = ( dtr(2) - dbl(2) ) / ny
    
  count = 0
  countinA = 0
  do i = 1,nx
     p(1) = dbl(1) + ( i - 0.5d0 ) * hx
     do j = 1,ny
        p(2) = dbl(2) + ( j - 0.5d0 ) * hy

        call polygon_contains_point_2d_3 ( nvb, put(1:2, vb(1:nvb)), p, inA )
           
        if ( .not. inA ) then
           call polygon_contains_point_2d_convex( 4, put(1:2, cuba(1:4)), p, inA )
        end if
           
        if ( .not. inA ) then
           call polygon_contains_point_2d_convex( 4, put(1:2, tf(1:4)), p, inA )
        end if

        if ( inA ) then
           k = 0
           insomeball = .false.
           do while ( .not. insomeball .and. k .lt. nballs )
              k = k + 1
              if ( ( x(2*k-1) - p(1) ) ** 2 + ( x(2*k) - p(2) ) ** 2 .le. x(2*nballs+1) ** 2 ) then
                 insomeball = .true.
              end if
           end do
           if ( insomeball ) then
              count = count + 1
           end if
           countinA = countinA + 1
        end if
     end do
  end do
	
  c = hx * hy * ( countinA - count )

end subroutine myevalc

! ******************************************************************
! ******************************************************************

subroutine myevaljac(n,x,ind,jcvar,jcval,jcnnz,lim,lmem,flag)

  use modamerica

  implicit none
  
  ! PARAMETERS
  integer, parameter :: nvb = 86

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcvar(lim)
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: jcval(lim)

  ! COMMON SCALARS
  integer :: nballs
  real(kind=8) :: h

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdata/ dbl,dtr,h,nballs

  ! LOCAL SCALARS
  logical :: inA,insomeball
  integer :: i,j,k,l,nclose,ntheta
  real(kind=8) :: htheta,theta

  ! LOCAL ARRAYS
  integer :: close(nballs),vb(nvb),cuba(4),tf(4)
  real(kind=8) :: p(2)
  
  ! Constant pi
  real(kind=8), parameter :: PI = acos( - 1.0d0 )
  
  ! Indices of the vertices (outer boundary)
  
  vb(1:nvb) = (/ 43, 80, 44, 78, 124, 45, 79, 46, 83, 35, 120, 56, 121, 33, 118, 104, 116,  &
                119, 105, 106, 117, 107, 86, 113, 114, 31, 92, 93, 30, 28, 27, 128, 59,     &
                130, 57, 58, 2, 69, 67, 66, 4, 3, 5, 68, 76, 48, 129, 47, 125, 17, 19, 122, &
                123, 24, 29, 91, 90, 32, 87, 115, 88, 103, 89, 85, 34, 94, 62, 98, 61, 60,  &
                96, 97, 36, 99, 100, 37, 38, 39, 101, 102, 40, 41, 42, 127, 81, 82 /)
                
  cuba(1:4) = (/ 49, 50, 51, 52 /)
  
  tf(1:4) = (/ 73, 77, 71, 72 /)

  flag = 0
  lmem = .false.
  
  if ( ind .ne. 1 ) then
     write(*,*) 'Wrong ind in evaljac.'
     flag = - 1
     return
  end if

  if ( lim .lt. n ) then
     write(*,*) 'lim smaller than n in evaljac.'
     flag = - 1
     lmem = .true.
     return
  end if
     
  jcnnz = n
  jcvar(1:n) = (/ (i,i=1,n) /)
  jcval(1:n) = 0.0d0
  
  if ( x(2*nballs+1) .eq. 0.0d0 ) then
     return
  end if
  
  ntheta = ceiling( ( 2.0d0 * PI * x(2*nballs+1) ) / h )
  htheta = 2.0d0 * PI * x(2*nballs+1) / ntheta

  do i = 1,nballs
     nclose = 0
     do j = 1,nballs
        if ( j .ne. i ) then
           if ( ( x(2*j-1) - x(2*i-1) ) ** 2 + ( x(2*j) - x(2*i) ) ** 2 .le. 4.0d0 * x(2*nballs+1) ** 2 ) then
              nclose = nclose + 1
              close(nclose) = j
           end if
        end if
     end do
     
     do k = 1,ntheta
        theta = ( k - 0.5d0 ) * htheta / x(2*nballs+1)
        p(1) = x(2*i-1) + x(2*nballs+1) * cos(theta)
        p(2) = x(2*i)   + x(2*nballs+1) * sin(theta)
        
        call polygon_contains_point_2d_3 ( nvb, put(1:2, vb(1:nvb)), p, inA )
           
        if ( .not. inA ) then
           call polygon_contains_point_2d_convex( 4, put(1:2, cuba(1:4)), p, inA )
        end if
           
        if ( .not. inA ) then
           call polygon_contains_point_2d_convex( 4, put(1:2, tf(1:4)), p, inA )
        end if

        if ( inA ) then
           l = 0
           insomeball = .false.
           do while ( .not. insomeball .and. l .lt. nclose )
              l = l + 1
              j = close(l)
              if ( ( x(2*j-1) - p(1) ) ** 2 + ( x(2*j) - p(2) ) ** 2 .lt. x(2*nballs+1) ** 2 ) then
                 insomeball = .true.
              end if
           end do

           if ( .not. insomeball ) then
              jcval(2*nballs+1) = jcval(2*nballs+1) + 1.0d0
              jcval(2*i-1)      = jcval(2*i-1) + cos(theta)
              jcval(2*i)        = jcval(2*i)   + sin(theta)
           end if
        end if
     end do
  end do

  jcval(1:n) = - htheta * jcval(1:n)
     
end subroutine myevaljac

! ******************************************************************
! ******************************************************************

subroutine myevalhc(n,x,ind,hcrow,hccol,hcval,hcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: ind,lim,n
  integer, intent(out) :: flag,hcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hccol(lim),hcrow(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: hcval(lim)

  flag = - 1

end subroutine myevalhc

! ******************************************************************
! ******************************************************************

subroutine myevalfc(n,x,f,m,c,flag)

  implicit none

  ! SCALAR ARGUMENTS
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  real(kind=8), intent(out) :: f

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: c(m)

  flag = - 1

end subroutine myevalfc

! ******************************************************************
! ******************************************************************

subroutine myevalgjac(n,x,g,m,jcfun,jcvar,jcval,jcnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,m,n
  integer, intent(out) :: flag,jcnnz

  ! ARRAY ARGUMENTS
  integer, intent(out) :: jcfun(lim),jcvar(lim)
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: g(n),jcval(lim)

  flag = - 1

end subroutine myevalgjac

! ******************************************************************
! ******************************************************************

subroutine myevalgjacp(n,x,g,m,p,q,work,gotj,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: gotj
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  character, intent(in) :: work

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(inout) :: p(m),q(n)
  real(kind=8), intent(out) :: g(n)

  flag = - 1

end subroutine myevalgjacp

! ******************************************************************
! ******************************************************************

subroutine myevalhl(n,x,m,lambda,sf,sc,hlrow,hlcol,hlval,hlnnz,lim,lmem,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(out) :: lmem
  integer, intent(in) :: lim,m,n
  integer, intent(out) :: flag,hlnnz
  real(kind=8), intent(in) :: sf

  ! ARRAY ARGUMENTS
  integer, intent(out) :: hlcol(lim),hlrow(lim)
  real(kind=8), intent(in) :: lambda(m),sc(m),x(n)
  real(kind=8), intent(out) :: hlval(lim)

  flag = - 1

end subroutine myevalhl

! ******************************************************************
! ******************************************************************

subroutine myevalhlp(n,x,m,lambda,sf,sc,p,hp,goth,flag)

  implicit none

  ! SCALAR ARGUMENTS
  logical, intent(inout) :: goth
  integer, intent(in) :: m,n
  integer, intent(out) :: flag
  real(kind=8), intent(in) :: sf

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: lambda(m),p(n),sc(m),x(n)
  real(kind=8), intent(out) :: hp(n)

  flag = - 1

end subroutine myevalhlp

! ******************************************************************
! ******************************************************************

subroutine drawsol(n,x)

  use modamerica

  implicit none
  
  ! PARAMETER
  integer, parameter :: nvb = 86

  ! SCALAR ARGUMENTS
  integer, intent(in) :: n

  ! ARRAY ARGUMENTS
  real(kind=8), intent(in) :: x(n)

  ! COMMON SCALARS
  integer :: nballs
  real(kind=8) :: h

  ! COMMON ARRAYS
  real(kind=8) :: dbl(2),dtr(2)

  ! COMMON BLOCKS
  common /pdata/ dbl,dtr,h,nballs

  ! LOCAL SCALARS
  integer :: i,k

  ! LOCAL ARRAYS
  integer :: vb(nvb),cuba(4),tf(4)
  character(len=80) :: filename
  
  ! Indices of the vertices (outer boundary)
  
  vb(1:nvb) = (/ 43, 80, 44, 78, 124, 45, 79, 46, 83, 35, 120, 56, 121, 33, 118, 104, 116,  &
                119, 105, 106, 117, 107, 86, 113, 114, 31, 92, 93, 30, 28, 27, 128, 59,     &
                130, 57, 58, 2, 69, 67, 66, 4, 3, 5, 68, 76, 48, 129, 47, 125, 17, 19, 122, &
                123, 24, 29, 91, 90, 32, 87, 115, 88, 103, 89, 85, 34, 94, 62, 98, 61, 60,  &
                96, 97, 36, 99, 100, 37, 38, 39, 101, 102, 40, 41, 42, 127, 81, 82 /)
                
  cuba(1:4) = (/ 49, 50, 51, 52 /)
  
  tf(1:4) = (/ 73, 77, 71, 72 /)

  write(filename,"(A17,I0,A3)") 'covering-america-',nballs,'.mp'
  
  open(unit=10,file=trim(filename))

  ! BEGINING
  write(10,10) nballs
  
  write(10,*) ' path p;'
  write(10,*) ' path s;'
  write(10,*) ' path t;'
  
  ! BALLS FILL
  do i = 1,nballs
     write(10,29) 2.0d0*x(2*nballs+1),x(2*i-1),x(2*i)
  end do

  ! REGION A
  write(10,*) ' p := '
  do k = 1,nvb
    write(10,20) put(1,vb(k)),put(2,vb(k))
  end do
  write(10,25)
  
  write(10,*) ' s := '
  do k = 1,4
    write(10,20) put(1,cuba(k)),put(2,cuba(k))
  end do
  write(10,26)
  
  write(10,*) ' t := '
  do k = 1,4
    write(10,20) put(1,tf(k)),put(2,tf(k))
  end do
  write(10,27)
  
  
  ! BALLS DRAW
  do i = 1,nballs
     write(10,30) 2.0d0*x(2*nballs+1),x(2*i-1),x(2*i)
  end do

  write(10,40) x(2*nballs+1)

  close(10)

  ! NON-EXECUTABLE STATEMENTS

10 format('prologues := 3;',/, &
          'outputtemplate := "%j-%c.mps";',/, &
          'input mpcolornames;',/, &
          'beginfig(',i3,');',/, &
          'u = 8.0 cm;')
20 format(' (',f10.3,'u,',f10.3,'u)--')
25 format(' cycle; fill p withcolor 0.9Dandelion; draw p withpen pencircle scaled 0.7;')
26 format(' cycle; fill s withcolor 0.9Dandelion; draw s withpen pencircle scaled 0.7;')
27 format(' cycle; fill t withcolor 0.9Dandelion; draw t withpen pencircle scaled 0.7;')
29 format('fill fullcircle scaled ',f20.10,'u shifted (',f20.10,'u,',f20.10,'u) withcolor 0.9white;')
30 format('draw fullcircle scaled ',f20.10,'u shifted (',f20.10,'u,',f20.10,'u) withcolor RoyalBlue withpen pencircle scaled 0.7;')
40 format('% r = ',f20.10,/,'endfig;',/,'end;') 

end subroutine drawsol
