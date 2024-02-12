program svdc
    implicit none
!-------------------------Insert dimension of A(m x n)-------------------------------!
    integer, parameter :: m = 20, n = 20
    integer :: i, j, rs
!------------------------------------------------------------------------------------!
    real :: a(m,n), u(m,m), v(n,n), s(m,n), usvt(m,n)
!-------------------------------Insert A---------------------------------------------!
    do i = 1,m
        do j = 1,n
            a(i,j) = (i-1)*m+j
        end do
    end do
!------------------------------------------------------------------------------------!
    call svd(a, m, n, u, v, s, rs)
    print *, "a"
    print *, a
    print *, "u"
    print *, u
    print *, "v"
    print *, v
    print *, "s"
    print *, s(1,1)
    print *, "U*S*Vt"
    call check(u,s,v,usvt,m,n)
    print *, usvt
    if(rs==1) then
        print *, "svd was successfully computed"
    else
        print *, "svd was not successfully computed"
        print *, "big matrix occur big error"
        print *, "please adjust the e value higher(line 273)"
    end if
!------------------------------------------------------------------------------------!
end program

subroutine svd(a,m,n,u,v,s,rs)
!   compute SVD
!------------------------------------------------------------------------------------!
! a : matrix for svd (m x n), u : left-singular vector v : right-singular vector
! s : singular value, rs : success or not
!------------------------------------------------------------------------------------!
    real, parameter :: e=0.00001
    real :: a(m,n), at(n,m), aat(m,m), ata(n,n), vc(n,n), uc(m,m)
    real :: u(m,m), v(n,n), s1(m,m), s2(n,n), usvt(m,n), vt(n,n), s(m,n)
    integer :: m, n, i, j, rs
!------------------------------------------------------------------------------------!
    ! Compute At
    call matt(a,at,m,n)
    ! Compute A*At
    call matm(a,at,aat,m,n,m)
    ! Compute At*A
    call matm(at,a,ata,n,m,n)
    ! Compute A`
    call eig(aat,u,s1,m)
    call eig(ata,v,s2,n)
    ! Compute sigma
    s = 0.0
    !! sort sigma
    call sort(u,s1,m)
    call sort(v,s2,n)
    !! arrange sigma
    if (m>n) then
        do i = 1,n
            s(i,i) = s2(i,i)
        end do
    else
        do i = 1,m
            s(i,i) = s1(i,i)
        end do
    end if
    !! find eigenvector
    call check(u,s,v,usvt,m,n)
    call ex(usvt,a,rs,m,n)
    do i=1,m
        call ex(usvt,a,rs,m,n)
        do j=1,n
            call ex(usvt,a,rs,m,n)
            if(rs == 1) then
                exit
            else
                vc = v
                call negc(vc,n,n,j)
                call check(u,s,vc,usvt,m,n)
                call ex(usvt,a,rs,m,n)
                if(rs == 1) then
                    v = vc
                    exit
                end if
            end if
        end do
        if(rs == 1) then
            exit
        else
            uc = u
            call negc(uc,m,m,i)
            call check(uc,s,v,usvt,m,n)
            call ex(usvt,a,rs,m,n)
            if(rs == 1) then
            u = uc
            exit
            end if
        end if
    end do
end subroutine

subroutine matt(a,at,m,n)
! find transpose of a matrix
!------------------------------------------------------------------------------------!
! a : matrix for transpose (m x n), at : transpose matrix of a (n x m)
!------------------------------------------------------------------------------------!
    real :: a(m,n), at(n,m)
    integer :: m, n
    do i = 1,m
      do j = 1,n
        at(j,i) = a(i,j)
      end do
    end do
end subroutine

subroutine matm(a,b,ab,m,n,l)
! compute matrix multiplication
!------------------------------------------------------------------------------------!
! a : matrix for matrix multiplication (left side) [m x n]
! b : matrix for matrix multiplication (right side) [n x l]
! ab : matrix for matrix multiplication of a and b [m x l]
!------------------------------------------------------------------------------------!
    real :: a(m,n), b(n,l), ab(m,l)
    integer :: m, n, l
    do i = 1,m
      do j = 1,l
        ab(i,j)=0
        do k = 1,n
          ab(i,j) = ab(i,j) + a(i,k)*b(k,j)
        end do
      end do
    end do
end subroutine

subroutine check(u,s,v,usvt,m,n)
! compute A (=U*S*Vt)
!------------------------------------------------------------------------------------!
! u : left singular vector (m x m) v : right singular vector (n x n)
! s : singular value (m x n) usvt : U*S*Vt
!------------------------------------------------------------------------------------!
    real :: u(m,m), v(n,n), vt(n,n), s(m,n), us(m,n), usvt(m,n)
    integer :: m, n
    call matt(v,vt,n,n)
    call matm(u,s,us,m,m,n)
    call matm(us,vt,usvt,m,n,n)
end subroutine

subroutine matcc(u,m,n,p,q)
! change each matrix column p, q
!------------------------------------------------------------------------------------!
! u : matrix for changing (m x n)
! p, q : number of column to change
!------------------------------------------------------------------------------------!
    real :: u(m,n),c(m,1)
    integer :: m, n, p, q

    do i =1,m
        c(i,1) = u(i,p)
    end do
    do i =1,m
        u(i,p) = u(i,q)
        u(i,q) = c(i,1)
    end do
end subroutine

subroutine eig(a,u,s,m)
! compute eigenvector by jacobi method
!------------------------------------------------------------------------------------!
! u : matrix for changing (m x n)
! p, q : number of column to change
!------------------------------------------------------------------------------------!
    real, parameter :: e=0.0000001, pi=3.14159265359
    real :: a(m,m), u(m,m), s(m,m), ac(m,m),r(m,m),rt(m,m),rta(m,m),uc(m,m)
    integer :: m, p, q, i, j, k, fa
    u = 0.0
    do i=1,m
        u(i,i) = 1.0
    end do
    ac = a
    do k = 1,10
        call ij(ac,fa,m,m)
        if (fa==1) then

            exit
        else
        do p = 1,m-1
            do q = p+1,m
                if (ac(p,p)==ac(q,q)) then
                    t = pi/4.
                else
                    t = 0.5 * atan(2.*ac(p,q)/(ac(p,p)-ac(q,q)))
                end if
                r = 0
                do i =1,m
                    r(i,i)=1
                end do
                r(p,p) = cos(t)
                r(p,q) = -sin(t)
                r(q,p) = sin(t)
                r(q,q) = cos(t)
                !print *, r
                call matt(r,rt,m,m)
                call matm(rt,ac,rta,m,m,m)
                call matm(rta,r,ac,m,m,m)
                call matm(u,r,uc,m,m,m)
                u = uc
                !print *, ac
                do i = 1,m
                    do j = 1,m
                        if(ac(i,j)**2<e) then
                            ac(i,j) = 0.0
                        end if
                    end do
                end do
            end do
        end do
        end if
    end do
    do i=1,m
        do j=1,m
            if(ac(i,j)<0) then
            ac(i,j)=0
            end if
        end do
    end do
    s = ac**(0.5)
end subroutine

subroutine negc(a,m,n,c)
! change the column c to negative value of column c
!------------------------------------------------------------------------------------!
! a : matrix for changing (m x n)
! c : number of column to change
!------------------------------------------------------------------------------------!
    real :: a(m,n)
    integer :: m, n, c, i
    do i=1,m
        a(i,c) = (-1.)*a(i,c)
    end do
end subroutine

subroutine sort(u,s,m)
! sort sigma (s1>s2>s3....>ss) with singular vector
!------------------------------------------------------------------------------------!
! u : singular vector (left or right) to sort with singular value (m x m)
! s : singular value (m x m)
!------------------------------------------------------------------------------------!
    real :: s(m,m), big, u(m,m)
    integer i, j, n, m
        do j = 1, m-1
            do i = 1, m-1
                if(s(i,i)<s(i+1,i+1)) then
                    big = s(i+1,i+1)
                    s(i+1,i+1) = s(i,i)
                    s(i,i) = big
                    call matcc(u,m,m,i,i+1)
                end if
            end do
        end do
end subroutine

subroutine ex(usvt,a,rs,m,n)
! verify U*S*Vt == original A
!------------------------------------------------------------------------------------!
! usvt : U*S*Vt (m x n), a : original vector (m x n)
! rs : switch of verifying (1:on, 0:off)
!------------------------------------------------------------------------------------!
    real, parameter :: e=0.01
    real :: usvt(m,n), a(m,n), s
    integer :: m, n, rs, i, n
    rs = 0
    s = 0
    do i= 1, m
        do j = 1, n
            s = s + (usvt(i,j)-a(i,j))**2
        end do
    end do
    if(abs(s)<(e*((m*n)**2))) then
        rs = 1
    end if
end subroutine

subroutine ij(a,fa,m,n)
! verify Diagonalization
!------------------------------------------------------------------------------------!
! a : vector for verifying (m x n)
! fa : switch of verifying (1:on, 0:off)
!------------------------------------------------------------------------------------!
    real :: a(m,n),s
    real, parameter :: e=0.0000001
    integer :: i,j,m,n, fa
    s = 0
    fa = 0
    do i = 1,m
        do j = 1,n
            if(i /= j) then
            s = s + a(i,j)**2
            end if
        end do
    end do
    if (s<e) then
        fa = 1
    end if
end subroutine
