module parameter
    implicit none
    complex, parameter :: z = cmplx(0.0,1.0)
    real*8, parameter :: SR3 = sqrt(3.d0)
    !lattice constant
    integer, parameter :: comid = 1 !1: MoS2
    integer, parameter :: halfNorbital = 3
    integer, parameter :: Norbital = 6
    integer, parameter :: nkamax = 60
    integer :: nkqcutoff
    !physics constant
    real*8, parameter :: pi = dacos(-1.d0)
    real*8, parameter :: pix2 = 2.d0*pi
    real*8, parameter :: kB = 8.617333262145e-5 !eV/K
    integer, parameter :: Temperate = 300 !K
    real*8, parameter :: Thermalbeta = 1.d0/(kB*Temperate)
    complex*8, parameter :: varepsilon = cmplx(2.5d0,0.)
    real*8, parameter :: varepsilon0 = 8.8541878128d-12/(1.6022d-19)*1.d-9 !e^2/(eV.nm)
    real*8, parameter :: qe = -1.d0 
    real*8, parameter :: hb = 0.658229d0 ![eV.fs]
    real*8, parameter :: me = 9.1094d0/1.6022d0  ! eV.fs^2/nm^2; E=mc^2;kg=J(s/m)^2=1/qe eV (s/m)^2
    real*8, parameter :: Calpha = qe**2/(2*varepsilon0*varepsilon)
    real*8, parameter :: kcutoff = 3.d0 + 1.d-7 ![1/nm]
    !Elecfield constant V/cm -> V/nm *1e7*1e-3 = 1e-4
    real*8, parameter :: E0 = 1.d6*1d-7![V/nm] ~ 1e7 [V/cm] !extreme value of electric field to calculate the exciton
    !real, parameter :: epsilonL = 1.6 !ev
    real*8, parameter :: td = 5.d0 !fsF
    !real, parameter :: w0 = epsilonL/hb
    real*8 :: w0, Egap
    real*8, parameter :: ell = 1.d0, CEP = 0.
    !real*8, parameter :: tmin = - 100.0, tmax = 500.d0, dt = 1.d0/50.d0 !fs
    real*8, parameter :: tmin = - 4.d0*(td/sqrt(2.*log(2.d0))), tmax = 1000.d0, dt = 1.d0/50.d0 !fs
    integer, parameter :: ntmax = int((tmax-tmin)/dt)
    real*8, parameter :: detuning = 10.d-3 !eV, 1e3 meV
    !approximation check
    logical, parameter :: usecoulomb = .true. !true: use Coulomb interaction, false: no Coulomb interaction
    logical, parameter :: useshieldcoulomb = .false. !true: use shielded Coulomb interaction, false: bare Coulomb interaction 
    logical, parameter :: usekcutoff = .true.
    logical, parameter :: usedephasing = .true. !T2 approximation
    logical, parameter :: usegga = .false. !true: use GGA, false: use LDA
    logical, parameter :: usegammaeecor = .false.
    logical, parameter :: usecircularlight = .true.
    integer, parameter :: userelax = 0 !T1 approximation
    !approximation parameter
    real*8, parameter :: Tco2 = 20.d0 !Minh thesis use 15.d0
    real, parameter :: T1depha = 1.
    real, parameter :: gamma = 65 !e-e correlation in T2 nm^3/fs
    !
    real :: nepsilon = 10000.

    double precision :: hermitianerror = nkamax*nkamax*1d-15
end module parameter
module para_MX2
    use parameter
    double precision :: a,e1,e2,t0,t1,t2,t11,t12,t22,r0,r1,r2,r11,r12,&
                        u0,u1,u2,u11,u12,u22,lambda
    save
 contains
    subroutine tb_parameters(material, use)
       implicit none
       integer :: material
       logical :: use
       if (use) then
           call tb_parametersgga(material)
       else
           call tb_parameterslda(material)
       end if
    end subroutine
    subroutine tb_parameterslda(material)
        implicit none
        integer :: material
        double precision, dimension(6) :: &
             atab =(/0.3129,0.3132,0.3254,0.3253,0.3472,0.3476/),&
             e1tab =(/0.820,0.905,0.715,0.860,0.574,0.675/),&
             e2tab =(/1.931,2.167,1.687,1.892,1.410,1.489/),&
             t0tab =(/-0.176,-0.175,-0.154,-0.152,-0.148,-0.124/),&
             t1tab =(/-0.101,-0.090,-0.134,-0.125,-0.173,-0.159/),&
             t2tab =(/0.531,0.611,0.437,0.508,0.333,0.362/),&
             t11tab=(/0.084,0.043,0.124,0.094,0.203,0.196/),&
             t12tab=(/0.169,0.181,0.119,0.129,0.186,0.101/),&
             t22tab=(/0.070,0.008,0.072,0.009,0.127,0.044/),&
             r0tab =(/0.070,0.075,0.048,0.044,0.007,-0.009/),&
             r1tab =(/-0.252,-0.282,-0.248,-0.278,-0.280,-0.25/),&
             r2tab =(/0.084,0.127,0.090,0.129,-0.067,0.129/),&
             r11tab=(/0.019,0.001,0.066,0.059,0.073,0.131/),&
             r12tab=(/0.093,0.114,0.045,0.058,0.081,-0.007/),&
             u0tab =(/-0.043,-0.063,-0.067,-0.090,-0.054,-0.086/),&
             u1tab =(/0.047,0.047,0.041,0.039,0.008,0.012/),&
             u2tab =(/0.005,0.004,0.005,0.001,0.037,-0.020/),&
             u11tab =(/0.304,0.374,0.327,0.392,0.145,0.361/),&
             u12tab =(/-0.192,-0.224,-0.194,-0.224,-0.078,-0.193/),&
             u22tab =(/-0.162,-0.177,-0.151,-0.165,0.035,-0.129/),&
             lambdatab=(/0.073,0.211,0.091,0.228,0.107,0.237/)
        a=atab(material)
        e1=e1tab(material)
        e2=e2tab(material)
        t0=t0tab(material)
        t1=t1tab(material)
        t2=t2tab(material)
        t11=t11tab(material)
        t12=t12tab(material)
        t22=t22tab(material)
        r0=r0tab(material)
        r1=r1tab(material)
        r2=r2tab(material)
        r11=r11tab(material)
        r12=r12tab(material)
        u0=u0tab(material)
        u1=u1tab(material)
        u2=u2tab(material)
        u11=u11tab(material)
        u12=u12tab(material)
        u22=u22tab(material)
        lambda=lambdatab(material)
    end subroutine


    subroutine tb_parametersgga(material)
        implicit none
        integer :: material
        double precision, dimension(6) :: &
             atab =(/0.3190,0.3191,0.3326,0.3325,0.3557,0.3560/),&
             e1tab =(/0.683,0.717,0.684,0.728,0.588,0.697/),&
             e2tab =(/1.707,1.916,1.546,1.655,1.303,1.380/),&
             t0tab =(/-0.146,-0.152,-0.146,-0.146,-0.226,-0.109/),&
             t1tab =(/-0.114,-0.097,-0.130,-0.124,-0.234,-0.164/),&
             t2tab =(/0.506,0.590,0.432,0.507,0.036,0.368/),&
             t11tab=(/0.085,0.047,0.144,0.117,0.400,0.204/),&
             t12tab=(/0.162,0.178,0.117,0.127,0.098,0.093/),&
             t22tab=(/0.073,0.016,0.075,0.015,0.017,0.038/),&
             r0tab =(/0.060,0.069,0.039,0.036,0.003,-0.015/),&
             r1tab =(/-0.236,-0.261,-0.209,-0.234,-0.025,-0.209/),&
             r2tab =(/0.067,0.107,0.069,0.107,-0.169,0.107/),&
             r11tab=(/0.016,-0.003,0.052,0.044,0.082,0.115/),&
             r12tab=(/0.087,0.109,0.060,0.075,0.051,0.009/),&
             u0tab =(/-0.038,-0.054,-0.042,-0.061,0.057,-0.066/),&
             u1tab =(/0.046,0.045,0.036,0.032,0.103,0.011/),&
             u2tab =(/0.001,0.002,0.008,0.007,0.187,-0.013/),&
             u11tab =(/0.266,0.325,0.272,0.329,-0.045,0.312/),&
             u12tab =(/-0.176,-0.206,-0.172,-0.202,-0.141,-0.177/),&
             u22tab =(/-0.150,-0.163,-0.150,-0.164,0.087,-0.132/),&
             lambdatab=(/0.073,0.211,0.091,0.228,0.107,0.237/)
        a=atab(material)
        e1=e1tab(material)
        e2=e2tab(material)
        t0=t0tab(material)
        t1=t1tab(material)
        t2=t2tab(material)
        t11=t11tab(material)
        t12=t12tab(material)
        t22=t22tab(material)
        r0=r0tab(material)
        r1=r1tab(material)
        r2=r2tab(material)
        r11=r11tab(material)
        r12=r12tab(material)
        u0=u0tab(material)
        u1=u1tab(material)
        u2=u2tab(material)
        u11=u11tab(material)
        u12=u12tab(material)
        u22=u22tab(material)
        lambda=lambdatab(material)
     end subroutine



    subroutine HTBTNN(H,kx,ky)
        implicit none
        real*8,intent(in)::kx,ky
        complex*16,intent(out)::H(6,6)
        complex*16::V1,V2,V12,H3(3,3),Lz(3,3)
        real*8::V0,V11,V22,alpha,beta,c2a,ca,cb,sa,sb,s2a,&
            c3a,c2b,c4a,s3a,s2b
        !real*8::e1,e2,t0,t1,t2,t11,t12,t22,r0,r1,r2,r11,r12,&
        !    u0,u1,u2,u11,u12,u22,lambda
        !common/tbparameter/e1,e2,t0,t1,t2,t11,t12,t22,r0,r1,r2,r11,r12,&
        !    u0,u1,u2,u11,u12,u22,lambda
        integer::m,n
        alpha=kx*a/2.d0;
        beta=ky*a*SR3/2.d0
        c2a=cos(2.d0*alpha); s2a=sin(2.d0*alpha)
        ca=cos(alpha); cb=cos(beta)
        sa=sin(alpha); sb=sin(beta)
        c3a=cos(3.d0*alpha); c2b=cos(2.d0*beta)
        c4a=cos(4.d0*alpha); s3a=sin(3.d0*alpha)
        s2b=sin(2.d0*beta)
        V0=e1+2.d0*t0*(2.d0*ca*cb+c2a)&
            +2.d0*r0*(2.d0*c3a*cb+c2b)&
            +2.d0*u0*(2.d0*c2a*c2b+c4a)
        V1=-2.d0*SR3*t2*sa*sb&
            +2.d0*(r1+r2)*s3a*sb&
            -2.d0*SR3*u2*s2a*s2b&
            +z*(2.d0*t1*sa*(2.d0*ca+cb)&
            +2.d0*(r1-r2)*s3a*cb&
            +2.d0*u1*s2a*(2.d0*c2a+c2b))
        V2=2.d0*t2*(c2a-ca*cb)&
            -2.d0*(r1+r2)*(c3a*cb-c2b)/SR3&
            +2.d0*u2*(c4a-c2a*c2b)&
            +z*(2.d0*SR3*t1*ca*sb&
            +2.d0*(r1-r2)*sb*(c3a+2.d0*cb)/SR3&
            +2.d0*SR3*u1*c2a*s2b)
        V11=e2+(t11+3.d0*t22)*ca*cb&
            +2.d0*t11*c2a+4.d0*r11*c3a*cb&
            +2.d0*(r11+SR3*r12)*c2b&
            +(u11+3.d0*u22)*c2a*c2b+2.d0*u11*c4a
        V12=SR3*(t22-t11)*sa*sb&
            +4.d0*r12*s3a*sb&
            +SR3*(u22-u11)*s2a*s2b&
            +z*(4.d0*t12*sa*(ca-cb)&
            +4.d0*u12*s2a*(c2a-c2b))
        V22=e2+(3.d0*t11+t22)*ca*cb&
            +2.d0*t22*c2a+2.d0*r11*(2.d0*c3a*cb+c2b)&
            +2.d0*r12*(4.d0*c3a*cb-c2b)/SR3&
            +(3.d0*u11+u22)*c2a*c2b+2.d0*u22*c4a
        
        H3=0.d0; Lz=0.d0
        !upper triangle of H3x3
        H3(1,1)=V0;
        H3(1,2)=V1;
        H3(1,3)=V2;
        !print *, 'H3(1,3)', H3(1,3)
        H3(2,2)=V11;
        H3(2,3)=V12;
        
        H3(3,3)=V22;
        !lower triangle of H3x3
        do m=1,2
            do n=m+1,3
                H3(n,m)=conjg(H3(m,n))
            end do
        end do
        
        Lz(2,3)=z*lambda;
        Lz(3,2)=-z*lambda;
        H=0.d0
        do m=1,3
            do n=1,3
                H(m,n)=H3(m,n)+Lz(m,n)
                H(m+3,n+3)=H3(m,n)-Lz(m,n)
            end do
        end do
    end subroutine
end module para_MX2


module prepare
    use para_MX2
    use parameter
    implicit none
    contains
    subroutine tbham(kx,ky,hamu,hamd)
        implicit none
        real*8 :: kx,ky,alpha,beta
        real*8 :: c2a,ca,cb,sa,sb,s2a,c3a,c2b,c4a,s3a,s2b
        complex*16, dimension(3,3) :: ham0,ham1,hamu,hamd
        alpha=kx*a/2
        beta=SR3*ky*a/2
        c2a=cos(2.d0*alpha); s2a=sin(2.d0*alpha)
        ca=cos(alpha); cb=cos(beta)
        sa=sin(alpha); sb=sin(beta)
        c3a=cos(3.d0*alpha); c2b=cos(2.d0*beta)
        c4a=cos(4.d0*alpha); s3a=sin(3.d0*alpha)
        s2b=sin(2.d0*beta)
        ham0(1,1)=e1+2*t0*(2*ca*cb+c2a)&
                 +2*r0*(2*c3a*cb+c2b)&
                 +2*u0*(2*c2a*c2b+c4a)
        ham0(1,2)=-2.d0*sqrt(3.d0)*t2*sa*sb&
                        +2.d0*(r1+r2)*s3a*sb&
                        -2.d0*sqrt(3.d0)*u2*s2a*s2b+z*(&
                        2.d0*t1*sa*(2.d0*ca+cb)&
                        +2.d0*(r1-r2)*s3a*cb&
                        +2.d0*u1*s2a*(2.*c2a+c2b))
        ham0(1,3)=2.d0*t2*(c2a-ca*cb)&
        -2.d0*(r1+r2)*(c3a*cb-c2b)/SR3&
        +2.d0*u2*(c4a-c2a*c2b)&
        +z*(2.d0*SR3*t1*ca*sb&
        +2.d0*(r1-r2)*sb*(c3a+2.d0*cb)/SR3&
        +2.d0*SR3*u1*c2a*s2b)
        !print*, 'ham13', ham0(1,3)
        ham0(2,1)=conjg(ham0(1,2))
        ham0(2,2)=e2+(t11+3*t22)*ca*cb+2*t11*c2a&
                 +4*r11*c3a*cb+2*(r11+SR3*r12)*c2b&
                 +(u11+3*u22)*c2a*c2b+2*u11*c4a
        ham0(2,3)=SR3*(t22-t11)*sa*sb&
                        +4*r12*s3a*sb&
                        +SR3*(u22-u11)*s2a*s2b+z*(&
                        4*t12*sa*(ca-cb)&
                        +4*u12*s2a*(c2a-c2b))
        ham0(3,1)=conjg(ham0(1,3))
        ham0(3,2)=conjg(ham0(2,3))
        ham0(3,3)=e2+(3*t11+t22)*ca*cb+2*t22*c2a&
                 +2*r11*(2*c3a*cb+c2b)&
                 +2*r12*(4*c3a*cb-c2b)/SR3&
                 +(3*u11+u22)*c2a*c2b+2*u22*c4a
        ham1(1,1)=0.
        ham1(1,2)=0.
        ham1(1,3)=0.
        ham1(2,1)=conjg(ham1(1,2))
        ham1(2,2)=0.
        ham1(2,3)=lambda*cmplx(0.,1.)
        ham1(3,1)=conjg(ham1(1,3))
        ham1(3,2)=conjg(ham1(2,3))
        ham1(3,3)=0.
        hamu=ham0+ham1
        hamd=ham0-ham1
     end subroutine
     subroutine gradient_H(kxs, kys, grad_H)
        implicit none
        double precision, intent(In) :: kxs, kys
        double COMPLEX, allocatable, dimension(:,:) :: hamu_plus, hamu_minus, hamd_plus, hamd_minus
        double COMPLEX, dimension(6,6,2), INTENT(OUT) :: grad_H
        real :: deltas = 5e-2
    allocate(hamd_plus(3,3), hamd_minus(3,3), hamu_plus(3,3), hamu_minus(3,3))
       call tbham(kxs + deltas, kys, hamu_plus , hamd_plus )
       call tbham(kxs - deltas, kys, hamu_minus, hamd_minus)
       grad_H(1:3,1:3,1) = (hamu_plus - hamu_minus) / (2 * deltas)
       grad_H(4:6,4:6,1) = (hamd_plus - hamd_minus) / (2 * deltas)
       call tbham(kxs, kys + deltas, hamu_plus , hamd_plus )
       call tbham(kxs, kys - deltas, hamu_minus, hamd_minus)
       grad_H(1:3,1:3,2) = (hamu_plus - hamu_minus) / (2 * deltas)
       grad_H(4:6,4:6,2) = (hamd_plus - hamd_minus) / (2 * deltas)
    deallocate(hamd_plus, hamd_minus, hamu_plus, hamu_minus)
        end subroutine gradient_H
end module prepare


module vector
    implicit none
    contains
    subroutine commu( As, Bs, comm)
        implicit none
        integer :: nk1, nk2
        double complex, dimension(:,:,:,:), intent(in) :: As, Bs
        double complex, dimension(size(As,1), size(As, 2), size(As,3), size(As,4)) :: comm
        double complex, dimension(size(As,1), size(As, 2)) :: mata, matb
        comm = 0.
        !$OMP PARALLEL DO PRIVATE(mata, matb)
        do nk2 = 1, size(As,4)
        do nk1 = 1, size(As,3)
            mata = matmul(As(:,:,nk1,nk2), Bs(:,:,nk1, nk2))
            matb = matmul(Bs(:,:,nk1,nk2), As(:,:,nk1, nk2))
            comm(:,:,nk1,nk2) = mata - matb
        end do
        end do
        !$OMP END PARALLEL DO
    end subroutine commu
    function deltakronechker(i,j)
        implicit none
        integer, intent(in) :: i, j
        integer :: deltakronechker
        deltakronechker = 0.
        if (i == j) then
            deltakronechker = 1.
            return
        end if
    end function deltakronechker
end module vector


module varia
    use parameter
    use para_MX2
    use prepare
    implicit none
    real*8, dimension(nkamax, nkamax, 2) :: grid
    double complex, allocatable, dimension(:,:,:,:,:) :: p, Xi
    double complex, allocatable, dimension(:,:,:,:) :: V
    double precision:: Elec(2), Alec(2), Atx, Aty, Etx, Ety, rk1, rk2
    double precision:: E(Norbital, nkamax, nkamax), ratio, cconst
    integer :: ctr
    integer , allocatable :: map(:,:)
    logical, allocatable, dimension(:,:,:) :: cutoff
    contains
    logical function isHermitian(As)
        implicit none
        double complex, dimension(:,:,:,:) :: As
        integer::i,j,m,n
        isHermitian=.true.
        do i=1, size(As, 4)
            do j=1, size(As, 3)
                do m=1, size(As, 2)
                    if (abs(imag(As( m, m, j, i))) .gt. 1.d-12) then
                        print *, As( m, m, j, i)
                        isHermitian=.false.
                        return
                    end if
                    do n= m+1, size(As,1)
                        if (abs(As( n, m, j, i) - conjg(As( m, n, j, i))) .gt. 1.d-12) then
                            print *, As( n, m, j, i), As( m, n, j, i)
                            isHermitian=.false.
                            return
                        end if
                    end do
                end do
            end do
        end do
    end function



    ! in k cutoff to reduce number of points




    subroutine variable()
        implicit none
        integer :: Deltas(Norbital, Norbital)

        real*8, dimension(2,2) :: B
        double complex, dimension(:,:), allocatable :: Vpkx, Vpky, hamu, hamd, ham
        double complex, allocatable :: Vp(:,:,:)
        double complex :: WORK(11)
        double precision :: kx,ky, RWORK(16), E_sorted(6), array(2), dS, dkx, dky, k1max, k2max, dk1, dk2
        double precision :: k1(nkamax), k2(nkamax)
        integer :: nk1, nk2, nu, mu, INFO
        allocate(ham(Norbital,Norbital),hamu(halfNorbital,halfNorbital),&
        hamd(halfNorbital,halfNorbital), &
        Vpkx(Norbital,Norbital), Vpky(Norbital,Norbital), Vp(Norbital,Norbital,2))
        allocate(p(Norbital,Norbital,nkamax,nkamax,2), Xi(Norbital,Norbital,nkamax,nkamax,2))
        allocate(V(Norbital,Norbital,nkamax,nkamax))
    
        Deltas = (0.,0.)
        V = (0.,0.)
        print *, 'Begin calculate variable'
        call tb_parameters(comid, usegga)
        B(1,1) =  2.*pi/a
        B(1,2) =  2.*pi/a
        B(2,1) =  2.*pi/(SR3*a)
        B(2,2) = -2.*pi/(SR3*a)
        print *, 'Grid calculate begin'
        

        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            array(1) = real(nk1 - 1)/real(nkamax - 1)
            array(2) = real(nk2 - 1)/real(nkamax - 1)
            grid( nk1, nk2,:) = matmul(B, array)
            grid( nk1, nk2, 1) = grid( nk1, nk2, 1) - 2.*pi/a
            grid( nk1, nk2, 2) = grid( nk1, nk2, 2) - 0.
        end do
        end do
    
        
        open(unit = 1, file = 'OUTPUT/grid.txt')
        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            write(1,*) grid( nk1, nk2, 1), grid( nk1, nk2, 2)
        end do
            write(1,*)
        end do
        close(1)
    


        k1max=(2.d0*pi)/(SR3*a)
        k2max=(2.d0*pi)/(SR3*a)
        dk1=2.d0*k1max/(Nkamax-1)
        dk2=2.d0*K2max/(Nkamax-1)
        do nk1=1,nkamax
            k1(nk1)=-k1max+(nk1-1)*dk1;
        end do
        do nk2=1,Nkamax
            k2(nk2)=-k2max+(nk2-1)*dk2;
        end do
        dkx=abs(SR3*(k1(1)+k2(1))/2.d0-SR3*(k1(2)+k2(1))/2.d0)
        dky=abs((k2(1)-k1(1))/2.d0-(k2(1)-k1(2))/2.d0)
        print*, "dkx", dkx
        print*, "dky", dky
        dS = dkx*dky
        ratio = dS/(pix2)**2.d0
        print *, 'ratio:', ratio
        print *, 'Initial Density of Electron in Valence Band:', ratio*2*nkamax**2
        cconst = qe**2/(2.d0*varepsilon0*varepsilon)*dkx*dky/(pix2)**2
        print*, 'Coulomb Interaction Constant:', cconst

        print*, 'Grid calculate done. Begin calculate E, V, p, Xi'
    
        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            E_sorted = 0.
            kx = grid(nk1,nk2,1)
            ky = grid(nk1,nk2,2)
            call tbham(kx, ky, hamu, hamd)
            ham = 0.
            ham(1:3, 1:3) = hamu
            ham(4:6, 4:6) = hamd
            call zheev('V', 'U', 6, ham, 6, E_sorted, WORK, 11, RWORK, INFO)
            E( :, nk1, nk2) = E_sorted
            V( :, :, nk1, nk2) = ham
        end do
        end do
    
    
        p = (0.,0.)

        Egap = minval(E( 3, :, :)) - maxval(E( 2, :, :))
        print *, "Energy Band Gap:", Egap, " (eV). Detuning:", detuning, " (eV)"
        w0 = (Egap + detuning)/hb
    
        open(unit = 1, file = 'OUTPUT/E.txt')
        open(unit = 2, file = 'OUTPUT/EBS.txt')
        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            write(1,*) grid( nk1, nk2, 1), grid( nk1, nk2, 2), (E( nu, nk1, nk2) , nu = 1,6)
            if (grid(nk1,nk2,2).eq. 0d0) then
                write(2,*) grid( nk1, nk2, 1), (E( nu, nk1, nk2) , nu = 1,6)
            end if
        end do
            write(1,*)
        end do
        close(1)
    
        !momentum calculation

        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            call gradient_H(grid(nk1,nk2,1), grid(nk1,nk2,2), Vp)
            p( :, :, nk1, nk2, 1) = me/hb*matmul(conjg(transpose(V(:,:,nk1,nk2))), matmul(Vp(:,:,1), V(:,:,nk1,nk2)))
            p( :, :, nk1, nk2, 2) = me/hb*matmul(conjg(transpose(V(:,:,nk1,nk2))), matmul(Vp(:,:,2), V(:,:,nk1,nk2)))
            do mu = 1, Norbital
            do nu = 1, Norbital
            if (abs(E( nu, nk1, nk2) - E( mu, nk1, nk2)) < 1e-2) then
                Deltas(nu, mu) = Deltas(nu, mu) + 1
            end if
            end do
            end do
        end do
        end do

        !approximation
        !p( :, 5, :, :, :) =0.; p( 5, :, :, :, :) =0.; p( :, 6, :, :, :) =0.; p( 6, :, :, :, :) =0.
        !p( :, 4, :, :, :) =0.; p( 4, :, :, :, :) =0.!; p( :, 3, :, :, :) =0.; p( 3, :, :, :, :) =0.
        !p( :, 1, :, :, :) =0.; p( 1, :, :, :, :) =0.; 
        !p( :, 2, :, :, :) =0.; p( 2, :, :, :, :) =0.
        if ((isHermitian(p(:,:,:,:,1))) .and. (isHermitian(p(:,:,:,:,2)))) then
            print *, "p is Hermitian"
        else
            print *, "Error: p is not Hermitian"
        end if
        call Hermitize(p(:,:,:,:,1))
        call Hermitize(p(:,:,:,:,2))
        
        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            do nu = 1, Norbital
            do mu = nu + 1, Norbital
        if (Deltas(nu,mu) == 0.) then
            Xi( nu, mu, nk1, nk2, 1) = - z * hb / me &
            * p( nu, mu, nk1, nk2, 1)/(E( nu, nk1, nk2) - E( mu, nk1, nk2))
            Xi( nu, mu, nk1, nk2, 2) = - z * hb / me &
            * p( nu, mu, nk1, nk2, 2)/(E( nu, nk1, nk2) - E( mu, nk1, nk2))
            Xi( mu, nu, nk1, nk2, 1) = conjg(Xi( nu, mu, nk1, nk2, 1))
            Xi( mu, nu, nk1, nk2, 2) = conjg(Xi( nu, mu, nk1, nk2, 2))
        end if
            end do
            end do
        end do
        end do
    
        print*, "Begin checking Hermitian matrix"
        if ((isHermitian(Xi(:,:,:,:,1))) .and. (isHermitian(Xi(:,:,:,:,2)))) then
            print *, "Xi is Hermitian"
        else
            print *, "Error: Xi is not Hermitian"
        end if
        call Hermitize(Xi(:,:,:,:,1))
        call Hermitize(Xi(:,:,:,:,2))

        !open(unit = 1, file = 'OUTPUT/p.txt')
        !do nk2 = 1, nkamax
        !do nk1 = 1, nkamax
        !    do nu = 1, Norbital
        !    do mu = 1, Norbital
        !write(1,'(2(I5,2x),6(E25.15, 2x))') mu, nu, grid( nk1, nk2, 1), grid( nk1, nk2, 2), &
        !real(p(mu, nu, nk1, nk2, 1)), imag(p(mu, nu, nk1, nk2, 1)), &
        !real(p(mu, nu, nk1, nk2, 2)), imag(p(mu, nu, nk1, nk2, 2))
        !    end do
        !    end do
        !end do
        !    write(1,*)
        !end do
        !close(1)

        !open(unit = 1, file = 'OUTPUT/Xi.txt')
        !do nk2 = 1, nkamax
        !do nk1 = 1, nkamax
        !    do nu = 1, Norbital
        !    do mu = 1, Norbital
        !write(1,'(2(I5,2x),6(E25.15, 2x))') mu, nu, grid( nk1, nk2, 1), grid( nk1, nk2, 2), &
        !real(xi(mu, nu, nk1, nk2, 1)), imag(xi(mu, nu, nk1, nk2, 1)), &
        !real(xi(mu, nu, nk1, nk2, 2)), imag(xi(mu, nu, nk1, nk2, 2))
        !    end do
        !    end do
        !end do
        !    write(1,*)
        !end do
        !close(1)


        print *, 'Calculate E, V, p, Xi done. Calculate integral ratio...'
        print *, 'Calcualte variable done'
    
        print *, 'Begin calculate map array'
        allocate(cutoff(nkamax, nkamax, 2))        
        cutoff = .false.
            ctr=1
            do nk1 = 1, nkamax
            do nk2 = 1, nkamax
                rk1 = sqrt((grid(nk1,nk2,1)-4.d0*pi/(3.d0*a))**2+grid(nk1,nk2,2)**2)
                rk2 = sqrt((grid(nk1,nk2,1)+4.d0*pi/(3.d0*a))**2+grid(nk1,nk2,2)**2)
                    if (rk1.le.kcutoff) then
                       cutoff(nk1,nk2,1) = .true.
                       ctr = ctr + 1
                    elseif (rk2.le.kcutoff) then
                       cutoff(nk1,nk2,2) = .true.
                       ctr = ctr +1
                    end if
            end do
            end do
        print *, 'Calculate map array done'
            ctr = ctr - 1 !ctr now is the number of k-cutoff points
            nkqcutoff = ctr
            print *, 'Number of k-cutoff points:', nkqcutoff
        allocate(map(nkqcutoff, 2))

        !mapping k in cutoff to map array
        ctr = 1
        do nk1 = 1, nkamax
        do nk2 = 1, nkamax
            if (cutoff(nk1,nk2,1)) then
                map(ctr, 1) = nk1
                map(ctr, 2) = nk2
                ctr = ctr + 1
            end if
        end do
        end do

        do nk1 = 1, nkamax
        do nk2 = 1, nkamax
            if (cutoff(nk1,nk2,2)) then
                map(ctr, 1) = nk1
                map(ctr, 2) = nk2
                ctr = ctr + 1
            end if
        end do
        end do

        if (ctr /= nkqcutoff + 1) stop 'Error: mapping k to map array'


        open(unit = 1, file = 'OUTPUT/gridcutoff.txt')
        do ctr = 1, nkqcutoff
            write(1,*) grid(map(ctr, 1), map(ctr, 2), 1), grid(map(ctr, 1), map(ctr, 2), 2)
        end do

        close(1)
        end subroutine variable
        
        !**********************************************************************
        subroutine Et(t,Ex,Ey)
            implicit none
            real*8,intent(in)::t
            real*8,intent(out)::Ex,Ey
            real*8::env
            env=E0*exp(-t**2/td**2)
            Ex=env*cos(w0*t)/sqrt(1.d0+ell**2)
            Ey=env*sin(w0*t)*ell/sqrt(1.d0+ell**2)
        end subroutine
        !**********************************************************************
        subroutine genAt()
            implicit none
            double precision :: time
            real*8 :: Ext, Eyt, Axt, Ayt
            integer :: nt
            Axt = 0. ; Ayt = 0.
            Ext = 0. ; Eyt = 0.
            open(unit = 1, file = 'OUTPUT/EAt.txt')
            print *, "w0", w0, "td", td, "E0", E0
            print *, "tmin", tmin, "tmax", tmax
            do nt = 1, ntmax
                time = tmin + (nt) * dt
                call Et(time, Ext, Eyt)
                Axt = Axt - Ext * dt
                Ayt = Ayt - Eyt * dt
                write(1,*) time, Ext, Eyt, Axt, Ayt
            end do
            close(1)
            !111 format(f10.4,2x, 5(E15.8,2x))
        end subroutine genAt

    !*********************************************************************
    subroutine Hermitize(ma)
        implicit none
        complex*16,intent(out)::ma(:,:,:,:)
        integer ::i,j,m,n
        do i=1,size(ma, 4)
        do j=1,size(ma, 3)
            do m=1,size(ma, 2)
            ma(m,m,j,i)=real(ma(m,m,j,i))
            do n=m+1,size(ma, 1)
                ma(m,n,i,j)=(ma(m,n,i,j)+conjg(ma(n,m,i,j)))/2.d0
                ma(n,m,i,j)=conjg(ma(m,n,i,j))
            end do
            end do
        end do
        end do
    end subroutine Hermitize
end module varia


module SBE
    use varia
    use vector
    implicit none
    complex*16, allocatable, dimension(:,:,:,:) :: Vcs, rho, Kappa
    real*8, allocatable, dimension(:,:) :: Jt, Pt, dent
    real*8, allocatable, dimension(:) :: dene, denh
    real*8 :: Atxd, Atyd
    contains
    subroutine calculateVcs
        implicit none
        integer :: i, j,m,n
        complex*16 :: s(6,6)
        allocate(Vcs(Nkqcutoff,Nkqcutoff,4,4))
        Vcs=0.d0
        if (usecoulomb .eqv. .false.) then
            print *, 'No Coulomb Interaction. Vcs = 0.'
            return
        end if
        !$omp parallel do private(S)
        do i=1,Nkqcutoff
        do j=1,Nkqcutoff
                if ((cutoff(map(i,1),map(i,2),1)&
                    .and.cutoff(map(j,1),map(j,2),1))&
                    .or.(cutoff(map(i,1),map(i,2),2)&
                    .and.cutoff(map(j,1),map(j,2),2))) then
                    S=matmul(transpose(conjg(V(:,:,map(j,1),map(j,2)))),&
                        V(:,:,map(i,1),map(i,2)))
                    do m=1,4
                        do n=1,4
                            Vcs(j,i,m,n)=S(m,n)
                        end do
                    end do
                end if
            end do
        end do
        !$omp end parallel do
        print *, 'Calculate Vcs done'
    end subroutine calculateVcs
    !subroutine calkappa
    !    implicit none
    !    integer::i,j,m,n
    !    complex*16::S(6,6)
    !    write(*,'(a)')"..Calculating kappa"
    !    write(*,'(a,i4)')"Number of k-points cutoff ",nkqcutoff
    !    allocate(kappa(Nkqcutoff,Nkqcutoff,4,4))
    !    kappa=0.d0
    !    !$omp parallel do private(S)
    !    do i=1,Nkqcutoff
    !        do j=1,Nkqcutoff
    !            if ((cutoff(map(i,1),map(i,2),1)&
    !                .and.cutoff(map(j,1),map(j,2),1))&
    !                .or.(cutoff(map(i,1),map(i,2),2)&
    !                .and.cutoff(map(j,1),map(j,2),2))) then
    !                S=matmul(transpose(conjg(V(:,:,map(j,1),map(j,2)))),&
    !                    V(:,:,map(i,1),map(i,2)))
    !                do m=1,4
    !                    do n=1,4
    !                        kappa(j,i,m,n)=S(m,n)
    !                    end do
    !                end do
    !            end if
    !        end do
    !    end do
    !    !$omp end parallel do
    !    write(*,'(a)')"done!"
    !end subroutine
    subroutine calHF(u, fh)
        implicit none
        complex*16, intent(in) :: u(Norbital, Norbital, nkamax, nkamax)
        complex*16, intent(out) :: fh(Norbital, Norbital, nkamax, nkamax)
        integer :: nk1, nk2, nu, mu, al, be
        double complex :: coeff2, coeff1, dene, kappa1
        double complex :: s1
        !cal density of electron on condution band
        dene = 0.d0; kappa1 = 0.d0
        if (useshieldcoulomb) then
            do nu = 3, Norbital 
                dene = dene + sum(u(nu, nu, :, :))
            end do
            dene = dene*ratio
            kappa1 = 2*me*qe**2/(varepsilon0*varepsilon*hb**2)*(1.d0 - exp(- hb**2*Thermalbeta*pi*dene/me))
        end if
        fh = 0.d0
        if (usecoulomb .eqv. .false.) return
        coeff1=qe**2/(2.d0*varepsilon0*varepsilon)*ratio
    	!print*,"COULOMB coeff", coeff1
        !stop
        !$OMP PARALLEL DO PRIVATE(s1,coeff2)
        do nk1 = 1, nkqcutoff
            do mu = 1, 4
            do nu = mu, 4
                s1 = 0.d0
        do nk2 = 1, nkqcutoff
            if ((cutoff(map(nk1,1),map(nk1,2),1)&
            .and.cutoff(map(nk2,1),map(nk2,2),1))&
            .or.(cutoff(map(nk1,1),map(nk1,2),2)&
            .and.cutoff(map(nk2,1),map(nk2,2),2))) then
                if (nk2.ne.nk1) then
                    coeff2 = sqrt((grid(map(nk2,1),map(nk2,2),1)&
                        -grid(map(nk1,1),map(nk1,2),1))**2&
                        +(grid(map(nk2,1),map(nk2,2),2)&
                        -grid(map(nk1,1),map(nk1,2),2))**2)
                    coeff2=1.d0/(coeff2)
                    do be=1,4
                    do al=1,4
                        !if ((al.le.2 .and. be .ge. 3).or.(al.ge.3 .and. be .le. 2)) then
            if (al.ne.be) then
                s1 = s1 - Vcs(nk2,nk1,al,nu)*Vcs(nk1,nk2,mu,be)*coeff2*u(be,al,map(nk2,1),map(nk2,2))
            elseif (al.gt.3) then
                s1 = s1 - Vcs(nk2,nk1,al,nu)*Vcs(nk1,nk2,mu,be)*coeff2*u(be,al,map(nk2,1),map(nk2,2))
            !else
            !    s1 = s1 + Vcs(nk2,nk1,al,nu)*Vcs(nk1,nk2,mu,be)*coeff2*(1-u(be,al,map(nk2,1),map(nk2,2)))
            end if
                    end do
                    end do
                end if
            end if
        end do
                fh( mu, nu, map(nk1,1), map(nk1,2)) = s1*cconst
                if (nu.ne.mu) fh( nu, mu, map(nk1,1), map(nk1,2))=conjg(fh( mu, nu, map(nk1,1), map(nk1,2)))
            end do
            end do
        end do
        !$OMP END PARALLEL DO
    end subroutine calHF

    subroutine calHFC(u,fh)
        implicit none
        complex*16,intent(in)::u(6,6,nkamax,nkamax)
        complex*16,intent(out)::fh(nkamax,nkamax,6,6)
        real*8::coeff1,coeff2, scrwn, dene
        integer::i,j,m,n,mm,nn
        complex*16::s1
        fh=0.d0
        coeff1=qe**2/(2.d0*varepsilon0)*ratio
        !$omp parallel do private(s1,coeff2)
        do i=1,nkqcutoff
            do m=1,4
            do n=m,4
                s1=0.d0
                do j=1,nkqcutoff
                    if ((cutoff(map(i,1),map(i,2),1)&
                        .and.cutoff(map(j,1),map(j,2),1))&
                        .or.(cutoff(map(i,1),map(i,2),2)&
                        .and.cutoff(map(j,1),map(j,2),2))) then
                    if (j.ne.i) then
                        coeff2=sqrt((grid(map(j,1),map(j,2),1)&
                            -grid(map(i,1),map(i,2),1))**2&
                            +(grid(map(j,1),map(j,2),2)&
                            -grid(map(i,1),map(i,2),2))**2)
                        coeff2=1.d0/(coeff2)
                        do mm=1,4
                        do nn=1,4
                            !----------------------------------------
                            !if (nn.eq.mm.and.mm.le.2) then
                            !	s1=s1-Vcs(j,i,mm,n)&
                            !		*Vcs(i,j,m,nn)&
                            !		*(u(map(j,1),map(j,2),nn,mm) - 1.d0)*coeff2
                            !elseif (nn.eq.mm.and.mm.ge.3) then
                            !	s1=s1-Vcs(j,i,mm,n)&
                            !		*Vcs(i,j,m,nn)&
                            !		*u(map(j,1),map(j,2),nn,mm)*coeff2
                            !else
                            if (nn.ne.mm) then
                                s1=s1-Vcs(j,i,mm,n)&
                                    *Vcs(i,j,m,nn)&
                                    *u(nn,mm,map(j,1),map(j,2))*coeff2
                            end if
                            !----------------------------------------
                        end do
                        end do
                    end if
                    end if
                end do
                fh(map(i,1),map(i,2),m,n)=s1*coeff1
                if (n.ne.m) fh(map(i,1),map(i,2),n,m)=conjg(fh(map(i,1),map(i,2),m,n))
            end do
            end do
        end do
        !$omp end parallel do
    end subroutine

    subroutine RHSVGsbe( ind, u, fr)
        implicit none
        complex*16, intent(in) :: u(Norbital, Norbital, nkamax, nkamax)
        complex*16, intent(out) :: fr(Norbital, Norbital, nkamax, nkamax)
        complex*16, allocatable, dimension(:,:,:,:) :: HF, coff, H, HC
        real*8 :: eAx, eAy, Tcd2, dene
        real*8, INTENT(in) :: ind
        integer :: nk1, nk2, nu, mu
        allocate(HF(Norbital, Norbital, nkamax, nkamax))
        allocate(H(Norbital, Norbital, nkamax, nkamax))
        allocate(HC(nkamax,nkamax,Norbital,Norbital))
        allocate(coff(Norbital, Norbital, nkamax, nkamax))
        dene = ratio*real(sum(u(3,3,:,:)) + sum(u(4,4,:,:)))
        if (usegammaeecor) then
            Tcd2 = 1./Tco2 + dene*gamma
        else
            Tcd2 = 1./Tco2
        end if
        fr = 0.d0 ; H = 0.d0; HF = 0.d0
        HC = 0.d0
        !print*, "calHF"
        eAx = qe /me * ( Atx - ind* Etx* dt)
        eAy = qe /me * ( Aty - ind* Ety* dt)
        !print*, maxval(real(- eAx*p(:,:,:,:,1) - eAy*p(:,:,:,:,2)))
        !$OMP PARALLEL DO
        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            do nu = 1, Norbital
                H(nu, nu, nk1, nk2) = H(nu, nu, nk1, nk2) + E(nu, nk1, nk2)
            end do
        end do
        end do
        !$OMP END PARALLEL DO
        H = H - eAx*p(:,:,:,:,1) - eAy*p(:,:,:,:,2)
        if (usecoulomb) then 
            call calHF(u, HF)
            H = H + HF
        end if
        !if (usecoulomb) then
        !    call calHFC(u, HC)
        !    do nu = 1, Norbital
        !    do mu = 1, norbital
        !        H(nu, mu, :, :) = H(nu, mu, :, :) + HC(:,:,nu,mu)
        !    end do
        !    end do
        !end if
        call commu(H, u, coff)
        !$OMP PARALLEL DO
        do nk2 = 1, nkamax
        do nk1 = 1, nkamax
            do nu = 1, Norbital
            do mu = nu, Norbital
                fr(nu, mu, nk1, nk2) = - z* coff(nu, mu, nk1, nk2)/hb
                if (nu /= mu) then
                    fr(nu, mu, nk1, nk2) = fr(nu, mu, nk1, nk2) - u(nu, mu, nk1, nk2)*Tcd2
                    fr(mu, nu, nk1, nk2) = conjg(fr(nu, mu, nk1, nk2))
                end if
            end do
            end do
        end do
        end do
        !$OMP END PARALLEL DO
        deallocate(HF, H, coff)
        deallocate(HC)
    end subroutine RHSVGsbe



    subroutine RK4sbeVG(rhos)
        implicit none
        double complex, intent(inout) :: rhos(Norbital, Norbital, nkamax, nkamax)
        double complex, allocatable, dimension(:,:,:,:) :: uu0, uu1, RHS
        
        allocate(uu0(Norbital, Norbital, nkamax, nkamax))
        allocate(uu1(Norbital, Norbital, nkamax, nkamax))
        allocate(RHS(Norbital, Norbital, nkamax, nkamax))
        uu0 = rhos
        call RHSVGsbe(0.d0, uu0, RHS) !RHS=k1
        rhos = rhos + RHS*dt/6.d0


        uu1 = uu0 + dt/2.d0 * RHS
        call RHSVGsbe(1/2.d0,uu1, RHS) !RHS=k2
        rhos = rhos + dt/3.d0 * RHS
        

        uu1 = uu0 + dt/2.d0 * RHS
        call RHSVGsbe(1/2.d0,uu1, RHS) !RHS=k3
        rhos = rhos + dt/3.d0 * RHS


        uu1 = uu0 + dt * RHS
        call RHSVGsbe(1.d0,uu1, RHS) !RHS=k4
        rhos = rhos + dt/6.d0 * RHS

        deallocate(uu0, uu1, RHS)
    end subroutine RK4sbeVG


    subroutine caljVG(u, j)
        implicit none
        double complex, intent(in) :: u(Norbital, Norbital, nkamax, nkamax)
        real*8, intent(out) :: j(2)
        real*8 :: As(2), s1, s2
        integer :: nk1, nk2, nu, mu, k
        As(1) = Atx*qe/me ; As(2) = Aty*qe/me; Jt = 0.
        do k = 1, 2
        s1 = 0. ; s2 = 0.
        !$OMP parallel do reduction(+:s1,s2)
        do nk1 = 1, nkamax
        do nk2 = 1, nkamax
            do mu = 1,6
                    s1=s1+real((p( mu, mu, nk1, nk2, k)/me - As(k))*u(mu, mu, nk1, nk2))
                do nu = mu +1,6
                    s2=s2+2.d0*real(p(mu, nu, nk1, nk2,k)/me*u(nu, mu, nk1, nk2))
                end do
            end do
        end do
        end do
        !$OMP end parallel do
        !check Hermitian
        j(k) = (s1+s2)*ratio
        end do
    end subroutine caljVG



    subroutine calP( u, j)
        implicit none
        double complex, intent(in) :: u(Norbital, Norbital, nkamax, nkamax)
        real*8, intent(out) :: j(2)
        real*8 :: s1, s2
        integer :: nk1, nk2, nu, mu
        j = 0.
        s1 = 0.; s2 = 0.
        !$OMP parallel do reduction(+:s1, s2)
        do nk1 = 1, nkamax
        do nk2 = 1, nkamax
            do nu = 1, Norbital
            do mu = nu+1,6
                s1=s1+2.d0*real(xi(mu, nu, nk1, nk2, 1)*u(nu, mu, nk1, nk2))
                s2=s2+2.d0*real(xi(mu, nu, nk1, nk2, 2)*u(nu, mu, nk1, nk2))
            end do
            end do
        end do
        end do
        !$OMP end parallel do
        j(1) = qe*s1*ratio
        j(2) = qe*s2*ratio
    end subroutine calP

    subroutine calden( rhos, dents, denes, denhs)
        implicit none
        double complex, intent(in) :: rhos(Norbital, Norbital, nkamax, nkamax)
        real*8, intent(out) :: dents(norbital), denes, denhs
        integer :: nk1, nk2, mu
        dent = 0.; dene = 0. ; denh = 0.
        !print*, maxval(real(rhos(3:6,3:6,:,:)))
        !check Hermitian
        do nk1 = 1, nkamax
        do nk2 = 1, nkamax
            do mu = 1, Norbital
                if (imag(rhos(mu, mu, nk1, nk2)) .gt. 0.) then
                    print*, "Error: rho diag is not Hermitian"
                    stop
                end if
            end do
        end do
        end do
        !$OMP parallel do reduction(+:dents, denhs)
        do nk1 = 1, nkamax
        do nk2 = 1, nkamax
            do mu = 1, Norbital
                dents(mu) = dents(mu) + real(rhos(mu, mu, nk1, nk2))
                if (mu .le. 2) then
                    denhs = denhs + real(1 - rhos(mu, mu, nk1, nk2))
                end if
            end do
        end do
        end do
        !$OMP END PARALLEL DO
        dents = dents*ratio
        denes = sum(dents(3:6))
        denhs = denhs*ratio
    end subroutine calden


    subroutine calAbsorp(j)
        implicit none
        real*8, intent(in)::j(0:ntmax,2)
        integer,parameter::Nw=2000
        real*8 ::wa,wb,w,dw,t,Ex,Ey, theta
        complex*16::eiwt,s1,s2,s3,s4,s5,s6
        integer::ii,tt
        write(*,'(a)')"..Calculating AbSpect"
        wa = (Egap - 1.d0)/hb
        wb = (Egap + 1.d0)/hb
        dw = (wb - wa)/Nw
        open(unit=10,file='OUTPUT/AbSpect.txt')
        do ii=0,Nw
            w=wa+ii*dw
            s1=0.d0; s2=0.d0; s3=0.d0; s4=0.d0; s5=0.d0; s6=0.d0 !s5 is P\vec{e} and s6 is \vec{E}
            do tt=0,ntmax
                t=tmin + (tt) * dt
                call Et(t,Ex,Ey)
                if (abs(Ex).le.1.d-16) then
                    theta = pi/2
                    else
                    theta = atan(Ey/Ex)
                end if
                eiwt=exp(z*w*t)
                s1=s1 + eiwt*j(tt,1)*dt
                s2=s2 + eiwt*j(tt,2)*dt
                s3=s3 + eiwt*Ex*dt
                s4=s4 + eiwt*Ey*dt
                s5=s5 + eiwt*(j(tt,1)*Ex + j(tt,2)*Ey)/sqrt(Ex**2 + Ey**2 + 1e-16)*dt
                s6=s6 + eiwt*(Ex*cos(theta) + Ey*sin(theta))*dt
            end do
            !pause
            write(10,'(20(E35.25, 2x))')w*hb,w*hb-Egap,real(s1),imag(s1),real(s2),imag(s2),&
                real(s3),imag(s3),real(s4),imag(s4),real(s1/s3),imag(s1/s3),&
                real(s2/s4),imag(s2/s4), real(s5), imag(s5), real(s6), imag(s6), real(s5/s6), imag(s5/s6)
        end do
        close(10)
    end subroutine



    
    subroutine writeinforonscreen(nt, start, finish)
        implicit none
        real*8, intent(in) :: start, finish
        integer, intent(in) :: nt
        write(*,*) "---------------------------------------------------------"
        write(*,*) "Time taken for simulation", 10, "fs:", finish - start, "seconds"
        write(*,*) "Current at", tmin + nt*dt, "fs:", real(Jt(nt, 1)), real(Jt(nt, 2))
        write(*,*) "Polarization at", tmin + nt*dt, "fs:", real(Pt(nt, 1)), real(Pt(nt, 2))
        write(*,*) 'Hole in Valence Band', real(denh(nt))
        write(*,*) 'Conduction Band', real(dene(nt))
    end subroutine




    subroutine calculateSBE
        use omp_lib
        implicit none
        integer :: nt, m
        double precision :: time
        real*8 :: start_time, end_time
        allocate(rho(norbital, norbital, nkamax, nkamax))
        allocate(Jt(ntmax, 2), Pt(ntmax, 2), dent(ntmax, Norbital))
        allocate(dene(ntmax), denh(ntmax))
        !initial condition
        rho = 0.
        rho(1, 1, :, :) = 1.
        rho(2, 2, :, :) = 1.
        Etx = 0. ; Ety = 0.
        Atx = 0. ; Aty = 0.

        open(unit = 1, file = 'OUTPUT/Current.txt')
        open(unit = 2, file = 'OUTPUT/Polar.txt')
        open(unit = 3, file = 'OUTPUT/Density.txt')

        !simulation
        do nt = 1, ntmax
            start_time = omp_get_wtime()

            
            time = tmin + (nt - 1) * dt
            call Et(time, Etx, Ety)
            call RK4sbeVG(rho)
            
            Atx = Atx - Etx * dt
            Aty = Aty - Ety * dt

            if (mod(tmin + nt*dt, 10.) == 0) then
                end_time = omp_get_wtime()
                call writeinforonscreen(nt, start_time, end_time)
                start_time = omp_get_wtime()
            end if

            call caljVG(rho, jt(nt,:))
            call calP(rho, Pt(nt,:))
            call calden(rho, dent(nt,:), dene(nt), denh(nt))
            
            write(1,*) time, real(jt(nt, 1)), real(jt(nt, 2))
            write(2,*) time, real(Pt(nt, 1)), real(Pt(nt, 2))
            write(3,'(F15.5,2x, 8(E15.5, 2x))') time, (real(dent(nt, m)), m = 1,6), real(dene(nt)), real(denh(nt))
        end do
        close(1)
        close(2)
        close(3)
    end subroutine calculateSBE
end module SBE


program main
    use parameter
    use varia
    use SBE
    implicit none
    call system('mkdir -p OUTPUT')
    call variable()
    print *, 'Band Gap', Egap
    call genAt()
    call calculateVcs()
    call calculateSBE()
    !allocate(Pt(ntmax,2))
    !Pt = 1.d0
    call calAbsorp(Pt)
end program main

