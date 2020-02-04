
!-------------Variables------------------! 


!-------------Routines------------------!
    module pulse
    use global_variables
    use domain
    use pulse_var
    use Functions
    implicit none
    private 
    save 
    public::Pulse_Static,One_Step_Pulse
    contains 
    subroutine Pulse_Static(z,Ac)  ! wistm k0 and Foc is normalized with kp 
    ! a0 is the normalized laser potential
    ! wist0 is the spot size
    ! k0 is the laser k in kp unit
    ! Foc is the focus plane measured in \xi. normally keep the same as z0
    ! z0 is the laser center measured in \xi
    ! Length the half length of laser 
    ! Types swicth the laser type
    ! Ac complex potential, is defined on half-integer gird
    ! As |Ac|^2
    implicit none
    real(RP)    ::  rhalf(1:Nr)
    real(RP)    ::  z
    real(RP)    ::  Wistz,CurvR,GouyPh,arg,argtot,phase
    integer     ::  k
    complex(Rp) ::  Ac(1:Nr)


    do k=1,Nr
        rhalf(k)=(k-1.5)*hr-Rmax-xcenter
    enddo

     select case(Types)
    
    case(0)     ! gaussian
    arg=exp(-2.77258872*((z-zcenter)/FWHM)**2);
     z=z-FocalPlane;
    RayLen=0.5d0*Wist0**2*k0
    Wistz=Wist0*sqrt(1+(z/RayLen)**2)
    GouyPh=atan(-z/RayLen)
    if (z==0) z=1d-8;
    CurvR=(-z)*(1+(RayLen/z)**2)
    !arg=(3.234*(z-zcenter)/FWHM-1.617)**2;

    case(1)     ! flat top pulse cosine drop

    if (abs(z-zcenter)<=FWHM*0.5d0) then
        arg=1d0
    else if(abs(z-zcenter)<=FWHM*0.5d0+rise)  then
        arg=cos((abs(z-zcenter)-FWHM*0.5d0)/rise*pi*0.5d0)
    else 
        arg=0d0
    end if

    case(2)     ! Trapezoidal

    if (abs(z-zcenter)<=FWHM) then
        arg=1d0
    else if(abs(z-zcenter)<=FWHM+rise)  then
        arg=1d0-(abs(z-zcenter)-FWHM)/rise
    else 
        arg=0d0
    end if

    end select

    do k=1,Nr
    argtot=(rhalf(k)/Wist0)**2
    !phase=k0*rhalf(k)**2/(2*CurvR)-GouyPh
    !Ac(k)=0.5*a0*(Wist0/Wistz)*exp(-argtot)*arg!*(cos(phase)+ci*sin(phase))  ! the definition is a0=a*exp()+c.c. so a=a0/2

    Ac(k)=0.5*a0*exp(-argtot)*arg!*(cos(phase)+ci*sin(phase))  ! the definition is a0=a*exp()+c.c. so a=a0/2
   ! if( abs(Ac(k))>0) write(*,*)   Ac(k)
    end do
    




    end subroutine Pulse_Static


    subroutine One_Step_Pulse(dd0,dd1,dd2,Chi,Aold,Anew,dz,dt)
    implicit none
    integer:: k,INFO
    real(RP) :: dd0,dt,dz,dzdt,hr2,delta
    real(RP) :: Chi(1:Nr),rhalf(1:Nr)
    complex(RP) :: Aold(1:Nr),Anew(1:Nr),dd1(1:Nr),dd2(1:Nr),dd3(1:Nr),source(1:Nr),laplace(1:Nr)
    complex(RP) ::aa(1:Nr-1),bb(1:Nr),cc(1:Nr-1),dd(1:Nr),ik0dt,beta
    dzdt=dz*dt
    ik0dt=ci*k0/dt
    hr2=1/hr**2



!---------laplace(A)---------
    laplace(1)=(Aold(2)-2d0*Aold(1))*hr2
    do k=2,Nr-1
    laplace(k)=(Aold(k-1)-2d0*Aold(k)+Aold(k+1))*hr2
    end do
!---------A_BC=0, tempeoray---------    
    laplace(Nr)=(Aold(Nr-1)-2d0*Aold(Nr))*hr2

!---------source---------
    do k=1,Nr
    source(k)=(Chi(k)*0.25+ik0dt-dd0/dzdt)*Aold(k)-laplace(k)*0.25+(-2d0*dd1(k)+0.5d0*dd2(k))/dzdt
    end do
!--------inversion coefficient

    if(Periodic) then
    !see Numerical Partial Differential Equations: Finite Difference Methods P. 225
    bb(1)=ik0dt-dd0/dzdt-Chi(1)*0.25-0.75*hr2
    cc(1)=0.25*hr2
    do k=2,Nr-1
        aa(k-1)=0.25*hr2
        bb(k)  =ik0dt-dd0/dzdt-Chi(k)*0.25-0.5*hr2
        cc(k)  =0.25*hr2
    end do
!---------A_BC=0, tempeoray---------  
    aa(Nr-1)=0.25*hr2
    bb(Nr)  =ik0dt-dd0/dzdt-Chi(Nr)*0.25-0.75*hr2
    dd(:)=source(:)
    call CGTSV(Nr, 1, aa, bb, cc, dd,Nr, INFO )
    beta=(dd(1)+dd(Nr))*0.25*hr2/(1+0.25*hr2*(1d0/bb(1)+1d0/bb(Nr)))
    dd(1)=dd(1)-beta/bb(1);dd(Nr)=dd(Nr)-beta/bb(Nr);
    
    else 

    bb(1)=ik0dt-dd0/dzdt-Chi(1)*0.25-0.5*hr2
    cc(1)=0.25*hr2
    do k=2,Nr-1
        aa(k-1)=0.25*hr2
        bb(k)  =ik0dt-dd0/dzdt-Chi(k)*0.25-0.5*hr2
        cc(k)  =0.25*hr2
    end do
!---------A_BC=0, tempeoray---------  
    aa(Nr-1)=0.25*hr2
    bb(Nr)  =ik0dt-dd0/dzdt-Chi(Nr)*0.25-0.5*hr2
    dd(:)=source(:)
    call CGTSV(Nr, 1, aa, bb, cc, dd,Nr, INFO )
    end if 

    Anew(:)=dd(:)
    end subroutine One_Step_Pulse


end module pulse

