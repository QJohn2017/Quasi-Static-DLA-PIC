module initialization

        use global_variables
        use pulse 
        use pulse_var
        use particles_var
        use fields_var
        use output_var
        use TestPart_var
        use Test_Particle
        use domain
        use InAndOut
        USE OneTstep

        implicit none
        private
        save
        public:: set_varibles
        contains

    subroutine set_varibles(start)
    implicit none
    real(RP):: start
    complex(Rp),allocatable:: Ac(:)
    integer:: n
    open(222,file='./QS.ini',status='old',action='read') 
    read(222,NML=Lasers)
    read(222,NML=Simulation)
    read(222,NML=TestParticles)
    close(222)
    write(*, *) "========Simulation Parameters========"
    write(*, *) "===                              ===="
    write(*,Lasers)
    write(*,Simulation)
    write(*,TestParticles)
    write(*, *) "===******************************===="


    if(MOD(Part_per_cell,2)==1)  then
    write(*, *) "=======Part_Per_Cell Must Even======="
    call EXIT(0)
    end if

    hr=Rmax/Nr
    dr=Hr/Part_per_cell
    Np=Nr*Part_per_cell+Part_per_cell/2  
    dz0=Zmax/Nz
    !NWake=Int(WakeLen/dz0)

    if (TestParts) then

    NTP_i=0
    do n=1,Num_of_Beam
    NTP_i(n)=floor(TP_len(n)*TP_perZ(n))*((1+floor(TP_Rad(n)*TP_perR(n)))*(floor(TP_Rad(n)*TP_perR(n))*N_per_peri(n))/2)
    end do
    NTP=sum(NTP_i)

    write(*, *) "=====Number of Test Beams   ==="
    write(*, *) Num_of_Beam
    write(*, *) "=====Number of TestParticles==="
    write(*, *) '(',NTP_i(1:Num_of_Beam),')'
    allocate(XXT(1:NTP));allocate(YYT(1:NTP));allocate(ZZT(1:NTP))
    allocate(PXT(1:NTP));allocate(PYT(1:NTP));allocate(PZT(1:NTP))
    allocate(WX_L(1:NTP));allocate(WZ_L(1:NTP));allocate(WX_W(1:NTP));allocate(WY_W(1:NTP));allocate(WZ_W(1:NTP))
    

    allocate(XXT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(YYT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(ZZT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(PXT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(PYT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(PZT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))

    end if 
    
    allocate(r(1:Np));allocate(V(1:Np));allocate(force(1:Np));allocate(mass(1:Np));allocate(gama(1:Np))
    allocate(psi(1:Nr+1));allocate(Er(1:Nr+1));allocate(Et(1:Nr+1));allocate(Bth(1:Nr+1))
    allocate(den(1:Nr+1));allocate(den_e(1:Nr+1));allocate(jr(1:Nr+1));allocate(jv(1:Nr+1));allocate(jz(1:Nr+1))
    

    allocate(r2d(1:Np,1:Nz+1));allocate(V2d(1:Np,1:Nz+1));allocate(force2d(1:Np,1:Nz+1));allocate(gama2d(1:Np,1:Nz+1))
    allocate(psi2d(1:Nr,1:Nz+1));allocate(Er2d(1:Nr,1:Nz+1));allocate(Et2d(1:Nr,1:Nz+1));allocate(Bth2d(1:Nr,1:Nz+1))
    allocate(den2d(1:Nr,1:Nz+1));allocate(den_e2d(1:Nr,1:Nz+1));allocate(jr2d(1:Nr,1:Nz+1));allocate(jz2d(1:Nr,1:Nz+1))
    allocate(Ac2d(1:Nr,1:Nz+1));allocate(Chi2d(1:Nr,1:Nz+1));allocate(Acold2d(1:Nr,1:Nz+1));allocate(Atemp2d(1:Nr,1:Nz+1))
    allocate(Ac(1:Nr));
    allocate(EL2d(0:Nr,0:Nz+1));allocate(BL2d(0:Nr,0:Nz+1)) 
    allocate(ELZ2d(0:Nr,0:Nz+1));allocate(ELZ2dold(0:Nr,0:Nz+1)) 
    allocate(EL2dold(0:Nr,0:Nz+1));allocate(BL2dold(0:Nr,0:Nz+1)) 
    allocate(Err2d(0:Nr,0:Nz+1));allocate(Ezz2d(0:Nr,0:Nz+1));allocate(Bthh2d(0:Nr,0:Nz+1))
    allocate(Err2dold(0:Nr,0:Nz+1));allocate(Ezz2dold(0:Nr,0:Nz+1));allocate(Bthh2dold(0:Nr,0:Nz+1))
 
    do n=1,Np
        r(n)=dr*(n-0.5d0)
        mass(n)=2d0*pi*r(n)*dr ! this actually is the 2-d grid size
    enddo 

    V=0d0
    gama=0.5d0*(2+V**2)
    r2d(:,1)=r;V2d=0d0;force2d(:,1)=0d0;gama2d(:,1)=gama
    psi2d=0d0;Er2d=0d0;Et2d=0d0;Bth2d=0d0;den2d=0d0;den_e2d=0d0;jr2d=0d0;jz2d=0d0
    BL2d=0d0;EL2d=0d0;Err2d=0d0;Ezz2d=0d0;Bthh2d=0d0
    BL2dold=0d0;EL2dold=0d0;Err2dold=0d0;Ezz2dold=0d0;Bthh2dold=0d0
    ELZ2d=0d0;ELZ2dold=0d0
  

    do n=1,Nz+1
    call Pulse_Static((n-1)*dz0,Ac)
    Ac2d(:,n)=Ac(:)
    end do



   a0 = 5.0d0
   Wist0 = 2.6
   k0 = 15.0
   FocalPlane = 10.61765
   zcenter = 10.61765
   FWHM = 2.0577160532


  ! pulse three
    do n=1,Nz+1
    call Pulse_Static((n-1)*dz0,Ac)
    Ac2d(:,n)=Ac2d(:,n)+Ac(:)
    end do



    call getlaserfield(Ac2d,Ac2d,Ac2d,dt0,dz0,k0)


!  test_particles
    if (TestParts) then
        
    WX_L=0d0;WZ_L=0d0;WX_W=0d0;WY_W=0d0;WZ_W=0d0

    XXT2d=0d0;YYT2d=0d0;ZZT2d=0d0;PXT2d=0d0;PYT2d=0d0;PZT2d=0d0
    call SeedTestParticles
    else 
    write(*, '(A,/)') "=====      No TestParticles     ====="
    end if

    call cpu_time(start)
    write(*, *) "===   All Varibles Are Defined   ===="
    write(*, *) "===            START             ===="
    write(*, *) "===******************************===="
    write(*,*) " "

    end subroutine set_varibles






end module initialization
