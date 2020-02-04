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
        use OneTstep

        implicit none
        private
        save
        public:: set_varibles
        contains

    subroutine set_varibles(start,Centri_old)
    implicit none
    real(RP):: start,Centri_old
    complex(Rp),allocatable:: Ac(:)
    integer:: n,NTPall
    complex(RP),allocatable:: BL(:)
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

    Np=(2*Nr+1)*Part_per_cell  
    
    Nr=Nr*2+2  ! planar Geometry

    dz0=Zmax/Nz
    !NWake=Int(WakeLen/dz0)
    NTP_i=0
    do n=1,Num_of_Beam
    NTP_i(n)=floor(TP_len(n)*TP_perZcell(n)/dz0)*floor(TP_Rad(n)*TP_perRcell(n)/hr)*2
    end do
    NTP=sum(NTP_i)
    
    if (TestParts) then
    write(*, *) "=====Number of Test Beams   ==="
    write(*, *) Num_of_Beam
    write(*, *) "=====Number of TestParticles==="
open(222,file='output/Num_of_Beam',status="replace",action="write")
    write(*, *) '(',NTP_i(1:Num_of_Beam),')'
    write(222, *) NTP_i(1:Num_of_Beam)
close(222)

    allocate(XXT(1:NTP));allocate(ZZT(1:NTP))
    allocate(PXT(1:NTP));allocate(PZT(1:NTP))
    allocate(WX_L(1:NTP));allocate(WZ_L(1:NTP))
    allocate(WX_W(1:NTP));allocate(WZ_W(1:NTP))
    allocate(beamidx(1:NTP));allocate(denT(1:NTP))

    allocate(ifpush(1:NTP))

    allocate(XXT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(ZZT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(PXT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    allocate(PZT2d(1:Dt_per_out*N_per_dt/Tp_Dt_per_out,1:NTP))
    end if


    allocate(r(1:Np));allocate(r0(1:Np));allocate(V(1:Np));allocate(force(1:Np));allocate(mass(1:Np));allocate(gama(1:Np))
   
    allocate(psi(1:Nr));allocate(Er(1:Nr));allocate(Et(1:Nr));allocate(Bth(1:Nr))
    allocate(den(1:Nr));allocate(den_e(1:Nr));allocate(jr(1:Nr));allocate(jv(1:Nr));allocate(jz(1:Nr))
    allocate(Ac(1:Nr));  allocate(BL(0:Nz+1));


    allocate(r2d(1:Np,1:Nz+1));allocate(V2d(1:Np,1:Nz+1));allocate(force2d(1:Np,1:Nz+1));allocate(gama2d(1:Np,1:Nz+1))
    
    allocate(denb2d(1:Nr,1:Nz+1));allocate(jzb2d(1:Nr,1:Nz+1));

    allocate(psi2d(1:Nr,1:Nz+1));allocate(Er2d(1:Nr,1:Nz+1));   allocate(Et2d(1:Nr,1:Nz+1));   allocate(Bth2d(1:Nr,1:Nz+1))
    allocate(den2d(1:Nr,1:Nz+1));allocate(den_e2d(1:Nr,1:Nz+1));allocate(jr2d(1:Nr,1:Nz+1));   allocate(jz2d(1:Nr,1:Nz+1))
    allocate(Ac2d(1:Nr,1:Nz+1)); allocate(Chi2d(1:Nr,1:Nz+1));  allocate(Acold2d(1:Nr,1:Nz+1));allocate(Atemp2d(1:Nr,1:Nz+1))
  
    allocate(EL2d(1:Nr,0:Nz+1));allocate(BL2d(1:Nr,0:Nz+1))
    allocate(ELZ2d(1:Nr,0:Nz+1));allocate(ELZ2dold(1:Nr,0:Nz+1)) 
    allocate(EL2dold(1:Nr,0:Nz+1));allocate(BL2dold(1:Nr,0:Nz+1)) 
    allocate(Err2d(1:Nr,0:Nz+1));allocate(Ezz2d(1:Nr,0:Nz+1));allocate(Bthh2d(1:Nr,0:Nz+1))
    allocate(Err2dold(1:Nr,0:Nz+1));allocate(Ezz2dold(1:Nr,0:Nz+1));allocate(Bthh2dold(1:Nr,0:Nz+1))
 

    do n=1,Np
        r(n)=-Rmax-Part_per_cell*dr/2+dr*(n-1)
    enddo 


    r0=r;
    mass=dr ! this actually is the 2-d grid size
    
    V=0d0
    gama=0.5d0*(2+V**2)
    r2d(:,1)=r;V2d=0d0;force2d(:,1)=0d0;gama2d(:,1)=gama
    psi2d=0d0;Er2d=0d0;Et2d=0d0;Bth2d=0d0;den2d=0d0;denb2d=0d0;jzb2d=0d0;den_e2d=0d0;jr2d=0d0;jz2d=0d0
    BL2d=0d0;EL2d=0d0;Err2d=0d0;Ezz2d=0d0;Bthh2d=0d0
    BL2dold=0d0;EL2dold=0d0;Err2dold=0d0;Ezz2dold=0d0;Bthh2dold=0d0
    ELZ2d=0d0;ELZ2dold=0d0
   




    if(Laser) then
    Ac2d=0d0
    Ac=0d0
    do n=1,Nz+1
    call Pulse_Static((n-1)*dz0,Ac)
    Ac2d(:,n)=Ac(:)
    end do


    do n=1,Nz
    BL2d(1:Nr,n)=ci*k0*Ac2d(:,n)-(Ac2d(:,n+1)-Ac2d(:,n))/dz0
    EL2d(1:Nr,n)=BL2d(1:Nr,n)
    end do 
    BL(0:Nz+1)=BL2d(Nr/2,:)
    Centri_old=LaserCentriod(BL)
    end if 


    if (TestParts) then
    WX_L=0d0;WX_W=0d0;WZ_W=0d0;WZ_L=0d0;
    XXT2d=0d0;ZZT2d=0d0;PXT2d=0d0;PZT2d=0d0
    call SeedTestParticles
    else 
    write(*, '(A,/)') "=====      No TestParticles     ====="
    end if



    if(Restart>0) then 

        if(Laser) then
        call ReadLaserPulse(Restart)
        do n=1,Nz
        BL2d(1:Nr,n)=ci*k0*Ac2d(:,n)-(Ac2d(:,n+1)-Ac2d(:,n))/dz0
        EL2d(1:Nr,n)=BL2d(1:Nr,n)
        end do 
        BL(0:Nz+1)=BL2d(Nr/2,:)
        Centri_old=LaserCentriod(BL)
        end if
        if(TestParts) then
        call ReadTestParts(Restart)
        end if

    write(*, *) "===          Restart             ===="
    write(*, *) "===******************************===="
    end if 
    



    call cpu_time(start)
    write(*, *) "===   All Varibles Are Defined   ===="
    write(*, *) "===            START             ===="
    write(*, *) "===******************************===="
    write(*,*) " "

    end subroutine set_varibles

    subroutine ReadTestParts(Restart)
    integer :: m,n,k,j,Restart, ncidtp,var_Tx,var_Tz,var_Tpx,var_Tpz,var_Wxl,var_Wzl,var_Wxw,var_Wzw
    character(18) ::filename
    real(8) ::ddz,ddr
    write(filename,'(I5)'), Restart

! read ini weight of particles
    call check(nf90_open("output/1TestParts.nc", nf90_nowrite, ncidtp))
    call check(nf90_inq_varid(ncidtp, "TPxx", var_Tx))
    call check(nf90_inq_varid(ncidtp, "TPzz", var_Tz))
    call check(nf90_get_var(ncidtp, var_Tx,  XXT(1:NTP_i(1)),start=(/1,1/),count=(/1,NTP_i(1)/)))
    call check(nf90_get_var(ncidtp, var_Tz,  ZZT(1:NTP_i(1)),start=(/1,1/),count=(/1,NTP_i(1)/)))
 j=0
do m=1,1
ddz=dz0/TP_perZcell(m)
ddr=hr/TP_perRcell(m)
 do n=1,floor(TP_Rad(m)*TP_perRcell(m)/hr)*2
    do k=1,floor(TP_len(m)*TP_perZcell(m)/dz0)
  j=j+1
  select case(Beam_Type(m))
  case(1)      ! rectangle
  denT(j)=Den_Beam(m)*polarity(m)
  case(2)    ! triangle
  denT(j)=Den_Beam(m)*polarity(m)*(ZZT(j)-TP_z0(m)+TP_len(m)*0.5)/TP_len(m)
  case(3)    ! rectangle with cos in front
  if  ((ZZT(j)-TP_z0(m)+TP_len(m)*0.5)>TP_PlateauBe(m)) then
  denT(j)=Den_Beam(m)*polarity(m)
  else 
  denT(j)=Den_Beam(m)*polarity(m)*cos((TP_PlateauBe(m)-(ZZT(j)-TP_z0(m)+TP_len(m)*0.5))/TP_PlateauBe(m)*pi/2d0)
  end if 
  case(4)    ! Reverse triangle
  denT(j)=Den_Beam(m)*polarity(m)*(-ZZT(j)+TP_z0(m)+TP_len(m)*0.5)/TP_len(m)
  case(5)    ! Longitudinal_Gaussian
  denT(j)=Den_Beam(m)*polarity(m)*exp(-2.77258872*(ZZT(j)-TP_z0(m))**2/TP_FWHM(m)**2)
  case(6)   !Gaussian
  denT(j)=Den_Beam(m)*polarity(m)*exp(-2.77258872*(ZZT(j)-TP_z0(m))**2/TP_FWHM(m)**2)*exp(-XXT(j)**2/TP_w0(m)**2)
  end select 
    end do
 end do
end do


    call check(nf90_open("output/"//TRIM(ADJUSTL(filename))//"TestParts.nc", nf90_nowrite, ncidtp))

    call check(nf90_inq_varid(ncidtp, "TPxx", var_Tx))
    call check(nf90_inq_varid(ncidtp, "TPzz", var_Tz))
    call check(nf90_inq_varid(ncidtp, "TPpz", var_Tpz))
    call check(nf90_inq_varid(ncidtp, "TPpx", var_Tpx))

    call check(nf90_inq_varid(ncidtp, "WXll", var_Wxl))
    call check(nf90_inq_varid(ncidtp, "WZll", var_Wzl))
    call check(nf90_inq_varid(ncidtp, "WZww", var_Wzw))
    call check(nf90_inq_varid(ncidtp, "WXww", var_Wxw))

    call check(nf90_get_var(ncidtp, var_Tx,  XXT(1:NTP_i(1)),start=(/Dt_per_out*N_per_dt/Tp_Dt_per_out,1/),count=(/1,NTP_i(1)/)))
    call check(nf90_get_var(ncidtp, var_Tz,  ZZT(1:NTP_i(1)),start=(/Dt_per_out*N_per_dt/Tp_Dt_per_out,1/),count=(/1,NTP_i(1)/)))
    call check(nf90_get_var(ncidtp, var_Tpx, PXT(1:NTP_i(1)),start=(/Dt_per_out*N_per_dt/Tp_Dt_per_out,1/),count=(/1,NTP_i(1)/)))
    call check(nf90_get_var(ncidtp, var_Tpz, PZT(1:NTP_i(1)),start=(/Dt_per_out*N_per_dt/Tp_Dt_per_out,1/),count=(/1,NTP_i(1)/)))

    call check(nf90_get_var(ncidtp, var_Wxl, WX_L(1:NTP_i(1)),start=(/1,1/),count=(/NTP_i(1),1/)))
    call check(nf90_get_var(ncidtp, var_Wzl, WZ_L(1:NTP_i(1)),start=(/1,1/),count=(/NTP_i(1),1/)))
    call check(nf90_get_var(ncidtp, var_Wxw, WX_W(1:NTP_i(1)),start=(/1,1/),count=(/NTP_i(1),1/)))
    call check(nf90_get_var(ncidtp, var_Wzw, WZ_W(1:NTP_i(1)),start=(/1,1/),count=(/NTP_i(1),1/)))

    call check( nf90_close(ncidtp)) 

    PZT(1:NTP_i(1))=-PZT(1:NTP_i(1))

    end subroutine ReadTestParts



    subroutine ReadLaserPulse(Restart)
    integer :: Restart, ncidf,var_AR,var_AI
    real(RP) ::Ac2dReadR(1:Nr,1:Nz+1),Ac2dReadI(1:Nr,1:Nz+1)
    character(18) ::filename
    write(filename,'(I5)'), Restart
    call check(nf90_open("output/"//TRIM(ADJUSTL(filename))//"Fields.nc", nf90_nowrite, ncidf))
    call check(nf90_inq_varid(ncidf, "Aaa_R", var_AR))
    call check(nf90_inq_varid(ncidf, "Aaa_I", var_AI))

    call check(nf90_get_var(ncidf, var_AR,  Ac2dReadR(1:Nr,1:Nz),start=(/1,1/),count=(/Nr,Nz/)))
    call check(nf90_get_var(ncidf, var_AI,  Ac2dReadI(1:Nr,1:Nz),start=(/1,1/),count=(/Nr,Nz/)))
    call check( nf90_close(ncidf)) 

    Ac2d=Ac2dReadR*0.5d0+ci*Ac2dReadI*0.5d0

    end subroutine ReadLaserPulse


end module initialization
