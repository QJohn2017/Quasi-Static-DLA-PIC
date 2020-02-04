module Test_Particle
use global_variables
use domain
use TestPart_var
use InAndOut

implicit none
private
save 
public::SeedTestParticles,PushTestParticlesNest,defile_TP,TestPartOut,TestParticlesdensity

        
contains

subroutine SeedTestParticles
implicit none
integer:: n,k,m,j
real(RP)::ddr,ddz,phi



j=0
ifpush=1



do m=1,Num_of_Beam 

ddz=dz0/TP_perZcell(m)
ddr=hr/TP_perRcell(m)

 do n=1,floor(TP_Rad(m)*TP_perRcell(m)/hr)*2
	do k=1,floor(TP_len(m)*TP_perZcell(m)/dz0)

  j=j+1
  beamidx(j)=m

  ZZT(j)=  TP_z0(m)-TP_len(m)*0.5d0+(k-1)*ddz
  XXT(j)= -TP_Rad(m)+ddr*0.5d0+(n-1)*ddr
  
  PXT(j)= r8_normal_01()*XSpread(m)
  PZT(j)= -sqrt(GamaT(m)**2-1d0-PXT(j)**2)
  

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

  case(7)   !trapezoidal
  denT(j)=Den_Beam(m)*polarity(m)*((ZZT(j)-TP_z0(m))/TP_len(m)*2d0*beamchirp(m)+1d0)

  end select 

	end do
 end do

end do

end subroutine SeedTestParticles





! uniform_normal_distribution_generator.
function r8_normal_01 ()
  implicit none
  real(8) :: r1,r2,r8_normal_01
  r1=rand()
  r2=rand()
  r8_normal_01= sqrt (-2.0d0*log(r1))*cos(2.0d0*pi*r2)
  return
end



subroutine PushTestParticlesNest(k0,nzc,npctp)

implicit none 
integer :: n,k,i1,i2,j1,j2,npctp,npp,indx,nzc
real(RP):: Exxp,Eyyp,Ezzp,Bxxp,Byyp,Erp,Btp,dr1,dr2,dz1,dz2,deltat,ddt,k0,TP_time
complex(RP)::ELtemp,BLtemp,ELZtemp
real(RP):: PZ_,PX_,PY_,ddd,dd1,dd2,Xtemp,Ytemp,Ztemp,gatemp,Rtemp,dentemp


ddt=dt0/N_per_dt

  
   do n=1,NTP
    
   j1=floor(ZZT(n)/dz0)+1

   if (j1/=nzc.or.ifpush(n)==0) cycle
   
   do npp=1,N_per_dt

   if (XXT(n)>=Rmax .or. XXT(n)<=-Rmax .or. ZZT(n)>=Zmax-dz0 .or. ZZT(n)<=0d0) cycle

   deltat=npp/N_per_dt
   Xtemp=XXT(n);Ztemp=ZZT(n)

   j1=floor(ZZT(n)/dz0)+1
   dz1=(ZZT(n)-dz0*(j1-1))/dz0
!--------------------------
   i1=floor((XXT(n)+Rmax)/hr)+1
   dr1=((XXT(n)+Rmax)-hr*(i1-1))/hr

   i2=floor((XXT(n)+Rmax+hr/2d0)/hr)+1
   dr2=((XXT(n)+Rmax+hr/2d0)-hr*(i2-1))/hr

   j2=floor((ZZT(n)-dz0/2d0)/dz0)+1
   dz2=((ZZT(n)-dz0/2d0)-dz0*(j2-1))/dz0


  ! write(*,*) n,XXT(n),ZZT(n),j1,i1,j2,i2
!---------------------------------------------
   ELtemp=0d0;BLtemp=0d0;ELZtemp=0d0

        !particle pusher--step-(1)
   gatemp=1d0/sqrt(1+PZT(n)**2+PXT(n)**2)

   XXT(n)=XXT(n)+0.5*ddt*PXT(n)*gatemp
   ZZT(n)=ZZT(n)+0.5*ddt*PZT(n)*gatemp+0.5*ddt 

   if (Laser) then
   ! interp laser field 
   !EL--- EX
   ELtemp =(EL2d(i2,j2)*(1-dr2)*(1-dz2)+EL2d(i2+1,j2)*(dr2)*(1-dz2)&
           +EL2d(i2,j2+1)*(1-dr2)*(dz2)+EL2d(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))
   ELtemp = ELtemp*deltat+(EL2dold(i2,j2)*(1-dr2)*(1-dz2)+EL2dold(i2+1,j2)*(dr2)*(1-dz2)&
           +EL2dold(i2,j2+1)*(1-dr2)*(dz2)+EL2dold(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))*(1-deltat)


   !BL---> BY
   BLtemp =(BL2d(i2,j2)*(1-dr2)*(1-dz2)+BL2d(i2+1,j2)*(dr2)*(1-dz2)&
           +BL2d(i2,j2+1)*(1-dr2)*(dz2)+BL2d(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))
   BLtemp = BLtemp*deltat+(BL2dold(i2,j2)*(1-dr2)*(1-dz2)+BL2dold(i2+1,j2)*(dr2)*(1-dz2)&
           +BL2dold(i2,j2+1)*(1-dr2)*(dz2)+BL2dold(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))*(1-deltat)

   ! interp field

   ELZtemp =(ELZ2d(i1,j1)*(1-dr1)*(1-dz1)+ELZ2d(i1+1,j1)*(dr1)*(1-dz1)&
           +ELZ2d(i1,j1+1)*(1-dr1)*(dz1)+ELZ2d(i1+1,j1+1)*(dr1)*(dz1))*exp(-ci*k0*ZZT(n))
   ELZtemp = ELZtemp*deltat+(ELZ2dold(i1,j1)*(1-dr1)*(1-dz1)+ELZ2dold(i1+1,j1)*(dr1)*(1-dz1)&
           +ELZ2dold(i1,j1+1)*(1-dr1)*(dz1)+ELZ2dold(i1+1,j1+1)*(dr1)*(dz1))*exp(-ci*k0*ZZT(n))*(1-deltat)

  

   if(.Not.LaserEZ) ELZtemp=0d0

   end if 

   !(Btheta for electron)
   Btp    =-(Bthh2d(i1,j1)*(1-dr1)*(1-dz1)+Bthh2d(i1+1,j1)*(dr1)*(1-dz1)&
            +Bthh2d(i1,j1+1)*(1-dr1)*(dz1)+Bthh2d(i1+1,j1+1)*(dr1)*(dz1))

   Btp    = Btp*deltat+(-(Bthh2dold(i1,j1)*(1-dr1)*(1-dz1)+Bthh2dold(i1+1,j1)*(dr1)*(1-dz1)&
            +Bthh2dold(i1,j1+1)*(1-dr1)*(dz1)+Bthh2dold(i1+1,j1+1)*(dr1)*(dz1)))*(1-deltat)


  !(E_r for electron)
   Erp    =Err2d(i1,j1)*(1-dr1)*(1-dz1)+Err2d(i1+1,j1)*(dr1)*(1-dz1)&
          +Err2d(i1,j1+1)*(1-dr1)*(dz1)+Err2d(i1+1,j1+1)*(dr1)*(dz1)

   Erp    =Erp*deltat+(Err2dold(i1,j1)*(1-dr1)*(1-dz1)+Err2dold(i1+1,j1)*(dr1)*(1-dz1)&
          +Err2dold(i1,j1+1)*(1-dr1)*(dz1)+Err2dold(i1+1,j1+1)*(dr1)*(dz1))*(1-deltat)+Btp 
   
   ! (E_xi for electron)
   Ezzp   =Ezz2d(i2,j1)*(1-dr2)*(1-dz1)+Ezz2d(i2+1,j1)*(dr2)*(1-dz1)&
          +Ezz2d(i2,j1+1)*(1-dr2)*(dz1)+Ezz2d(i2+1,j1+1)*(dr2)*(dz1)
   
   Ezzp   =Ezzp*deltat+(Ezz2dold(i2,j1)*(1-dr2)*(1-dz1)+Ezz2dold(i2+1,j1)*(dr2)*(1-dz1)&
          +Ezz2dold(i2,j1+1)*(1-dr2)*(dz1)+Ezz2dold(i2+1,j1+1)*(dr2)*(dz1))*(1-deltat)

   
   
   Exxp   = (Erp  +real(ELtemp)*2) *polarity(beamidx(n))
   Byyp   =( Btp  +real(BLtemp)*2) *polarity(beamidx(n))
   Ezzp   = (Ezzp +real(ELZtemp)*2)*polarity(beamidx(n))

    !--------------------------   
   
        !particle pusher--step-(2) 
         PX_=PXT(n)+0.5*ddt*Exxp
         PZ_=PZT(n)+0.5*ddt*Ezzp
         

        select case(pushtype)

        case(0)
        
         ddd=Byyp*ddt*0.5d0/sqrt(1+PX_**2+PZ_**2)
         dd1=2*ddd/(1+ddd**2)

         PXT(n)=PX_-dd1*ddd*PX_+dd1*PZ_
         PZT(n)=-dd1*PX_+PZ_-dd1*ddd*PZ_
        case(1)

         ddd=0.5*ddt/sqrt(1+PZ_**2+PX_**2)
         dd1=2*ddd/(1+ddd**2*(Byyp**2))
         dd2=dd1*ddd
         !particle pusher--step-(3)
         PXT(n)=(1-Byyp**2*dd2)*PX_+Byyp*dd1*PZ_
         PZT(n)=(-Byyp*dd1*PX_)+(1-Byyp**2*dd2)*PZ_
         
        end select

         !particle pusher--step-(4)
         PXT(n)=PXT(n)+0.5*ddt*Exxp
         PZT(n)=PZT(n)+0.5*ddt*Ezzp
         

         gatemp=1d0/sqrt(1+PZT(n)**2+PXT(n)**2)
   		   XXT(n)=XXT(n)+0.5*ddt*PXT(n)*gatemp
         ZZT(n)=ZZT(n)+0.5*ddt*PZT(n)*gatemp+0.5*ddt 

         WX_L(n)=WX_L(n)+(XXT(n)-Xtemp)*(real(ELtemp)*2)
         WZ_L(n)=WZ_L(n)+(ZZT(n)-Ztemp-ddt)*(real(ELZtemp)*2)
         WX_W(n)=WX_W(n)+(XXT(n)-Xtemp)*(Exxp-real(ELtemp)*2)
         WZ_W(n)=WZ_W(n)+(ZZT(n)-Ztemp-ddt)*(Ezzp-real(ELZtemp)*2)

        
        if (npp==N_per_dt) ifpush(n)=0

        indx=(npctp*N_per_dt+npp)/Tp_Dt_per_out
        
        if(Tp_Dt_per_out==1) then
        XXT2d(indx,n)=XXT(n)
        ZZT2d(indx,n)=ZZT(n);
        PXT2d(indx,n)=PXT(n);
        PZT2d(indx,n)=PZT(n);

        else

        if(indx<Dt_per_out*N_per_dt/Tp_Dt_per_out) then
        XXT2d(indx+1,n)=XXT(n)
        ZZT2d(indx+1,n)=ZZT(n);
        PXT2d(indx+1,n)=PXT(n);
        PZT2d(indx+1,n)=PZT(n);
        end if

        end if 
  
   end do

 end do


end subroutine PushTestParticlesNest


  subroutine TestParticlesdensity(nzc)
  implicit none 
  integer :: n,i2,j1,nzc
  real(RP):: dr2,dz1,dentemp,gatemp
  denb2d(:,nzc+1)=0d0;Jzb2d(:,nzc+1)=0d0
   do n=1,NTP
       if (XXT(n)>=Rmax .or. XXT(n)<=-Rmax .or. ZZT(n)>=Zmax-dz0 .or. ZZT(n)<=0d0) cycle
       j1=floor(ZZT(n)/dz0)+1
       if (j1<nzc .or.j1>nzc+1) cycle
        dentemp=1d0/TP_perZcell(beamidx(n))/TP_perRcell(beamidx(n))*denT(n)
        !dentemp=denT(n)

        i2=floor((XXT(n)+Rmax+hr/2d0)/hr)+1
        dr2=((XXT(n)+Rmax+hr/2d0)-hr*(i2-1))/hr

        dz1=(ZZT(n)-dz0*(j1-1))/dz0
        
        gatemp=1d0/sqrt(1+PZT(n)**2+PXT(n)**2)
        if(j1==nzc+1) then
        denb2d(i2,j1)    =denb2d(i2,j1)+(1-dr2)*(1-dz1)*dentemp
        denb2d(i2+1,j1)  =denb2d(i2+1,j1)+(dr2)*(1-dz1)*dentemp
        jzb2d(i2,j1)    =jzb2d(i2,j1)-(1-dr2)*(1-dz1)*dentemp*PZT(n)*gatemp
        jzb2d(i2+1,j1)  =jzb2d(i2+1,j1)-(dr2)*(1-dz1)*dentemp*PZT(n)*gatemp
        else 
        denb2d(i2,j1+1)  =denb2d(i2,j1+1)+(1-dr2)*(dz1)*dentemp
        denb2d(i2+1,j1+1)=denb2d(i2+1,j1+1)+(dr2)*(dz1)*dentemp
        jzb2d(i2,j1+1)  =jzb2d(i2,j1+1)-(1-dr2)*(dz1)*dentemp*PZT(n)*gatemp
        jzb2d(i2+1,j1+1)=jzb2d(i2+1,j1+1)-(dr2)*(dz1)*dentemp*PZT(n)*gatemp
        end if 
  end do 


end subroutine TestParticlesdensity

subroutine defile_TP(outct)
implicit none
integer :: outct
character(12):: filename
write(filename,'(I5)'), outct
call check(nf90_create("output/"//TRIM(ADJUSTL(filename))//"TestParts.nc", NF90_CLOBBER, ncidtp) )
call check( nf90_def_dim(ncidtp, "time", Dt_per_out*N_per_dt/Tp_Dt_per_out,ztime_dim) )
call check( nf90_def_dim(ncidtp, "Ntp", Ntp, p_dim) )
dims=(/ztime_dim,p_dim/)
call check( nf90_def_var(ncidtp, "TPxx", NF90_REAL, dims, var_Tx) )
call check( nf90_def_var(ncidtp, "TPzz", NF90_REAL, dims, var_Tz) )
call check( nf90_def_var(ncidtp, "TPpx", NF90_REAL, dims, var_Tpx)) 
call check( nf90_def_var(ncidtp, "TPpz", NF90_REAL, dims, var_Tpz)) 

call check( nf90_def_dim(ncidtp, "LastT", 1,ztime_dim) )
dims=(/p_dim,ztime_dim/)
call check( nf90_def_var(ncidtp, "WXll", NF90_REAL, dims, var_Wxl)) 
call check( nf90_def_var(ncidtp, "WZll", NF90_REAL, dims, var_Wzl)) 
call check( nf90_def_var(ncidtp, "WXww", NF90_REAL, dims, var_Wxw)) 
call check( nf90_def_var(ncidtp, "WZww", NF90_REAL, dims, var_Wzw)) 
call check( nf90_enddef(ncidtp) )

end subroutine defile_TP


subroutine TestPartOut
  implicit none 
    call check(nf90_put_var(ncidtp,var_Tx,XXT2D))
    call check(nf90_put_var(ncidtp,var_Tz,ZZT2D))
    call check(nf90_put_var(ncidtp,var_Tpx,PXT2D))
    call check(nf90_put_var(ncidtp,var_Tpz,-PZT2D))

    call check(nf90_put_var(ncidtp,var_Wxl,WX_L))
    call check(nf90_put_var(ncidtp,var_Wzl,WZ_L))
    call check(nf90_put_var(ncidtp,var_Wxw,WX_W))
    call check(nf90_put_var(ncidtp,var_Wzw,WZ_W))

    call check( nf90_close(ncidtp)) 
end subroutine TestPartOut





end module Test_Particle
