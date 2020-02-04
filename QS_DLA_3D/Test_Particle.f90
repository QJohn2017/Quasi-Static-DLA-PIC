module Test_Particle
use global_variables
use domain
use TestPart_var
use InAndOut

implicit none
private
save 
public::SeedTestParticles,PushTestParticles,defile_TP,TestPartOut

        
contains

subroutine SeedTestParticles
implicit none
integer:: n,k,l,j,m
real(RP)::ddr,ddz,phi
j=0

do m=1,Num_of_Beam

ddz=1d0/TP_perZ(m);ddr=1d0/TP_perR(m)

 do n=1,floor(TP_Rad(m)*TP_perR(m))
  do l=1,n*N_per_peri(m)
	do k=1,floor(TP_len(m)*TP_perZ(m))
    j=j+1
    phi=2*pi/N_per_peri(m)/n*(l-1)
	XXT(j)=(ddr*0.5+(n-1)*ddr)*cos(phi)
  YYT(j)=(ddr*0.5+(n-1)*ddr)*sin(phi)
	ZZT(j)=TP_z0(m)-TP_len(m)*0.5+(k-1)*ddz
  PZT(j)=-sqrt(GamaT(m)**2-1d0)
      end do
	end do
 end do

end do

PXT(:)=0d0;PYT(:)=0d0;
	
end subroutine SeedTestParticles



subroutine PushTestParticles(k0,TP_time,npp,npctp)

implicit none 
integer :: n,k,i1,i2,j1,j2,outct,npp,indx,npctp
real(RP):: Exxp,Eyyp,Ezzp,Bxxp,Byyp,Erp,Btp,dr1,dr2,dz1,dz2,ddt,deltat,k0,TP_time
complex(RP)::ELtemp,BLtemp,ELZtemp
real(RP):: PZ_,PX_,PY_,ddd,dddx,dddy,dd1,dd2,dd1x,dd1y,Xtemp,Ytemp,Ztemp,gatemp,Rtemp

ddt=dt0/N_per_dt
deltat=npp/N_per_dt

Err2d(0,:) =Err2d(1,:); Ezz2d(0,:) =Ezz2d(1,:);  
Bthh2d(0,:)=Bthh2d(1,:); EL2d(0,:)  =EL2d(1,:);   
ELZ2d(0,:)  =ELZ2d(1,:);   
BL2d(0,:)  =BL2d(1,:);  

   do n=1,NTP

    Xtemp=XXT(n);Ytemp=YYT(n);Ztemp=ZZT(n)

   !particle pusher--step-(1)
    gatemp=1d0/sqrt(1+PZT(n)**2+PYT(n)**2+PXT(n)**2)

    XXT(n)=XXT(n)+0.5*ddt*PXT(n)*gatemp
    YYT(n)=YYT(n)+0.5*ddt*PYT(n)*gatemp
    ZZT(n)=ZZT(n)+0.5*ddt*PZT(n)*gatemp+0.5*ddt 
 
    Rtemp=sqrt(XXT(n)**2+YYT(n)**2)

   if (Rtemp>=Rmax-hr .or. ZZT(n)>=Zmax-dz0 .or. ZZT(n)<=0d0) cycle


!--------------------------
   i1=floor(Rtemp/hr)+1
   dr1=(Rtemp-hr*(i1-1))/hr

   j1=floor(ZZT(n)/dz0)+1
   dz1=(ZZT(n)-dz0*(j1-1))/dz0

   i2=floor((Rtemp-hr/2d0)/hr)+1
   dr2=((Rtemp-hr/2d0)-hr*(i2-1))/hr


   j2=floor((ZZT(n)-dz0/2d0)/dz0)+1
   dz2=((ZZT(n)-dz0/2d0)-dz0*(j2-1))/dz0
!---------------------------------------------


   ! interp laser field 
   !EL--- EX
   ELtemp =(EL2d(i2,j2)*(1-dr2)*(1-dz2)+EL2d(i2+1,j2)*(dr2)*(1-dz2)&
           +EL2d(i2,j2+1)*(1-dr2)*(dz2)+EL2d(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))
   
   ELtemp = ELtemp*deltat+(EL2dold(i2,j2)*(1-dr2)*(1-dz2)+EL2dold(i2+1,j2)*(dr2)*(1-dz2)&
           +EL2dold(i2,j2+1)*(1-dr2)*(dz2)+EL2dold(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))*(1-deltat)

   !BL---> BY
   BLtemp =(BL2d(i2,j2)*(1-dr2)*(1-dz2)+BL2d(i2+1,j2)*(dr2)*(1-dz2)&
           +BL2d(i2,j2+1)*(1-dr2)*(dz2)+BL2d(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))
   ! time interp field
   BLtemp = BLtemp*deltat+(BL2dold(i2,j2)*(1-dr2)*(1-dz2)+BL2dold(i2+1,j2)*(dr2)*(1-dz2)&
           +BL2dold(i2,j2+1)*(1-dr2)*(dz2)+BL2dold(i2+1,j2+1)*(dr2)*(dz2))*exp(-ci*k0*ZZT(n))*(1-deltat)

   !EL--- EZ

   ELZtemp =(ELZ2d(i1,j1)*(1-dr1)*(1-dz1)+ELZ2d(i1+1,j1)*(dr1)*(1-dz1)&
           +ELZ2d(i1,j1+1)*(1-dr1)*(dz1)+ELZ2d(i1+1,j1+1)*(dr1)*(dz1))*exp(-ci*k0*ZZT(n))
   
   ELZtemp =ELZtemp*deltat+(ELZ2dold(i1,j1)*(1-dr1)*(1-dz1)+ELZ2dold(i1+1,j1)*(dr1)*(1-dz1)&
           +ELZ2dold(i1,j1+1)*(1-dr1)*(dz1)+ELZ2dold(i1+1,j1+1)*(dr1)*(dz1))*exp(-ci*k0*ZZT(n))*(1-deltat)

   if(.Not.LaserEZ) ELZtemp=0d0
   
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

   
   Exxp=Erp*XXT(n)/Rtemp+real(ELtemp)*2;
   Eyyp=Erp*YYT(n)/Rtemp
   Bxxp=-Btp*YYT(n)/Rtemp; 
   Byyp=Btp*XXT(n)/Rtemp+real(BLtemp)*2
   Ezzp=Ezzp+real(ELZtemp)*2*XXT(n)/Rtemp

    !--------------------------   
   
        !particle pusher--step-(2) 
         PX_=PXT(n)+0.5*ddt*Exxp
         PY_=PYT(n)+0.5*ddt*Eyyp
         PZ_=PZT(n)+0.5*ddt*Ezzp
         

        select case(pushtype)
        case(0)

        dddx=Bxxp*ddt*0.5d0/sqrt(1+PX_**2+PY_**2+PZ_**2)
        dddy=Byyp*ddt*0.5d0/sqrt(1+PX_**2+PY_**2+PZ_**2)
        
        dd1x=2*dddx/(1+dddx**2+dddy**2)
        dd1y=2*dddy/(1+dddx**2+dddy**2)

        PXT(n)=(1-dd1y*dddy)*PX_+dd1y*dddx*PY_+dd1y*PZ_
        PYT(n)=dd1x*dddy*PX_+(1-dd1x*dddx)*PY_-dd1x*PZ_
        PZT(n)=-dd1y*PX_+dd1x*PY_-(dd1x*dddx+dd1y*dddy-1)*PZ_

        case(1)
         ddd=0.5*ddt/sqrt(1+PZ_**2+PY_**2+PX_**2)
         dd1=2*ddd/(1+ddd**2*(Bxxp**2+Byyp**2))
         dd2=dd1*ddd
         !particle pusher--step-(3)
         PXT(n)=(1-Byyp**2*dd2)*PX_+Bxxp*Byyp*dd2*PY_+Byyp*dd1*PZ_
         PYT(n)=(Bxxp*Byyp*dd2*PX_)+(1-Bxxp**2*dd2)*PY_-Bxxp*dd1*PZ_
         PZT(n)=-(Byyp*dd1*PX_)+Bxxp*dd1*PY_+(1+(-Bxxp**2-Byyp**2)*dd2)*PZ_
        end select


         !particle pusher--step-(4)
         PXT(n)=PXT(n)+0.5*ddt*Exxp
         PYT(n)=PYT(n)+0.5*ddt*Eyyp
         PZT(n)=PZT(n)+0.5*ddt*Ezzp
         

         !particle pusher--step-(5)

         gatemp=1d0/sqrt(1+PZT(n)**2+PYT(n)**2+PXT(n)**2)
   		   XXT(n)=XXT(n)+0.5*ddt*PXT(n)*gatemp
         YYT(n)=YYT(n)+0.5*ddt*PYT(n)*gatemp
         ZZT(n)=ZZT(n)+0.5*ddt*PZT(n)*gatemp+0.5*ddt 



        
        ! work done of the fields
        WX_L(n)=WX_L(n)+(XXT(n)-Xtemp)*(real(ELtemp)*2)
        WZ_L(n)=WZ_L(n)+(ZZT(n)-Ztemp-ddt)*(real(ELZtemp)*(XXT(n)+Xtemp)/Rtemp)
        WX_W(n)=WX_W(n)+(XXT(n)-Xtemp)*(Exxp-real(ELtemp)*2)
        WY_W(n)=WY_W(n)+(YYT(n)-Ytemp)*(Eyyp)
        WZ_W(n)=WZ_W(n)+(ZZT(n)-Ztemp-ddt)*(Ezzp-(real(ELZtemp)*(XXT(n)+Xtemp)/Rtemp))



        indx=(npctp*N_per_dt+npp)/Tp_Dt_per_out
        if(Tp_Dt_per_out==1) then
        XXT2d(indx,n)=XXT(n)
        YYT2d(indx,n)=YYT(n)
        ZZT2d(indx,n)=ZZT(n);
        PXT2d(indx,n)=PXT(n);
        PYT2d(indx,n)=PYT(n);
        PZT2d(indx,n)=PZT(n);
        else
        if(indx<Dt_per_out*N_per_dt/Tp_Dt_per_out) then
        XXT2d(indx+1,n)=XXT(n)
        YYT2d(indx+1,n)=YYT(n)
        ZZT2d(indx+1,n)=ZZT(n);
        PXT2d(indx+1,n)=PXT(n);
        PYT2d(indx+1,n)=PYT(n);
        PZT2d(indx+1,n)=PZT(n);
        end if
        end if 
       

   end do

 TP_time=TP_time+ddt

end subroutine PushTestParticles

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
call check( nf90_def_var(ncidtp, "TPyy", NF90_REAL, dims, var_Ty) )
call check( nf90_def_var(ncidtp, "TPzz", NF90_REAL, dims, var_Tz) )
call check( nf90_def_var(ncidtp, "TPpx", NF90_REAL, dims, var_Tpx) ) 
call check( nf90_def_var(ncidtp, "TPpy", NF90_REAL, dims, var_Tpy) ) 
call check( nf90_def_var(ncidtp, "TPpz", NF90_REAL, dims, var_Tpz) ) 
call check( nf90_def_dim(ncidtp, "LastT", 1,ztime_dim) )
dims=(/p_dim,ztime_dim/)
call check( nf90_def_var(ncidtp, "WXll", NF90_REAL, dims, var_Wxl) ) 
call check( nf90_def_var(ncidtp, "WZll", NF90_REAL, dims, var_Wzl)) 
call check( nf90_def_var(ncidtp, "WXww", NF90_REAL, dims, var_Wxw) ) 
call check( nf90_def_var(ncidtp, "WYww", NF90_REAL, dims, var_Wyw) ) 
call check( nf90_def_var(ncidtp, "WZww", NF90_REAL, dims, var_Wzw) ) 
call check( nf90_enddef(ncidtp) )

end subroutine defile_TP


subroutine TestPartOut
	implicit none 
	integer :: npc
	  call check(nf90_put_var(ncidtp,var_Tx, XXT2D))
    call check(nf90_put_var(ncidtp,var_Ty, YYT2D))
    call check(nf90_put_var(ncidtp,var_Tz, ZZT2D))
    call check(nf90_put_var(ncidtp,var_Tpx,PXT2D))
    call check(nf90_put_var(ncidtp,var_Tpy,PYT2D))
    call check(nf90_put_var(ncidtp,var_Tpz,-PZT2D))

    call check(nf90_put_var(ncidtp,var_Wxl,WX_L))
    call check(nf90_put_var(ncidtp,var_Wzl,WZ_L))
    call check(nf90_put_var(ncidtp,var_Wxw,WX_W))
    call check(nf90_put_var(ncidtp,var_Wyw,WY_W))
    call check(nf90_put_var(ncidtp,var_Wzw,WZ_W))
    call check(nf90_close(ncidtp))
    
end subroutine TestPartOut




end module Test_Particle
