
        module  OneZstep
        use global_variables
        use wake_f
        use magnetic_f
        use electron_pusher
        use particles_var
        use fields_var
        use domain
        use Functions
   
        implicit none
        private
        save 
        public::one_step_z,resetplasma

        
        contains

        subroutine one_step_z(ztime,Ac,denb1d,jzb1d,Chi,nzc,f_cross,time) ! chi is the only output here because I need only chi for laser advance.
        implicit none 
        integer :: k,f_cross,nzc,n
        complex(Rp) :: Ac(1:Nr)
        real(RP):: v1(1:Np),r1(1:Np)
        real(RP):: f1(1:Np),fv1(1:Np),f2(1:Np),fv2(1:Np),Chi(1:Nr)
        real(RP):: As(1:Nr),ponderF(1:Nr),rhof(1:Nr),denb1d(1:Nr),jzb1d(1:Nr)
        real(RP):: ztime,dz,dentemp,time
        
        do k=1,Nr
        As(k)=real(Ac(k)*conjg(Ac(k)))
        end do

        dentemp=1d0
        if (time<=rampBegin) then
        dentemp=0d0
        else if (time<=rampBegin+ramp_len) then
        dentemp=(time-rampBegin)/ramp_len
        else 
        dentemp=1d0
        end if 


        dz=dz0
        

  	do  n=1,Np
        if(abs(V(n))*dz0>10d0*hr)  V(n)=sign(10d0*hr/dz0,V(n))
        end do 

        if (maxval(abs(V))*dz0>hr) then
        dz=hr/maxval(abs(V))/Push_per_cell*1d0
        if ((ztime+dz)>(nzc+1)*dz0)  dz=(nzc+1)*dz0-ztime
        end if 



        if (Push_Order==1) r=r+V*dz*0.5d0
            
        call electron_density_r(r,mass,dentemp,den)
        call electron_current_jr(r,V,mass,dentemp,jr)
        call wake_fields(den,denb1d,jzb1d,Er,Psi,dentemp)
        call wake_field_Et(jr,Et)
        call electron_current_jv(r,V,mass,dentemp,jv)
        call electron_current_jz(jv,den,As,Psi,jz,ponderF)
        call force_rhoF(Psi,den,Et,Er,jv,jr,As,ponderF,rhoF)
        call magnetic_field(jv,jzb1d,rhoF,jz,psi,den,Bth) 
        call electron_push(r0,r,V,gama,fv1,f1,Er,Et,Bth,ponderF,Psi,As)
        
        if (Push_Order==1) then
            V=V+f1*dz
            r=r+V*dz*0.5d0
        else 
            r1=r+fv1*dz*0.5d0
            V1=V+ f1*dz*0.5d0
        end if 


        if (Push_Order/=1) then
        call electron_density_r(r1,mass,dentemp,den)
        call electron_current_jr(r1,V1,mass,dentemp,jr)
        call wake_fields(den,denb1d,jzb1d,Er,Psi,dentemp)
        call wake_field_Et(jr,Et)
        call electron_current_jv(r1,V1,mass,dentemp,jv)
        call electron_current_jz(jv,den,As,Psi,jz,ponderF)
        call force_rhoF(Psi,den,Et,Er,jv,jr,As,ponderF,rhoF)
        call magnetic_field(jv,jzb1d,rhoF,jz,psi,den,Bth) 
        call electron_push(r0,r1,V1,gama,fv2,f2,Er,Et,Bth,ponderF,Psi,As)
       
        r=r+fv2*dz
        V=V+f2*dz
        end if 


        call electron_density_e(jv,psi,den,As,den_e)
        call Susceptibility(den,Psi,Chi)
        force(:)=f2
        ztime=ztime+dz
        nzc=(ztime+1e-8)/dz0
       ! write(*,*) ztime,nzc
    
        ! if(nzc==244)  then
        ! open(999,file='output/debug')
        ! write(999,*) rhoF
        ! end if 

        end subroutine one_step_z



        subroutine resetplasma(ztime,Chi0,nzc,nzc0,f_cross)
        implicit none
        integer:: n,nzc,nzc0,f_cross
        real(RP) :: ztime,Chi0(1:Nr)
        ztime=0d0

        do n=1,Np
        r(n)=-Rmax-Part_per_cell*dr/2+dr*(n-1)
        enddo 

        V=0d0
        gama=1d0
        psi=0d0;Er=0d0;Et=0d0;Bth=0d0;den=0d0;den_e=0d0;jr=0d0;jz=0d0
        nzc=0
        nzc0=0
        f_cross=0
        Chi0(:)=1d0
        end subroutine resetplasma





  


end module OneZstep
