
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

        subroutine one_step_z(ztime,Ac,Chi,nzc,f_cross,time) ! chi is the only output here because I need only chi for laser advance.
        implicit none 
        integer :: k,f_cross,nzc,n
        complex(Rp) :: Ac(1:Nr)
        real(RP):: v1(1:Np),r1(1:Np)
        real(RP):: f1(1:Np),fv1(1:Np),f2(1:Np),fv2(1:Np),Chi(1:Nr+1)
        real(RP):: As(1:Nr),ponderF(1:Nr+1),rhof(1:Nr+1)
        real(RP):: ztime,dz,dentemp,time
        do k=1,Nr
        As(k)=real(Ac(k)*conjg(Ac(k)))
        end do

        dz=dz0
       
        dentemp=1d0
        if (time<=rampBegin) then
            dentemp=0d0
        else if (time<=rampBegin+ramp_len) then
        dentemp=(time-rampBegin)/ramp_len
        else 
        dentemp=1d0
        end if 
        

        if (abs(minval(V))*dz0>hr .and. minval(V)<0) then
        dz=hr/abs(minval(V))/Push_per_cell*1d0
        f_cross=1
        if ((ztime+dz)>(nzc+1)*dz0)  dz=(nzc+1)*dz0-ztime
        end if 

        if (f_cross==1) then
         if (max(FastestInDomain(V,r),abs(minval(V)))*dz0>hr) then
           dz=hr/max(FastestInDomain(V,r),abs(minval(V)))/Push_per_cell*1d0
           if ((ztime+dz)>(nzc+1)*dz0)  dz=(nzc+1)*dz0-ztime
         else 
            dz=(nzc+1)*dz0-ztime
            f_cross=0
         end if 
        end if   


        if (Push_Order==1) then
        r=r+V*dz*0.5d0
            do n=1,Np
            if(r(n)<=hr/2d0.and.V(n)<0)then
                r(n)=abs(r(n))
                V(n)=abs(V(n))
            endif 
            enddo
        end if 

        call electron_density_r(r,mass,dentemp,den)
        call electron_current_jr(r,V,mass,dentemp,jr)
        call wake_fields(den,Er,Psi,dentemp)
        call wake_field_Et(jr,Et)
        call electron_current_jv(r,V,mass,dentemp,jv)
        call electron_current_jz(jv,den,As,Psi,jz,ponderF)
        call force_rhoF(Psi,den,Et,Er,jv,jr,As,ponderF,rhoF)
        call magnetic_field(jv,rhoF,jz,psi,den,Bth) 
        call electron_push(r,V,gama,fv1,f1,Er,Et,Bth,ponderF,Psi,As)
        
        if (Push_Order==1) then

            V=V+f1*dz
            r=r+V*dz*0.5d0
            
            do n=1,Np
            if(r(n)<=hr/2d0.and.V(n)<0)then
                r(n)=abs(r(n))
                V(n)=abs(V(n))
            endif 
            enddo

        else 
        
            r1=r+fv1*dz*0.5d0
            V1=V+ f1*dz*0.5d0
       ! end do
!
            do n=1,Np
            if(r1(n)<=hr/2d0.and.V1(n)<0)then
                r1(n)=abs(r1(n))
                V1(n)=abs(V1(n))  
            endif

            if(V1(n)>hr/dz) V1(n)=hr/dz

            enddo

        end if 

    

        if (Push_Order/=1) then
        call electron_density_r(r1,mass,dentemp,den)
        call electron_current_jr(r1,V1,mass,dentemp,jr)
        call wake_fields(den,Er,Psi,dentemp)
        call wake_field_Et(jr,Et)
        call electron_current_jv(r1,V1,mass,dentemp,jv)
        call electron_current_jz(jv,den,As,Psi,jz,ponderF)
        call force_rhoF(Psi,den,Et,Er,jv,jr,As,ponderF,rhoF)
        call magnetic_field(jv,rhof,jz,psi,den,Bth)        
        call electron_push(r1,V1,gama,fv2,f2,Er,Et,Bth,ponderF,Psi,As)

        r=r+fv2*dz
        V=V+f2*dz

        do n=1,Np
        if(r(n)<=hr/2d0.and.V(n)<0)then
                r(n)=abs(r(n))
                V(n)=abs(V(n))
            endif 
        enddo
        
        end if 


        call electron_density_e(jv,psi,den,As,den_e)
        call Susceptibility(den,Psi,Chi)
        force(:)=f2
        ztime=ztime+dz
        nzc=(ztime+1e-8)/dz0

        end subroutine one_step_z



        subroutine resetplasma(ztime,Chi0,nzc,nzc0,f_cross)
        implicit none
        integer:: n,nzc,nzc0,f_cross
        real(RP) :: ztime,Chi0(1:Nr)
        ztime=0d0
        do n=1,Np
            r(n)=dr*(n-0.5d0)
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