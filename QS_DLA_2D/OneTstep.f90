


        module  OneTstep
        use global_variables
        use domain
        use InAndOut
        use netcdf
        use pulse
        use OneZstep
        use fields_var
        use particles_var
        use TestPart_var
        use output_var
        use Test_Particle
        implicit none
        private
        save 
        public::Get_Chi,AdvancePulse,one_step_defile,one_step_output,getlaserfield,shiftdomain,LaserCentriod

        contains  


            subroutine Get_Chi(Acom2d,time,npctp)
            implicit none 
            real(RP):: ztime,delta,Chi(1:Nr),time,denb1d(1:Nr),jzb1d(1:Nr)
            complex(RP) :: Acom2d(1:Nr,1:Nz+1),Ac(1:Nr)
            integer :: outflag,npctp
            integer :: nzc,nzc0,f_cross,npp
            real(RP):: dd0,dt
           complex(RP) :: dd1(1:Nr),dd2(1:Nr),dd3(1:Nr)
           complex(RP) :: Aold(1:Nr),Anew(1:Nr)

            call resetplasma(ztime,Chi,nzc,nzc0,f_cross)
            Chi2d(1:Nr,1)=Chi(1:Nr);Ac(:)=Acom2d(:,1);
            denb1d(:)=denb2d(:,1);jzb1d(:)=jzb2d(:,1);
            dd0=1.5d0;dd1=0d0;dd2=0d0
            
            if (TestParts)   ifpush=1

            do while (nzc<Nz)

                  call one_step_z(ztime,Ac,denb1d,jzb1d,Chi,nzc,f_cross,time)!Advance Trajectory
                
                  !write(*,*) ztime
                  if (nzc>nzc0) then
                   nzc0=nzc  

                   Chi2d(1:Nr,nzc+1)=Chi(1:Nr)
                   denb1d(:)=denb2d(:,nzc+1);Jzb1d(:)=jzb2d(:,nzc+1)

    !-----------------------------------
                       r2d(:,nzc+1)=r;V2d(:,nzc+1)=V;force2d(:,nzc+1)=force;gama2d(:,nzc+1)=gama

                       psi2d(:,nzc+1)=Psi(1:Nr);Er2d(:,nzc+1)=Er(1:Nr);Et2d(:,nzc+1)=Et(1:Nr);Bth2d(:,nzc+1)=Bth(1:Nr)
                       den2d(:,nzc+1)=den(1:Nr);den_e2d(:,nzc+1)=den_e(1:Nr);jr2d(:,nzc+1)=jr(1:Nr);jz2d(:,nzc+1)=jz(1:Nr)
                       
                       Err2d(1:Nr,nzc+1)=Er(1:Nr);Ezz2d(1:Nr,nzc+1)=Et(1:Nr);Bthh2d(1:Nr,nzc+1)=Bth(1:Nr)
    

                      if(laser) then
                        Aold(:)=Acom2d(:,nzc+1)
                        call One_Step_Pulse(dd0,dd1,dd2,Chi2d(:,nzc+1),Aold,Anew,dz0,dt0)
                        dd2(:)=dd1(:)
                        dd1(:)=Anew(:)-Aold(:)
                        Acom2d(:,nzc+1)=Anew(:)
                        call getlaserfield(Acold2d,Acom2d,dt0,dz0,k0,nzc)
                      end if 
                        Ac(:)=Acom2d(:,nzc+1);


                      if (TestParts) then 
                      call PushTestParticlesNest(k0,nzc,npctp)
                      if(Interaction) call TestParticlesdensity(nzc)
                      end if
                      
    !-----------------------------------
                  else 
                  delta=(ztime-nzc*dz0)/dz0
                  Ac(:)=Acom2d(:,nzc+1)*delta+Acom2d(:,nzc)*(1-delta); 
                  denb1d(:)=denb2d(:,nzc+1)*delta+denb2d(:,nzc)*(1-delta); 
                  jzb1d(:)= jzb2d(:,nzc+1)*delta +jzb2d(:,nzc)*(1-delta); 
                  endif 

            enddo


            end subroutine Get_Chi

            subroutine getlaserfield(Ac2d_0,Ac2d_dt,dt,dz,k0,nzc)
            use TestPart_var
            implicit none
            integer::i,nzc
            real(RP) :: dt,dz,k0
            complex(RP) :: Ac2d_0(1:Nr,1:Nz+1),Ac2d_hdt(1:Nr,2),Ac2d_dt(1:Nr,1:Nz+1)  

            Ac2d_hdt=0.5d0*(Ac2d_0(:,nzc:nzc+1)+Ac2d_dt(:,nzc:nzc+1))
            BL2d(:,nzc)=ci*k0*(Ac2d_hdt(:,1)+Ac2d_hdt(:,2))*0.5-(Ac2d_hdt(:,2)-Ac2d_hdt(:,1))/dz
            EL2d(:,nzc)=BL2d(1:Nr,nzc)-0.5*(Ac2d_dt(:,nzc+1)+Ac2d_dt(:,nzc)-Ac2d_0(:,nzc+1)-Ac2d_0(:,nzc))/dt
            ELZ2d(1:Nr-1,nzc+1)=(Ac2d_hdt(2:Nr,2)-Ac2d_hdt(1:Nr-1,2))/hr
  
            end subroutine getlaserfield


           subroutine shiftdomain(Centri_old,Ac2d,shift)
           implicit none
           real(RP):: Centri_New,Centri_old
           complex(RP)::Ac2d(1:Nr,1:Nz+1)
           integer:: shift,n
           complex(RP):: BL(0:Nz+1)
           shift=0
           BL(0:Nz+1)=BL2d(Nr/2,:)
           Centri_New=LaserCentriod(BL)
           write(*,'(A,F10.5,A,F10.5)') '==Centriod==',Centri_New,';  =old=',Centri_old
           if (Centri_New-Centri_old>=1) then
           shift=floor(Centri_New-Centri_old)

           do n=1,Nz+1-shift
           Ac2d(:,n)=Ac2d(:,n+shift)
           Err2dold(:,n)=Err2dold(:,n+shift);Ezz2dold(:,n)=Ezz2dold(:,n+shift)
           Bthh2dold(:,n)=Bthh2dold(:,n+shift);EL2dold(:,n)=EL2dold(:,n+shift)
           BL2dold(:,n)=BL2dold(:,n+shift);ELZ2dold(:,n)=ELZ2dold(:,n+shift) 
           end do
           Ac2d(:,Nz+2-shift:Nz+1)=0d0
           Err2dold(:,Nz+2-shift:Nz+1)=0d0;Ezz2dold(:,Nz+2-shift:Nz+1)=0d0
           Bthh2dold(:,Nz+2-shift:Nz+1)=0d0;EL2dold(:,Nz+2-shift:Nz+1)=0d0
           BL2dold(:,Nz+2-shift:Nz+1)=0d0;ELZ2dold(:,Nz+2-shift:Nz+1)=0d0
  
           if(TestParts) ZZT=ZZT-shift*dz0
           
           Centri_old=Centri_New-shift
           end if


           
           end subroutine shiftdomain

                    
           real(RP) function LaserCentriod(BL)
           implicit none
           complex(RP):: BL(0:Nz+1)
           real(RP)   :: Int(0:Nz+1),ENR
           integer :: n
           EnR=0d0
            Int=Abs(BL)**2
            do n=0,Nz+1
            EnR=EnR+n*Int(n)
            end do
            LaserCentriod=EnR/sum(Int)
            return
           end function LaserCentriod

       		subroutine AdvancePulse(Acom2d,dt)
            use fields_var
        	implicit none 
        	real(RP):: dd0,dt
        	complex(RP) :: dd1(1:Nr),dd2(1:Nr),dd3(1:Nr)
    		complex(RP) :: Aold(1:Nr),Anew(1:Nr),Acom2d(1:Nr,1:Nz+1)
    		integer :: k
            dd0=1.5d0;dd1=0d0;dd2=0d0
            !Chi2d=1d0
                do k=1,Nz
                 Aold(:)=Acom2d(:,k)
                 call One_Step_Pulse(dd0,dd1,dd2,Chi2d(:,k),Aold,Anew,dz0,dt)
                 dd2(:)=dd1(:)
                 dd1(:)=Anew(:)-Aold(:)
                 Acom2d(:,k)=Anew(:)
                enddo

        	end subroutine AdvancePulse


            subroutine one_step_defile(ntc)
            implicit none
            integer :: ntc
            character(12):: filename

            write(filename,'(I5)'), ntc
            call check(nf90_create("output/"//TRIM(ADJUSTL(filename))//"Fields.nc", NF90_CLOBBER, ncidf) )
          !  call check(nf90_create("output/"//TRIM(ADJUSTL(filename))//"Particle.nc", NF90_CLOBBER, ncidp) )
 
            call check( nf90_def_dim(ncidf, "ztime", Nz, ztime_dim) )
            call check( nf90_def_dim(ncidf, "Rgird", Nr, r_dim) )
           ! call check( nf90_def_dim(ncidp, "ztime", Nz, ztime_dim) )
           ! call check( nf90_def_dim(ncidp, "Part",  Np, p_dim) )
            dims=(/r_dim,ztime_dim/)
            call check( nf90_def_var(ncidf, "Den_e", NF90_REAL, dims, var_de) )
            call check( nf90_def_var(ncidf, "Jjj_r", NF90_REAL, dims, var_jr) )
            call check( nf90_def_var(ncidf, "Den_r", NF90_REAL, dims, var_dr) )
            call check( nf90_def_var(ncidf, "Bth_a", NF90_REAL, dims, var_bt) )
            call check( nf90_def_var(ncidf, "Psi_i", NF90_REAL, dims, var_ps) )
            call check( nf90_def_var(ncidf, "Jjj_z", NF90_REAL, dims, var_jz) )
            call check( nf90_def_var(ncidf, "Eee_z", NF90_REAL, dims, var_ez) )
            call check( nf90_def_var(ncidf, "Eee_r", NF90_REAL, dims, var_er) )

            if(Laser) then
            call check( nf90_def_var(ncidf, "Aaa_R", NF90_REAL, dims, var_AR) )
            call check( nf90_def_var(ncidf, "Aaa_I", NF90_REAL, dims, var_AI) )
            call check( nf90_def_var(ncidf, "Byl_R", NF90_REAL, dims, var_BR) )
            call check( nf90_def_var(ncidf, "Byl_I", NF90_REAL, dims, var_BI) )
            call check( nf90_def_var(ncidf, "Exl_R", NF90_REAL, dims, var_ELR) )
            call check( nf90_def_var(ncidf, "Exl_I", NF90_REAL, dims, var_ELI) )
            call check( nf90_def_var(ncidf, "Ezl_R", NF90_REAL, dims, var_EZLR) )
            call check( nf90_def_var(ncidf, "Ezl_I", NF90_REAL, dims, var_EZLI) )
            end if

           if(Interaction) then
           call check( nf90_def_var(ncidf, "Den_b", NF90_REAL, dims, var_deb) )
           call check( nf90_def_var(ncidf, "Jjz_b", NF90_REAL, dims, var_jzb) )
           end if 

            call check(nf90_put_att(ncidf,var_de,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_jr,"Grid","Int") )
            call check(nf90_put_att(ncidf,var_dr,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_bt,"Grid","Int") )
            call check(nf90_put_att(ncidf,var_ps,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_jz,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_ez,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_er,"Grid","Int") )
            if(Laser) then
            call check(nf90_put_att(ncidf,var_AR,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_AI,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_BR,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_BI,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_ELR,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_ELI,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_EZLR,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_EZLI,"Grid","half-z,half-r") )
            end if 

           ! dims=(/p_dim,ztime_dim/)
           ! call check( nf90_def_var(ncidp, "Vvv_r", NF90_REAL, dims, var_vr) )
           ! call check( nf90_def_var(ncidp, "Rrr_r", NF90_REAL, dims, var_rr) )
            !call check( nf90_def_var(ncidp, "Fff_r", NF90_REAL, dims, var_fr) ) 
           ! call check( nf90_def_var(ncidp, "Gam_a", NF90_REAL, dims, var_ga) ) 

            call check( nf90_enddef(ncidf) )
           ! call check( nf90_enddef(ncidp) )

            end subroutine one_step_defile


        subroutine one_step_output
        implicit none
        integer :: n
        complex(RP) :: temp(1:Nr,0:Nz+1)

        do n=1,Nz
        ELZ2d(:,n)=(ELZ2d(:,n+1)+ELZ2d(:,n))*0.5d0
        end do
        do n=1,Nr-1
        temp(n+1,:)=(ELZ2d(n,:)+ELZ2d(n+1,:))*0.5d0
        end do
        ELZ2d=temp

        call check(nf90_put_var(ncidf,var_de,den_e2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_jr,jr2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_dr,den2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_bt,Bth2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_ps,Psi2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_jz,jz2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_ez,Et2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_er,Er2d(1:Nr,1:Nz)))
        if(Laser) then
        call check(nf90_put_var(ncidf,var_AR,real(Ac2d(1:Nr,1:Nz))*2d0))
        call check(nf90_put_var(ncidf,var_AI,AImag(Ac2d(1:Nr,1:Nz))*2d0))
        call check(nf90_put_var(ncidf,var_BR,real(BL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_BI,AImag(BL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_ELR,real(EL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_ELI,AImag(EL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_EZLR,real(ELZ2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_EZLI,AImag(ELZ2d(1:Nr,1:Nz))))    
        end if 

         if(Interaction) then
        call check(nf90_put_var(ncidf,var_deb,denb2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_jzb, jzb2d(1:Nr,1:Nz)))
         end if 

        !call check(nf90_put_var(ncidp,var_vr,V2d(1:Np,1:Nz)))
        !call check(nf90_put_var(ncidp,var_rr,r2d(1:Np,1:Nz)))
        !call check(nf90_put_var(ncidp,var_fr,force2d(1:Np,1:Nz)))
       ! call check(nf90_put_var(ncidp,var_ga,gama2d(1:Np,1:Nz)))
        call check( nf90_close(ncidf))
        !call check( nf90_close(ncidp)) 
        end subroutine one_step_output
           



        end   module  OneTstep
