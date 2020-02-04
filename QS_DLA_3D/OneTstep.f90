      


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
        implicit none
        private
        save 
        public::Get_Chi,AdvancePulse,one_step_defile,one_step_output,getlaserfield

        contains  


            subroutine Get_Chi(Acom2d,outflag,time)
            implicit none 
            real(RP):: ztime,delta,Chi(1:Nr+1),time
            complex(RP) :: Acom2d(1:Nr,1:Nz+1),Ac(1:Nr)
            integer :: outflag
            integer :: nzc,nzc0,f_cross

            call resetplasma(ztime,Chi,nzc,nzc0,f_cross)
            Chi2d(1:Nr,1)=Chi(1:Nr)
            Ac(:)=Acom2d(:,1);
            
            if(outflag==0) then
            Err2dold=Err2d; Ezz2dold=Ezz2d;Bthh2dold=Bthh2d
            end if

            do while (nzc<Nz)
                  call one_step_z(ztime,Ac,Chi,nzc,f_cross,time)!Advance Trajectory
                  if (nzc>nzc0) then
                   nzc0=nzc
                   Ac(:)=Acom2d(:,nzc+1)
                   Chi2d(1:Nr,nzc+1)=Chi(1:Nr)
                       if(outflag==1)  then 
                       r2d(:,nzc+1)=r;V2d(:,nzc+1)=V;force2d(:,nzc+1)=force;gama2d(:,nzc+1)=gama
                       psi2d(:,nzc+1)=Psi(1:Nr);Er2d(:,nzc+1)=Er(1:Nr);Et2d(:,nzc+1)=Et(1:Nr);Bth2d(:,nzc+1)=Bth(1:Nr)
                       den2d(:,nzc+1)=den(1:Nr);den_e2d(:,nzc+1)=den_e(1:Nr);jr2d(:,nzc+1)=jr(1:Nr);jz2d(:,nzc+1)=jz(1:Nr)
                       else 
                       Err2d(1:Nr,nzc+1)=Er(1:Nr);Ezz2d(1:Nr,nzc+1)=Et(1:Nr)
                       Bthh2d(1:Nr,nzc+1)=Bth(1:Nr)
                       end if
                  else 
                   delta=(ztime-nzc*dz0)/dz0
                   Ac(:)=Acom2d(:,nzc+1)*delta+Acom2d(:,nzc)*(1-delta); 
                  endif 
            enddo
            end subroutine Get_Chi

            subroutine getlaserfield(Ac2d_0,Ac2d_hdt,Ac2d_dt,dt,dz,k0)
            use TestPart_var
            implicit none
            integer::i
            real(RP) :: dt,dz,k0
            complex(RP) :: Ac2d_0(1:Nr,1:Nz+1),Ac2d_hdt(1:Nr,1:Nz+1),Ac2d_dt(1:Nr,1:Nz+1)
            EL2dold=EL2d;BL2dold=BL2d;ELZ2dold=ELZ2d;
            
            Ac2d_hdt=(Ac2d_0+Ac2d_dt)*0.5d0


            BL2d(1:Nr,1:Nz)=ci*k0*(Ac2d_hdt(:,1:Nz)+Ac2d_hdt(:,2:Nz+1))*0.5-(Ac2d_hdt(:,2:Nz+1)-Ac2d_hdt(:,1:Nz))/dz
            EL2d(1:Nr,1:Nz)=BL2d(1:Nr,1:Nz)-0.5*(Ac2d_dt(:,2:Nz+1)+Ac2d_dt(:,1:Nz)-Ac2d_0(:,2:Nz+1)-Ac2d_0(:,1:Nz))/dt
        
            ELZ2d(2:Nr,1:Nz+1)=(Ac2d_hdt(2:Nr,:)-Ac2d_hdt(1:Nr-1,:))/hr
            

            end subroutine getlaserfield



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
           ! call check(nf90_create("output/"//TRIM(ADJUSTL(filename))//"Particle.nc", NF90_CLOBBER, ncidp) )
 
            call check( nf90_def_dim(ncidf, "ztime", Nz, ztime_dim) )
            call check( nf90_def_dim(ncidf, "Rgird", Nr, r_dim) )
            !call check( nf90_def_dim(ncidp, "ztime", Nz, ztime_dim) )
            !call check( nf90_def_dim(ncidp, "Part",  Np, p_dim) )
            dims=(/r_dim,ztime_dim/)
            call check( nf90_def_var(ncidf, "Den_e", NF90_REAL, dims, var_de) )
            !call check( nf90_def_var(ncidf, "Jjj_r", NF90_REAL, dims, var_jr) )
            !call check( nf90_def_var(ncidf, "Den_r", NF90_REAL, dims, var_dr) )
            call check( nf90_def_var(ncidf, "Bth_a", NF90_REAL, dims, var_bt) )
            call check( nf90_def_var(ncidf, "Psi_i", NF90_REAL, dims, var_ps) )
            !call check( nf90_def_var(ncidf, "Jjj_z", NF90_REAL, dims, var_jz) )
            call check( nf90_def_var(ncidf, "Eee_z", NF90_REAL, dims, var_ez) )
            call check( nf90_def_var(ncidf, "Eee_r", NF90_REAL, dims, var_er) )
            call check( nf90_def_var(ncidf, "Aaa_R", NF90_REAL, dims, var_AR) )
            call check( nf90_def_var(ncidf, "Aaa_I", NF90_REAL, dims, var_AI) )
            call check( nf90_def_var(ncidf, "Byl_R", NF90_REAL, dims, var_BR) )
            call check( nf90_def_var(ncidf, "Byl_I", NF90_REAL, dims, var_BI) )
            call check( nf90_def_var(ncidf, "Exl_R", NF90_REAL, dims, var_ELR) )
            call check( nf90_def_var(ncidf, "Exl_I", NF90_REAL, dims, var_ELI) )
            call check( nf90_def_var(ncidf, "Ezl_R", NF90_REAL, dims, var_EZLR) )
            call check( nf90_def_var(ncidf, "Ezl_I", NF90_REAL, dims, var_EZLI) )

            call check(nf90_put_att(ncidf,var_de,"Grid","Half") )
            !call check(nf90_put_att(ncidf,var_jr,"Grid","Int") )
            !call check(nf90_put_att(ncidf,var_dr,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_bt,"Grid","Int") )
            call check(nf90_put_att(ncidf,var_ps,"Grid","Half") )
            !call check(nf90_put_att(ncidf,var_jz,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_ez,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_er,"Grid","Int") )
            call check(nf90_put_att(ncidf,var_AR,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_AI,"Grid","Half") )
            call check(nf90_put_att(ncidf,var_BR,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_BI,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_ELR,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_ELI,"Grid","half-z,half-r") )
            call check(nf90_put_att(ncidf,var_EZLR,"Grid","full-z,full-r") )
            call check(nf90_put_att(ncidf,var_EZLI,"Grid","full-z,full-r") )
   

            !dims=(/p_dim,ztime_dim/)
            !call check( nf90_def_var(ncidp, "Vvv_r", NF90_REAL, dims, var_vr) )
            !call check( nf90_def_var(ncidp, "Rrr_r", NF90_REAL, dims, var_rr) )
            !call check( nf90_def_var(ncidp, "Fff_r", NF90_REAL, dims, var_fr) ) 
           ! call check( nf90_def_var(ncidp, "Gam_a", NF90_REAL, dims, var_ga) ) 

            call check( nf90_enddef(ncidf) )
            !call check( nf90_enddef(ncidp) )

            end subroutine one_step_defile


        subroutine one_step_output
        implicit none
        call check(nf90_put_var(ncidf,var_de,den_e2d(1:Nr,1:Nz)))
        !call check(nf90_put_var(ncidf,var_jr,jr2d(1:Nr,1:Nz)))
        !call check(nf90_put_var(ncidf,var_dr,den2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_bt,Bth2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_ps,Psi2d(1:Nr,1:Nz)))
        !call check(nf90_put_var(ncidf,var_jz,jz2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_ez,Et2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_er,Er2d(1:Nr,1:Nz)))
        call check(nf90_put_var(ncidf,var_AR,real(Ac2d(1:Nr,1:Nz))*2d0))
        call check(nf90_put_var(ncidf,var_AI,AImag(Ac2d(1:Nr,1:Nz))*2d0))
        call check(nf90_put_var(ncidf,var_BR,real(BL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_BI,AImag(BL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_ELR,real(EL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_ELI,AImag(EL2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_EZLR,real(ELZ2d(1:Nr,1:Nz))))
        call check(nf90_put_var(ncidf,var_EZLI,AImag(ELZ2d(1:Nr,1:Nz))))     
        !call check(nf90_put_var(ncidp,var_vr,V2d(1:Np,1:Nz)))
        !call check(nf90_put_var(ncidp,var_rr,r2d(1:Np,1:Nz)))
        !call check(nf90_put_var(ncidp,var_fr,force2d(1:Np,1:Nz)))
        !call check(nf90_put_var(ncidp,var_ga,gama2d(1:Np,1:Nz)))
        call check( nf90_close(ncidf))
        !call check( nf90_close(ncidp)) 
        end subroutine one_step_output
           



        end   module  OneTstep