    PROGRAM NEW_WAKE
      USE OneTstep
      USE initialization
      USE global_variables
      USE domain
      use pulse_var
      use TestPart_var
      use Test_Particle

      
      implicit none

      real(RP):: start,finish
      real(RP):: time,TP_time
      integer :: ntc,outflag,outct,npctp,npp,npc


    time=0d0  ! time
    TP_time=0d0
    ntc=0
    npc=0
    npctp=0
   
    outct=0

!------------------------------------------------------------!
!---------------------- Initialization-----------------------!
!------------------------------------------------------------!
call set_varibles(start)
!----------------------------------------------------------------!
!----------------------     T-loop    ---------------------------!
!----------------------------------------------------------------!


do while(time<=Tmax)

write(*,'(A,F10.5)') '====time====',time
      
       Acold2d=Ac2d;Atemp2d=Ac2d
  !(1) =============== 
      call Get_Chi(Ac2d,1,time)  ! get chi at t
      
      if(ntc/Dt_per_out==outct) then
      call one_step_defile(outct) 
      call one_step_output ! output Plasma and laser same time


      if(outct>0 .and.TestParts) then
      call defile_TP(outct)
      call TestPartOut
      end if 

      outct=outct+1;npctp=0
      call cpu_time(finish)
      write(*,'(A,I3,A,F8.4,A,F8.4,A)') " Output:=",outct-1,";===Time:=",Time,".(ckp)===,",Time/RayLen,'.Rayleigh Len'
      write(*,'(A,F7.3,A,F8.4,A)') " Progress:=====",(Time/Tmax*100)*1d0,'%.===Last==',(finish-start)/60,'mins=========='
      write(*,*) ''
      end if
  
  !(2) ===============   
      call AdvancePulse(Ac2d,dt0*0.5d0) ! A(t)-> A(t+0.5)
  !(3) ===============    
      call Get_Chi(Ac2d,0,time)  ! get chi at t+0.5

  !(4) ===============     
      call AdvancePulse(Acold2d,dt0) ! A(t)-> A(t+1)  
      
      call getlaserfield(Atemp2d,Ac2d,Acold2d,dt0,dz0,k0)
      if (ntc==0) then 
      Err2dold=Err2d; Ezz2dold=Ezz2d;Bthh2dold=Bthh2d
      EL2dold=EL2d;BL2dold=BL2d;ELZ2dold=ELZ2d
      end if

      Ac2d=Acold2d
      ntc=ntc+1
      time=time+dt0

     
     !  ======Test_Particle=========  
     
    if (TestParts) then 
      do npp=1,N_per_dt 
      call PushTestParticles(k0,TP_time,npp,npctp)
      end do
    npctp=npctp+1
    end if 


    end do
    
    !call one_step_defile(ntc)  ! final output
    write(*,'(A)') " ==================Done======================"
    STOP


  END PROGRAM NEW_WAKE
