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
      real(RP):: time,TP_time,Centri_old
      integer :: ntc,outflag,outct,npp,npc,npctp,shift


    time=0d0  ! time
    ntc=0
    npctp=0
    outct=0




!------------------------------------------------------------!
!---------------------- Initialization-----------------------!
!------------------------------------------------------------!
call set_varibles(start,Centri_old)
!----------------------------------------------------------------!
!----------------------     T-loop    ---------------------------!
!----------------------------------------------------------------!
open(222,file='output/DomainShift',status="replace",action="write")




    if(Restart>0) then
    time=Restart*dt0*Dt_per_out
    ntc=Restart*Dt_per_out
    outct=Restart


    end if 


do while(time<=Tmax)
 

      if(mod(ntc*1d0,Dt_per_out/10d0)==0) write(*,'(A,F10.5)') '/////====time====////',time


      if(Laser) Acold2d=Ac2d
      !(1) =============== 
      Err2dold=Err2d; Ezz2dold=Ezz2d;Bthh2dold=Bthh2d
      EL2dold=EL2d; BL2dold=BL2d; ELZ2dold=ELZ2d


      call Get_Chi(Ac2d,time,npctp)  ! get chi at t
      if (TestParts) npctp=npctp+1

      !---------output--------------!
      if(ntc/Dt_per_out==outct) then
      call one_step_defile(outct) 
      call one_step_output ! output Plasma and laser same time

        if(outct>Restart .and.TestParts) then
        call defile_TP(outct)
        call TestPartOut
        end if 

      outct=outct+1;npctp=0
      call cpu_time(finish)
      write(*,'(A,I3,A,F10.6,A,F10.6,A)') " Output:=",outct-1,";===Time:=",Time,".(ckp)===,",Time/RayLen,'.Rayleigh Len'
      write(*,'(A,F9.5,A,F10.6,A)') " Progress:=====",(Time/Tmax*100)*1d0,'%.===Last==',(finish-start)/60,'mins=========='
      write(*,*) ''
      end if



      ntc=ntc+1

      time=time+dt0

   

    if(Laser.and.Domain_Shift) then

    call shiftdomain(Centri_old,Ac2d,shift)
    if(shift>0) then
    write(222,*) time,shift
    write(*,*) " =========Shift Laser Forward=================",shift
    flush(222)
    end if 

    end if


    end do
    close(222)
    !call one_step_defile(ntc)  ! final output
    write(*,'(A)') " ==================Done======================"

    STOP


  END PROGRAM NEW_WAKE
