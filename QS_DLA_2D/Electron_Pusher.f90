!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------

    module electron_pusher
    use global_variables
    use domain
    implicit none
    private 
    save 
    public:: electron_push
    		
    	
!  F,S       jr(3),       Bth(3)          Er(3)  ************************   r=2*hr
! v  Ez(2)      jz(2),            den(2), psi(2) ------------------------   r=h1+hr/2
!  F,S       jr(2),       Bth(2)          Er(2)  ************************   r=hr
! v  Ez(1)      jz(1),            den(1), psi(1) ------------------------   r=hr/2
!  F,S       jr(1),       Bth(1),         Er(1)  ************************   r=0
      	

    contains


          

    SUBROUTINE electron_push(r0,r,V,gama,fv,f,Er,Et,Bth,ponderF,Psi,As)
    implicit none
    real(RP) ::delta,Er_p,Bth_p,Et_p,psi_p,As_p,ponder_p
    real(RP) :: r0(1:Np),r(1:Np),V(1:Np),fv(1:Np),f(1:Np),gama(1:Np)
    real(RP) :: Er(1:Nr),Et(1:Nr),Bth(1:Nr),ponderF(1:Nr),Psi(1:Nr),As(1:Nr)
    integer:: i,n

        fv(:)=v(:)
        f=0d0

     do n=1,Np

        if(r(n)>=Rmax-hr .or. r(n)<=-Rmax+hr) then
        f(n)=0d0; fv(n)=0d0;v(n)=0d0;
        cycle
        end if

        i=floor((r(n)+Rmax)/hr)
        delta=(r(n)+Rmax-hr*i)/hr
        Er_p=Er(i+1)*(1-delta)+Er(i+2)*delta
        Bth_p=Bth(i+1)*(1-delta)+Bth(i+2)*delta
        ponder_p=ponderF(i+1)*(1-delta)+ponderF(i+2)*delta


        i=floor((r(n)+hr/2d0+Rmax)/hr)+1
        delta= ((r(n)+hr/2d0+Rmax)-hr*(i-1))/hr
        psi_p=psi(i)*(1-delta)+psi(i+1)*delta
        Et_p=Et(i)*(1-delta)+Et(i+1)*delta
        As_p=As(i)*(1-delta)+As(i+1)*delta

                   
        !if(psi_p<-0.9)psi_p=-0.9
        gama(n)=(0.5d0/(1+psi_p))*(1+V(n)**2*(1+psi_p)**2+2*As_p+(1+psi_p)**2)
  !      fv(n)=v(n)
    !f(n)=-(gama(n)-0*(1+psi_p)*V(n)**2)*Er_p/(1+psi_p)**2+V(n)/(1+psi_p)*Ez_p           
        f(n)=((gama(n)*Er_p-ponder_p)/(1+psi_p)-(V(n)*Er_p+Et_p)*V(n))/(1+psi_p)-Bth_p/(1+psi_p)

         

    enddo

    end subroutine electron_push


   
 !****************fastest particle in domain*************tianhong  
         



    end module electron_pusher
