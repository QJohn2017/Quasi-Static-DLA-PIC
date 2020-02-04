    
module magnetic_f
    use global_variables
    use domain
    use Functions
    implicit none
    private 
    save 
    public:: electron_current_jz,electron_current_jv,force_rhoF,magnetic_field


contains
! this program uses V instead of pr: V=pr/(1+psi)
! tianhong Feb-02-2017 No problem No need to check or change!!!
    subroutine electron_current_jz(jv,den,As,Psi,jz,ponderF)   !???????
! density current jr is defined at integer points r=hr*(i)
! jz=n*vz=n_r/(1-vz)*vz=n_r*pz/(1+psi)=n_r*0.5/(1+psi)**2*( 1+V**2*(1+psi)**2-(1+psi)**2)
! jz=0.5*jv+0.5*n_r*(1-(1+psi)**2)/(1+psi)**2
! ponder F is on integer grid
implicit none
real(RP) :: jv(1:Nr+1),den(1:Nr+1),As(1:Nr+1),Psi(1:Nr+1),jz(1:Nr+1),ponderF(1:Nr+1)
integer:: k

jz=0d0
ponderF=0d0

do k=1,Nr
    jz(k)=0.5d0*jv(k)+0.5d0*den(k)*(1+2*As(k)-(1+psi(k))**2)/(1+psi(k))**2
    if (k==1) then
      ponderF(1)=0
    else
      ponderF(k)=(As(k)-As(k-1))/hr
  end if 
  enddo
end subroutine electron_current_jz


!  F         jr(3),       Bth(3)          Er(3)  ************************   r=2*hr
! v  Ez(2)      jz(2),            den(2), psi(2) ------------------------   r=h1+hr/2
!  F         jr(2),       Bth(2)          Er(2)  ************************   r=hr
! v  Ez(1)      jz(1),            den(1), psi(1) ------------------------   r=hr/2
!  F         jr(1),       Bth(1),         Er(1)  ************************   r=0

! tianhong Feb-02-2017 No problem No need to check or change!!!
! change something see old
subroutine electron_current_jv(r,V,mass,dentemp,jv)
! ensity current jv is defined at half integer points r=hr*(i)
implicit none
real(RP) :: r(1:Np),V(1:Np),mass(1:Np),jv(1:Nr+1),masstemp(1:Np)
real(RP) :: delta,dentemp
integer:: i,n,k

jv=0d0
masstemp=mass*dentemp
do n=1,Np
    if(r(n)>Rmax+hr/2d0) cycle
    if(r(n)<hr/2d0)then
        jv(1)=jv(1)+masstemp(n)*V(n)**2*r(n)/(hr/2d0)
            !old jv(1)=jv(1)+mass(n)*V(n)**2/(hr/2)
        endif
        if(r(n)>=hr/2d0)then
            i=floor((r(n)-hr/2)/hr)+1
            delta=((r(n)-hr/2)-hr*(i-1))/hr
            jv(i)=jv(i)+ masstemp(n)*(1-delta)*V(n)**2
            jv(i+1)=jv(i+1)+masstemp(n)*delta*V(n)**2
        endif
        enddo
        do k=1,Nr
            jv(k)=jv(k)/(2d0*pi*(hr*k-hr/2d0)*hr)
            enddo

        end subroutine electron_current_jv
        
!  F,S       jr(3),       Bth(3)          Er(3)  ************************   r=2*hr
! v  Ez(2)      jz(2),            den(2), psi(2) ------------------------   r=h1+hr/2
!  F,S       jr(2),       Bth(2)          Er(2)  ************************   r=hr
! v  Ez(1)      jz(1),            den(1), psi(1) ------------------------   r=hr/2
!  F,S       jr(1),       Bth(1),         Er(1)  ************************   r=0

! tianhong Feb-03-2017 No problem No need to check or change!!!
! I changed something, see old 
subroutine force_rhoF(Psi,den,Et,Er,jv,jr,As,ponderF,rhoF)
! ere F is \rho*F, it  is defined at integer points r=hr*i
implicit none
real(RP) :: Psi(1:Nr+1),den(1:Nr+1),Et(1:Nr+1),As(1:Nr+1),Er(1:Nr+1),jr(1:Nr+1),jv(1:Nr+1),rhoF(1:Nr+1),ponderF(1:Nr+1)
real(RP) :: psi_a,jv_a,Et_a,den_a,As_a
integer:: k
rhoF=0d0
do k=1,Nr           
    if(k==1)then
        psi_a=1d0+0.5d0*(psi(k)+psi(k))
        jv_a=0d0
            ! old       jv_a=0.5d0*(jv(k)+jv(k))
            den_a=0.5d0*(den(k)+den(k))
            Et_a=0.5d0*(Et(k)+Et(k))
            As_a=As(1)

            !rhoF(k)=1/psi_a**2*(0.5d0/psi_a*(den_a+jv_a*psi_a**2+den_a*psi_a**2)-jv_a*psi_a**2)*Er(k)
            rhoF(k)=1d0/psi_a**2*(0.5d0/psi_a*(den_a*(1+2*As_a+psi_a**2)+jv_a*psi_a**2)-jv_a*psi_a)*Er(k)
            rhoF(k)=rhoF(k)-jr(k)/psi_a*Et_a-(den_a/psi_a**2)*ponderF(1)
        endif

        if(k>1)then
            psi_a=1d0+0.5d0*(psi(k)+psi(k-1))
            jv_a=0.5d0*(jv(k)+jv(k-1))
            den_a=0.5d0*(den(k)+den(k-1))
            Et_a=0.5d0*(Et(k)+Et(k-1))
            As_a=0.5d0*(As(k)+As(k-1))
            !rhoF(k)=1/psi_a**2*(0.5d0/psi_a*(den_a+jv_a*psi_a**2+den_a*psi_a**2)-jv_a*psi_a**2)*Er(k)
            rhoF(k)=1d0/psi_a**2*(0.5d0/psi_a*(den_a*(1+2*As_a+psi_a**2)+jv_a*psi_a**2)-jv_a*psi_a)*Er(k)
            rhoF(k)=rhoF(k)-jr(k)/psi_a*Et_a-(den_a/psi_a**2)*ponderF(k)
        endif
        enddo    
    end subroutine force_rhoF

!  F,S       jr(3),       Bth(3)          Er(3)  ************************   r=2*hr
! v  Ez(2)      jz(2),            den(2), psi(2) ------------------------   r=h1+hr/2
!  F,S       jr(2),       Bth(2)          Er(2)  ************************   r=hr
! v  Ez(1)      jz(1),            den(1), psi(1) ------------------------   r=hr/2
!  F,S       jr(1),       Bth(1),         Er(1)  ************************   r=0


! tianhong Feb-03-2017 No problem No need to check or change!!!      
subroutine source_for_B(jv,rhof,jz,S)
! ere F is \rho*F, it  is defined at integer points r=hr*i
implicit none
real(RP) ::rint(1:Nr),rhalf(1:Nr),jv(1:Nr+1),rhof(1:Nr+1),jz(1:Nr+1),S(1:Nr)
integer:: k

s=0d0
do k=1,Nr
    rint(k)=hr*k-hr
    rhalf(k)=hr*k-hr/2d0
    enddo


    do k=1,Nr           
        if(k==1)then
            S(k)=0d0
        endif

        if(k>1)then
            S(k)=(rhalf(k)*jv(k)-rhalf(k-1)*jv(k-1))/(rint(k)*hr)-rhoF(k)-(jz(k)-jz(k-1))/hr
        endif
        enddo
        
    end subroutine source_for_B


    
!  F         jr(3),       Bth(3)          Er(3)  ************************   r=2*hr
! v  Ez(2)      jz(2),            den(2), psi(2) ------------------------   r=h1+hr/2
!  F         jr(2),       Bth(2)          Er(2)  ************************   r=hr
! v  Ez(1)      jz(1),            den(1), psi(1) ------------------------   r=hr/2
!  F         jr(1),       Bth(1),         Er(1)  ************************   r=0

    ! tianhong Feb-03-2017 No problem No need to check or change!!!
    subroutine magnetic_field(jv,rhof,jz,psi,den,Bth)	
        implicit none
        real(RP) ::rint(1:Nr),rhalf(1:Nr),jv(1:Nr+1),rhof(1:Nr+1),jz(1:Nr+1),S(1:Nr),psi(1:Nr+1),den(1:Nr+1),Bth(1:Nr+1)
        real(RP) :: a(1:Nr-1),b(1:Nr),c(1:Nr-1),d(1:Nr)
        integer:: k,info
        real(RP) :: psi_a,den_a
        bth=0d0
        do k=1,Nr
            rint(k)=hr*k-hr
            rhalf(k)=hr*k-hr/2d0
            enddo

            call source_for_B(jv,rhof,jz,S)

            b(1)=1d0; c(1)=0d0; d(1)=0d0     
            do k=2,Nr-1
                psi_a=1d0+0.5d0*(psi(k)+psi(k-1))
                den_a=0.5d0*(den(k)+den(k-1))
                a(k-1)=rint(k-1)/rhalf(k-1)/hr**2;
                b(k)=-rint(k)*(1/rhalf(k)+1/rhalf(k-1))/hr**2-den_a/psi_a
                c(k)=rint(k+1)/rhalf(k)/hr**2
                d(k)=S(k)
                enddo              
                k=Nr
                psi_a=1d0+0.5d0*(psi(k)+psi(k-1))
                den_a=0.5d0*(den(k)+den(k-1))
                a(k-1)=rint(k-1)/(rhalf(k-1))/hr**2;
                b(k)=-rint(k)*(1/rhalf(k)+1/rhalf(k-1))/hr**2-den_a/psi_a
                d(k)=S(k)
                call DGTSV( Nr, 1, a, b, c, d,Nr, INFO )
                Bth(1:Nr)=d
                Bth(Nr+1)=0d0   
                
            end subroutine magnetic_field

        end module magnetic_f
