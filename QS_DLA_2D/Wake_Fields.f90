    module wake_f
    use global_variables
    use domain
    implicit none
    private 
    save 
    public::  electron_density_r,electron_current_jr,wake_fields,wake_field_Et,electron_density_e,Susceptibility


    contains

! this program uses V instead of pr: V=pr/(1+psi)
! tianhong Feb-02-2017 No problem No need to check or change!!!
    subroutine electron_density_r(r,mass,dentemp,den)     
! ensity defined at half-integer points r=hr*(i-1/2)
    implicit none
    real(RP) :: r(1:Np),mass(1:Np),den(1:Nr),masstemp(1:Np)
    real(RP) :: delta,dentemp
    integer:: i,n,k
    den=0d0
    masstemp=mass*dentemp
    do n=1,Np
        !---bebug
    if(r(n)>=Rmax+hr/2 .or. r(n)<=-Rmax-hr/2) cycle
            i=floor((r(n)+Rmax+hr/2d0)/hr)+1
    if(i<1 .or. i>Nr-1) cycle
            delta= ((r(n)+Rmax+hr/2d0)-hr*(i-1))/hr
            den(i)=den(i)+masstemp(n)*(1-delta)
            den(i+1)=den(i+1)+masstemp(n)*delta
    enddo
     
    den=den/hr
    
    den(1)=dentemp;den(Nr)=dentemp;

    end subroutine electron_density_r
    	
!  F         jr(3),       Bth(3)          Er(3)  ************************   r=2*hr
! v  Ez(2)      jz(2),            den(2), psi(2) ------------------------   r=h1+hr/2
!  F         jr(2),       Bth(2)          Er(2)  ************************   r=hr
! v  Ez(1)      jz(1),            den(1), psi(1) ------------------------   r=hr/2
!  F         jr(1),       Bth(1),         Er(1)  ************************   r=0
           
    subroutine electron_density_e(jv,psi,den,As,den_e)
! ensity defined at half-integer points r=hr*(i-1/2)
    implicit none
    real(RP) :: jv(1:Nr),den(1:Nr),Psi(1:Nr),Den_e(1:Nr),As(1:Nr)
    real(RP) :: delta
    integer:: i,n,k
    den_e=0d0
    do k=1,Nr
        den_e(k)=0.5d0*jv(k)+0.5d0*den(k)*(1+2*As(k)+(1+psi(k))**2)/(1+psi(k))**2
    enddo
! rint*,den(35)
! top
    end subroutine electron_density_e

    subroutine wake_fields(den,denb1d,jzb1d,Er,Psi,dentemp)
    implicit none
    real(RP) :: rho(1:Nr),den(1:Nr),Er(1:Nr),Psi(1:Nr),Er2(1:Nr),Psi2(1:Nr),denb1d(1:Nr),jzb1d(1:Nr)
    real(RP) :: SS
    real(RP) :: delta,dentemp
    integer :: i,n,k

    !====tianhong====!       
            
    Er=0d0
    Er2=0d0
    Psi2=0d0
    Psi=0d0
    rho=den-dentemp!+denb1d-jzb1d

    do k=1,Nr
        !====tianhong====!  
        if (k==1) then 
        Er(1)=rho(1)*hr
        else 
        Er(k)=Er(k-1)+rho(k)*hr
        end if 
    enddo

    do k=Nr-1,1,-1   ! will end at 1 not 0
    psi(k)=psi(k+1)-Er(k)*hr
    enddo

    do k=Nr-1,1,-1
    Er2(k)=Er2(k+1)-rho(k+1)*hr
    enddo
 
    do k=1,Nr-1  ! will end at 1 not 0
    psi2(k+1)=psi2(k)+Er2(k)*hr
    enddo

    Er=(Er+Er2)*0.5d0
    Psi=(Psi+Psi2)*0.5d0
                    
! rint*,den(35)
! stop
    end subroutine wake_fields
! tianhong Feb-02-2017 No problem No need to check or change!!!
! change something see old
! And I need to check how to decomposite the velocity-related quantity. Is that self-consistent ?
    subroutine electron_current_jr(r,V,mass,dentemp,jr)
! ensity current jr is defined at integer points r=hr*(i)
    implicit none
    real(RP) :: r(1:Np),V(1:Np),mass(1:Np),jr(1:Nr),masstemp(1:Np)
    real(RP) :: delta,dentemp
    integer:: i,n,k

    jr=0d0     ! jr(1) is on axis
    masstemp=mass*dentemp
    do n=1,Np

        if(r(n)>Rmax .or. r(n)<-Rmax) cycle
        
            i=floor((r(n)+Rmax)/hr)
            delta=(r(n)+Rmax-hr*i)/hr
            jr(i+1)=jr(i+1)+masstemp(n)*V(n)*(1-delta)
            jr(i+2)=jr(i+2)+masstemp(n)*V(n)*delta
    enddo

    jr=jr/hr



    end subroutine electron_current_jr
! tianhong Feb-02-2017 No problem No need to check or change!!!

    subroutine wake_field_Et(jr,Et)
! t =E_{\xi} is defined at half-integer points r=hr*i-hr/2;

! t/dr=-jr,  jr positive
    use global_variables
    implicit none
    real(RP) :: jr(1:Nr+1),Et(1:Nr),Et2(1:Nr)
    real(RP) :: SS
    real(RP) :: delta
    integer:: i,n,k

    Et=0d0
    Et2=0d0

    do k=Nr-1,1,-1
        Et(k)=Et(k+1)+jr(k)*hr
    enddo

    do k=1,Nr-1
        Et2(k+1)=Et2(k)-jr(k)*hr
    enddo

    Et=0.5d0*(Et+Et2)

    end subroutine wake_field_Et


    subroutine  Susceptibility(den,Psi,Chi)
    implicit none
    real(RP) :: den(1:Nr),Psi(1:Nr),Chi(1:Nr)
    integer :: i,n,k
    
    Chi(:)=den(:)/(1+Psi(:))
        
    end subroutine Susceptibility



    end module wake_f

!     ************************************************