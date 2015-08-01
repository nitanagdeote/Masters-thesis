!-----------------------------------------------------
        subroutine REGULAR(ak,X,P,IP,AJ,LMAX,R)
! depends
!        convert_R2S(X,W)
!        associated_legendre( P,IP,n,W(2))
! slatec
!        dbesj(akr,0.5d0,LMAX,AJ,NBES)
!-----------------------------------------------------
	
        parameter (PI =3.14159265359d0, HALF_PI=1.57079632679d0)
        parameter (TWOPI=6.28318530718d0,FOURPI=12.5663706144d0)
        real*8 AK,X(3),P(0:LMAX-1)
        integer IP(0:LMAX-1)
        real*8  AJ(0:LMAX-1)
        complex*16 R(0:LMAX-1, -LMAX+1:LMAX-1) 

        real*8 AKR,W(3)

        call convert_R2S(X,W)
        AKR= W(1)*AK

        call dbesj(akr,0.5d0,LMAX,AJ,NBES)
        AJ = sqrt(HALF_PI/akr)*AJ   

        R(0,0) = AJ(0)/sqrt(FOURPI)
        do n=1,LMAX-1
         call associated_legendre( P,IP,n,W(2))
        do m= -n,n,1
         R(n,m) = AJ(n)
     &    * P(abs(m))
     &    * exp(cmplx(0.0d0,W(3)*m))
     &    * (-1.d0)**m/sqrt(TWOPI)
        enddo
        enddo

        return
        end

