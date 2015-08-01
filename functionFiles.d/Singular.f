!-----------------------------------------------------
        subroutine SINGULAR(AK,X,P,IP,AJ,AY,LMAX,S)
!depends:
!        associated_legendre( P,IP,n,W(2))
!        convert_R2S(X,W)
! slatec:
!        dbesj(akr,0.5d0,LMAX,AJ,NBES)
!        dbesy(akr,0.5d0,LMAX,AY,NBES)
!-----------------------------------------------------
        parameter (PI =3.14159265359d0, HALF_PI=1.57079632679d0)
        parameter (TWOPI=6.28318530718d0,FOURPI=12.5663706144d0)
        real*8 AK,X(3), P (0:LMAX-1)
        integer IP(0:LMAX-1)
        real*8  AJ(0:LMAX-1), AY(0:LMAX-1)
        complex*16 S(0:LMAX-1, -LMAX+1:LMAX-1) 

        real*8 AKR,W(3)

        call convert_R2S(X,W)
        AKR= W(1)*AK

        call dbesj(akr,0.5d0,LMAX,AJ,NBES)
        call dbesy(akr,0.5d0,LMAX,AY,NBES)

        AJ = sqrt(HALF_PI/akr)*AJ   
        AY = sqrt(HALF_PI/akr)*AY   
        S(0,0) = cmplx(AJ(0),AY(0))/sqrt(FOURPI)
        do n=1,LMAX-1
         call associated_legendre( P,IP,n,W(2))
        do m= -n,n,1
           S(n,m) = cmplx(AJ(n),AY(n))
     1    *P(abs(m))
     2    *exp(cmplx(0.0d0,W(3)*m))
     3    *(-1.0d0)**m/sqrt(TWOPI)
        enddo
        enddo

        return
        end
