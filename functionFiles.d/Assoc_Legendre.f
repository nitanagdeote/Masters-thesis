!------------------------------------------------------------        
        subroutine associated_legendre( P,IP, norder,theta)
! depends
! slatec
!       DXLEGF(1.d0*norder,0,0,norder,tmp,ID,P,IP,IERROR)
!------------------------------------------------------------        
        parameter (PI =3.14159265359d0,HALF_PI=1.57079632679d0)
        parameter (ID=4)
        real*8 theta,P(0:0)
        integer IP(0:0)

        real*8 darg,tmp
        
!==========================
        if (theta .eq.0.0d0) then
                tmp=1.0d-6 
        elseif (theta.eq.PI) then
                tmp = PI*(1.d0 - 1.0d-6)
        elseif (theta.eq.HALF_PI) then
                tmp =HALF_PI*(1.d0 - 1.0d-6)
        else
                tmp = theta
        endif
!===========================
        IF(tmp.lt.HALF_PI) then
         call DXLEGF(1.d0*norder,0,0,norder,tmp,ID,P,IP,IERROR)
        ELSE
         tmp= PI-tmp
         call DXLEGF(1.d0*norder,0,0,norder,tmp,ID,P,IP,IERROR)
                do m=0,norder
                        P(m)=P(m)*(-1)**(norder+M)
                enddo
        ENDIF

        return 
        end

