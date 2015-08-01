       parameter LMAX=2
      parameter (PI =3.14159265359d0, HALF_PI=1.57079632679d0)
        parameter (TWOPI=6.28318530718d0,FOURPI=12.5663706144d0)
        real*8 AK,X(3),P(0:LMAX-1)
        integer IP(0:LMAX-1)
        real*8  AJ(0:LMAX-1)
        complex*16 S(0:LMAX-1, -LMAX+1:LMAX-1) 

        real*8 AKR,W(3)
         ak   = 2
         X(1) = 1
         x(2) = 1
         x(3) = 1
         P  = 1
         IP = 1
         AY = 1
      call singular(ak,X,P,IP,AJ,AY,LMAX,S)
      write(*,*)S

      end
