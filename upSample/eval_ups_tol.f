      SUBROUTINE eval_ups_tol(nterms, pterms, hmpole, hqmpole, res) 
      implicit real *8 (a-h,o-z)
      
      integer :: nterms, pterms
      complex *16 :: ctmp
      real *8 :: res
      complex *16 :: hmpole(0:pterms, -pterms:pterms)
      complex *16 :: hqmpole(0:nterms, -nterms:nterms)
            
      sum1 = 0d0
      sum2 = 0d0
      do n = 0,pterms
        do m = -n,n
          ctmp = hmpole(n,m)-hqmpole(n,m)
         ! write(44,*) "n,m,hmpole(n,m)-hqmpole(n,m)",n,m,ctmp
          sum1 = sum1 + dble(ctmp*dconjg(ctmp))
          sum2 = sum2 + dble(hqmpole(n,m)*dconjg(hqmpole(n,m)))
        enddo
      enddo
      res = dsqrt(sum1/sum2)
      write(*,*) "res=", res
      END SUBROUTINE


      SUBROUTINE eval_ups_tol_real(nterms, pterms, hmpole, hqmpole, res) 
        implicit real *8 (a-h,o-z)
        
        integer :: nterms, pterms
        complex *16 :: ctmp
        real *8 :: res
        complex *16 :: hmpole(0:pterms, 0:pterms)
        complex *16 :: hqmpole(0:nterms, 0:nterms)
              
        sum1 = 0d0
        sum2 = 0d0
        do n = 0,pterms

          do m = 0,n
            ctmp = hmpole(n,m)-hqmpole(n,m)
         !   write(44,*) "n,m,hmpole(n,m)-hqmpole(n,m)",n,m,ctmp
            sum1 = sum1 + dble(ctmp*dconjg(ctmp))
            sum2 = sum2 + dble(hqmpole(n,m)*dconjg(hqmpole(n,m)))
          enddo
        enddo
        res = dsqrt(sum1/sum2)
        write(*,*) "res=", res
      !  write(44,*) "res=", res
        END SUBROUTINE