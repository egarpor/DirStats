      Subroutine lscv_dir_vmf(data_dir, h, Cq, Cq2, n, nq, lenh, kk, B, CV)
      implicit none

      ! Arguments
      Integer n, nq, kk, lenh
      Double Precision data_dir(n, nq + 1), B(kk)
      Double Precision h(lenh), Cq(lenh), Cq2(lenh)
      Double Precision CV(lenh)

      ! Local variables
      Double Precision CV1(n, n, lenh), CV2(n, n, lenh)
      Double Precision aux1, aux2, aux4, aux5, ALPHA
      Integer hi, i, j, NCALC, IZE
      Real, parameter :: Pi = 3.1415926536

      ! For calling RIBESL
      ALPHA = 0.5 * (nq - 1) - (kk - 1)
      IZE = 1

      ! Set CV1 and CV2 to zero
      CV1 = 0
      CV2 = 0

      ! Computation of CV1 and CV2
      do i = 2, n
        do j = 1, i - 1

           aux1 = exp(dot_product(data_dir(i, :), data_dir(j, :)))
           aux2 = Sqrt(Sum((data_dir(i, :) + data_dir(j, :)) ** 2))

           do hi = 1, lenh

                ! CV1
                CV1(i, j, hi) = aux1 ** (h(hi) ** (-2.0))

                ! CV2
                aux4 = aux2 / (h(hi) ** 2)
                call RIBESL(aux4, ALPHA, kk, IZE, B, NCALC)
                aux5 = ((2 * Pi) ** (0.5 * (nq + 1))) * B(kk)
                CV2(i, j, hi) = aux5 / (aux4 ** (0.5 * (nq - 1)))

           end do

        end do
      end do

      CV = (4 * Cq) / (n * (n - 1)) * Sum(Sum(CV1, dim = 1), dim = 1)
      CV = CV - (Cq ** 2) / (n * Cq2)
      CV = CV - (2 * Cq ** 2) / (n ** 2) * Sum(Sum(CV2, dim = 1), dim = 1)
      CV = -CV

      Return
      End
