      Subroutine lscv_dir_vmf(data_dir, h, Cq, Cq2, n, nq, lenh, kk, B, CV)
      implicit none

      ! Arguments
      Integer n, nq, kk, lenh
      Double Precision data_dir(n, nq + 1), B(kk)
      Double Precision h(lenh), Cq(lenh), Cq2(lenh)
      Double Precision CV(lenh)

      ! Local variables
      Double Precision sumCV1(lenh), sumCV2(lenh)
      Double Precision dotij, aux2, aux4, aux5, ALPHA
      Integer hi, i, j, NCALC, IZE
      Double Precision, parameter :: Pi = 3.141592653589793d0

      ! For calling RIBESL
      ALPHA = 0.5 * (nq - 1) - (kk - 1)
      IZE = 1

      ! Accumulate the pairwise sums directly instead of storing the full
      ! n x n x lenh arrays (memory O(lenh)).
      sumCV1 = 0
      sumCV2 = 0

      ! Computation of the CV1 and CV2 sums
      do i = 2, n
        do j = 1, i - 1

           dotij = dot_product(data_dir(i, :), data_dir(j, :))
           aux2 = Sqrt(Sum((data_dir(i, :) + data_dir(j, :)) ** 2))

           do hi = 1, lenh

                ! CV1: exp(dotij / h^2), avoiding the exp-then-power round trip
                sumCV1(hi) = sumCV1(hi) + exp(dotij / (h(hi) ** 2))

                ! CV2
                aux4 = aux2 / (h(hi) ** 2)
                call RIBESL(aux4, ALPHA, kk, IZE, B, NCALC)
                aux5 = ((2 * Pi) ** (0.5 * (nq + 1))) * B(kk)
                sumCV2(hi) = sumCV2(hi) + aux5 / (aux4 ** (0.5 * (nq - 1)))

           end do

        end do
      end do

      CV = (4 * Cq) / (n * (n - 1)) * sumCV1
      CV = CV - (Cq ** 2) / (n * Cq2)
      CV = CV - (2 * Cq ** 2) / (n ** 2) * sumCV2
      CV = -CV

      Return
      End
