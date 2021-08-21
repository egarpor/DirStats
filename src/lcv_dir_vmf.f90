      Subroutine lcv_dir_vmf(data_dir, h, ch, n, nq, lenh, CV)
      implicit none

      ! Arguments
      Integer n, nq, lenh
      Double Precision data_dir(n, nq + 1)
      Double Precision h(lenh), ch(lenh)
      Double Precision CV(lenh)

      ! Local variables
      Double Precision mat_L(n, n)
      Double Precision aux1, aux2(n, n)
      Integer hi, i
      Real, parameter :: Pi = 3.1415926536

      ! Computation of kernel matrix
      mat_L = exp(-1 + matmul(data_dir, transpose(data_dir)))
      do i = 1, n
        mat_L(i, i)=0
      end do

      ! Computation of CV
      do hi = 1, lenh

        aux1 = ch(hi) / (n - 1)
        aux2 = aux1 * (mat_L ** (h(hi) ** (-2.0)))
        CV(hi) = -Sum(log(Sum(aux2, dim = 1))) / n

      end do

      Return
      End
