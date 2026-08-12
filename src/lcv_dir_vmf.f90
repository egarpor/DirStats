      Subroutine lcv_dir_vmf(data_dir, h, ch, n, nq, lenh, CV)
      implicit none

      ! Arguments
      Integer n, nq, lenh
      Double Precision data_dir(n, nq + 1)
      Double Precision h(lenh), ch(lenh)
      Double Precision CV(lenh)

      ! Local variables
      Double Precision M(n, n)
      Double Precision aux1, aux2(n, n)
      Integer hi, i

      ! Log-kernel matrix M = log(exp(-1 + <x_i, x_j>)) to avoid recomputing
      M = matmul(data_dir, transpose(data_dir)) - 1

      ! Computation of CV
      do hi = 1, lenh

        aux1 = ch(hi) / (n - 1)
        aux2 = exp(M / (h(hi) ** 2))

        ! Leave-one-out: drop the self-contributions
        do i = 1, n
          aux2(i, i) = 0
        end do

        CV(hi) = -Sum(log(aux1 * Sum(aux2, dim = 1))) / n

      end do

      Return
      End
