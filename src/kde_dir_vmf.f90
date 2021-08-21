      Subroutine kde_dir_vmf(x, data_dir, h, ch, n, q, nx, kde)
      implicit none

      ! Arguments
      Integer n, q, nx
      Double Precision x(nx, q + 1), data_dir(n, q + 1)
      Double Precision h, ch
      Double Precision kde(nx)

      ! Local variables
      Integer i
      Double Precision hh

      ! Square bandwidth
      hh = h * h

      ! Computation of kernel estimator
      do i = 1, nx
        kde(i) = sum(exp(-(1 - matmul(data_dir, x(i, :))) / hh))
      end do

      ! Normalize
      kde = (ch / n) * kde

      Return
      End
