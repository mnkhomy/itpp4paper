module all_var
    integer(4) N_full, Nt, Nx, Ny, Nz, NxNy, NxNyNz, Neq, NxC, NyC, NzC
    integer(4) vec_step, save_step, save_step_planes, nnorm, save_raw_mnmq, shape, last_saved_plane
    real(8) T_end, ht, pi, tet_zero, phi_zero
    real(8) Lx, hx, Ly, hy, Lz, hz, Rbig, hx2, hy2, hz2, Lbig, K3big
    real(8) x0, y0, z0, n0x, n0y, n0z, rvi, rvi2, wn, wn2, wq2, wE
    real(8) bdn_a, bdn_b, bdn_c, bdq_a, bdq_b, bdq_c
    integer(4) bdn_func, bdq_func

    real(8) k1, k2, k3, k23, k4, k24, Dpar, Dper, Da, mu1, mu3, mu13, gam, A, B, C
    real(8) Rsph, Rsph_in, Rcur, dR, Rsph2, R2, qh, Q0, Qroot, Qroot2, kap2, gtet, gphi, eps
    real(8) h_small

    real(8), allocatable :: mvecni(:,:,:,:), mQ(:,:,:,:), mvecE(:,:,:,:)
    real(8) m_tmp1(1), m_tmp2(1)
    character(4) :: ss, ss_shape, bdn_func_ss, bdq_func_ss

    integer(4) scheme

    end module all_var

