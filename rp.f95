module rp_mod
    use all_var
    use dif_oper
    use init_mod
    implicit none
    
    real(8), allocatable, private :: grad_Q(:,:,:,:), rot_n(:,:,:,:), div_n(:,:,:,:)
    real(8), allocatable, private :: skobka_1(:,:,:,:), skobka_2(:,:,:,:), skobka_3(:,:,:,:)
    real(8), allocatable, private :: grad_nx(:,:,:,:), grad_ny(:,:,:,:), grad_nz(:,:,:,:)
    !real(8), allocatable, private :: vec_comp(:,:,:)
    real(8), allocatable, private :: mvecni_next(:,:,:,:), mQ_next(:,:,:,:)
    real(8), allocatable, private :: dndx(:,:,:,:), dndy(:,:,:,:), dndz(:,:,:,:)
    real(8), allocatable, private :: skobka_P(:,:,:,:)

    
    contains 
    
    subroutine init_mass()
        allocate(grad_Q(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(rot_n(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(div_n(1, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(skobka_1(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(skobka_2(1, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(skobka_3(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(skobka_P(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(grad_nx(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(grad_ny(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(grad_nz(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(dndx(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(dndy(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        allocate(dndz(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        dndx = 0d0
        dndy = 0d0
        dndz = 0d0
        !allocate(vec_comp(0:Nx+1, 0:Ny+1, 0:Nz+1))
        if (scheme == 2) then
            allocate(mvecni_next(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
            allocate(mQ_next(1, 0:Nx+1, 0:Ny+1, 0:Nz+1))
        endif
    end subroutine init_mass
    
    subroutine boundary_U(Ui, xl, xr, yl, yr, zl, zr, edge)
        implicit none
        real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1)
        real(8) xl, xr, yl, yr, zl, zr
        integer(4) ix, iy, iz, edge
    
        do iy = edge, Ny+1-edge
            do iz = edge, Nz+1-edge
                Ui(1, edge, iy, iz) = xl
                Ui(1, Nx+1-edge, iy, iz) = xr
            end do
        end do
    
        do ix = edge, Nx+1-edge
            do iz = edge, Nz+1-edge
                Ui(1, ix, edge, iz) = yl
                Ui(1, ix, Ny+1-edge, iz) = yr
            end do
        end do
    
        do ix = edge, Nx+1-edge
            do iy = edge, Ny+1-edge
                Ui(1, ix, iy, edge) = zl
                Ui(1, ix, iy, Nz+1-edge) = zr
            end do
        end do
    
    end subroutine boundary_U
    
    subroutine bound_Utrig(Ui, edge, estim)
        use all_var
        implicit none
        real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1)
        integer(4) edge
        real(8) estim
        external estim
    
        integer(4) ix, iy, iz
        real(8) xx, x1, x2, yy, y1, y2, zz, z1, z2
    
        x1 = (edge-1)*hx - x0
        x2 = (Nx+1-edge-1)*hx - x0
        do iy = edge, Ny+1-edge
            yy = (iy-1)*hy - y0
            do iz = edge, Nz+1-edge
                zz = (iz-1)*hz - z0
                Ui(1, edge, iy, iz) = estim(x1, yy, zz)
                Ui(1, Nx+1-edge, iy, iz) = estim(x2, yy, zz)
            end do
        end do
    
        y1 = (edge-1)*hy - y0
        y2 = (Ny+1-edge-1)*hy - y0
        do ix = edge, Nx+1-edge
            xx = (ix-1)*hx - x0
            do iz = edge, Nz+1-edge
                zz = (iz-1)*hz - z0
                Ui(1, ix, edge, iz) = estim(xx, y1, zz)
                Ui(1, ix, Ny+1-edge, iz) = estim(xx, y2, zz)
            end do
        end do
    
        z1 = (edge-1)*hz - z0
        z2 = (Nz+1-edge-1)*hz - z0
        do ix = edge, Nx+1-edge
            xx = (ix-1)*hx - x0
            do iy = edge, Ny+1-edge
                yy = (iy-1)*hy - y0
                Ui(1, ix, iy, edge) = estim(xx, yy, z1)
                Ui(1, ix, iy, Nz+1-edge) = estim(xx, yy, z2)
            end do
        end do
    
    end subroutine bound_Utrig
    
    subroutine rp(N, t, U, dU_dt)
        use all_var
        implicit none
        integer(4) N
        real(8) t
        real(8) U(N), dU_dt(N)        
    
        integer(4) ix, iy, iz, edge_do, edge_difoper, edge_bound
        real(8) x, y, z, Qi, Qi2, Uq, lambd  !hx2, hy2, hz2, d2U_dx2, d2U_dy2, d2U_dz2
        real(8) grad_Qi(3), vec_ni(3), rot_ni(3), grad_div_n(3), grad_skobka_2(3), rot_skobka_3(3), gradQ_cross_rotn(3), Funcn(3)
        real(8) vec_E(3)
        real(8) n_dot_grad_Q, n_dot_rot_n, div_skobka_1, n_dot_grad_div_n, n_cross_rot_n(3), n_cross_grad_Q(3)
        real(8) n_cross_rot_n2, div_ni, div_ni2, div_P, Lagr, n_dot_E, n_dot_E2 !, mod_ni !, rvi
        real(8) dnidnk, dn2x, dn2y, dn2z, absq
        !real(4) time_rp

        edge_bound = 0
        if (scheme == 1) call line2tensor(U, mvecni, mQ, nnorm)
   
        if (shape == 1) then !сфера
            bdn_a = 0d0; bdn_b = 1d0; bdn_c = 0d0;
            bdq_a = 0d0; bdq_b = 1d0; bdq_c = 0d0;
            bdq_func = 0
            if (bdq_func == 0) then
                call bound_all_abc(mvecni(1,:,:,:), bdn_a, bdn_b, bdn_c, edge_bound, zero_all)
                call bound_all_abc(mvecni(2,:,:,:), bdn_a, bdn_b, bdn_c, edge_bound, zero_all)
                call bound_all_abc(mvecni(3,:,:,:), bdn_a, bdn_b, bdn_c, edge_bound, zero_all)
                call bound_all_abc(mQ, bdq_a, bdq_b, bdq_c, edge_bound, zero_all)
            elseif (bdq_func == 100) then
                call bound_all_abc(mvecni(1,:,:,:), bdn_a, bdn_b, bdn_c, edge_bound, init_nx)
                call bound_all_abc(mvecni(2,:,:,:), bdn_a, bdn_b, bdn_c, edge_bound, init_ny)
                call bound_all_abc(mvecni(3,:,:,:), bdn_a, bdn_b, bdn_c, edge_bound, init_nz)
                call bound_all_abc(mQ, bdq_a, bdq_b, bdq_c, edge_bound, zero_all)
            endif
        end if

        if (shape == 2) then !полусфера
            call bound_all_abc(mvecni(1,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
            call bound_all_abc(mvecni(2,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
            call bound_all_abc(mvecni(3,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
            call bound_all_abc(mQ, 0d0, 1d0, 0d0, edge_bound, zero_all)

            call bound_abc_zmin(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, one_all)
            call bound_abc_zmin(mQ, 1d0, 0d0, 0d0, edge_bound, one_all)
        endif

        if (shape == 3) then !цилиндр
            call bound_all_abc(mvecni(1,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
            call bound_all_abc(mvecni(2,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
            call bound_all_abc(mvecni(3,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
            call bound_all_abc(mQ, 0d0, 1d0, 0d0, edge_bound, zero_all)

            call bound_abc_xmin(mvecni(1,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nx)
            call bound_abc_xmin(mvecni(2,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_ny)
            call bound_abc_xmin(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nz)

            call bound_abc_xmax(mvecni(1,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nx)
            call bound_abc_xmax(mvecni(2,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_ny)
            call bound_abc_xmax(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nz)

            call bound_abc_ymin(mvecni(1,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nx)
            call bound_abc_ymin(mvecni(2,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_ny)
            call bound_abc_ymin(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nz)

            call bound_abc_ymax(mvecni(1,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nx)
            call bound_abc_ymax(mvecni(2,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_ny)
            call bound_abc_ymax(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, init_nz)
        endif

        !call bound_abc_zmax(mvecni(1,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
        !call bound_abc_zmin(mvecni(1,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
        !call bound_abc_zmax(mvecni(2,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
        !call bound_abc_zmin(mvecni(2,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
        !call bound_abc_zmax(mvecni(3,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)
        !call bound_abc_zmin(mvecni(3,:,:,:), 0d0, 1d0, 0d0, edge_bound, zero_all)

        ! call bound_abc_zmax(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, zero_all)
        ! call bound_abc_zmin(mvecni(3,:,:,:), 1d0, 0d0, 0d0, edge_bound, one_all)
        ! call bound_abc_zmin(mQ, 1d0, 0d0, 0d0, edge_bound, one_all)

        edge_do = 0
        edge_difoper = 0
    
    !$omp parallel shared(skobka_1, skobka_2, skobka_3, skobka_P, grad_Q, rot_n, div_n, dndx, dndy, dndz)
    !$omp do private(ix, iy, iz, x, y, z, Qi, Qi2, grad_Qi, vec_ni, n_dot_grad_Q, rot_ni, n_dot_rot_n, n_cross_grad_Q, n_cross_rot_n, div_ni)
    
        do iz = edge_do, Nz+1-edge_do
            z = (iz-1)*hz - z0
            do iy = edge_do, Ny+1-edge_do
                y = (iy-1)*hy - y0
                do ix = edge_do, Nx+1-edge_do
                    x = (ix-1)*hx - x0
    !                rvi = dsqrt(x**2 + y**2 + z**2)
    !                if (rvi < Rbig) then
                        vec_ni(1) = mvecni(1, ix, iy, iz)
                        vec_ni(2) = mvecni(2, ix, iy, iz)
                        vec_ni(3) = mvecni(3, ix, iy, iz)
                        Qi = mQ(1, ix, iy, iz)
                        Qi2 = Qi**2
    
        !               !grad_Q
                        grad_Qi = fgrad(mQ, ix, iy, iz, edge_difoper)
                        grad_Q(1, ix, iy, iz) = grad_Qi(1)
                        grad_Q(2, ix, iy, iz) = grad_Qi(2)
                        grad_Q(3, ix, iy, iz) = grad_Qi(3)
    
                        !n_dot_grad_Q
                        n_dot_grad_Q = dot_product(vec_ni, grad_Qi)
    
                        !skobka_1
                        skobka_1(1, ix, iy, iz) = Dper * grad_Qi(1) + Da * vec_ni(1) * n_dot_grad_Q
                        skobka_1(2, ix, iy, iz) = Dper * grad_Qi(2) + Da * vec_ni(2) * n_dot_grad_Q
                        skobka_1(3, ix, iy, iz) = Dper * grad_Qi(3) + Da * vec_ni(3) * n_dot_grad_Q
    
                        !div_n
                        div_n(1, ix, iy, iz) = fdiv(mvecni, ix, iy, iz, edge_difoper)
                        div_ni = div_n(1, ix, iy, iz)
    
                        !rot_n
                        rot_ni = frot(mvecni, ix, iy, iz, edge_difoper)
                        rot_n(1, ix, iy, iz) = rot_ni(1)
                        rot_n(2, ix, iy, iz) = rot_ni(2)
                        rot_n(3, ix, iy, iz) = rot_ni(3)
    
                        !n_dot_rot_n
                        n_dot_rot_n = dot_product(vec_ni, rot_ni)

                        !skobka_2
                        !в соответствии с уточненными уравнениями от 2023-03-07: (k1-k24)
                        skobka_2(1, ix, iy, iz) = (k1-k24) * Qi2 * div_n(1, ix, iy, iz) + mu1 * Qi * n_dot_grad_Q 
    
                        !skobka_3
                        !в соответствии с уточненными уравнениями от 2023-03-07: (k3+k24)
                        skobka_3(1, ix, iy, iz) = Qi2 * (k23 * vec_ni(1) * n_dot_rot_n + (k3 + k24) * rot_ni(1))
                        skobka_3(2, ix, iy, iz) = Qi2 * (k23 * vec_ni(2) * n_dot_rot_n + (k3 + k24) * rot_ni(2))
                        skobka_3(3, ix, iy, iz) = Qi2 * (k23 * vec_ni(3) * n_dot_rot_n + (k3 + k24) * rot_ni(3))
    
                        !n_cross_grad_Q
                        n_cross_grad_Q = fcross(mvecni, grad_Q, ix, iy, iz)
    
                        skobka_3(1, ix, iy, iz) = skobka_3(1, ix, iy, iz) + mu3*Qi*n_cross_grad_Q(1) + 2*k2*Qi2*qh*vec_ni(1)
                        skobka_3(2, ix, iy, iz) = skobka_3(2, ix, iy, iz) + mu3*Qi*n_cross_grad_Q(2) + 2*k2*Qi2*qh*vec_ni(2)
                        skobka_3(3, ix, iy, iz) = skobka_3(3, ix, iy, iz) + mu3*Qi*n_cross_grad_Q(3) + 2*k2*Qi2*qh*vec_ni(3)

                        n_cross_rot_n = fcross(mvecni, rot_n, ix, iy, iz)
                        skobka_P(1, ix, iy, iz) = mu1 * vec_ni(1) * div_ni - mu3 * n_cross_rot_n(1)
                        skobka_P(2, ix, iy, iz) = mu1 * vec_ni(2) * div_ni - mu3 * n_cross_rot_n(2)
                        skobka_P(3, ix, iy, iz) = mu1 * vec_ni(3) * div_ni - mu3 * n_cross_rot_n(3)

                        !function fdUdx(Ui, ic, ix, iy, iz, edge) result (fdUdx_)
                        !fdUdx(mvecni, 1, ix, iy, iz, edge_difoper)
                        dndx(1, ix, iy, iz) = Qi2 * fdUdx(mvecni, 1, ix, iy, iz, edge_difoper) !dnx_dx
                        dndx(2, ix, iy, iz) = Qi2 * fdUdx(mvecni, 2, ix, iy, iz, edge_difoper) !dny_dx
                        dndx(3, ix, iy, iz) = Qi2 * fdUdx(mvecni, 3, ix, iy, iz, edge_difoper) !dnz_dx

                        dndy(1, ix, iy, iz) = Qi2 * fdUdy(mvecni, 1, ix, iy, iz, edge_difoper) !dnx_dy
                        dndy(2, ix, iy, iz) = Qi2 * fdUdy(mvecni, 2, ix, iy, iz, edge_difoper) !dny_dy
                        dndy(3, ix, iy, iz) = Qi2 * fdUdy(mvecni, 3, ix, iy, iz, edge_difoper) !dnz_dy

                        dndz(1, ix, iy, iz) = Qi2 * fdUdz(mvecni, 1, ix, iy, iz, edge_difoper) !dnx_dz
                        dndz(2, ix, iy, iz) = Qi2 * fdUdz(mvecni, 2, ix, iy, iz, edge_difoper) !dny_dz
                        dndz(3, ix, iy, iz) = Qi2 * fdUdz(mvecni, 3, ix, iy, iz, edge_difoper) !dnz_dx

    !                end if
                end do
            end do
        end do
    
    !$omp enddo
    !$omp endparallel
    
        edge_difoper = 0 !0
        edge_do = 1
    
    !$omp parallel shared(dU_dt, mvecni_next, mQ_next)
    !$omp do private(ix, iy, iz, x, y, z, Qi, Qi2, vec_ni, Uq, grad_div_n, grad_skobka_2, vec_E, &
    !$omp              rot_skobka_3, gradQ_cross_rotn, n_dot_rot_n, n_dot_grad_Q, n_dot_E, div_ni, dnidnk, &
    !$omp              Funcn, lambd, div_skobka_1, n_dot_grad_div_n, n_cross_rot_n, n_cross_rot_n2, div_ni2, div_P, Lagr, &
    !$omp              dn2x, dn2y, dn2z)
        
        
        do iz = edge_do, Nz+1-edge_do
            z = (iz-1)*hz - z0
            do iy = edge_do, Ny+1-edge_do
                y = (iy-1)*hy - y0
                do ix = edge_do, Nx+1-edge_do
                    x = (ix-1)*hx - x0
    !                rvi = dsqrt(x**2 + y**2 + z**2)
    !                if (rvi < Rbig) then
    
                        vec_ni(1) = mvecni(1, ix, iy, iz)
                        vec_ni(2) = mvecni(2, ix, iy, iz)
                        vec_ni(3) = mvecni(3, ix, iy, iz)
    
                        vec_E(1) = mvecE(1, ix, iy, iz)
                        vec_E(2) = mvecE(2, ix, iy, iz)
                        vec_E(3) = mvecE(3, ix, iy, iz)
    
                        Qi = mQ(1, ix, iy, iz)
                        Qi2 = Qi**2
                        Uq = A*Qi - B*Qi2 + C*Qi2*Qi
    
                        !grad_div_n
                        grad_div_n = fgrad(div_n, ix, iy, iz, edge_difoper)
    
                        !grad_skobka_2
                        grad_skobka_2 = fgrad(skobka_2, ix, iy, iz, edge_difoper)
    
                        !rot_skobka_3
                        rot_skobka_3 = frot(skobka_3, ix, iy, iz, edge_difoper)
    
                        !gradQ_cross_rotn
                        gradQ_cross_rotn = fcross(grad_Q, rot_n, ix, iy, iz)
    
                        !n_dot_rot_n
                        n_dot_rot_n = fdot(mvecni, rot_n, ix, iy, iz)
    
                        !n_dot_grad_Q
                        n_dot_grad_Q = fdot(mvecni, grad_Q, ix, iy, iz)
    
                        !n_dot_E
                        n_dot_E = dot_product(vec_ni, vec_E)
                        n_dot_E2 = n_dot_E**2
    
                        div_ni = div_n(1, ix, iy, iz)

                        !function fdUdx(Ui, ic, ix, iy, iz, edge) result (fdUdx_)
                        !fdUdx(mvecni, 1, ix, iy, iz, edge_difoper)
                        dn2x = fdUdx(dndx, 1, ix, iy, iz, edge_difoper)
                        dn2x = dn2x + fdUdy(dndx, 2, ix, iy, iz, edge_difoper)
                        dn2x = dn2x + fdUdz(dndx, 3, ix, iy, iz, edge_difoper)

                        dn2y = fdUdx(dndy, 1, ix, iy, iz, edge_difoper)
                        dn2y = dn2y + fdUdy(dndy, 2, ix, iy, iz, edge_difoper)
                        dn2y = dn2y + fdUdz(dndy, 3, ix, iy, iz, edge_difoper)

                        dn2z = fdUdx(dndz, 1, ix, iy, iz, edge_difoper)
                        dn2z = dn2z + fdUdy(dndz, 2, ix, iy, iz, edge_difoper)
                        dn2z = dn2z + fdUdz(dndz, 3, ix, iy, iz, edge_difoper)

                        h_small = 1d0
                        
                        Funcn(1) = grad_skobka_2(1) - rot_skobka_3(1)
                        Funcn(1) = Funcn(1) - k23 * Qi2 * n_dot_rot_n * rot_n(1, ix, iy, iz)
                        Funcn(1) = Funcn(1) - Da * grad_Q(1, ix, iy, iz) * n_dot_grad_Q
                        Funcn(1) = Funcn(1) - mu1 * Qi * grad_Q(1, ix, iy, iz) * div_ni
                        Funcn(1) = Funcn(1) + mu3 * Qi * gradQ_cross_rotn(1)
                        Funcn(1) = Funcn(1) + kap2 * vec_E(1) * n_dot_E + k24 * dn2x * h_small
                        
                        Funcn(2) = grad_skobka_2(2) - rot_skobka_3(2)
                        Funcn(2) = Funcn(2) - k23 * Qi2 * n_dot_rot_n * rot_n(2, ix, iy, iz)
                        Funcn(2) = Funcn(2) - Da * grad_Q(2, ix, iy, iz) * n_dot_grad_Q
                        Funcn(2) = Funcn(2) - mu1 * Qi * grad_Q(2, ix, iy, iz) * div_ni
                        Funcn(2) = Funcn(2) + mu3 * Qi * gradQ_cross_rotn(2)
                        Funcn(2) = Funcn(2) + kap2 * vec_E(2) * n_dot_E + k24 * dn2y * h_small
                        
                        Funcn(3) = grad_skobka_2(3) - rot_skobka_3(3)
                        Funcn(3) = Funcn(3) - k23 * Qi2 * n_dot_rot_n * rot_n(3, ix, iy, iz)
                        Funcn(3) = Funcn(3) - Da * grad_Q(3, ix, iy, iz) * n_dot_grad_Q
                        Funcn(3) = Funcn(3) - mu1 * Qi * grad_Q(3, ix, iy, iz) * div_ni
                        Funcn(3) = Funcn(3) + mu3 * Qi * gradQ_cross_rotn(3)
                        Funcn(3) = Funcn(3) + kap2 * vec_E(3) * n_dot_E + k24 * dn2z * h_small
                        
                        lambd = dot_product(vec_ni, Funcn) 
                        
                        !div_skobka_1
                        div_skobka_1 = fdiv(skobka_1, ix, iy, iz, edge_difoper)
    
                        !n_dot_grad_div_n
                        n_dot_grad_div_n = dot_product(vec_ni, grad_div_n)
    
                        !n_cross_rot_n
                        n_cross_rot_n = fcross(mvecni, rot_n, ix, iy, iz)
                        n_cross_rot_n2 = dot_product(n_cross_rot_n, n_cross_rot_n)
    
                        !div_P
                        div_ni2 = div_ni**2
                        div_P = fdiv(skobka_P, ix, iy, iz, edge_difoper) !mu13 * (div_ni2 + n_dot_grad_div_n)
    
                        !поверхностные члены
                        dnidnk =          fdUdx(mvecni, 1, ix, iy, iz, edge_difoper)  * fdUdx(mvecni, 1, ix, iy, iz, edge_difoper)
                        dnidnk = dnidnk + fdUdy(mvecni, 1, ix, iy, iz, edge_difoper)  * fdUdx(mvecni, 2, ix, iy, iz, edge_difoper)
                        dnidnk = dnidnk + fdUdz(mvecni, 1, ix, iy, iz, edge_difoper)  * fdUdx(mvecni, 3, ix, iy, iz, edge_difoper)
    
                        dnidnk = dnidnk + fdUdx(mvecni, 2, ix, iy, iz, edge_difoper)  * fdUdy(mvecni, 1, ix, iy, iz, edge_difoper)
                        dnidnk = dnidnk + fdUdy(mvecni, 2, ix, iy, iz, edge_difoper)  * fdUdy(mvecni, 2, ix, iy, iz, edge_difoper)
                        dnidnk = dnidnk + fdUdz(mvecni, 2, ix, iy, iz, edge_difoper)  * fdUdy(mvecni, 3, ix, iy, iz, edge_difoper)
    
                        dnidnk = dnidnk + fdUdx(mvecni, 3, ix, iy, iz, edge_difoper)  * fdUdz(mvecni, 1, ix, iy, iz, edge_difoper)
                        dnidnk = dnidnk + fdUdy(mvecni, 3, ix, iy, iz, edge_difoper)  * fdUdz(mvecni, 2, ix, iy, iz, edge_difoper)
                        dnidnk = dnidnk + fdUdz(mvecni, 3, ix, iy, iz, edge_difoper)  * fdUdz(mvecni, 3, ix, iy, iz, edge_difoper)
    
                        Lagr = k1 * div_ni2 + k2 * (n_dot_rot_n + qh)**2 + k3 * n_cross_rot_n2 + k24*(dnidnk - div_ni2)

                        !mvecni = tempn
                        !mQ = tempq
                        !Qi2 = Qi2 + eps

                        absq = mQ(1,ix,iy,iz) + ((div_skobka_1 + Qi * (div_P - Lagr) - Uq) / gam)*ht
                            
                        if (dabs(absq) > 1d-8) then
                            mvecni_next(1, ix, iy, iz) = mvecni(1,ix,iy,iz) + ((Funcn(1) - lambd * vec_ni(1)) / (gam))*ht
                            mvecni_next(2, ix, iy, iz) = mvecni(2,ix,iy,iz) + ((Funcn(2) - lambd * vec_ni(2)) / (gam))*ht
                            mvecni_next(3, ix, iy, iz) = mvecni(3,ix,iy,iz) + ((Funcn(3) - lambd * vec_ni(3)) / (gam))*ht
                        else
                            mvecni_next(1, ix, iy, iz) = mvecni(1,ix,iy,iz)
                            mvecni_next(2, ix, iy, iz) = mvecni(2,ix,iy,iz)
                            mvecni_next(3, ix, iy, iz) = mvecni(3,ix,iy,iz)
                        endif
                        mQ_next(1, ix, iy, iz) = absq
    
    !                else
    !
    !                    tet = dacos(z/rvi) * (1 - dexp(-rvi2/wn2))
    !                    phi = datan2(y,x) * (1 - dexp(-rvi2/wn2))
    !
    !                    mvecni_next(1, ix, iy, iz) = dsin(tet) * dcos(phi)       !nx
    !                    mvecni_next(2, ix, iy, iz) = dsin(tet) * dsin(phi)       !ny
    !                    mvecni_next(3, ix, iy, iz) = dcos(tet)       !nz
    !                    mQ_next(1, ix, iy, iz) = 1d0 - dexp(-rvi2/wq2)       !Q
    
    !                    dU_dt(get_i(1, ix, iy, iz)) = 0d0
    !                    dU_dt(get_i(2, ix, iy, iz)) = 0d0
    !                    dU_dt(get_i(3, ix, iy, iz)) = 0d0
    !                    dU_dt(get_i(4, ix, iy, iz)) = 0d0
    
    !                end if
                end do
            end do
        end do
    
    !$omp enddo
    !$omp endparallel
    
        mvecni = mvecni_next
        mQ = mQ_next
        if (nnorm >= 1) then 
            call normalize(mvecni, mQ)
        endif
 
    end subroutine rp
end module rp_mod
    