module init_mod

contains

function init_nx(rx, ry, rz) result (init_nx_)
    use all_var
    real(8) rx, ry, rz
    real(8) an(3)
    real(8) init_nx_

    an = init_n(rx, ry, rz)
    init_nx_ = an(1)    !nx

 end function init_nx


 function init_ny(rx, ry, rz) result (init_ny_)
    use all_var
    real(8) rx, ry, rz
    real(8) an(3)
    real(8) init_ny_

    an = init_n(rx, ry, rz) 
    init_ny_ = an(2)    !ny

 end function init_ny


 function init_nz(rx, ry, rz) result (init_nz_)
    use all_var
    real(8) rx, ry, rz
    real(8) an(3)
    real(8) init_nz_

    an = init_n(rx, ry, rz)
    init_nz_ = an(3)    !nz

end function init_nz

function init_n(rx, ry, rz) result (n_)
    use all_var
    real(8) rx, ry, rz
    real(8) n_(3)

    real(8) r_, dr_, tet_, phi_, r0_, phi0_, tet0_

    if (shape == 1) then !сфера
        tet_zero = 90d0
        phi_zero = 0d0

        tet0_ = tet_zero * pi / 180d0
        phi0_ = phi_zero * pi / 180d0 !pi/3d0 !80*pi/180 !3d0 !0.5d0 !pi/4d0
        !r0_ = 0.8*(Lx/2d0)    !Для кубика все r0_ равны

        r_ = dsqrt(rx**2 + ry**2 + rz**2)
        !dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)
        tet_ = dacos(rz/(r_+eps)) + tet0_ !pi/2d0 !dacos(rz/(r_+eps)) + tet0_ !pi / 2d0  !dacos(rz/(r_+eps)) + tet0_ !pi /2d0 !dacos(rz/(r_+eps)) + tet0_ !(pi/2d0)*(rz+0.5d0)  !dacos(rz/(r_+1d-5)) !* (1 - dexp(-rvi2/wn2))
        phi_ = datan2(ry,rx) + phi0_ !datan2(ry,rx) + phi0_ !* (0.5*erf((r_-r0_)/wn)+0.5) !dexp(-(jrvi-r0)**2/wn2)  !* (1 - dexp(-rvi2/wn2))

        n_(1) = dsin(tet_) * dcos(phi_) !nx
        n_(2) = dsin(tet_) * dsin(phi_) !ny
        n_(3) = dcos(tet_) !nz
    
    elseif (shape == 2) then !полусфера
        phi0_ = pi/2d0
        tet0_ = pi/2d0
        r0_ = 0.8*(Lx/2d0)    !Для кубика все r0_ равны

        r_ = dsqrt(rx**2 + ry**2 + rz**2)
        dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)
        tet_ = dacos(rz/(r_+eps)) + tet0_ !(pi/2d0)*(rz+0.5d0)  !dacos(rz/(r_+1d-5)) !* (1 - dexp(-rvi2/wn2))
        phi_ = datan2(ry,rx) !+ phi0_ !datan2(ry,rx) + phi0_ !* (0.5*erf((r_-r0_)/wn)+0.5) !dexp(-(jrvi-r0)**2/wn2)  !* (1 - dexp(-rvi2/wn2))

        n_(1) = dsin(tet_) * dcos(phi_) !nx
        n_(2) = dsin(tet_) * dsin(phi_) !ny
        n_(3) = dcos(tet_) !nz

    elseif (shape == 3) then !цилиндр
        r_ = dsqrt(rx**2 + ry**2)
        phi0_ = 0d0 !pi/2d0
        tet0_ = pi/2d0
        phi_ = datan2(ry,rx) + phi0_
        tet_ = tet0_
        n_(1) = dsin(tet_) * dcos(phi_)  !0d0 !rx/(r_+eps) !nx
        n_(2) = dsin(tet_) * dsin(phi_) !0d0 !ry/(r_+eps) !ny
        n_(3) = dcos(tet_) !1d0 !nz

    endif

end function init_n


function init_q(rx, ry, rz) result (init_q_)
    use all_var
    real(8) rx, ry, rz
    integer(4) forma
    real(8) init_q_

    real(8) r_, r2_, tet_, phi_

    if (shape == 1) then !сфера
        r2_ = rx**2 + ry**2 + rz**2
        !r_ = dsqrt(r2_)
        !tet_ = dacos(rz/(r_+eps)) !* (1 - dexp(-rvi2/wn2))
        !phi_ = datan2(ry,rx) !* (1 - dexp(-rvi2/wn2))

        init_q_ = Qroot !* (1 - dexp(-r2_/wq2)) !*Qroot !Q
    
    elseif (shape == 2) then !полусфера
        r2_ = rx**2 + ry**2 + rz**2
        !init_q_ = (1 - dexp(-r2_/wq2))*Qroot !Q
        init_q_ = Qroot
    
    elseif (shape == 3) then !цилиндр
        r2_ = rx**2 + ry**2
        !init_q_ = (1 - dexp(-r2_/wq2))*Qroot !Q
        init_q_ = Qroot

    endif

end function init_q


function init_E(rx, ry, rz, forma) result (e_)
    use all_var
    real(8) rx, ry, rz
    integer(4) forma    
    real(8) e_(3)

    real(8) r_, r0_, erf_

    if (forma == 1) then !сфера
        r0_ = 0.8*(Lx/2d0)    !Для кубика все r0_ равны
        r_ = dsqrt(rx**2 + ry**2 + rz**2)
        !dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)
        erf_ = 0.5d0*erf((r_-r0_)/wE)+0.5d0

        e_(1) = erf_ * (rx/(r_+eps)) !Ex
        e_(2) = erf_ * (ry/(r_+eps)) !Ey
        e_(3) = erf_ * (rz/(r_+eps)) !Ez
    
    elseif (forma == 2) then !полусфера
        r0_ = 0.8*(Lx/2d0)    !Для кубика все r0_ равны
        r_ = dsqrt(rx**2 + ry**2 + rz**2)
        !dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)
        erf_ = 0.5d0*erf((r_-r0_)/wE)+0.5d0

        e_(1) = erf_ * (rx/(r_+eps)) !Ex
        e_(2) = erf_ * (ry/(r_+eps)) !Ey
        e_(3) = erf_ * (rz/(r_+eps)) !Ez

    elseif (forma == 3) then !цилиндр
        r0_ = 0.8*(Lx/2d0)
        r_ = dsqrt(rx**2 + ry**2)

        erf_ = 0.5d0*erf((r_-r0_)/wE)+0.5d0
        e_(1) = erf_ * (rx/(r_+eps)) !Ex
        e_(2) = erf_ * (ry/(r_+eps)) !Ey
        e_(3) = 0d0 !Ez
    endif

end function init_E


subroutine init()
    use all_var
    use dif_oper
    implicit none
    

    integer(4) ix, iy, iz, i
    integer(4) N_saves_per_sec
    real(8) x, y, z
    real(8) nl(3), el(3)
    real(4) time_init

    !real(8) mod_ni

    time_init = secnds(0e0)

    pi = 4d0*datan(1d0)
    eps = 1d-8

    Neq = 4

    ! Можно сохранять в файл не каждый элемент массива mvecni,
    ! для удобства построения графиков-векторов
    vec_step = 1

    shape = 1           !1 - сфера, 2 - полусфера, 3 - цилиндр

    if (shape == 1) then
        Nx = 31
        Ny = Nx
        Nz = Nx

        Lx = 1d0
        Ly = 1d0
        Lz = 1d0

        x0 = Lx / 2d0
        y0 = Ly / 2d0
        z0 = Lz / 2d0

    elseif (shape == 2) then !полусфера
        Nx = 31
        Ny = Nx
        Nz = Nx/2 + 1

        Lx = 1d0
        Ly = 1d0
        Lz = 0.5d0

        x0 = Lx / 2d0
        y0 = Ly / 2d0
        z0 = 0d0

    elseif (shape == 3) then
        Nx = 101
        Ny = Nx
        Nz = 301

        Lx = 1d0
        Ly = 1d0
        Lz = 3d0 !1d0

        x0 = Lx / 2d0
        y0 = Ly / 2d0
        z0 = Lz / 2d0

    endif

    NxNy = Nx*Ny
    NxNyNz = NxNy*Nz
    N_full = Neq*NxNyNz !50000 !*1000

    
    hx = Lx / (Nx-1)
    hx2 = hx**2
    hy = Ly / (Ny-1)
    hy2 = hy**2
    hz = Lz / (Nz-1)
    hz2 = hz**2
    wn = 0.05d0
    wn2 = wn**2
    wq2 = (0.05d0)**2

    wE = 2d0*dsqrt((hx**2+hy**2+hz**2)/3d0)
    kap2 = -100d0 !35d0 !-35d0
    Rsph_in = 0.8*(Lx/2d0)
    Rsph = Rsph_in + wE

    gam = 1d0

    k1 = 0.75d0 !0.75d0 !1d0 !2d0
    k2 = 0.39d0 !0.39d0 !0.39d0 !0.1d0 !0.39d0 !3d-1
    k3 = 1.0d0 !1d0 !3d0
    k4 = 0d0 !0d0 !1d0
    k23 = k2-k3
    k24 = -0.075d0 ! k2 / 2d0 !0d0 !-0.075d0 !0d0 !0.2d0 !-0.075d0

    qh = 0d0 !2*pi*Lx / Rsph !0d0 !pi

!    mu1 = 0.05d0
!    mu3 = 0.1d0
    mu1 = 0.6d0 !0.05d0
    mu3 = 0d0 !0.8d0 !0.7d0
    mu13 = mu1 + mu3

    !Блинов
    A = -0.13d6 !-9.42d3 !-9.42d3 !-3d2
    B = 1.6d6 !1.16d5 !1.16d5 !1.2d2
    C = 3.9d6 !2.83d5 !2.83d5 !2d3

    !Kos
    !A = -0.172d6
    !B = 2.13d6
    !C = 3.73d6 

    Lbig = 1d-6 * 0.3d0

    K3big = 13.8d-12

    ! A = -0.5d2 
    ! B = 0d0
    ! C = 0.5d2

    Qroot = (B + dsqrt(B**2 - 4*A*C)) / (2*C)
    Qroot2 = Qroot**2

    A = (Qroot2 * A * Lbig**2) / K3big !-9.42d3 !-3d2
    B = (Qroot2 * B * Lbig**2) / K3big !1.16d5 !1.2d2
    C = (Qroot2 * C * Lbig**2) / K3big !2.83d5 !2d3

    !A = -0.5d2
    !B = 0d0
    !C = 0.5d2

    Dpar = 1d0 !1d0 !2d0
    Dper = 0.5d0 * Dpar
    Da = Dpar - Dper

    T_end = 5d0

    N_saves_per_sec = 200

    scheme = 2         !1 - неявная; 2 - явная;
    nnorm = 1           !делаем нормировку nx, ny, nz
    save_raw_mnmq = 0   !сохраняем массивы n и q по времени

    ht = 0.2d0*min(hx**2, hy**2, hz**2) / max(k1, k2, k3, Dpar, Dper, k3+k24, k1-k24)  !ht=1d-5 !T_end/(Nt-1)!0.001d0
    Nt = nint(T_end / ht) !Nt=1d5+1 !nint(T_end/ht)

    save_step = nint((Nt/T_end)/N_saves_per_sec)
    save_step_planes = save_step * 2
    
    print *, 'scheme', scheme
    print *, 'hx2', hx**2
    print *, 'ht', ht
    print *, 'T_end', T_end
    print *, 'Nt', Nt
    print *, 'save_step', save_step
    print *, 'save_step_planes', save_step_planes
    
    allocate(mvecE(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
    mvecE = 0d0

    do iz = 1, Nz
        z = (iz-1)*hz - z0
        do iy = 1, Ny
            y = (iy-1)*hy - y0
            do ix = 1, Nx
                x = (ix-1)*hx - x0

                el = init_E(x, y, z, shape)
                mvecE(1, ix, iy, iz) = el(1)
                mvecE(2, ix, iy, iz) = el(2)
                mvecE(3, ix, iy, iz) = el(3)
            end do
        end do
    end do

    allocate(mvecni(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1))
    allocate(mQ(1, 0:Nx+1, 0:Ny+1, 0:Nz+1))

    do ix = 0, Nx+1
        do iy = 0, Ny+1
            do iz = 0, Nz+1
                mvecni(1, ix, iy, iz) = 0d0
                mvecni(2, ix, iy, iz) = 0d0
                mvecni(3, ix, iy, iz) = 0d0
                mQ(1, ix, iy, iz) = 0d0
            end do
        end do
    end do

    do iz = 0, Nz+1
        z = (iz-1)*hz - z0
        do iy = 0, Ny+1
            y = (iy-1)*hy - y0
            do ix = 0, Nx+1
                x = (ix-1)*hx - x0

                nl = init_n(x, y, z)
                mvecni(1, ix, iy, iz) = nl(1)           !nx
                mvecni(2, ix, iy, iz) = nl(2)           !ny
                mvecni(3, ix, iy, iz) = nl(3)           !nz
                mQ(1, ix, iy, iz) = init_q(x, y, z)     !Q
            end do
        end do
    end do

    time_init = secnds(time_init)
    print *, "end init...it takes...", time_init, "...seconds"

end subroutine init
end module init_mod
