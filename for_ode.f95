module for_ode_mod
    contains

    subroutine save_nprj_xy(filename, radius)
        use all_var
        implicit none
        character(101) :: filename
        real(8) radius

        integer(4) ix_, iy_, iz_, miniz_
        real(8) x_, y_, z_, mindr_, minz_, dR_, Rcur_, nox_, noy_, noz_, mod_ni_
        !проекции компонентов вектора с верхней полусферы на плоскость xy
        !write(ch,'("./results/nprj_xy_out.txt")')
        !read (ch,'(a)') header
        open(58, file=filename)
        do iy_ = 1, Ny, vec_step
            y_ = (iy_-1)*hy - y0
            do ix_ = 1, Nx, vec_step
                x_ = (ix_-1)*hx - x0
                mindr_ = 1d3
                miniz_ = -1
                do iz_ = int(Nz/2)+1, Nz, vec_step
                    z_ = (iz_-1)*hz - z0            
                    Rcur_ = dsqrt(x_**2 + y_**2 + z_**2)
                    dR_ = radius - Rcur_
                    if ((dR_>0) .AND. (dR_ < mindr_)) then
                        mindr_ = dR_
                        minz_ = z_
                        miniz_ = iz_
                    end if    
                end do
                if (miniz_ >= 0) then
                    nox_ = mvecni(1, ix_, iy_, miniz_)
                    noy_ = mvecni(2, ix_, iy_, miniz_)
                    noz_ = mvecni(3, ix_, iy_, miniz_)
                    mod_ni_ = dsqrt(nox_**2 + noy_**2)
                    write(58,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                                x_, y_, minz_, nox_, noy_, noz_, mod_ni_
                else
                    nox_ = 0d0
                    noy_ = 0d0
                    noz_ = 0d0
                    minz_ = -1d0
                    mod_ni_ = 0d0
                    write(58,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                                x_, y_, minz_, nox_, noy_, noz_, mod_ni_
                end if
            end do
        end do
        close(58)
    end subroutine save_nprj_xy


    subroutine save_nprj_xz(filename, radius)
        use all_var
        implicit none
        character(101) :: filename
        real(8) radius

        integer(4) ix_, iy_, iz_, miniy_
        real(8) x_, y_, z_, mindr_, miny_, dR_, Rcur_, nox_, noy_, noz_, mod_ni_    
        !проекции компонентов вектора с верхней полусферы на плоскость xz
        !write(ch,'("./results/nprj_xz_out.txt")')
        !read (ch,'(a)') header
        open(58, file=filename)
        do iz_ = 1, Nz, vec_step
            z_ = (iz_-1)*hz - z0
            do ix_ = 1, Nx, vec_step
                x_ = (ix_-1)*hx - x0
                mindr_ = 1d3
                miniy_ = -1
                do iy_ = int(Ny/2)+1, Ny, vec_step
                    y_ = (iy_-1)*hy - y0            
                    Rcur_ = dsqrt(x_**2 + y_**2 + z_**2)
                    dR_ = radius - Rcur_
                    if ((dR_>0) .AND. (dR_ < mindr_)) then
                        mindr_ = dR_
                        miny_ = y_
                        miniy_ = iy_
                    end if    
                end do
                if (miniy_ >= 0) then
                    nox_ = mvecni(1, ix_, miniy_, iz_)
                    noy_ = mvecni(2, ix_, miniy_, iz_)
                    noz_ = mvecni(3, ix_, miniy_, iz_)
                    mod_ni_ = dsqrt(nox_**2 + noz_**2)
                    write(58,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                                x_, miny_, z_, nox_, noy_, noz_, mod_ni_
                else
                    nox_ = 0d0
                    noy_ = 0d0
                    noz_ = 0d0
                    miny_ = -1d0
                    mod_ni_ = 0d0
                    write(58,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                                x_, miny_, z_, nox_, noy_, noz_, mod_ni_
                end if
            end do
        end do
        close(58)
    end subroutine save_nprj_xz

    subroutine save_nvec_xy_z0(filename, nka)
        use all_var
        implicit none
        character(101) :: filename
        real(8) nka(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)

        integer(4) ix, iy, iz
        real(8) x_, y_, nox_, noy_

        open(16, file=filename)
        
        iz = int(Nz/2)+1
        do ix = 1, Nx, vec_step
            x_ = (ix-1)*hx - x0
            do iy = 1, Ny, vec_step
                y_ = (iy-1)*hy - y0
                nox_ = nka(1, ix, iy, iz)
                noy_ = nka(2, ix, iy, iz)
                !noz = mvecni(3, ix, iy, iz)
                write(16,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') x_, y_, nox_, noy_, dsqrt(nox_**2 + noy_**2)
            end do
        end do

        close(16)

    end subroutine save_nvec_xy_z0

    subroutine save_Qplane_xy_z0(filename, qushka)
        use all_var
        implicit none
        character(101) :: filename
        real(8) qushka(1, 0:Nx+1, 0:Ny+1, 0:Nz+1)

        integer(4) ix, iy, iz
        real(8) x_, y_
        open(16, file=filename)

        iz = int(Nz/2)+1
        do ix = 1, Nx
            x_ = (ix-1)*hx - x0
            do iy = 1, Ny
                y_ = (iy-1)*hy - y0
                write(16,'(F0.10, 1x, F0.10, 1x, F0.10)') x_, y_, qushka(1, ix, iy, iz)
            end do
            write(16, *)
        end do

        close(16)
    
    end subroutine save_Qplane_xy_z0

    subroutine save_transmit(filename, nka, Lcur)
        use all_var
        use dif_oper
        implicit none
        character(101) :: filename
        real(8) nka(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
        real(8) Lcur

        integer(4) ix, iy
        real(8) x_, y_
        real(8) transmit_(Nx, Ny)

        !subroutine find_transmittance(tshka, nka, Lmax)
        call find_transmittance(transmit_, nka, Lcur)

        open(16, file=filename)
        do ix = 1, Nx
            x_ = (ix-1)*hx - x0
            do iy = 1, Ny
                y_ = (iy-1)*hy - y0
                write(16,'(F0.10, 1x, F0.10, 1x, F0.10)') x_, y_, transmit_(ix, iy)
            end do
        write(16, *)
        end do
        close(16)

    end subroutine save_transmit


    subroutine save_transmit_y(filename, nka, Lcur)
        use all_var
        use dif_oper
        implicit none
        character(101) :: filename
        real(8) nka(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
        real(8) Lcur

        integer(4) ix, iz
        real(8) x_, z_
        real(8) transmit_(Nx, Nz)

        !subroutine find_transmittance(tshka, nka, Lmax)
        call find_transmittance_y(transmit_, nka, Lcur)

        open(16, file=filename)
        do ix = 1, Nx
            x_ = (ix-1)*hx - x0
            do iz = 1, Nz
                z_ = (iz-1)*hz - z0
                write(16,'(F0.10, 1x, F0.10, 1x, F0.10)') x_, z_, transmit_(ix, iz)
            end do
        write(16, *)
        end do
        close(16)

    end subroutine save_transmit_y


    subroutine solve()
        use all_var
        use init_mod
        use rp_mod
        use dif_oper
        implicit none
    
        integer(4) it, i, ix, iy, iz, i1, i2
        integer(4) ixv1, iyv1, izv1, ixv2, iyv2, izv2, ixv3, iyv3, izv3
        real(8) xv1, yv1, zv1, xv2, yv2, zv2, xv3, yv3, zv3
        real(8) t, t1, t2, ht_loc, x, y, z, nox, noy, noz, nop, cp, mod_ni, proc
        real(8),allocatable :: mphi_x(:), mdphi_x(:), transmit(:,:)
        real(4) time1
        real(8) minx, miny, minz, mindr
        integer(4) minix, miniy, miniz
    
        character(101) :: header, ch
    
        call init()
        print *,'hello'
    
        open(11, file="./results/nx_axisX_Y0_Z0_start.txt")
        iz = int(Nz/2)+1
        iy = int(Ny/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            write(11,'(F0.10, 1x, F0.10)') x, mvecni(1, ix, iy, iz)
        end do
        close(11)
    
        open(11, file="./results/Q_axisX_Y0_Z0_start.txt")
        iz = int(Nz/2)+1
        iy = int(Ny/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            write(11,'(F0.10, 1x, F0.10)') x, mQ(1, ix, iy, iz)
        end do
        close(11)
    
        open(11, file="./results/nx_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecni(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/nvec_XY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx, vec_step
            x = (ix-1)*hx - x0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') x, y, nox, noy, dsqrt(nox**2 + noy**2)
            end do
        end do
        close(11)
    
        open(11, file="./results/nvec_XZ_Y0_start.txt")
        iy = int(Ny/2)+1
        do ix = 1, Nx, vec_step
            x = (ix-1)*hx - x0
            do iz = 1, Nz, vec_step
                z = (iz-1)*hz - z0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') x, z, nox, noz, dsqrt(nox**2 + noz**2)
            end do
        end do
        close(11)
    
        open(11, file="./results/nvec_YZ_X0_start.txt")
        ix = int(Nx/2)+1
        do iz = 1, Nz, vec_step
            z = (iz-1)*hz - z0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') y, z, noy, noz, dsqrt(noy**2 + noz**2)
            end do
        end do
        close(11)
    
        !проекции компонентов вектора n на наклонную плоскость...
        !...работает только для кубика
        !open(11, file="./results/nvec_X(Y=-Z)_start.txt")
        !do iy = 1, Ny, vec_step
        !    y = (iy-1)*hy - y0
        !    z = -y
        !    iz = nint((z + z0)/hz + 1)
        !    do ix = 1, Nx, vec_step
        !        x = (ix-1)*hx - x0
        !        nox = mvecni(1, ix, iy, iz)
        !        noy = mvecni(2, ix, iy, iz)
        !        noz = mvecni(3, ix, iy, iz)
        !        nop = noy*(1d0/dsqrt(2d0)) - noz*(1d0/dsqrt(2d0))
        !        cp = dsqrt((y+y0)**2+(y+y0)**2) - dsqrt(2d0)/2d0
        !        write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') cp, x, nop, nox, dsqrt(nox**2 + nop**2)
        !    end do
        !end do
        !close(11)

        open(11, file="./results/nvec_XYZ_start.txt")
        iz = int(Nz/2)+1
        z = (iz-1)*hz - z0
        do ix = 1, Nx, vec_step
            x = (ix-1)*hx - x0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                x, y, z, nox, noy, noz, dsqrt(nox**2 + noy**2)
            end do
        end do
        ix = int(Nx/2)+1
        x = (ix-1)*hx - x0
        do iz = 1, Nz, vec_step
            z = (iz-1)*hz - z0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                x, y, z, nox, noy, noz, dsqrt(noz**2 + noy**2)
            end do
        end do
        iy = int(Ny/2)+1
        y = (iy-1)*hy - y0
        do iz = 1, Nz, vec_step
            z = (iz-1)*hz - z0
            do ix = 1, Nx, vec_step
                x = (ix-1)*hx - x0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                x, y, z, nox, noy, noz, dsqrt(noz**2 + nox**2)
            end do
        end do
        close(11)
    
        open(11, file="./results/ny_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecni(2, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/nz_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecni(3, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/Q_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mQ(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/Ex_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecE(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)

        open(11, file="./results/Ex_planeXZ_Y0_start.txt")
        iy = int(Ny/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iz = 1, Nz
                z = (iz-1)*hz - z0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, z, mvecE(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)        
    
        open(11, file="./results/Ey_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecE(2, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/Ez_planeXY_Z0_start.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecE(3, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        gtet = 0d0
        gphi = 0d0
        !subroutine n2trig(rnx, rny, rnz, rtet, rphi)
        allocate(mphi_x(Nx))
        allocate(mdphi_x(Nx))
        mphi_x(:) = 0d0
        mdphi_x(:) = 0d0

        iy = int(Ny/2)+1
        iz = int(Nz/2)+1
        do ix = 1, Nx
            call n2trig(mvecni(1, ix, iy, iz), mvecni(2, ix, iy, iz), mvecni(3, ix, iy, iz), gtet, gphi)
            mphi_x(ix) = gphi
        end do
    
        !grad_(3) = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
        open(20, file="./results/n2trig_start.txt")
        ix = 1
        x = (ix-1)*hx - x0
        mdphi_x(ix) = (-3d0*mphi_x(ix) + 4d0*mphi_x(ix+1) - mphi_x(ix+2)) / (2d0*hx)
        write(20,'(F0.10, 1x, F0.10, 1x, F0.10)') x, mphi_x(ix), mdphi_x(ix)
        do ix = 2, Nx-1
            x = (ix-1)*hx - x0
            !grad_(3) = (Ui(1, ix, iy, iz+1) - Ui(1, ix, iy, iz-1)) / (2d0*hz)
            mdphi_x(ix) = (mphi_x(ix+1) - mphi_x(ix-1)) / (2d0*hx)
            write(20 ,'(F0.10, 1x, F0.10, 1x, F0.10)') x, mphi_x(ix), mdphi_x(ix)
        end do
        ix = Nx
        x = (ix-1)*hx - x0
        !grad_(3) = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
        mdphi_x(ix) = (3d0*mphi_x(ix) - 4d0*mphi_x(ix-1) + mphi_x(ix-2)) / (2d0*hx)
        write(20,'(F0.10, 1x, F0.10, 1x, F0.10)') x, mphi_x(ix), mdphi_x(ix)
        close(20)
    
        xv1 = Lx/4d0
        yv1 = 0d0
        zv1 = 0d0
        ixv1 = int((xv1 + x0)/hx) + 1
        iyv1 = int((yv1 + y0)/hy) + 1
        izv1 = int((zv1 + z0)/hz) + 1
    
        xv2 = Lx/4d0
        yv2 = Ly/4d0
        zv2 = 0d0
        ixv2 = int((xv2 + x0)/hx) + 1
        iyv2 = int((yv2 + y0)/hy) + 1
        izv2 = int((zv2 + z0)/hz) + 1
    
        xv3 = Lx/4d0
        yv3 = Ly/4d0
        zv3 = Lz/4d0
        ixv3 = int((xv3 + x0)/hx) + 1
        iyv3 = int((yv3 + y0)/hy) + 1
        izv3 = int((zv3 + z0)/hz) + 1

        !stop
        write(ch,'("./results/nprj_xy_start.txt")')
        read (ch,'(a)') header
        call save_nprj_xy(header, Rsph)

        open(11, file="./results/min_max_val_nx(t).txt")
        open(12, file="./results/min_max_val_ny(t).txt")
        open(14, file="./results/min_max_val_nz(t).txt")
        open(15, file="./results/min_max_val_Q(t).txt")
        open(17, file="./results/nq_pv1.txt")
        open(18, file="./results/nq_pv2.txt")
        open(19, file="./results/nq_pv3.txt")

        call init_mass()
        allocate(transmit(Nx, Ny))
    
        time1 = secnds(0e0)
        
        ht_loc = ht

        if (scheme == 1) then
            Nt = Nt-1
            ss = 'gamd'
        elseif (scheme == 2) then
            ss = 'expl'
        endif

        if (shape == 1) then
            ss_shape = 'sph'
        elseif (shape == 2) then
            ss_shape = 'hsp'
        elseif (shape == 3) then
            ss_shape = 'cyl'
        endif

        do it = 1, Nt
            t1=ht*(it-1)
            t2=ht*(it)

            if (mod((it-1), save_step) == 0) then
                write(11, '(F0.10, 1x, F0.10, 1x, F0.10)') t1, minval(mvecni(1,:,:,:)), maxval(mvecni(1,:,:,:))   !nx
                write(12, '(F0.10, 1x, F0.10, 1x, F0.10)') t1, minval(mvecni(2,:,:,:)), maxval(mvecni(2,:,:,:))   !ny
                write(14, '(F0.10, 1x, F0.10, 1x, F0.10)') t1, minval(mvecni(3,:,:,:)), maxval(mvecni(3,:,:,:))   !nz
                write(15, '(F0.10, 1x, F0.10, 1x, F0.10)') t1, minval(mQ(1,:,:,:)), maxval(mQ(1,:,:,:))           !Q

                write(17, '(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') t1, mvecni(1, ixv1, iyv1, izv1), &
                                                                            mvecni(2, ixv1, iyv1, izv1), &
                                                                            mvecni(3, ixv1, iyv1, izv1), &
                                                                                mQ(1, ixv1, iyv1, izv1)
        
                write(18, '(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') t1, mvecni(1, ixv2, iyv2, izv2), &
                                                                            mvecni(2, ixv2, iyv2, izv2), &
                                                                            mvecni(3, ixv2, iyv2, izv2), &
                                                                            mQ(1, ixv2, iyv2, izv2)
        
                write(19, '(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') t1, mvecni(1, ixv3, iyv3, izv3), &
                                                                            mvecni(2, ixv3, iyv3, izv3), &
                                                                            mvecni(3, ixv3, iyv3, izv3), &
                                                                            mQ(1, ixv3, iyv3, izv3)
            end if

            if (mod((it-1), save_step_planes) == 0) then
                last_saved_plane = it-1

                write(ch,'("./results/nprj_xy_out/nprj_xy_",I0.8,".txt")') it-1
                read (ch,'(a)') header
                call save_nprj_xy(header, Rsph)

                write(ch,'("./results/nprj_xz_out/nprj_xz_",I0.8,".txt")') it-1
                read (ch,'(a)') header
                call save_nprj_xz(header, Rsph)

                write(ch,'("./results/nprj_xy_in/nprj_xy_",I0.8,".txt")') it-1
                read (ch,'(a)') header
                call save_nprj_xy(header, Rsph_in)

                write(ch,'("./results/nprj_xz_in/nprj_xz_",I0.8,".txt")') it-1
                read (ch,'(a)') header
                call save_nprj_xz(header, Rsph_in)


                !write(ch,'("./results/npic/nvec_XY_Z0_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_nvec_xy_z0(header, mvecni)
        
                !write(ch,'("./results/qpic/Q_planeXY_Z0_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_Qplane_xy_z0(header, mQ)

                !write(ch,'("./results/transmit_000/transmit_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_transmit(header, mvecni, 0d0)

                !write(ch,'("./results/transmit_050/transmit_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_transmit(header, mvecni, 0.50d0)

                !write(ch,'("./results/transmit_045/transmit_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_transmit(header, mvecni, 0.46d0)

                !write(ch,'("./results/transmit_040/transmit_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_transmit(header, mvecni, 0.40d0)

                !write(ch,'("./results/transmit_025/transmit_",I0.8,".txt")') it-1
                !read (ch,'(a)') header
                !call save_transmit(header, mvecni, 0.25d0)

                ! write(ch,'("./results/transmit_031/transmit_",I0.8,".txt")') it-1
                ! read (ch,'(a)') header
                ! call save_transmit(header, mvecni, 0.31d0)

                ! write(ch,'("./results/transmit_y000/transmit_",I0.8,".txt")') it-1
                ! read (ch,'(a)') header
                ! call save_transmit_y(header, mvecni, 0d0)
                
                if (save_raw_mnmq >= 1) then
                    write(ch,'("./results/mnmq/mn_",I0.8,".bin")') it-1
                    read (ch,'(a)') header
                    open(31,file=header, form='unformatted')
                    write(31) 3, Nx,Ny,Nz
                    write(31) Lx,Ly,Lz
                    write(31) mvecni
                    close(31)

                    write(ch,'("./results/mnmq/mq_",I0.8,".bin")') it-1
                    read (ch,'(a)') header
                    open(31,file=header, form='unformatted')
                    write(31) 1, Nx,Ny,Nz
                    write(31) Lx,Ly,Lz
                    write(31) mQ
                    close(31)
                end if

            end if                                                

            call rp(1,1d0,m_tmp1,m_tmp2)

            if (mod(it, save_step) == 0) then
                time1 = secnds(time1)
                proc = (real(it)/real(Nt))*100
                write(*, '(a4"/"a4": "F0.2"% "I0.0"/"I0.0" <time "es8.2">  <nx "F0.2">  <ny "F0.2">  <nz "F0.2">  <Q: "F0.2" to "F0.2">  <tps "F0.2">")') ss, ss_shape, &
                                proc, it, Nt, t1, &
                                maxval(dabs(mvecni(1,:,:,:))), maxval(dabs(mvecni(2,:,:,:))), &
                                maxval(dabs(mvecni(3,:,:,:))), minval(dabs(mQ(1,:,:,:))), maxval(dabs(mQ(1,:,:,:))), time1
                time1 = secnds(0e0)
            end if
        enddo
      
        write(11, '(F0.10, 1x, F0.10, 1x, F0.10)') t2, minval(mvecni(1,:,:,:)), maxval(mvecni(1,:,:,:))                         !nx
        write(12, '(F0.10, 1x, F0.10, 1x, F0.10)') t2, minval(mvecni(2,:,:,:)), maxval(mvecni(2,:,:,:))   !ny
        write(14, '(F0.10, 1x, F0.10, 1x, F0.10)') t2, minval(mvecni(3,:,:,:)), maxval(mvecni(3,:,:,:))    !nz
        write(15, '(F0.10, 1x, F0.10, 1x, F0.10)') t2, minval(mQ(1,:,:,:)), maxval(mQ(1,:,:,:))    !Q
    
        write(17, '(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') t2, mvecni(1, ixv1, iyv1, izv1), &
                                                                    mvecni(2, ixv1, iyv1, izv1), &
                                                                    mvecni(3, ixv1, iyv1, izv1), &
                                                                    mQ(1, ixv1, iyv1, izv1)
    
        write(18, '(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') t2, mvecni(1, ixv2, iyv2, izv2), &
                                                                    mvecni(2, ixv2, iyv2, izv2), &
                                                                    mvecni(3, ixv2, iyv2, izv2), &
                                                                    mQ(1, ixv2, iyv2, izv2)
    
        write(19, '(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') t2, mvecni(1, ixv3, iyv3, izv3), &
                                                                    mvecni(2, ixv3, iyv3, izv3), &
                                                                    mvecni(3, ixv3, iyv3, izv3), &
                                                                    mQ(1, ixv3, iyv3, izv3)
    
        close(11)
        close(12)
        close(14)
        close(15)
        close(17)
        close(18)
        close(19)

        open(31,file='./nka.bin',form='unformatted')
        write(31) 3, Nx,Ny,Nz
        write(31) Lx,Ly,Lz
        write(31) mvecni
        close(31)

        open(31,file='./qushka.bin',form='unformatted')
        write(31) 1, Nx,Ny,Nz
        write(31) Lx,Ly,Lz
        write(31) mQ
        close(31)

        print *,"nka_qushka_saved_ok"

        gtet = 0d0
        gphi = 0d0
        !subroutine n2trig(rnx, rny, rnz, rtet, rphi)
        mphi_x(:) = 0d0
        mdphi_x(:) = 0d0
        iy = int(Ny/2)+1
        iz = int(Nz/2)+1
    
        do ix = 1, Nx
            call n2trig(mvecni(1, ix, iy, iz), mvecni(2, ix, iy, iz), mvecni(3, ix, iy, iz), gtet, gphi)
            mphi_x(ix) = gphi
        end do
    
        !grad_(3) = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
        open(20, file="./results/n2trig_end.txt")
        ix = 1
        x = (ix-1)*hx - x0
        mdphi_x(ix) = (-3d0*mphi_x(ix) + 4d0*mphi_x(ix+1) - mphi_x(ix+2)) / (2d0*hx)
        write(20,'(F0.10, 1x, F0.10, 1x, F0.10)') x, mphi_x(ix), mdphi_x(ix)
        do ix = 2, Nx-1
            x = (ix-1)*hx - x0
            !grad_(3) = (Ui(1, ix, iy, iz+1) - Ui(1, ix, iy, iz-1)) / (2d0*hz)
            mdphi_x(ix) = (mphi_x(ix+1) - mphi_x(ix-1)) / (2d0*hx)
            write(20 ,'(F0.10, 1x, F0.10, 1x, F0.10)') x, mphi_x(ix), mdphi_x(ix)
        end do
        ix = Nx
        x = (ix-1)*hx - x0
        !grad_(3) = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
        mdphi_x(ix) = (3d0*mphi_x(ix) - 4d0*mphi_x(ix-1) + mphi_x(ix-2)) / (2d0*hx)
        write(20,'(F0.10, 1x, F0.10, 1x, F0.10)') x, mphi_x(ix), mdphi_x(ix)
        close(20)
        
        !print *,"nka_qushka_start_saving..."

        !subroutine find_transmittance(tshka, nka)
        !call find_transmittance_y(transmit, mvecni, 0d0)

        !subroutine save_transmit_y(filename, nka, Lcur)
        !write(ch,'("./results/transmit_end000.txt")')
        !read (ch,'(a)') header
        !call save_transmit_y(header, mvecni, 0d0)
        !call save_transmit(header, mvecni, 0d0)

        !write(ch,'("./results/transmit_end025.txt")')
        !read (ch,'(a)') header
        !call save_transmit_y(header, mvecni, 0.25d0)
        !call save_transmit(header, mvecni, 0.25d0)

        !write(ch,'("./results/transmit_end050.txt")')
        !read (ch,'(a)') header
        !call save_transmit_y(header, mvecni, 0.50d0)
        !call save_transmit(header, mvecni, 0.50d0)

        !print *,"save_transmint_y_done"
        !print *,"transmit_saved_ok"

        open(11, file="./results/nx_axisX_Y0_Z0_end.txt")
        iz = int(Nz/2)+1
        iy = int(Ny/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            write(11,'(F0.10, 1x, F0.10)') x, mvecni(1, ix, iy, iz)
        end do
        close(11)
    
        open(11, file="./results/Q_axisX_Y0_Z0_end.txt")
        iz = int(Nz/2)+1
        iy = int(Ny/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            write(11,'(F0.10, 1x, F0.10)') x, mQ(1, ix, iy, iz)
        end do
        close(11)

        open(11, file="./results/Q_axisZ_X0_Y0_end.txt")
        ix = int(Nx/2)+1
        iy = int(Ny/2)+1
        do iz = 1, Nz
            z = (iz-1)*hz - z0
            write(11,'(F0.10, 1x, F0.10)') z, mQ(1, ix, iy, iz)
        end do
        close(11)        
    
        open(11, file="./results/nx_planeXY_Z0_end.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecni(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
  
       open(11, file="./results/nvec_XY_Z0_end.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx, vec_step
            x = (ix-1)*hx - x0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                !noz = a_U(get_i(3, ix, iy, iz))
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') x, y, nox, noy, dsqrt(nox**2 + noy**2)
            end do
        end do
        close(11)
    
        open(11, file="./results/nvec_XZ_Y0_end.txt")
        iy = int(Ny/2)+1
        do ix = 1, Nx, vec_step
            x = (ix-1)*hx - x0
            do iz = 1, Nz, vec_step
                z = (iz-1)*hz - z0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') x, z, nox, noz, dsqrt(nox**2 + noz**2)
            end do
        end do
        close(11)
    
        open(11, file="./results/nvec_YZ_X0_end.txt")
        ix = int(Nx/2)+1
        do iz = 1, Nz, vec_step
            z = (iz-1)*hz - z0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') y, z, noy, noz, dsqrt(noy**2 + noz**2)
            end do
        end do
        close(11)

        !проекции компонентов вектора с верхней полусферы на плоскость xy
        write(ch,'("./results/nprj_xy_out.txt")')
        read (ch,'(a)') header
        call save_nprj_xy(header, Rsph)
        ! open(11, file=header)
        ! do iy = 1, Ny, vec_step
        !     y = (iy-1)*hy - y0
        !     do ix = 1, Nx, vec_step
        !         x = (ix-1)*hx - x0
        !         mindr = 1d3
        !         miniz = -1
        !         do iz = int(Nz/2)+1, Nz, vec_step
        !             z = (iz-1)*hz - z0            
        !             Rcur = dsqrt(x**2 + y**2 + z**2)
        !             dR = Rsph - Rcur
        !             if ((dR>0) .AND. (dR < mindr)) then
        !                 mindr = dR
        !                 minz = z
        !                 miniz = iz
        !             end if    
        !         end do
        !         if (miniz >= 0) then
        !             nox = mvecni(1, ix, iy, miniz)
        !             noy = mvecni(2, ix, iy, miniz)
        !             noz = mvecni(3, ix, iy, miniz)
        !             mod_ni = dsqrt(nox**2 + noy**2)
        !             write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
        !                         x, y, minz, nox, noy, noz, mod_ni
        !         else
        !             nox = 0d0
        !             noy = 0d0
        !             noz = 0d0
        !             minz = -1d0
        !             mod_ni = 0d0
        !             write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
        !                         x, y, minz, nox, noy, noz, mod_ni
        !         end if
        !     end do
        ! end do
        ! close(11)

        !проекции компонентов вектора с верхней полусферы на плоскость xz
        write(ch,'("./results/nprj_xz_out.txt")')
        read (ch,'(a)') header
        call save_nprj_xz(header, Rsph)
        ! open(11, file=header)
        ! do iz = 1, Nz, vec_step
        !     z = (iz-1)*hz - z0
        !     do ix = 1, Nx, vec_step
        !         x = (ix-1)*hx - x0
        !         mindr = 1d3
        !         miniy = -1
        !         do iy = int(Ny/2)+1, Ny, vec_step
        !             y = (iy-1)*hy - y0            
        !             Rcur = dsqrt(x**2 + y**2 + z**2)
        !             dR = Rsph - Rcur
        !             if ((dR>0) .AND. (dR < mindr)) then
        !                 mindr = dR
        !                 miny = y
        !                 miniy = iy
        !             end if    
        !         end do
        !         if (miniy >= 0) then
        !             nox = mvecni(1, ix, miniy, iz)
        !             noy = mvecni(2, ix, miniy, iz)
        !             noz = mvecni(3, ix, miniy, iz)
        !             mod_ni = dsqrt(nox**2 + noz**2)
        !             write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
        !                         x, miny, z, nox, noy, noz, mod_ni
        !         else
        !             nox = 0d0
        !             noy = 0d0
        !             noz = 0d0
        !             miny = -1d0
        !             mod_ni = 0d0
        !             write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
        !                         x, miny, z, nox, noy, noz, mod_ni
        !         end if
        !     end do
        ! end do
        ! close(11)

        !проекции компонентов вектора n на наклонную плоскость...
        !...работает только для кубика
        !open(11, file="./results/nvec_X(Y=-Z)_end.txt")
        !do iy = 1, Ny, vec_step
        !    y = (iy-1)*hy - y0
        !    z = -y
        !    iz = nint((z + z0)/hz + 1)
        !    do ix = 1, Nx, vec_step
        !        x = (ix-1)*hx - x0
        !        nox = mvecni(1, ix, iy, iz)
        !        noy = mvecni(2, ix, iy, iz)
        !        noz = mvecni(3, ix, iy, iz)
        !        nop = noy*(1d0/dsqrt(2d0)) - noz*(1d0/dsqrt(2d0))
        !        cp = dsqrt((y+y0)**2+(y+y0)**2) - dsqrt(2d0)/2d0
        !        write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') cp, x, nop, nox, dsqrt(nox**2 + nop**2)
        !    end do
        !end do
        !close(11)
 
        ! open(11, file="./results/nvec_sphere_end.txt")
        ! iz = int(Nz/2)+1
        ! z = (iz-1)*hz - z0
        ! do ix = 1, Nx, vec_step
        !     x = (ix-1)*hx - x0
        !     do iy = 1, Ny, vec_step
        !         y = (iy-1)*hy - y0
        !         do iz = 1, Nz, vec_step
        !             z = (iz-1)*hz - z0            
        !             nox = mvecni(1, ix, iy, iz)
        !             noy = mvecni(2, ix, iy, iz)
        !             noz = mvecni(3, ix, iy, iz)
        !             mod_ni = dsqrt(nox**2 + noy**2 + noz**2)
        !             Rcur = dsqrt(x**2 + y**2 + z**2)
        !             dR = Rsph - Rcur
        !             if ((dR <= hx) .AND. (dR > 0d0)) then
        !                 write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
        !                     x, y, z, nox, noy, noz, mod_ni
        !             end if
        !         end do
        !     end do
        ! end do
        ! close(11)

        open(11, file="./results/nvec_XYZ_end.txt")
        iz = int(Nz/2)+1
        z = (iz-1)*hz - z0
        do ix = 1, Nx, vec_step
            x = (ix-1)*hx - x0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                x, y, z, nox, noy, noz, dsqrt(nox**2 + noy**2)
            end do
        end do
        ix = int(Nx/2)+1
        x = (ix-1)*hx - x0
        do iz = 1, Nz, vec_step
            z = (iz-1)*hz - z0
            do iy = 1, Ny, vec_step
                y = (iy-1)*hy - y0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                x, y, z, nox, noy, noz, dsqrt(noz**2 + noy**2)
            end do
        end do
        iy = int(Ny/2)+1
        y = (iy-1)*hy - y0
        do iz = 1, Nz, vec_step
            z = (iz-1)*hz - z0
            do ix = 1, Nx, vec_step
                x = (ix-1)*hx - x0
                nox = mvecni(1, ix, iy, iz)
                noy = mvecni(2, ix, iy, iz)
                noz = mvecni(3, ix, iy, iz)
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10, 1x, F0.10)') &
                x, y, z, nox, noy, noz, dsqrt(noz**2 + nox**2)
            end do
        end do
        close(11)
    
        open(11, file="./results/ny_planeXY_Z0_end.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecni(2, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/nz_planeXY_Z0_end.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mvecni(3, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        open(11, file="./results/Q_planeXY_Z0_end.txt")
        iz = int(Nz/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iy = 1, Ny
                y = (iy-1)*hy - y0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, y, mQ(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)

        open(11, file="./results/Q_planeXZ_Y0_end.txt")
        iy = int(Ny/2)+1
        do ix = 1, Nx
            x = (ix-1)*hx - x0
            do iz = 1, Nz
                z = (iz-1)*hz - z0
                write(11,'(F0.10, 1x, F0.10, 1x, F0.10)') x, z, mQ(1, ix, iy, iz)
            end do
            write(11, *)
        end do
        close(11)
    
        print *,"end_ode"
    
    end subroutine solve
    
    end module for_ode_mod
    