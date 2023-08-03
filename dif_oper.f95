module dif_oper

    use all_var
    implicit none

contains


!subroutine save_to_file(filename)
!
!end subroutine save_to_file

function get_i(eq_n, jx, jy, jz) result (get_i_)
    use all_var
    implicit none
    integer(4) eq_n, jx, jy, jz
    integer(4) get_i_

    get_i_ = NxNyNz*(eq_n-1) + NxNy*(jz-1) + Nx*(jy-1) + jx

end function get_i

subroutine find_transmittance_y(tshka, nka, Lmax)
    use all_var
    implicit none
    real(8) tshka(Nx, Nz)
    real(8) nka(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) Lmax

    real(8) psi_(0:Nx+1, 0:Nz+1)
    real(8) e_par, e_perp, n_par, n_perp, e_a, hii, ni, ni2, n1, nn, Ldl, teti, phii
    integer(4) ix, iy, iz, Nmax

    hii = pi / 2
    n_par = 1.6
    e_par = n_par**2
    n_perp = 1.52
    e_perp = n_perp**2
    e_a = e_par - e_perp

    psi_ = 0d0
    Ldl = Lbig / 0.5d-6 !Ly / 0.5 !20d0 !Lz/liambda

    Nmax = nint((Lmax+y0)/hy) + 1 !верхний диапазон интегрирования
    print *,"Nmax:", Nmax

    do ix = 1, Nx
        do iz = 1, Nz
            do iy = 1, Nmax
                !subroutine n2trig(rnx, rny, rnz, rtet, rphi)
                call n2trig_y(nka(1, ix, iy, iz), nka(2, ix, iy, iz), nka(3, ix, iy, iz), teti, phii)
                ni2 = (e_par * e_perp) / (e_perp + e_a *(dcos(teti))**2)
                ni = dsqrt(ni2)
                !интеграл по методу трапеций
                if (iy == 1) then
                    n1 = ni
                elseif (iy == Nmax) then
                    nn = ni
                else
                    psi_(ix, iz) = psi_(ix, iz) + ni*hy
                endif
            end do
            psi_(ix, iz) = psi_(ix, iz) + 0.5*(n1+nn)*hy
        end do
    end do    

    psi_ = psi_*2*pi*Ldl - 2*pi*Ldl*n_perp
    print *,"psi_ calculated"

    do ix = 1, Nx
        do iz = 1, Nz
            call n2trig_y(nka(1, ix, Nmax, iz), nka(2, ix, Nmax, iz), nka(3, ix, Nmax, iz), teti, phii)
            tshka(ix, iz) = (dsin(2*phii))**2 * (dsin(psi_(ix, iz)))**2
        end do
    end do
    print *,"tshka calculated"
end subroutine find_transmittance_y


subroutine find_transmittance(tshka, nka, Lmax)
    use all_var
    implicit none
    real(8) tshka(Nx, Ny)
    real(8) nka(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) Lmax

    real(8) psi_(0:NX+1, 0:Ny+1)
    real(8) e_par, e_perp, n_par, n_perp, e_a, hii, ni, ni2, n1, nn, Lzdl, teti, phii
    integer(4) ix, iy, iz, Nzmax

    hii = pi / 2
    n_par = 1.6
    e_par = n_par**2
    n_perp = 1.52
    e_perp = n_perp**2
    e_a = e_par - e_perp

    psi_ = 0d0
    Lzdl = Lz / 0.5 !20d0 !Lz/liambda

    Nzmax = nint((Lmax+z0)/hz) + 1 !Nz/2 + 1 !верхний диапазон интегрирования

    do ix = 1, Nx
        do iy = 1, Ny
            do iz = 1, Nzmax
                !subroutine n2trig(rnx, rny, rnz, rtet, rphi)
                call n2trig(nka(1, ix, iy, iz), nka(2, ix, iy, iz), nka(3, ix, iy, iz), teti, phii)
                ni2 = (e_par * e_perp) / (e_perp + e_a * dcos(teti)**2)
                ni = dsqrt(ni2)
                !интеграл по методу трапеций
                if (iz == 1) then
                    n1 = ni
                elseif (iz == Nzmax) then
                    nn = ni
                else
                    psi_(ix, iy) = psi_(ix, iy) + ni*hz
                endif
            end do
            psi_(ix, iy) = psi_(ix, iy) + 0.5*(n1+nn)*hz
        end do
    end do    

    psi_ = psi_*2*pi*Lzdl - 2*pi*Lzdl*n_perp

    do ix = 1, Nx
        do iy = 1, Ny
            call n2trig(nka(1, ix, iy, Nzmax), nka(2, ix, iy, Nzmax), nka(3, ix, iy, Nzmax), teti, phii)
            tshka(ix, iy) = dcos(hii)**2 - dsin(2*phii) * dsin(2*(phii-hii)) * dsin(psi_(ix, iy)/2d0)**2
        end do
    end do

end subroutine find_transmittance


subroutine linear_norm(flatU)
    use all_var
    implicit none
    real(8) flatU(N_full)
    
    integer(4) ix, iy, iz
    real(8) mod_ni_, mod_q_

    do iz = 1, Nz
        do iy = 1, Ny
            do ix = 1, Nx
                mod_ni_ = dsqrt(flatU(get_i(1, ix, iy, iz))**2 + flatU(get_i(2, ix, iy, iz))**2 + flatU(get_i(3, ix, iy, iz))**2)
                if (mod_ni_ > 1d0) then
                    flatU(get_i(1, ix, iy, iz)) = flatU(get_i(1, ix, iy, iz)) / mod_ni_
                    flatU(get_i(2, ix, iy, iz)) = flatU(get_i(2, ix, iy, iz)) / mod_ni_
                    flatU(get_i(3, ix, iy, iz)) = flatU(get_i(3, ix, iy, iz)) / mod_ni_
                endif

            !    mod_q_ = dabs(qushka(1, ix, iy, iz))
            !    if (mod_q_ > 1d0) then
            !        qushka(1, ix, iy, iz) = qushka(1, ix, iy, iz) / mod_q_
            !    endif
            end do
        end do
    end do            
end subroutine linear_norm

subroutine normalize(nka_, qushka_)
    use all_var
    implicit none
    real(8) nka_(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) qushka_(1, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    
    integer(4) ix, iy, iz
    real(8) mod_ni_, mod_q_

    do iz = 1, Nz
        do iy = 1, Ny
            do ix = 1, Nx
                mod_ni_ = dsqrt(nka_(1, ix, iy, iz)**2 + nka_(2, ix, iy, iz)**2 + nka_(3, ix, iy, iz)**2)
                if (mod_ni_ > 1d0) then
                    nka_(1, ix, iy, iz) = nka_(1, ix, iy, iz) / mod_ni_
                    nka_(2, ix, iy, iz) = nka_(2, ix, iy, iz) / mod_ni_
                    nka_(3, ix, iy, iz) = nka_(3, ix, iy, iz) / mod_ni_
                endif

            !    mod_q_ = dabs(qushka(1, ix, iy, iz))
            !    if (mod_q_ > 1d0) then
            !        qushka(1, ix, iy, iz) = qushka(1, ix, iy, iz) / mod_q_
            !    endif
            end do
        end do
    end do            
end subroutine normalize


subroutine line2tensor(uline, nka, qushka, do_norm)
    use all_var
    implicit none
    real(8) uline(N_full)
    real(8) nka(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) qushka(1, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) do_norm

    integer(4) ix, iy, iz

    do iz = 1, Nz
        do iy = 1, Ny
            do ix = 1, Nx
                nka(1, ix, iy, iz) = uline(get_i(1, ix, iy, iz))
                nka(2, ix, iy, iz) = uline(get_i(2, ix, iy, iz))
                nka(3, ix, iy, iz) = uline(get_i(3, ix, iy, iz))
                qushka(1, ix, iy, iz) = uline(get_i(4, ix, iy, iz))
            end do
        end do
    end do            

    if (do_norm >= 1) then
        call normalize(nka, qushka)
    endif

end subroutine line2tensor


subroutine bound_all_abc(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy, iz
    real(8) xx, x1, x2, yy, y1, y2, zz, z1, z2
    real(8) u2, u3, un1, un2


    x1 = (edge-1)*hx - x0
    x2 = (Nx+1-edge-1)*hx - x0
    do iy = edge, Ny+1-edge
        yy = (iy-1)*hy - y0
        do iz = edge, Nz+1-edge
            zz = (iz-1)*hz - z0
            u2 = Ui(1, edge +1, iy, iz)
            u3 = Ui(1, edge +2, iy, iz)
            !левая s := -(4*BB*h*u2 - BB*h*u3 - 2*F*h^2 - 4*CC*u2 + 2*CC*u3) / (2*AA*h^2 - 3*BB*h + 2*CC)
            Ui(1, edge, iy, iz) = -(4*bb*hx*u2 - bb*hx*u3 - 2*estim(x1, yy, zz)*hx2 - 4*cc*u2 + 2*cc*u3)
            Ui(1, edge, iy, iz) = Ui(1, edge, iy, iz) / (2*aa*hx2 - 3*bb*hx + 2*cc)
            !правая s := -(BB*h*u1 - 4*BB*h*u2 - 2*F*h^2 + 2*CC*u1 - 4*CC*u2) / (2*AA*h^2 + 3*BB*h + 2*CC)
            un1 = Ui(1, Nx+1-edge -2, iy, iz)
            un2 = Ui(1, Nx+1-edge -1, iy, iz)
            Ui(1, Nx+1-edge, iy, iz) = -(bb*hx*un1 - 4*bb*hx*un2 - 2*estim(x2, yy, zz)*hx2 + 2*cc*un1 - 4*cc*un2)
            Ui(1, Nx+1-edge, iy, iz) = Ui(1, Nx+1-edge, iy, iz) / (2*aa*hx2 + 3*bb*hx + 2*cc)
        end do
    end do

    y1 = (edge-1)*hy - y0
    y2 = (Ny+1-edge-1)*hy - y0
    do ix = edge, Nx+1-edge
        xx = (ix-1)*hx - x0
        do iz = edge, Nz+1-edge
            zz = (iz-1)*hz - z0
            u2 = Ui(1, ix, edge +1, iz)
            u3 = Ui(1, ix, edge +2, iz)
            Ui(1, ix, edge, iz) = -(4*bb*hy*u2 - bb*hy*u3 - 2*estim(xx, y1, zz)*hy2 - 4*cc*u2 + 2*cc*u3)
            Ui(1, ix, edge, iz) = Ui(1, ix, edge, iz) / (2*aa*hy2 - 3*bb*hy + 2*cc)
            un1 = Ui(1, ix, Ny+1-edge -2, iz)
            un2 = Ui(1, ix, Ny+1-edge -1, iz)
            Ui(1, ix, Ny+1-edge, iz) = -(bb*hy*un1 - 4*bb*hy*un2 - 2*estim(xx, y2, zz)*hy2 + 2*cc*un1 - 4*cc*un2)
            Ui(1, ix, Ny+1-edge, iz) = Ui(1, ix, Ny+1-edge, iz) / (2*aa*hy2 + 3*bb*hy + 2*cc)
        end do
    end do

    z1 = (edge-1)*hz - z0
    z2 = (Nz+1-edge-1)*hz - z0
    do ix = edge, Nx+1-edge
        xx = (ix-1)*hx - x0
        do iy = edge, Ny+1-edge
            yy = (iy-1)*hy - y0
            u2 = Ui(1, ix, iy, edge +1)
            u3 = Ui(1, ix, iy, edge +2)
            Ui(1, ix, iy, edge) = -(4*bb*hz*u2 - bb*hz*u3 - 2*estim(xx, yy, z1)*hz2 - 4*cc*u2 + 2*cc*u3)
            Ui(1, ix, iy, edge) = Ui(1, ix, iy, edge) / (2*aa*hz2 - 3*bb*hz + 2*cc)
            un1 = Ui(1, ix, iy, Nz+1-edge -2)
            un2 = Ui(1, ix, iy, Nz+1-edge -1)
            Ui(1, ix, iy, Nz+1-edge) = -(bb*hz*un1 - 4*bb*hz*un2 - 2*estim(xx, yy, z2)*hz2 + 2*cc*un1 - 4*cc*un2)
            Ui(1, ix, iy, Nz+1-edge) = Ui(1, ix, iy, Nz+1-edge) / (2*aa*hz2 + 3*bb*hz + 2*cc)
        end do
    end do

end subroutine bound_all_abc


subroutine bound_abc_zmax(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy
    real(8) xx, yy, zmax
    real(8) un1, un2

    zmax = (Nz+1-edge-1)*hz - z0
    do ix = edge, Nx+1-edge
        xx = (ix-1)*hx - x0
        do iy = edge, Ny+1-edge
            yy = (iy-1)*hy - y0
            un1 = Ui(1, ix, iy, Nz+1-edge -2)
            un2 = Ui(1, ix, iy, Nz+1-edge -1)
            Ui(1, ix, iy, Nz+1-edge) = -(bb*hz*un1 - 4*bb*hz*un2 - 2*estim(xx, yy, zmax)*hz2 + 2*cc*un1 - 4*cc*un2)
            Ui(1, ix, iy, Nz+1-edge) = Ui(1, ix, iy, Nz+1-edge) / (2*aa*hz2 + 3*bb*hz + 2*cc)
        end do
    end do
end subroutine bound_abc_zmax


subroutine bound_abc_zmin(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy
    real(8) xx, yy, zmin
    real(8) u2, u3

    zmin = (edge-1)*hz - z0
    do ix = edge, Nx+1-edge
        xx = (ix-1)*hx - x0
        do iy = edge, Ny+1-edge
            yy = (iy-1)*hy - y0
            u2 = Ui(1, ix, iy, edge +1)
            u3 = Ui(1, ix, iy, edge +2)
            Ui(1, ix, iy, edge) = -(4*bb*hz*u2 - bb*hz*u3 - 2*estim(xx, yy, zmin)*hz2 - 4*cc*u2 + 2*cc*u3)
            Ui(1, ix, iy, edge) = Ui(1, ix, iy, edge) / (2*aa*hz2 - 3*bb*hz + 2*cc)
        end do
    end do
end subroutine bound_abc_zmin


subroutine bound_abc_xmax(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy, iz
    real(8) xx, yy, zz, xmax
    real(8) un1, un2

    xmax = (Nx+1-edge-1)*hx - x0
    do iz = edge, Nz+1-edge
        zz = (iz-1)*hz - z0
        do iy = edge, Ny+1-edge
            yy = (iy-1)*hy - y0
            un1 = Ui(1, Nx+1-edge -2, iy, iz)
            un2 = Ui(1, Nx+1-edge -1, iy, iz)
            Ui(1, Nx+1-edge, iy, iz) = -(bb*hx*un1 - 4*bb*hx*un2 - 2*estim(xmax, yy, zz)*hx2 + 2*cc*un1 - 4*cc*un2)
            Ui(1, Nx+1-edge, iy, iz) = Ui(1, Nx+1-edge, iy, iz) / (2*aa*hx2 + 3*bb*hx + 2*cc)
        end do
    end do
end subroutine bound_abc_xmax


subroutine bound_abc_xmin(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy, iz
    real(8) xx, yy, zz, xmin
    real(8) u2, u3

    xmin = (edge-1)*hx - x0
    do iz = edge, Nz+1-edge
        zz = (iz-1)*hz - z0
        do iy = edge, Ny+1-edge
            yy = (iy-1)*hy - y0
            u2 = Ui(1, edge +1, iy, iz)
            u3 = Ui(1, edge +2, iy, iz)
            Ui(1, edge, iy, iz) = -(4*bb*hx*u2 - bb*hx*u3 - 2*estim(xmin, yy, zz)*hx2 - 4*cc*u2 + 2*cc*u3)
            Ui(1, edge, iy, iz) = Ui(1, edge, iy, iz) / (2*aa*hx2 - 3*bb*hx + 2*cc)
        end do
    end do
end subroutine bound_abc_xmin

subroutine bound_abc_ymax(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy, iz
    real(8) xx, yy, zz, ymax
    real(8) un1, un2

    ymax = (Ny+1-edge-1)*hy - y0
    do iz = edge, Nz+1-edge
        zz = (iz-1)*hz - z0
        do ix = edge, Nx+1-edge
            xx = (ix-1)*hx - x0
            un1 = Ui(1, ix, Ny+1-edge -2, iz)
            un2 = Ui(1, ix, Ny+1-edge -1, iz)
            Ui(1, ix, Ny+1-edge, iz) = -(bb*hy*un1 - 4*bb*hy*un2 - 2*estim(xx, ymax, zz)*hy2 + 2*cc*un1 - 4*cc*un2)
            Ui(1, ix, Ny+1-edge, iz) = Ui(1, ix, Ny+1-edge, iz) / (2*aa*hy2 + 3*bb*hy + 2*cc)
        end do
    end do
end subroutine bound_abc_ymax


subroutine bound_abc_ymin(Ui, aa, bb, cc, edge, estim)
    use all_var
    implicit none
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1), aa, bb, cc
    integer(4) edge
    real(8) estim
    external estim

    integer(4) ix, iy, iz
    real(8) xx, yy, zz, ymin
    real(8) u2, u3

    ymin = (edge-1)*hy - y0
    do iz = edge, Nz+1-edge
        zz = (iz-1)*hz - z0
        do ix = edge, Nx+1-edge
            xx = (ix-1)*hx - x0
            u2 = Ui(1, ix, edge +1, iz)
            u3 = Ui(1, ix, edge +2, iz)
            Ui(1, ix, edge, iz) = -(4*bb*hy*u2 - bb*hy*u3 - 2*estim(xx, ymin, zz)*hy2 - 4*cc*u2 + 2*cc*u3)
            Ui(1, ix, edge, iz) = Ui(1, ix, edge, iz) / (2*aa*hy2 - 3*bb*hy + 2*cc)
        end do
    end do
end subroutine bound_abc_ymin

function zero_all(rx, ry, rz) result (zero_all_)
    use all_var
    real(8) rx, ry, rz
    real(8) zero_all_

    zero_all_ = 0d0

end function zero_all


function one_all(rx, ry, rz) result (one_all_)
    use all_var
    real(8) rx, ry, rz
    real(8) one_all_

    one_all_ = 1d0

end function one_all


function qroot_all(rx, ry, rz) result (qroot_all_)
    use all_var
    real(8) rx, ry, rz
    real(8) qroot_all_

    qroot_all_ = Qroot

end function qroot_all



subroutine n2trig_y(rnx, rny, rnz, rtet, rphi)
    !n2trig_y углы отсчитываются от y
    real(8) rnx, rny, rnz, rtet, rphi

    rtet = dacos(rny)
    rphi = datan2(rnx, rnz)

end subroutine n2trig_y


subroutine n2trig(rnx, rny, rnz, rtet, rphi)
    !n2trig_z углы отсчитываются от z
    real(8) rnx, rny, rnz, rtet, rphi

    rtet = dacos(rnz)
    rphi = datan2(rny, rnx)

end subroutine n2trig


subroutine trig2n(rtet, rphi, rnx, rny, rnz)
    real(8) rtet, rphi, rnx, rny, rnz

    rnx = dsin(rtet) * dcos(rphi) !nx
    rny = dsin(rtet) * dsin(rphi) !ny
    rnz = dcos(rtet) !nz
end subroutine trig2n



! function init_Ex(rx, ry, rz) result (init_Ex_)
!     use all_var
!     real(8) rx, ry, rz
!     real(8) init_Ex_

!     real(8) r_, dr_, r0_

!     r0_ = 0.8*(Lx/2)    !Для кубика все r0_ равны
!     r_ = dsqrt(rx**2 + ry**2 + rz**2)
!     dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)

!     init_Ex_ = (0.5*erf((r_-r0_)/wE)+0.5) * (rx/(r_+1d-5)) !Ex

! end function init_Ex


! function init_Ey(rx, ry, rz) result (init_Ey_)
!     use all_var
!     real(8) rx, ry, rz
!     real(8) init_Ey_

!     real(8) r_, dr_, r0_

!     r0_ = 0.8*(Lx/2)    !Для кубика все r0_ равны
!     r_ = dsqrt(rx**2 + ry**2 + rz**2)
!     dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)

!     init_Ey_ = (0.5*erf((r_-r0_)/wE)+0.5) * (ry/(r_+1d-5)) !Ey

! end function init_Ey


! function init_Ez(rx, ry, rz) result (init_Ez_)
!     use all_var
!     real(8) rx, ry, rz
!     real(8) init_Ez_

!     real(8) r_, dr_, r0_

!     r0_ = 0.8*(Lx/2)    !Для кубика все r0_ равны
!     r_ = dsqrt(rx**2 + ry**2 + rz**2)
!     dr_ = dsqrt((rx-r0_)**2 + (ry-r0_)**2 + (rz-r0_)**2)

!     init_Ez_ = (0.5*erf((r_-r0_)/wE)+0.5) * (rz/(r_+1d-5)) !Ez

! end function init_Ez


function get_di(vc, jx, jy, jz) result (get_di_)
    integer(4) jx, jy, jz, vc
    integer(4) get_di_

    get_di_ = NxNyNz*(vc-1) + NxNy*(jz-1) + Nx*(jy-1) + jx

end function get_di


function fdot(v1, v2, ix, iy, iz) result(dot_) !на вход должны подаваться массивы(срезы)/вектора длины 3 * NxNyNz
    real(8) v1(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1), v2(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ix, iy, iz
    real(8) dot_
    real(8) dot_x, dot_y, dot_z

    dot_x = v1(1, ix, iy, iz) * v2(1, ix, iy, iz)
    dot_y = v1(2, ix, iy, iz) * v2(2, ix, iy, iz)
    dot_z = v1(3, ix, iy, iz) * v2(3, ix, iy, iz)
    dot_ = dot_x + dot_y + dot_z

end function fdot


function fcross(v1, v2, ix, iy, iz) result(cross_) !на вход должны подаваться массивы(срезы)/вектора длины 3 * NxNyNz
    real(8) v1(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1), v2(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ix, iy, iz
    real(8) cross_(3)

    cross_(1) = v1(2, ix, iy, iz) * v2(3, ix, iy, iz) - v1(3, ix, iy, iz) * v2(2, ix, iy, iz)
    cross_(2) = v1(3, ix, iy, iz) * v2(1, ix, iy, iz) - v1(1, ix, iy, iz) * v2(3, ix, iy, iz)
    cross_(3) = v1(1, ix, iy, iz) * v2(2, ix, iy, iz) - v1(2, ix, iy, iz) * v2(1, ix, iy, iz)

end function fcross

function fdUdx(Ui, ic, ix, iy, iz, edge) result (fdUdx_)
    real(8) Ui(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ic, ix, iy, iz, edge
    real(8) fdUdx_, Ui0, Ui1, Ui2

    if (ix==edge) then
        Ui0 = Ui(ic, ix, iy, iz)
        Ui1 = Ui(ic, ix+1, iy, iz)
        Ui2 = Ui(ic, ix+2, iy, iz)
        fdUdx_ = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hx)
    elseif (ix==(Nx+1-edge)) then
        Ui0 = Ui(ic, ix, iy, iz)
        Ui1 = Ui(ic, ix-1, iy, iz)
        Ui2 = Ui(ic, ix-2, iy, iz)
        fdUdx_ = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hx)
    else
        fdUdx_ = (Ui(ic, ix+1, iy, iz) - Ui(ic, ix-1, iy, iz)) / (2d0*hx)
    end if

end function fdUdx

!function fdUdx(Ui, ix, iy, iz, edge) result (fdUdx_)
!    real(8) Ui(0:Nx+1, 0:Ny+1, 0:Nz+1)
!    integer(4) ix, iy, iz, edge
!    real(8) fdUdx_, Ui0, Ui1, Ui2
!
!    if (ix==edge) then
!        Ui0 = Ui(ix, iy, iz)
!        Ui1 = Ui(ix+1, iy, iz)
!        Ui2 = Ui(ix+2, iy, iz)
!        fdUdx_ = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hx)
!    elseif (ix==(Nx+1-edge)) then
!        Ui0 = Ui(ix, iy, iz)
!        Ui1 = Ui(ix-1, iy, iz)
!        Ui2 = Ui(ix-2, iy, iz)
!        fdUdx_ = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hx)
!    else
!        fdUdx_ = (Ui(ix+1, iy, iz) - Ui(ix-1, iy, iz)) / (2d0*hx)
!    end if
!
!end function fdUdx


function fdUdy(Ui, ic, ix, iy, iz, edge) result (fdUdy_)
    real(8) Ui(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ic, ix, iy, iz, edge
    real(8) fdUdy_, Ui0, Ui1, Ui2

    if (iy==edge) then
        Ui0 = Ui(ic, ix, iy, iz)
        Ui1 = Ui(ic, ix, iy+1, iz)
        Ui2 = Ui(ic, ix, iy+2, iz)
        fdUdy_ = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hy)
    elseif (iy==(Ny+1-edge)) then
        Ui0 = Ui(ic, ix, iy, iz)
        Ui1 = Ui(ic, ix, iy-1, iz)
        Ui2 = Ui(ic, ix, iy-2, iz)
        fdUdy_ = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hy)
    else
        fdUdy_ = (Ui(ic, ix, iy+1, iz) - Ui(ic, ix, iy-1, iz)) / (2d0*hy)
    end if

end function fdUdy


function fdUdz(Ui, ic, ix, iy, iz, edge) result (fdUdz_)
    real(8) Ui(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ic, ix, iy, iz, edge
    real(8) fdUdz_, Ui0, Ui1, Ui2

    if (iz==edge) then
        Ui0 = Ui(ic, ix, iy, iz)
        Ui1 = Ui(ic, ix, iy, iz+1)
        Ui2 = Ui(ic, ix, iy, iz+2)
        fdUdz_ = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
    elseif (iz==(Nz+1-edge)) then
        Ui0 = Ui(ic, ix, iy, iz)
        Ui1 = Ui(ic, ix, iy, iz-1)
        Ui2 = Ui(ic, ix, iy, iz-2)
        fdUdz_ = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
    else
        fdUdz_ = (Ui(ic, ix, iy, iz+1) - Ui(ic, ix, iy, iz-1)) / (2d0*hz)
    end if

end function fdUdz


function fgrad(Ui, ix, iy, iz, edge) result (grad_) !на вход должен подаваться массив размера (1, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) Ui(1, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ix, iy, iz, edge
    real(8) grad_(3)
    real(8) Ui0, Ui1, Ui2

    if (ix==edge) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix+1, iy, iz)
        Ui2 = Ui(1, ix+2, iy, iz)
        grad_(1) = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hx)
    elseif (ix==(Nx+1-edge)) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix-1, iy, iz)
        Ui2 = Ui(1, ix-2, iy, iz)
        grad_(1) = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hx)
    else
        grad_(1) = (Ui(1, ix+1, iy, iz) - Ui(1, ix-1, iy, iz)) / (2d0*hx)
    end if

    if (iy==edge) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy+1, iz)
        Ui2 = Ui(1, ix, iy+2, iz)
        grad_(2) = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hy)
    elseif (iy==(Ny+1-edge)) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy-1, iz)
        Ui2 = Ui(1, ix, iy-2, iz)
        grad_(2) = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hy)
    else
        grad_(2) = (Ui(1, ix, iy+1, iz) - Ui(1, ix, iy-1, iz)) / (2d0*hy)
    end if

    if (iz==edge) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy, iz+1)
        Ui2 = Ui(1, ix, iy, iz+2)
        grad_(3) = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
    elseif (iz==(Nz+1-edge)) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy, iz-1)
        Ui2 = Ui(1, ix, iy, iz-2)
        grad_(3) = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
    else
        grad_(3) = (Ui(1, ix, iy, iz+1) - Ui(1, ix, iy, iz-1)) / (2d0*hz)
    end if

!    grad_(1) = fdUdx(Ui(:,:,:), ix, iy, iz, edge)
!    grad_(2) = fdUdy(Ui(:,:,:), ix, iy, iz, edge)
!    grad_(3) = fdUdz(Ui(:,:,:), ix, iy, iz, edge)

end function fgrad


function fdiv(Ui, ix, iy, iz, edge) result (div_) !на вход должен подаваться массив размера (3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) Ui(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ix, iy, iz, edge
    real(8) div_
    real(8) divX, divY, divZ
    real(8) Ui0, Ui1, Ui2

    if (ix==edge) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix+1, iy, iz)
        Ui2 = Ui(1, ix+2, iy, iz)
        divX = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hx)
    elseif (ix==(Nx+1-edge)) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix-1, iy, iz)
        Ui2 = Ui(1, ix-2, iy, iz)
        divX = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hx)
    else
        divX = (Ui(1, ix+1, iy, iz) - Ui(1, ix-1, iy, iz)) / (2d0*hx)
    end if

    if (iy==edge) then
        Ui0 = Ui(2, ix, iy, iz)
        Ui1 = Ui(2, ix, iy+1, iz)
        Ui2 = Ui(2, ix, iy+2, iz)
        divY = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hy)
    elseif (iy==(Ny+1-edge)) then
        Ui0 = Ui(2, ix, iy, iz)
        Ui1 = Ui(2, ix, iy-1, iz)
        Ui2 = Ui(2, ix, iy-2, iz)
        divY = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hy)
    else
        divY = (Ui(2, ix, iy+1, iz) - Ui(2, ix, iy-1, iz)) / (2d0*hy)
    end if

   if (iz==edge) then
        Ui0 = Ui(3, ix, iy, iz)
        Ui1 = Ui(3, ix, iy, iz+1)
        Ui2 = Ui(3, ix, iy, iz+2)
        divZ = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
    elseif (iz==(Nz+1-edge)) then
        Ui0 = Ui(3, ix, iy, iz)
        Ui1 = Ui(3, ix, iy, iz-1)
        Ui2 = Ui(3, ix, iy, iz-2)
        divZ = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
    else
        divZ = (Ui(3, ix, iy, iz+1) - Ui(3, ix, iy, iz-1)) / (2d0*hz)
    end if

!    divX = fdUdx(Ui(1,:,:,:), ix, iy, iz, edge)
!    divY = fdUdy(Ui(2,:,:,:), ix, iy, iz, edge)
!    divZ = fdUdz(Ui(3,:,:,:), ix, iy, iz, edge)
    div_ = divX + divY + divZ

end function fdiv


function frot(Ui, ix, iy, iz, edge) result(rot_) !на вход должен подаваться массив размера (3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    real(8) Ui(1:3, 0:Nx+1, 0:Ny+1, 0:Nz+1)
    integer(4) ix, iy, iz, edge
    real(8) rot_(3)
    real(8) dU3dy, dU2dz, dU3dx, dU1dz, dU2dx, dU1dy
    real(8) Ui0, Ui1, Ui2

    !dU3dy = fdUdy(Ui(3,:,:,:), ix, iy, iz, edge) !(U(get_di(3, ix, iy+1, iz)) - U(get_di(3, ix, iy-1, iz))) / (2 * hy)
    if (iy==edge) then
        Ui0 = Ui(3, ix, iy, iz)
        Ui1 = Ui(3, ix, iy+1, iz)
        Ui2 = Ui(3, ix, iy+2, iz)
        dU3dy = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hy)
    elseif (iy==(Ny+1-edge)) then
        Ui0 = Ui(3, ix, iy, iz)
        Ui1 = Ui(3, ix, iy-1, iz)
        Ui2 = Ui(3, ix, iy-2, iz)
        dU3dy = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hy)
    else
        dU3dy = (Ui(3, ix, iy+1, iz) - Ui(3, ix, iy-1, iz)) / (2d0*hy)
    end if

    !dU2dz = fdUdz(Ui(2,:,:,:), ix, iy, iz, edge) !(U(get_di(2, ix, iy, iz+1)) - U(get_di(2, ix, iy, iz-1))) / (2 * hz)
    if (iz==edge) then
        Ui0 = Ui(2, ix, iy, iz)
        Ui1 = Ui(2, ix, iy, iz+1)
        Ui2 = Ui(2, ix, iy, iz+2)
        dU2dz = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
    elseif (iz==(Nz+1-edge)) then
        Ui0 = Ui(2, ix, iy, iz)
        Ui1 = Ui(2, ix, iy, iz-1)
        Ui2 = Ui(2, ix, iy, iz-2)
        dU2dz = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
    else
        dU2dz = (Ui(2, ix, iy, iz+1) - Ui(2, ix, iy, iz-1)) / (2d0*hz)
    end if

    !dU3dx = fdUdx(Ui(3,:,:,:), ix, iy, iz, edge) !(U(get_di(3, ix+1, iy, iz)) - U(get_di(3, ix-1, iy, iz))) / (2 * hx)
    if (ix==edge) then
        Ui0 = Ui(3, ix, iy, iz)
        Ui1 = Ui(3, ix+1, iy, iz)
        Ui2 = Ui(3, ix+2, iy, iz)
        dU3dx = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hx)
    elseif (ix==(Nx+1-edge)) then
        Ui0 = Ui(3, ix, iy, iz)
        Ui1 = Ui(3, ix-1, iy, iz)
        Ui2 = Ui(3, ix-2, iy, iz)
        dU3dx = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hx)
    else
        dU3dx = (Ui(3, ix+1, iy, iz) - Ui(3, ix-1, iy, iz)) / (2d0*hx)
    end if

    !dU1dz = fdUdz(Ui(1,:,:,:), ix, iy, iz, edge) !(U(get_di(1, ix, iy, iz+1)) - U(get_di(1, ix, iy, iz-1))) / (2 * hz)
    if (iz==edge) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy, iz+1)
        Ui2 = Ui(1, ix, iy, iz+2)
        dU1dz = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hz)
    elseif (iz==(Nz+1-edge)) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy, iz-1)
        Ui2 = Ui(1, ix, iy, iz-2)
        dU1dz = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hz)
    else
        dU1dz = (Ui(1, ix, iy, iz+1) - Ui(1, ix, iy, iz-1)) / (2d0*hz)
    end if

    !dU2dx = fdUdx(Ui(2,:,:,:), ix, iy, iz, edge) !(U(get_di(2, ix+1, iy, iz)) - U(get_di(2, ix-1, iy, iz))) / (2 * hx)
    if (ix==edge) then
        Ui0 = Ui(2, ix, iy, iz)
        Ui1 = Ui(2, ix+1, iy, iz)
        Ui2 = Ui(2, ix+2, iy, iz)
        dU2dx = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hx)
    elseif (ix==(Nx+1-edge)) then
        Ui0 = Ui(2, ix, iy, iz)
        Ui1 = Ui(2, ix-1, iy, iz)
        Ui2 = Ui(2, ix-2, iy, iz)
        dU2dx = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hx)
    else
        dU2dx = (Ui(2, ix+1, iy, iz) - Ui(2, ix-1, iy, iz)) / (2d0*hx)
    end if

    !dU1dy = fdUdy(Ui(1,:,:,:), ix, iy, iz, edge) !(U(get_di(1, ix, iy+1, iz)) - U(get_di(1, ix, iy-1, iz))) / (2 * hy)
    if (iy==edge) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy+1, iz)
        Ui2 = Ui(1, ix, iy+2, iz)
        dU1dy = (-3d0*Ui0 + 4d0*Ui1 - Ui2) / (2d0*hy)
    elseif (iy==(Ny+1-edge)) then
        Ui0 = Ui(1, ix, iy, iz)
        Ui1 = Ui(1, ix, iy-1, iz)
        Ui2 = Ui(1, ix, iy-2, iz)
        dU1dy = (3d0*Ui0 - 4d0* Ui1 + Ui2) / (2d0*hy)
    else
        dU1dy = (Ui(1, ix, iy+1, iz) - Ui(1, ix, iy-1, iz)) / (2d0*hy)
    end if

    rot_(1) = dU3dy - dU2dz
    rot_(2) = dU1dz - dU3dx
    rot_(3) = dU2dx - dU1dy

end function frot


end module dif_oper
