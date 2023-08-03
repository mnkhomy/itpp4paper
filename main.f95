program main
    use for_ode_mod
    use omp_lib
    use all_var

    real(4) total_time
    integer :: nthreads

    nthreads = 4
    total_time = secnds(0e0)
    
    call omp_set_num_threads(nthreads)
    !$omp parallel
    print *,"job done!"
    !$omp end parallel

    call solve()

    total_time = secnds(total_time)

    bdn_func_ss = "NDND"
    bdq_func_ss = "zero"
    if (bdn_func == 0) then
        bdn_func_ss = "zero"
    elseif (bdn_func == 100) then
        bdn_func_ss = "init"
    endif

    open(11, file='./results/log.txt')
        write(11, '("scheme:"a4" / shape:"a4)') ss, ss_shape
        write(11, '("nthreads used: "I0.0)') nthreads
        write(11, '("total time taken: "F0.1" min")') total_time/60
        write(11, '("Nx: "I0.0," Ny: "I0.0, " Nz: "I0.0)') Nx, Ny, Nz
        write(11, '("T_end: "F0.1)') T_end
        write(11, '("ht: "F0.10)') ht
        write(11, '("save_step_planes: "I0.0)') save_step_planes
        write(11, '("Last_saved_plane: "I0.0)') last_saved_plane
        write(11, '(" ")')
        write(11, '("A = "es10.2)') A
        write(11, '("B = "es10.2)') B
        write(11, '("C = "es10.2)') C
        write(11, '("Lbig = "es10.2)') Lbig
        write(11, '("Qroot = "F0.4)') Qroot
        write(11, '("Dpar = "F0.4)') Dpar
        write(11, '("Dper = "F0.4)') Dper
        write(11, '("Da = "F0.4)') Da
        write(11, '("mu1 = "F0.4)') mu1
        write(11, '("mu3 = "F0.4)') mu3
        write(11, '("k1 = "F0.4)') k1
        write(11, '("k2 = "F0.4)') k2
        write(11, '("k3 = "F0.4)') k3
        write(11, '("k4 = "F0.4)') k4        
        write(11, '("k23 = "F0.4)') k23
        write(11, '("k24 = "F0.4)') k24
        write(11, '("qh = "F0.4)') qh
        write(11, '("kap2 = "F0.4)') kap2
        write(11, '("phi0 = "F0.1)') phi_zero
        write(11, '("tet0 = "F0.1)') tet_zero
        write(11, '("bounds:")')
        write(11, '("n:("F0.1,"; ", F0.1,"; " F0.1") "a4)') bdn_a, bdn_b, bdn_c, bdn_func_ss
        write(11, '("q:("F0.1,"; ", F0.1,"; " F0.1") "a4)') bdq_a, bdq_b, bdq_c, bdq_func_ss
        write(11, '("Rsph = "F0.4)') Rsph
    close(11)

    open(11, file='./results/rsph.txt')
        write(11, '(F0.4)') Rsph
    close(11)

end program main

