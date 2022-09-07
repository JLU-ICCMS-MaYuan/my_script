! Module bondcrt defined in file ./FPP/BondCrt.fpp

subroutine f90wrap_storecfr(clycellv, atompsdr, atomid, cfr, n0, n1, n2, n3, n4)
    use bondcrt, only: storecfr
    implicit none
    
    real(8), intent(in), dimension(3,n0) :: clycellv
    real(8), intent(in), dimension(n1,n2) :: atompsdr
    integer, intent(in), dimension(n3) :: atomid
    real(8), intent(inout), dimension(n4) :: cfr
    integer :: n0
    !f2py intent(hide), depend(clycellv) :: n0 = shape(clycellv,1)
    integer :: n1
    !f2py intent(hide), depend(atompsdr) :: n1 = shape(atompsdr,0)
    integer :: n2
    !f2py intent(hide), depend(atompsdr) :: n2 = shape(atompsdr,1)
    integer :: n3
    !f2py intent(hide), depend(atomid) :: n3 = shape(atomid,0)
    integer :: n4
    !f2py intent(hide), depend(cfr) :: n4 = shape(cfr,0)
    call storecfr(clycellv=clycellv, atompsdr=atompsdr, atomid=atomid, cfr=cfr)
end subroutine f90wrap_storecfr

subroutine f90wrap_calstrdis(cfr1, cfr2, spgit1, spgit2, ret_calstrdis, ne, n0, n1)
    use bondcrt, only: calstrdis
    implicit none
    
    real(8), dimension(n0) :: cfr1
    real(8), dimension(n1) :: cfr2
    integer :: spgit1
    integer :: spgit2
    real(8), intent(out) :: ret_calstrdis
    integer :: ne
    integer :: n0
    !f2py intent(hide), depend(cfr1) :: n0 = shape(cfr1,0)
    integer :: n1
    !f2py intent(hide), depend(cfr2) :: n1 = shape(cfr2,0)
    ret_calstrdis = calstrdis(cfr1=cfr1, cfr2=cfr2, spgit1=spgit1, spgit2=spgit2, ne=ne)
end subroutine f90wrap_calstrdis

subroutine f90wrap_funcwr2(ret_funcwr2, x)
    use bondcrt, only: funcwr2
    implicit none
    
    real(8), intent(out) :: ret_funcwr2
    real(8) :: x
    ret_funcwr2 = funcwr2(x=x)
end subroutine f90wrap_funcwr2

subroutine f90wrap_funcwr3(ret_funcwr3, x)
    use bondcrt, only: funcwr3
    implicit none
    
    real(8), intent(out) :: ret_funcwr3
    real(8) :: x
    ret_funcwr3 = funcwr3(x=x)
end subroutine f90wrap_funcwr3

subroutine f90wrap_gasf(wf, b, gr, n0)
    use bondcrt, only: gasf
    implicit none
    
    real(8) :: wf
    real(8) :: b
    real(8), dimension(n0) :: gr
    integer :: n0
    !f2py intent(hide), depend(gr) :: n0 = shape(gr,0)
    call gasf(wf=wf, b=b, gr=gr)
end subroutine f90wrap_gasf

subroutine f90wrap_calcbond(cellv, atompsdr, atomid, n0, n1, n2, n3)
    use bondcrt, only: calcbond
    implicit none
    
    real(8), dimension(3,n0) :: cellv
    real(8), dimension(n1,n2) :: atompsdr
    integer, dimension(n3) :: atomid
    integer :: n0
    !f2py intent(hide), depend(cellv) :: n0 = shape(cellv,1)
    integer :: n1
    !f2py intent(hide), depend(atompsdr) :: n1 = shape(atompsdr,0)
    integer :: n2
    !f2py intent(hide), depend(atompsdr) :: n2 = shape(atompsdr,1)
    integer :: n3
    !f2py intent(hide), depend(atomid) :: n3 = shape(atomid,0)
    call calcbond(cellv=cellv, atompsdr=atompsdr, atomid=atomid)
end subroutine f90wrap_calcbond

subroutine f90wrap_vectorlength(ret_vectorlength, vc)
    use bondcrt, only: vectorlength
    implicit none
    
    real(8), intent(out) :: ret_vectorlength
    real(8), dimension(3) :: vc
    ret_vectorlength = vectorlength(vc=vc)
end subroutine f90wrap_vectorlength

subroutine f90wrap_recipvector(ret_recipvector, lat, n0, n1, n2)
    use bondcrt, only: recipvector
    implicit none
    
    real(8), intent(out), dimension(3,n0) :: ret_recipvector
    real(8), intent(in), dimension(n1,n2) :: lat
    integer :: n0
    integer :: n1
    !f2py intent(hide), depend(lat) :: n1 = shape(lat,0)
    integer :: n2
    !f2py intent(hide), depend(lat) :: n2 = shape(lat,1)
    ret_recipvector = recipvector(lat=lat)
end subroutine f90wrap_recipvector

subroutine f90wrap_volume(ret_volume, lat, n0, n1)
    use bondcrt, only: volume
    implicit none
    
    real(8), intent(out) :: ret_volume
    real(8), intent(in), dimension(n0,n1) :: lat
    integer :: n0
    !f2py intent(hide), depend(lat) :: n0 = shape(lat,0)
    integer :: n1
    !f2py intent(hide), depend(lat) :: n1 = shape(lat,1)
    ret_volume = volume(lat=lat)
end subroutine f90wrap_volume

subroutine f90wrap_crossp(va, ret_crossp, vb)
    use bondcrt, only: crossp
    implicit none
    
    real(8), intent(in), dimension(3) :: va
    real(8), dimension(3), intent(out) :: ret_crossp
    real(8), intent(in), dimension(3) :: vb
    ret_crossp = crossp(va=va, vb=vb)
end subroutine f90wrap_crossp

subroutine f90wrap_scc(xi, yi, ret_scc, npt, n0, n1)
    use bondcrt, only: scc
    implicit none
    
    real(8), dimension(n0) :: xi
    real(8), dimension(n1) :: yi
    real(8), intent(out) :: ret_scc
    integer(4) :: npt
    integer :: n0
    !f2py intent(hide), depend(xi) :: n0 = shape(xi,0)
    integer :: n1
    !f2py intent(hide), depend(yi) :: n1 = shape(yi,0)
    ret_scc = scc(xi=xi, yi=yi, npt=npt)
end subroutine f90wrap_scc

subroutine f90wrap_bondcrt__get__sn(f90wrap_sn)
    use bondcrt, only: bondcrt_sn => sn
    implicit none
    integer, intent(out) :: f90wrap_sn
    
    f90wrap_sn = bondcrt_sn
end subroutine f90wrap_bondcrt__get__sn

subroutine f90wrap_bondcrt__get__rcut(f90wrap_rcut)
    use bondcrt, only: bondcrt_rcut => rcut
    implicit none
    real(8), intent(out) :: f90wrap_rcut
    
    f90wrap_rcut = bondcrt_rcut
end subroutine f90wrap_bondcrt__get__rcut

subroutine f90wrap_bondcrt__set__rcut(f90wrap_rcut)
    use bondcrt, only: bondcrt_rcut => rcut
    implicit none
    real(8), intent(in) :: f90wrap_rcut
    
    bondcrt_rcut = f90wrap_rcut
end subroutine f90wrap_bondcrt__set__rcut

subroutine f90wrap_bondcrt__get__dr(f90wrap_dr)
    use bondcrt, only: bondcrt_dr => dr
    implicit none
    real(8), intent(out) :: f90wrap_dr
    
    f90wrap_dr = bondcrt_dr
end subroutine f90wrap_bondcrt__get__dr

subroutine f90wrap_bondcrt__set__dr(f90wrap_dr)
    use bondcrt, only: bondcrt_dr => dr
    implicit none
    real(8), intent(in) :: f90wrap_dr
    
    bondcrt_dr = f90wrap_dr
end subroutine f90wrap_bondcrt__set__dr

subroutine f90wrap_bondcrt__get__dr2(f90wrap_dr2)
    use bondcrt, only: bondcrt_dr2 => dr2
    implicit none
    real(8), intent(out) :: f90wrap_dr2
    
    f90wrap_dr2 = bondcrt_dr2
end subroutine f90wrap_bondcrt__get__dr2

subroutine f90wrap_bondcrt__set__dr2(f90wrap_dr2)
    use bondcrt, only: bondcrt_dr2 => dr2
    implicit none
    real(8), intent(in) :: f90wrap_dr2
    
    bondcrt_dr2 = f90wrap_dr2
end subroutine f90wrap_bondcrt__set__dr2

subroutine f90wrap_bondcrt__get__outr(f90wrap_outr)
    use bondcrt, only: bondcrt_outr => outr
    implicit none
    real(8), intent(out) :: f90wrap_outr
    
    f90wrap_outr = bondcrt_outr
end subroutine f90wrap_bondcrt__get__outr

subroutine f90wrap_bondcrt__set__outr(f90wrap_outr)
    use bondcrt, only: bondcrt_outr => outr
    implicit none
    real(8), intent(in) :: f90wrap_outr
    
    bondcrt_outr = f90wrap_outr
end subroutine f90wrap_bondcrt__set__outr

subroutine f90wrap_bondcrt__get__inr(f90wrap_inr)
    use bondcrt, only: bondcrt_inr => inr
    implicit none
    real(8), intent(out) :: f90wrap_inr
    
    f90wrap_inr = bondcrt_inr
end subroutine f90wrap_bondcrt__get__inr

subroutine f90wrap_bondcrt__set__inr(f90wrap_inr)
    use bondcrt, only: bondcrt_inr => inr
    implicit none
    real(8), intent(in) :: f90wrap_inr
    
    bondcrt_inr = f90wrap_inr
end subroutine f90wrap_bondcrt__set__inr

subroutine f90wrap_bondcrt__get__pi(f90wrap_pi)
    use bondcrt, only: bondcrt_pi => pi
    implicit none
    real(8), intent(out) :: f90wrap_pi
    
    f90wrap_pi = bondcrt_pi
end subroutine f90wrap_bondcrt__get__pi

subroutine f90wrap_bondcrt__set__pi(f90wrap_pi)
    use bondcrt, only: bondcrt_pi => pi
    implicit none
    real(8), intent(in) :: f90wrap_pi
    
    bondcrt_pi = f90wrap_pi
end subroutine f90wrap_bondcrt__set__pi

subroutine f90wrap_bondcrt__get__delta(f90wrap_delta)
    use bondcrt, only: bondcrt_delta => delta
    implicit none
    real(8), intent(out) :: f90wrap_delta
    
    f90wrap_delta = bondcrt_delta
end subroutine f90wrap_bondcrt__get__delta

subroutine f90wrap_bondcrt__set__delta(f90wrap_delta)
    use bondcrt, only: bondcrt_delta => delta
    implicit none
    real(8), intent(in) :: f90wrap_delta
    
    bondcrt_delta = f90wrap_delta
end subroutine f90wrap_bondcrt__set__delta

subroutine f90wrap_bondcrt__get__deltar(f90wrap_deltar)
    use bondcrt, only: bondcrt_deltar => deltar
    implicit none
    real(8), intent(out) :: f90wrap_deltar
    
    f90wrap_deltar = bondcrt_deltar
end subroutine f90wrap_bondcrt__get__deltar

subroutine f90wrap_bondcrt__set__deltar(f90wrap_deltar)
    use bondcrt, only: bondcrt_deltar => deltar
    implicit none
    real(8), intent(in) :: f90wrap_deltar
    
    bondcrt_deltar = f90wrap_deltar
end subroutine f90wrap_bondcrt__set__deltar

subroutine f90wrap_bondcrt__array__mtr(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_mtr => mtr
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(bondcrt_mtr)
    dloc = loc(bondcrt_mtr)
end subroutine f90wrap_bondcrt__array__mtr

subroutine f90wrap_bondcrt__array__tpfr(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_tpfr => tpfr
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    dshape(1:1) = shape(bondcrt_tpfr)
    dloc = loc(bondcrt_tpfr)
end subroutine f90wrap_bondcrt__array__tpfr

subroutine f90wrap_bondcrt__array__minbd(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_minbd => minbd
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 1
    dtype = 12
    if (allocated(bondcrt_minbd)) then
        dshape(1:1) = shape(bondcrt_minbd)
        dloc = loc(bondcrt_minbd)
    else
        dloc = 0
    end if
end subroutine f90wrap_bondcrt__array__minbd

subroutine f90wrap_bondcrt__array__wfscc(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_wfscc => wfscc
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(bondcrt_wfscc)) then
        dshape(1:2) = shape(bondcrt_wfscc)
        dloc = loc(bondcrt_wfscc)
    else
        dloc = 0
    end if
end subroutine f90wrap_bondcrt__array__wfscc

subroutine f90wrap_bondcrt__array__wfscc2(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_wfscc2 => wfscc2
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(bondcrt_wfscc2)) then
        dshape(1:2) = shape(bondcrt_wfscc2)
        dloc = loc(bondcrt_wfscc2)
    else
        dloc = 0
    end if
end subroutine f90wrap_bondcrt__array__wfscc2

subroutine f90wrap_bondcrt__array__bondinf(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_bondinf => bondinf
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 12
    if (allocated(bondcrt_bondinf)) then
        dshape(1:2) = shape(bondcrt_bondinf)
        dloc = loc(bondcrt_bondinf)
    else
        dloc = 0
    end if
end subroutine f90wrap_bondcrt__array__bondinf

subroutine f90wrap_bondcrt__array__ibondinf(dummy_this, nd, dtype, dshape, dloc)
    use bondcrt, only: bondcrt_ibondinf => ibondinf
    implicit none
    integer, intent(in) :: dummy_this(2)
    integer, intent(out) :: nd
    integer, intent(out) :: dtype
    integer, dimension(10), intent(out) :: dshape
    integer*8, intent(out) :: dloc
    
    nd = 2
    dtype = 5
    if (allocated(bondcrt_ibondinf)) then
        dshape(1:2) = shape(bondcrt_ibondinf)
        dloc = loc(bondcrt_ibondinf)
    else
        dloc = 0
    end if
end subroutine f90wrap_bondcrt__array__ibondinf

subroutine f90wrap_bondcrt__get__etypenub(f90wrap_etypenub)
    use bondcrt, only: bondcrt_etypenub => etypenub
    implicit none
    integer, intent(out) :: f90wrap_etypenub
    
    f90wrap_etypenub = bondcrt_etypenub
end subroutine f90wrap_bondcrt__get__etypenub

subroutine f90wrap_bondcrt__set__etypenub(f90wrap_etypenub)
    use bondcrt, only: bondcrt_etypenub => etypenub
    implicit none
    integer, intent(in) :: f90wrap_etypenub
    
    bondcrt_etypenub = f90wrap_etypenub
end subroutine f90wrap_bondcrt__set__etypenub

subroutine f90wrap_bondcrt__get__atomnub(f90wrap_atomnub)
    use bondcrt, only: bondcrt_atomnub => atomnub
    implicit none
    integer, intent(out) :: f90wrap_atomnub
    
    f90wrap_atomnub = bondcrt_atomnub
end subroutine f90wrap_bondcrt__get__atomnub

subroutine f90wrap_bondcrt__set__atomnub(f90wrap_atomnub)
    use bondcrt, only: bondcrt_atomnub => atomnub
    implicit none
    integer, intent(in) :: f90wrap_atomnub
    
    bondcrt_atomnub = f90wrap_atomnub
end subroutine f90wrap_bondcrt__set__atomnub

subroutine f90wrap_bondcrt__get__ifr(f90wrap_ifr)
    use bondcrt, only: bondcrt_ifr => ifr
    implicit none
    integer, intent(out) :: f90wrap_ifr
    
    f90wrap_ifr = bondcrt_ifr
end subroutine f90wrap_bondcrt__get__ifr

subroutine f90wrap_bondcrt__set__ifr(f90wrap_ifr)
    use bondcrt, only: bondcrt_ifr => ifr
    implicit none
    integer, intent(in) :: f90wrap_ifr
    
    bondcrt_ifr = f90wrap_ifr
end subroutine f90wrap_bondcrt__set__ifr

subroutine f90wrap_bondcrt__get__ctnt(f90wrap_ctnt)
    use bondcrt, only: bondcrt_ctnt => ctnt
    implicit none
    integer, intent(out) :: f90wrap_ctnt
    
    f90wrap_ctnt = bondcrt_ctnt
end subroutine f90wrap_bondcrt__get__ctnt

subroutine f90wrap_bondcrt__set__ctnt(f90wrap_ctnt)
    use bondcrt, only: bondcrt_ctnt => ctnt
    implicit none
    integer, intent(in) :: f90wrap_ctnt
    
    bondcrt_ctnt = f90wrap_ctnt
end subroutine f90wrap_bondcrt__set__ctnt

subroutine f90wrap_bondcrt__get__lsetup(f90wrap_lsetup)
    use bondcrt, only: bondcrt_lsetup => lsetup
    implicit none
    logical, intent(out) :: f90wrap_lsetup
    
    f90wrap_lsetup = bondcrt_lsetup
end subroutine f90wrap_bondcrt__get__lsetup

subroutine f90wrap_bondcrt__set__lsetup(f90wrap_lsetup)
    use bondcrt, only: bondcrt_lsetup => lsetup
    implicit none
    logical, intent(in) :: f90wrap_lsetup
    
    bondcrt_lsetup = f90wrap_lsetup
end subroutine f90wrap_bondcrt__set__lsetup

subroutine f90wrap_bondcrt__get__l2d(f90wrap_l2d)
    use bondcrt, only: bondcrt_l2d => l2d
    implicit none
    logical, intent(out) :: f90wrap_l2d
    
    f90wrap_l2d = bondcrt_l2d
end subroutine f90wrap_bondcrt__get__l2d

subroutine f90wrap_bondcrt__set__l2d(f90wrap_l2d)
    use bondcrt, only: bondcrt_l2d => l2d
    implicit none
    logical, intent(in) :: f90wrap_l2d
    
    bondcrt_l2d = f90wrap_l2d
end subroutine f90wrap_bondcrt__set__l2d

subroutine f90wrap_bondcrt__get__lct(f90wrap_lct)
    use bondcrt, only: bondcrt_lct => lct
    implicit none
    logical, intent(out) :: f90wrap_lct
    
    f90wrap_lct = bondcrt_lct
end subroutine f90wrap_bondcrt__get__lct

subroutine f90wrap_bondcrt__set__lct(f90wrap_lct)
    use bondcrt, only: bondcrt_lct => lct
    implicit none
    logical, intent(in) :: f90wrap_lct
    
    bondcrt_lct = f90wrap_lct
end subroutine f90wrap_bondcrt__set__lct

! End of module bondcrt defined in file ./FPP/BondCrt.fpp

