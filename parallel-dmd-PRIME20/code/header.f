
      module GLOBAL

      integer, parameter :: k12=selected_int_kind(12)
      integer, parameter :: k1=selected_int_kind(1)
      real drandm,dtime
      integer iflag,time
      external dtime,drandm,time,srand

      character*4 fname_digits,fname_digits_pre
      common fname_digits,fname_digits_pre

!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      integer(kind=k1) ev_code((nop1+nop2),(nop1+nop2))
      integer ncim1,ncai,ncaj,nnjp1
!LR: Added a third species numbeads variable. I could not replace these with noptotal variables due to FORTRAN limitations
!LR: Added hp3 to hold hydrophobic data for species 3
      integer hp(numbeads1+numbeads2),hp1(numbeads1),hp2(numbeads2)!, hp3(numbeads3)
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      integer res(0:360,0:360),chnnum((nop1+nop2)),coltype((nop1+nop2)+3),extra_repuls((nop1+nop2),4)
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      integer bptnr((nop1+nop2)),identity((nop1+nop2)),nptnr((nop1+nop2)+3)
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
	integer npt((nop1+nop2)+1),nb(maxnbs*(nop1+nop2)+1),na_npt((nop1+nop2)),npt_dn((nop1+nop2)+1),dnnab(maxnbs*(nop1+nop2)+1),nnabdn((nop1+nop2))
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      integer tlinks((nop1+nop2)+3),tlinks2((nop1+nop2)+3),bin(0:numbin+1),nbin
      integer, dimension(:), allocatable::cell,wrap_map
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
!LR: Added firstside3 to hold side-chain data for species 3
      integer clinks((nop1+nop2)),num_cell,fside1(chnln1),fside2(chnln2)!,fside3(chnln3)
      integer map(n_nab_cell)
#ifdef equil
      integer(kind=k12) quarter, hb_sum(4),hh_sum(4)
      integer qt_count(4)
#endif
      integer(kind=k12) coll,ncoll
!      integer*8 coll,ncoll     
      logical success
!LR: Added a third species chnln variable. I could not replace these with noptotal variables due to FORTRAN limitations
      real analysis_boxl,boxl,boxl_orig,t,bdln(chnln1+chnln2),analysis_setemp,setemp,ev_param(3,50)
!LR: Here, I did not add a chnln3 variable, because these are bond-only arrays, and the nanoparticle cannot have bonds.
      real bl_rn(chnln1+chnln2),bl_rc(chnln1+chnln2),hdelr,rlsq(50),sv(6,(nop1+nop2))
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      real old_rx((nop1+nop2)),old_ry((nop1+nop2)),old_rz((nop1+nop2))
      real sigma(28),welldia(28),epsilon(28)
      real sigma_sq(28,28),sigma_2b(28,28),welldia_sq(28,28),ep_sqrt(28,28),sig_max_all,shlddia_sq(28,28)
      real tim((nop1+nop2)+3),tfalse,pi,interval,t_fact,interval_max,sortsize,tbin_off,n_forced
      real width, half
!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
      real ep(28,28),bds(28,28),wel(28,28),bm((nop1+nop2)) !! by mookyung
!LR: Here, I did not add a chnln3 variable, because these are bond-only arrays, and the nanoparticle cannot have bonds.
      real del_bdln(chnln1+chnln2),del_blrn(chnln1+chnln2),del_blrc(chnln1+chnln2) !! by mookyung
      real sqz610(5,28) !! by mookyung
      integer noptotal

!LR: Added a third species nop variable. I could not replace these with noptotal variables due to FORTRAN limitations
!LR: Note - I am not sure who originally coded the -glycine flag, but I found it to be extremely buggy and incomplete.
!LR: Edits I made within -glycine flagged code do not constitute my usage of that feature, I was just trying to be as complete with my edits as possible.
#ifdef glycine 
      real bdln_dummy((nop1+nop2))
#endif

      real shder_dist1,shder_dist2,shder_dist3,shder_dist4

      end module global

