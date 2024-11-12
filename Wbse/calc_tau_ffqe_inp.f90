!
! Copyright (C) 2015-2024 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Ngoc Linh Nguyen, Victor Yu
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau_ffqe_inp()
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE pwcom,                ONLY : isk,npw,ngk
  USE wavefunctions,        ONLY : evc
  USE io_global,            ONLY : stdout
  USE westcom,              ONLY : lrwfc,iuwfc,ev,dvg,n_pdep_eigen_to_use,npwqx,nbnd_occ,l_pdep,&
                                 & spin_channel,l_bse
  USE lsda_mod,             ONLY : nspin
  USE pdep_db,              ONLY : pdep_db_read
  USE mp,                   ONLY : mp_bcast
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE wbse_dv,              ONLY : wbse_dv_setup
  USE buffers,              ONLY : get_buffer
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : pert,kpt_pool
  USE qbox_interface,       ONLY : init_qbox,finalize_qbox
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : using_evc,using_evc_d
#endif
  !
  IMPLICIT NONE
  !
  ! Workspace
  !
  INTEGER :: iks,current_spin
  INTEGER :: nbndval
  LOGICAL :: spin_resolve
  !
  CALL wbse_dv_setup(.FALSE.)
  !
  DO iks = 1,kpt_pool%nloc
     !
     current_spin = isk(iks)
     npw = ngk(iks)
     nbndval = nbnd_occ(iks)
     !
     IF(kpt_pool%nloc > 1) THEN
        IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
        CALL mp_bcast(evc,0,inter_image_comm)
        !
#if defined(__CUDA)
        CALL using_evc(2)
        CALL using_evc_d(0)
#endif
     ENDIF
     !
     CALL calc_tau_ffkc(current_spin,nbndval)
     CALL calc_tau_ffqe(current_spin,nbndval)
     !
  ENDDO
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau_ffkc(current_spin,nbndval)
  !-----------------------------------------------------------------------
  !
  ! Compute and store screened exchange integrals tau in case of BSE, or unscreened exchange
  ! integrals tau_u in case of hybrid TDDFT
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE io_push,              ONLY : io_push_title
  USE io_global,            ONLY : stdout
  USE types_coulomb,        ONLY : pot3D
  USE westcom,              ONLY : ev,dvg,wbse_init_calculation,wbse_init_save_dir,l_bse,l_pdep,&
                                 & chi_kernel,l_local_repr,overlap_thr,n_trunc_bands
  USE fft_base,             ONLY : dffts, dfftp
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE noncollin_module,     ONLY : npol
  USE gvect,                ONLY : ngm
  USE pwcom,                ONLY : npw,npwx,lsda, nbnd
  USE pdep_io,              ONLY : pdep_merge_and_write_G
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : bandpair
  USE mp,                   ONLY : mp_barrier,mp_sum
  USE lsda_mod,             ONLY : nspin
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_invfft_gamma
  USE mp_global,            ONLY : intra_image_comm,my_image_id,inter_bgrp_comm,nbgrp,&
                                 & intra_bgrp_comm,me_bgrp
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wbse_dv,              ONLY : wbse_dv_setup,wbse_dv_of_drho
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc_work=>evc_d,psic=>psic_d
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin,nbndval
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,ibnd_g,jbnd_g,tmp_size
  INTEGER :: il1,ig1,ir,do_idx,is
  INTEGER :: iu,ig,ierr,ip,ip_g
  INTEGER :: nbnd_do
  INTEGER :: dffts_nnr
  REAL(DP) :: uPi
  REAL(DP) :: factor
  REAL(DP) :: ovl_value
  REAL(DP), ALLOCATABLE :: ovl_matrix(:,:)
  INTEGER, ALLOCATABLE :: idx_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_loc(:,:)
  !
  COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux1_r(:,:),aux2_g(:,:), aux3(:)
  REAL(DP), ALLOCATABLE :: aux_r_test(:,:)
  COMPLEX(DP), ALLOCATABLE :: aux_g_test(:,:)
  !$acc declare device_resident(aux_r)
  !
  CHARACTER(LEN=:), ALLOCATABLE :: fname
  CHARACTER :: slabel
  CHARACTER(LEN=6) :: ilabel,jlabel
  CHARACTER(LEN=6) :: flabel
  !
  LOGICAL :: l_xcchi,l_skip,l_restart
  !
  TYPE(bar_type) :: barra
  !
  nbnd_do = nbnd
  !
  dffts_nnr = dffts%nnr
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  tmp_size = nbnd_do*nbnd_do
  ALLOCATE(idx_matrix(tmp_size,2))
  ALLOCATE(ovl_matrix(nbnd_do,nbnd_do))
  !
  ovl_matrix(:,:) = 0._DP
  idx_matrix(:,:) = 0
  !
  IF(l_local_repr) THEN
     !
     ALLOCATE(evc_loc(npwx,nbnd_do))
     !$acc enter data create(evc_loc)
     !
     CALL wbse_localization(current_spin,n_trunc_bands+1,nbndval,evc_loc,ovl_matrix,l_restart)
     !
  ENDIF
  !
  ! compute idx_matrix
  !
  do_idx = 0
  !
  DO jbnd = 1,nbnd_do
     DO ibnd = 1,nbnd_do
        !
        ovl_value = ovl_matrix(ibnd,jbnd)
        !
        IF(ovl_value >= overlap_thr) THEN
           IF(jbnd >= ibnd) THEN
              do_idx = do_idx+1
              !
              idx_matrix(do_idx,1) = ibnd+n_trunc_bands
              idx_matrix(do_idx,2) = jbnd+n_trunc_bands
           ENDIF
        ENDIF
        !
     ENDDO
  ENDDO
  !
  DEALLOCATE(ovl_matrix)
  !
  ! initialize the paralellization
  !
  bandpair = idistribute()
  CALL bandpair%init(do_idx,'i','n_pairs',.TRUE.)
  !
  ALLOCATE(aux_r(dffts%nnr))
  !$acc enter data create(tau,aux1_g) copyin(dvg)
  !
  CALL io_push_title('Applying Vc kernel to orbital density for KC screening factors')
  !
  CALL start_bar_type(barra,'tau',bandpair%nlocx)
  !
  ! parallel loop
  !
  DO il1 = 1,bandpair%nlocx
     !
     ig1  = bandpair%l2g(il1) ! global index of n_total
     !
     ibnd_g = idx_matrix(ig1,1)
     jbnd_g = idx_matrix(ig1,2)
     ibnd = ibnd_g-n_trunc_bands
     jbnd = jbnd_g-n_trunc_bands
     !
     l_skip = .FALSE.
     !
     IF (ibnd_g .NE. jbnd_g) l_skip = .TRUE.
     !
     IF(ig1 < 1 .OR. ig1 > do_idx) l_skip = .TRUE.
     !
     IF(.NOT. l_skip) THEN
        !
        IF(l_local_repr) THEN
           !$acc host_data use_device(evc_loc)
           CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
           !$acc end host_data
        ELSE
           CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd_g),evc_work(:,jbnd_g),psic,'Wave')
        ENDIF
        !
        !$acc parallel loop present(aux_r)
        DO ir = 1,dffts_nnr
           aux_r(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
        ENDDO
        !$acc end parallel
        !
        ALLOCATE(aux1_r(dffts%nnr,nspin))
        aux1_r(:,:) = (0._DP,0._DP)
        aux1_r(:,current_spin) = aux_r(:)
        !
        !CALL wbse_dv_of_drho_kc(current_spin, aux1_r,.FALSE.,.FALSE.)
        CALL wbse_dv_of_drho(aux1_r,.TRUE.,.FALSE.)
        !
        uPi=0.D0
        DO ir = 1, dfftp%nnr
           uPi = uPi + ( aux_r(ir) * DBLE(aux1_r(ir, current_spin)) )
        ENDDO
        uPi = uPi *omega/(dfftp%nr1*dfftp%nr2*dfftp%nr3)
        CALL mp_sum( uPi , intra_bgrp_comm )
        !
        WRITE(stdout,*) "Orbital Unrelaxed Koopmans uPi=", jbnd_g, current_spin, uPi
        !
        ALLOCATE( aux2_g(dfftp%ngm, nspin ))
        aux2_g(:,:) = (0.0_DP, 0.0_DP)
        !
        ALLOCATE( aux3( dfftp%nnr ) )
        DO is = 1, nspin
           aux3(:) = (0.0_DP, 0.0_DP)
           aux3(:) = CMPLX(DBLE(aux1_r(:, is)),0.0_DP)
           CALL fwfft('Rho', aux3, dfftp)
           aux2_g(:, is) = aux3(dfftp%nl(:))
        ENDDO
        DEALLOCATE( aux3 )
        !
        ! write vc_rho to disk
        !
        WRITE(ilabel,'(i6.6)') ibnd_g
        WRITE(jlabel,'(i6.6)') jbnd_g
        WRITE(slabel,'(i1)') current_spin
        !
        fname = 'PertpotKC'//ilabel//'_'//jlabel//'_'//slabel
        CALL write_rho_xml_ff(aux2_g, fname)
        !
        DEALLOCATE(aux1_r)
        DEALLOCATE(aux2_g)
        !
     ENDIF
     !
     !CALL update_bar_type(barra,'tau',1)
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'tau')
  !
  !$acc exit data delete(dvg,tau,aux1_g)
  DEALLOCATE(aux_r)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
  DEALLOCATE(idx_matrix)
  !
  IF(l_local_repr) THEN
     !$acc exit data delete(evc_loc)
     DEALLOCATE(evc_loc)
  ENDIF
  !
END SUBROUTINE
!
!
!-----------------------------------------------------------------------
SUBROUTINE calc_tau_ffqe(current_spin,nbndval)
  !-----------------------------------------------------------------------
  !
  ! Compute and store screened exchange integrals tau in case of BSE, or unscreened exchange
  ! integrals tau_u in case of hybrid TDDFT
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega
  USE io_push,              ONLY : io_push_title
  USE io_global,             ONLY : stdout
  USE types_coulomb,        ONLY : pot3D
  USE westcom,              ONLY : ev,dvg,wbse_init_calculation,wbse_init_save_dir,l_bse,l_pdep,&
                                 & chi_kernel,l_local_repr,overlap_thr,n_trunc_bands
  USE fft_base,             ONLY : dffts
      USE fft_base,              ONLY : dfftp
    USE fft_interfaces,        ONLY : fwfft,invfft
  USE noncollin_module,     ONLY : npol
  USE gvect,                 ONLY : ngm
  USE pwcom,                ONLY : npw,npwx,lsda, nbnd
  USE pdep_io,              ONLY : pdep_merge_and_write_G
  USE class_idistribute,    ONLY : idistribute
  USE distribution_center,  ONLY : bandpair
  USE mp,                   ONLY : mp_barrier,mp_sum
  USE lsda_mod,             ONLY : nspin
  USE fft_at_gamma,         ONLY : single_fwfft_gamma,single_invfft_gamma,double_invfft_gamma
  USE mp_global,            ONLY : intra_image_comm,my_image_id,inter_bgrp_comm,nbgrp,&
                                 & intra_bgrp_comm,me_bgrp
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE wbse_dv,              ONLY : wbse_dv_setup,wbse_dv_of_drho
#if defined(__CUDA)
  USE wavefunctions_gpum,   ONLY : evc_work=>evc_d,psic=>psic_d
  USE west_gpu,             ONLY : allocate_gpu,deallocate_gpu
#else
  USE wavefunctions,        ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin,nbndval
  !
  ! Workspace
  !
  INTEGER :: ibnd,jbnd,ibnd_g,jbnd_g,tmp_size
  INTEGER :: il1,ig1,ir,do_idx,is
  INTEGER :: iu,ig,ierr,ip,ip_g
  INTEGER :: nbnd_do
  INTEGER :: dffts_nnr
  REAL(DP) :: uPi
  REAL(DP) :: factor
  REAL(DP) :: ovl_value
  REAL(DP), ALLOCATABLE :: ovl_matrix(:,:)
  INTEGER, ALLOCATABLE :: idx_matrix(:,:)
  COMPLEX(DP), ALLOCATABLE :: evc_loc(:,:)
  COMPLEX(DP), ALLOCATABLE :: tau(:)
  !
  COMPLEX(DP), ALLOCATABLE :: aux_r(:),aux1_r(:,:),aux2_g(:,:), aux3(:)
  !$acc declare device_resident(aux_r)
  !
  CHARACTER(LEN=:), ALLOCATABLE :: fname
  CHARACTER :: slabel
  CHARACTER(LEN=6) :: ilabel,jlabel
  CHARACTER(LEN=6) :: flabel
  !
  LOGICAL :: l_xcchi,l_skip,l_restart
  !
  TYPE(bar_type) :: barra
  !
  nbnd_do = nbndval-n_trunc_bands
  !
  dffts_nnr = dffts%nnr
  !
#if defined(__CUDA)
  CALL allocate_gpu()
#endif
  !
  tmp_size = nbnd_do*nbnd_do
  ALLOCATE(idx_matrix(tmp_size,2))
  ALLOCATE(ovl_matrix(nbnd_do,nbnd_do))
  !
  ovl_matrix(:,:) = 0._DP
  idx_matrix(:,:) = 0
  !
  IF(l_local_repr) THEN
     !
     ALLOCATE(evc_loc(npwx,nbnd_do))
     !$acc enter data create(evc_loc)
     !
     CALL wbse_localization(current_spin,n_trunc_bands+1,nbndval,evc_loc,ovl_matrix,l_restart)
     !
  ENDIF
  !
  ! compute idx_matrix
  !
  do_idx = 0
  !
  DO jbnd = 1,nbnd_do
     DO ibnd = 1,nbnd_do
        !
        ovl_value = ovl_matrix(ibnd,jbnd)
        !
        IF(ovl_value >= overlap_thr) THEN
           IF(jbnd >= ibnd) THEN
              do_idx = do_idx+1
              !
              idx_matrix(do_idx,1) = ibnd+n_trunc_bands
              idx_matrix(do_idx,2) = jbnd+n_trunc_bands
           ENDIF
        ENDIF
        !
     ENDDO
  ENDDO
  !
  DEALLOCATE(ovl_matrix)
  !
  ! initialize the paralellization
  !
  bandpair = idistribute()
  CALL bandpair%init(do_idx,'i','n_pairs',.TRUE.)
  !
  ALLOCATE(aux_r(dffts%nnr))
  !$acc enter data create(tau,aux1_g) copyin(dvg)
  !
  CALL io_push_title('Applying Vc kernel to orbital density for BSE screening')
  !
  CALL start_bar_type(barra,'tau',bandpair%nlocx)
  !
  ! parallel loop
  !
  DO il1 = 1,bandpair%nlocx
     !
     ig1  = bandpair%l2g(il1) ! global index of n_total
     !
     ibnd_g = idx_matrix(ig1,1)
     jbnd_g = idx_matrix(ig1,2)
     ibnd = ibnd_g-n_trunc_bands
     jbnd = jbnd_g-n_trunc_bands
     !
     l_skip = .FALSE.
     !
     IF(ig1 < 1 .OR. ig1 > do_idx) l_skip = .TRUE.
     !
     IF(.NOT. l_skip) THEN
        !
        IF(l_local_repr) THEN
           !$acc host_data use_device(evc_loc)
           CALL double_invfft_gamma(dffts,npw,npwx,evc_loc(:,ibnd),evc_loc(:,jbnd),psic,'Wave')
           !$acc end host_data
        ELSE
           CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibnd_g),evc_work(:,jbnd_g),psic,'Wave')
        ENDIF
        !
        !$acc parallel loop present(aux_r)
        DO ir = 1,dffts_nnr
           aux_r(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
        ENDDO
        !$acc end parallel
        !
        ALLOCATE(aux1_r(dffts%nnr,nspin))
        aux1_r(:,:) = (0._DP,0._DP)
        aux1_r(:,current_spin) = aux_r(:)
        !
        ! aux1_r = vc*aux1_r()
        !
        CALL wbse_dv_of_drho(aux1_r,.TRUE.,.FALSE.)
        !
        ALLOCATE( aux2_g( ngm, nspin ))
        ALLOCATE( aux3( dfftp%nnr ) )
        DO is = 1, nspin
           aux3(:) = CMPLX(dble(aux1_r(:,is)),0.D0,kind=dp)
           !$acc host_data use_device(dvscf)
           CALL fwfft('Rho', aux3(:), dfftp)
           !$acc end host_data
           aux2_g(:,is) = aux3(dfftp%nl(:))
        ENDDO
        DEALLOCATE( aux3 )
        !
        ! write vc_rho to disk
        !
        WRITE(ilabel,'(i6.6)') ibnd_g
        WRITE(jlabel,'(i6.6)') jbnd_g
        WRITE(slabel,'(i1)') current_spin
        !
        fname = 'PertpotBSE'//ilabel//'_'//jlabel//'_'//slabel
        CALL write_rho_xml_ff(aux2_g, fname)
        !
        DEALLOCATE(aux2_g)
        DEALLOCATE(aux1_r)
        !
     ENDIF
     !
     CALL update_bar_type(barra,'tau',1)
     !
  ENDDO
  !
  CALL stop_bar_type(barra,'tau')
  !
  !$acc exit data delete(dvg,tau,aux1_g)
  DEALLOCATE(aux_r)
  !
#if defined(__CUDA)
  CALL deallocate_gpu()
#endif
  !
  DEALLOCATE(idx_matrix)
  !
  IF(l_local_repr) THEN
     !$acc exit data delete(evc_loc)
     DEALLOCATE(evc_loc)
  ENDIF
  !
END SUBROUTINE
