!
! Copyright (C) 2015-2023 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Yu Jin
!
!-----------------------------------------------------------------------
SUBROUTINE hybrid_kernel_term1_slow(current_spin, evc1, hybrid_kd1, sf)
  !-----------------------------------------------------------------------
  !
  ! \sum_{v'} (\int v_c \phi_{v'} \phi_{v}) a_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma
  USE mp_global,             ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,nbgrp,my_bgrp_id
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : kpt_pool,band_group
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  COMPLEX(DP), INTENT(INOUT) :: hybrid_kd1(npwx,band_group%nlocx)
  !
  ! Workspace
  !
  INTEGER :: lbnd, ibnd, ibndp, jbnd, jbndp, ir, ig, iks_do
  INTEGER :: nbndval, nbnd_do, current_spin_ikq, ikq
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: aux_hybrid1(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: aux_hybrid1
#endif
  COMPLEX(DP), ALLOCATABLE :: caux(:), gaux(:), raux(:)
  !$acc declare device_resident(caux,gaux,raux)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_hybrid1(npwx,nbndval0x-n_trunc_bands))
  !$acc enter data create(aux_hybrid1)
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,iks_do)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        !$acc kernels present(raux)
        raux(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        DO lbnd = 1, nbnd_do ! index to be summed
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           ! product of evc and evc
           !
           CALL double_invfft_gamma(dffts,npw,npwx,evc_work(:,ibndp),evc_work(:,jbndp),psic,'Wave')
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
           ENDDO
           !$acc end parallel
           !
           ! Apply the bare Coulomb potential
           !
           !$acc host_data use_device(caux,gaux)
           CALL single_fwfft_gamma(dffts,npw,npwx,caux,gaux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(gaux,pot3D)
           DO ig = 1, npw
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(gaux,evc1,caux)
           CALL double_invfft_gamma(dffts,npw,npwx,gaux,evc1(:,lbnd,ikq),caux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop present(raux)
           DO ir = 1, dffts_nnr
              raux(ir) = raux(ir)+psic(ir)
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        !$acc host_data use_device(raux,aux_hybrid1)
        CALL single_fwfft_gamma(dffts,npw,npwx,raux,aux_hybrid1(:,jbnd),'Wave')
        !$acc end host_data
        !
     ENDDO
     !
     !$acc update host(aux_hybrid1)
     CALL mp_sum(aux_hybrid1, inter_bgrp_comm)
     !$acc update device(aux_hybrid1)
     !
     !$acc parallel loop collapse(2) present(hybrid_kd1,aux_hybrid1)
     DO lbnd = 1, nbnd_do
        DO ig = 1, npw
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           hybrid_kd1(ig,lbnd) = hybrid_kd1(ig,lbnd) - aux_hybrid1(ig,ibnd) * exxalfa
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  !$acc exit data delete(aux_hybrid1)
  DEALLOCATE(aux_hybrid1)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE hybrid_kernel_term2(current_spin, evc1, hybrid_kd2, sf)
  !-----------------------------------------------------------------------
  !
  ! \sum_{v'} (\int v_c a_{v'} \phi_{v}) \phi_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma
  USE mp_global,             ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,nbgrp,my_bgrp_id
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : kpt_pool,band_group
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  COMPLEX(DP), INTENT(INOUT) :: hybrid_kd2(npwx,band_group%nlocx)
  !
  ! Workspace
  !
  INTEGER :: lbnd, ibnd, ibndp, jbnd, jbndp, ir, ig
  INTEGER :: current_spin_ikq, ikq, nbndval, nbnd_do
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: aux_hybrid2(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: aux_hybrid2
#endif
  COMPLEX(DP), ALLOCATABLE :: caux(:), gaux(:), raux(:)
  !$acc declare device_resident(caux,gaux,raux)
  !
  IF(sf) CALL errore('hybrid_kernel_gamma_term2', 'spin-flip is not supported', 1)
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_hybrid2(npwx,nbndval0x-n_trunc_bands))
  !$acc enter data create(aux_hybrid2)
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     nbndval = nbnd_occ(ikq)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        !$acc kernels present(raux)
        raux(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        DO lbnd = 1, nbnd_do ! index to be summed
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           IF(ibndp > nbndval) CYCLE
           !
           ! product of evc1 and evc
           !
           !$acc host_data use_device(evc1)
           CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,lbnd,ikq),evc_work(:,jbndp),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
           ENDDO
           !$acc end parallel
           !
           ! Apply the bare Coulomb potential
           !
           !$acc host_data use_device(caux,gaux)
           CALL single_fwfft_gamma(dffts,npw,npwx,caux,gaux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(gaux,pot3D)
           DO ig = 1, npw
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(gaux,caux)
           CALL double_invfft_gamma(dffts,npw,npwx,gaux,evc_work(:,ibndp),caux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop present(raux)
           DO ir = 1, dffts_nnr
              raux(ir) = raux(ir)+psic(ir)
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        !$acc host_data use_device(raux,aux_hybrid2)
        CALL single_fwfft_gamma(dffts,npw,npwx,raux,aux_hybrid2(:,jbnd),'Wave')
        !$acc end host_data
        !
     ENDDO
     !
     !$acc update host(aux_hybrid2)
     CALL mp_sum(aux_hybrid2, inter_bgrp_comm)
     !$acc update device(aux_hybrid2)
     !
     !$acc parallel loop collapse(2) present(hybrid_kd2,aux_hybrid2)
     DO lbnd = 1, nbnd_do
        DO ig = 1, npw
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           hybrid_kd2(ig,lbnd) = hybrid_kd2(ig,lbnd) - aux_hybrid2(ig,ibnd) * exxalfa
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  !$acc exit data delete(aux_hybrid2)
  DEALLOCATE(aux_hybrid2)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE hybrid_kernel_term3(current_spin, evc1, hybrid_kd3, sf)
  !-----------------------------------------------------------------------
  !
  ! \sum_{v'} (\int v_c a_{v'} \phi_{v}) a_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma
  USE mp_global,             ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,nbgrp,my_bgrp_id
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : kpt_pool,band_group
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  COMPLEX(DP), INTENT(INOUT) :: hybrid_kd3(npwx,band_group%nlocx)
  !
  ! Workspace
  !
  INTEGER :: lbnd, ibnd, jbnd, jbndp, ir, ig, iks_do, nbnd_do
  INTEGER :: current_spin_ikq, ikq, nbndval, flnbndval
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: aux_hybrid3(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: aux_hybrid3
#endif
  COMPLEX(DP), ALLOCATABLE :: caux(:), gaux(:), raux(:)
  !$acc declare device_resident(caux,gaux,raux)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_hybrid3(npwx,nbndval0x-n_trunc_bands))
  !$acc enter data create(aux_hybrid3)
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(ikq)
     flnbndval = nbnd_occ(iks_do)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= flnbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        jbndp = jbnd + n_trunc_bands
        !
        !$acc kernels present(raux)
        raux(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        DO lbnd = 1, nbnd_do
           !
           ! product of evc1 and evc
           !
           !$acc host_data use_device(evc1)
           CALL double_invfft_gamma(dffts,npw,npwx,evc1(:,lbnd,ikq),evc_work(:,jbndp),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
           ENDDO
           !$acc end parallel
           !
           ! Apply the bare Coulomb potential
           !
           !$acc host_data use_device(caux,gaux)
           CALL single_fwfft_gamma(dffts,npw,npwx,caux,gaux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(gaux,pot3D)
           DO ig = 1, npw
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(gaux,evc1,caux)
           CALL double_invfft_gamma(dffts,npw,npwx,gaux,evc1(:,lbnd,ikq),caux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop present(raux)
           DO ir = 1, dffts_nnr
              raux(ir) = raux(ir)+psic(ir)
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        !$acc host_data use_device(raux,aux_hybrid3)
        CALL single_fwfft_gamma(dffts,npw,npwx,raux,aux_hybrid3(:,jbnd),'Wave')
        !$acc end host_data
        !
     ENDDO
     !
     !$acc update host(aux_hybrid3)
     CALL mp_sum(aux_hybrid3, inter_bgrp_comm)
     !$acc update device(aux_hybrid3)
     !
     ! Recompute nbnd_do for the current spin channel
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     !$acc parallel loop collapse(2) present(hybrid_kd3,aux_hybrid3)
     DO lbnd = 1, nbnd_do
        DO ig = 1, npw
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           hybrid_kd3(ig,lbnd) = hybrid_kd3(ig,lbnd) - aux_hybrid3(ig,ibnd) * exxalfa
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  !$acc exit data delete(aux_hybrid3)
  DEALLOCATE(aux_hybrid3)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux)
  !
END SUBROUTINE
!
!-----------------------------------------------------------------------
SUBROUTINE hybrid_kernel_term4(current_spin, evc1, hybrid_kd4, sf)
  !-----------------------------------------------------------------------
  !
  ! \sum_{v'} (\int v_c a_{v'} a_{v}) \phi_{v'}
  !
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : omega
  USE fft_base,              ONLY : dffts
  USE types_coulomb,         ONLY : pot3D
  USE mp,                    ONLY : mp_sum,mp_bcast
  USE fft_at_gamma,          ONLY : single_fwfft_gamma,double_invfft_gamma
  USE mp_global,             ONLY : inter_image_comm,my_image_id,inter_bgrp_comm,nbgrp,my_bgrp_id
  USE pwcom,                 ONLY : npw,npwx,isk,ngk
  USE westcom,               ONLY : nbnd_occ,iuwfc,lrwfc,nbndval0x,n_trunc_bands
  USE exx,                   ONLY : exxalfa
  USE buffers,               ONLY : get_buffer
  USE distribution_center,   ONLY : kpt_pool,band_group
#if defined(__CUDA)
  USE wavefunctions_gpum,    ONLY : using_evc,using_evc_d,evc_work=>evc_d,psic=>psic_d
  USE wavefunctions,         ONLY : evc_host=>evc
#else
  USE wavefunctions,         ONLY : evc_work=>evc,psic
#endif
  !
  IMPLICIT NONE
  !
  ! I/O
  !
  INTEGER, INTENT(IN) :: current_spin
  COMPLEX(DP), INTENT(IN) :: evc1(npwx,band_group%nlocx,kpt_pool%nloc)
  LOGICAL, INTENT(IN) :: sf
  COMPLEX(DP), INTENT(INOUT) :: hybrid_kd4(npwx,band_group%nlocx)
  !
  ! Workspace
  !
  INTEGER :: lbnd, ibnd, ibndp, jbnd, ir, ig, iks_do
  INTEGER :: current_spin_ikq, ikq, nbndval, nbnd_do
  INTEGER :: dffts_nnr
  COMPLEX(DP), ALLOCATABLE :: aux_hybrid4(:,:), aux_evc1(:,:)
#if defined(__CUDA)
  ATTRIBUTES(PINNED) :: aux_hybrid4, aux_evc1
#endif
  COMPLEX(DP), ALLOCATABLE :: caux(:), gaux(:), raux(:)
  !$acc declare device_resident(caux,gaux,raux)
  INTEGER, PARAMETER :: flks(2) = [2,1]
  !
  dffts_nnr = dffts%nnr
  !
  ALLOCATE(aux_evc1(npwx,nbndval0x-n_trunc_bands))
  ALLOCATE(aux_hybrid4(npwx,nbndval0x-n_trunc_bands))
  !$acc enter data create(aux_evc1,aux_hybrid4)
  ALLOCATE(caux(dffts%nnr))
  ALLOCATE(gaux(npwx))
  ALLOCATE(raux(dffts%nnr))
  !
  DO ikq = 1,kpt_pool%nloc
     !
     current_spin_ikq = isk(ikq)
     IF(current_spin_ikq /= current_spin) CYCLE
     !
     IF(sf) THEN
        iks_do = flks(ikq)
     ELSE
        iks_do = ikq
     ENDIF
     !
     nbndval = nbnd_occ(ikq)
     !
     nbnd_do = 0
     DO lbnd = 1,band_group%nloc
        ibnd = band_group%l2g(lbnd)+n_trunc_bands
        IF(ibnd > n_trunc_bands .AND. ibnd <= nbndval) nbnd_do = nbnd_do+1
     ENDDO
     !
     ! ... Number of G vectors for PW expansion of wfs at k
     !
     npw = ngk(ikq)
     !
     ! ... read in GS wavefunctions ikq
     !
     IF(kpt_pool%nloc > 1) THEN
#if defined(__CUDA)
        IF(my_image_id == 0) CALL get_buffer(evc_host,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_host,0,inter_image_comm)
        !
        CALL using_evc(2)
        CALL using_evc_d(0)
#else
        IF(my_image_id == 0) CALL get_buffer(evc_work,lrwfc,iuwfc,ikq)
        CALL mp_bcast(evc_work,0,inter_image_comm)
#endif
     ENDIF
     !
     ! Collect evc1 from all bgrps
     !
     !$acc kernels present(aux_evc1)
     aux_evc1(:,:) = (0._DP,0._DP)
     !$acc end kernels
     !
     !$acc parallel loop collapse(2) present(aux_evc1,evc1)
     DO lbnd = 1, nbnd_do
        DO ig = 1, npw
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           aux_evc1(ig,ibnd) = evc1(ig,lbnd,iks_do)
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
     !$acc update host(aux_evc1)
     CALL mp_sum(aux_evc1,inter_bgrp_comm)
     !$acc update device(aux_evc1)
     !
     DO jbnd = 1, nbndval - n_trunc_bands ! index to be left
        !
        !$acc kernels present(raux)
        raux(:) = (0._DP,0._DP)
        !$acc end kernels
        !
        DO lbnd = 1, nbnd_do ! index to be summed
           !
           ibnd = band_group%l2g(lbnd)
           ibndp = ibnd + n_trunc_bands
           !
           ! product of evc1 and evc1
           !
           !$acc host_data use_device(aux_evc1,evc1)
           CALL double_invfft_gamma(dffts,npw,npwx,aux_evc1(:,jbnd),evc1(:,lbnd,iks_do),psic,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              caux(ir) = CMPLX(REAL(psic(ir),KIND=DP)*AIMAG(psic(ir))/omega,KIND=DP)
           ENDDO
           !$acc end parallel
           !
           ! Apply the bare Coulomb potential
           !
           !$acc host_data use_device(caux,gaux)
           CALL single_fwfft_gamma(dffts,npw,npwx,caux,gaux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(gaux,pot3D)
           DO ig = 1, npw
              gaux(ig) = gaux(ig) * (pot3D%sqvc(ig)**2)
           ENDDO
           !$acc end parallel
           !
           !$acc host_data use_device(gaux,caux)
           CALL double_invfft_gamma(dffts,npw,npwx,gaux,evc_work(:,ibndp),caux,'Wave')
           !$acc end host_data
           !
           !$acc parallel loop present(caux)
           DO ir = 1, dffts_nnr
              psic(ir) = CMPLX(REAL(caux(ir),KIND=DP)*AIMAG(caux(ir)),KIND=DP)
           ENDDO
           !$acc end parallel
           !
           !$acc parallel loop present(raux)
           DO ir = 1, dffts_nnr
              raux(ir) = raux(ir)+psic(ir)
           ENDDO
           !$acc end parallel
           !
        ENDDO
        !
        !$acc host_data use_device(raux,aux_hybrid4)
        CALL single_fwfft_gamma(dffts,npw,npwx,raux,aux_hybrid4(:,jbnd),'Wave')
        !$acc end host_data
        !
     ENDDO
     !
     !$acc update host(aux_hybrid4)
     CALL mp_sum(aux_hybrid4, inter_bgrp_comm)
     !$acc update device(aux_hybrid4)
     !
     !$acc parallel loop collapse(2) present(hybrid_kd4,aux_hybrid4)
     DO lbnd = 1, nbnd_do
        DO ig = 1, npw
           !
           ! ibnd = band_group%l2g(lbnd)
           !
           ibnd = nbgrp*(lbnd-1)+my_bgrp_id+1
           !
           hybrid_kd4(ig,lbnd) = hybrid_kd4(ig,lbnd) - aux_hybrid4(ig,ibnd) * exxalfa
           !
        ENDDO
     ENDDO
     !$acc end parallel
     !
  ENDDO
  !
  !$acc exit data delete(aux_evc1,aux_hybrid4)
  DEALLOCATE(aux_evc1)
  DEALLOCATE(aux_hybrid4)
  DEALLOCATE(caux)
  DEALLOCATE(gaux)
  DEALLOCATE(raux)
  !
END SUBROUTINE
