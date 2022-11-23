!
! Copyright (C) 2001-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE wbse_lanczos_diago()
  !---------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE lsda_mod,             ONLY : nspin
  USE pwcom,                ONLY : npw,npwx,ngk,nks,isk,current_spin
  USE westcom,              ONLY : nbnd_occ,lrwfc,iuwfc,nbnd_occ,wbse_calculation,d0psi,ipol_input,&
                                 & n_lanczos,beta_store,zeta_store,nbndval0x,l_bse_calculation,&
                                 & size_index_matrix_lz,n_steps_write_restart
  USE lanczos_db,           ONLY : lanczos_d0psi_read,lanczos_d0psi_write,lanczos_evcs_write,&
                                 & lanczos_evcs_read
  USE lanczos_restart,      ONLY : lanczos_restart_write,lanczos_restart_read,lanczos_postpro_write
  USE mp_global,            ONLY : my_image_id,inter_image_comm
  USE mp,                   ONLY : mp_bcast,mp_barrier
  USE wavefunctions,        ONLY : evc
  USE buffers,              ONLY : get_buffer
  USE distribution_center,  ONLY : aband,bandpair
  USE class_idistribute,    ONLY : idistribute
  USE io_push,              ONLY : io_push_title
  USE bar,                  ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  !
  IMPLICIT NONE
  !
  ! Local variables
  !
  LOGICAL :: l_from_scratch
  INTEGER :: ip,iip,pol_index,nipol_input
  INTEGER :: iter
  INTEGER :: iks,is,nbndval
  INTEGER :: ilan_restart,ilan_stopped,ipol_restart,ipol_stopped
  INTEGER :: size_index_matrix
  INTEGER, PARAMETER :: n_ipol = 3
  INTEGER, ALLOCATABLE :: pol_index_input(:)
  CHARACTER(LEN=3), ALLOCATABLE :: pol_label_input(:)
  REAL(DP) :: beta(nspin),gamma(nspin)
  COMPLEX(DP) :: zeta(nspin),wbse_dot_out(nspin)
  COMPLEX(DP), ALLOCATABLE :: evc1(:,:,:),evc1_old(:,:,:),evc1_new(:,:,:)
  TYPE(bar_type) :: barra
  !
  ! ... DISTRIBUTE lanczos
  !
  aband = idistribute()
  !
  CALL aband%init(nbndval0x,'i','nbndval',.TRUE.)
  !
  ! ... DISTRIBUTE bse_kernel
  !
  size_index_matrix = MAXVAL(size_index_matrix_lz)
  !
  IF(l_bse_calculation) THEN
     bandpair = idistribute()
     CALL bandpair%init(size_index_matrix,'i','n_pairs',.TRUE.)
  ENDIF
  !
  ! Main Lanzcos program
  !
  IF(n_lanczos < 1) CALL errore('wbse_lanczos_diago','n_lanczos must be > 0',1)
  !
  SELECT CASE(ipol_input)
  CASE('XX','xx')
     nipol_input = 1
     ALLOCATE(pol_index_input(1))
     ALLOCATE(pol_label_input(1))
     pol_index_input(1) = 1
     pol_label_input(1) = 'XX'
  CASE('YY','yy')
     nipol_input = 1
     ALLOCATE(pol_index_input(1))
     ALLOCATE(pol_label_input(1))
     pol_index_input(1) = 2
     pol_label_input(1) = 'YY'
  CASE('ZZ','zz')
     nipol_input = 1
     ALLOCATE(pol_index_input(1))
     ALLOCATE(pol_label_input(1))
     pol_index_input(1) = 3
     pol_label_input(1) = 'ZZ'
  CASE('XYZ','xyz')
     nipol_input = 3
     ALLOCATE(pol_index_input(3))
     ALLOCATE(pol_label_input(3))
     pol_index_input(1) = 1
     pol_label_input(1) = 'XX'
     pol_index_input(2) = 2
     pol_label_input(2) = 'YY'
     pol_index_input(3) = 3
     pol_label_input(3) = 'ZZ'
  CASE DEFAULT
     CALL errore('wbse_lanczos_diago','wrong ipol_input',1)
  END SELECT
  !
  ALLOCATE(d0psi(npwx,nbndval0x,nks,n_ipol))
  ALLOCATE(beta_store(nipol_input,n_lanczos,nspin))
  ALLOCATE(zeta_store(nipol_input,n_ipol,n_lanczos,nspin))
  !
  beta_store(:,:,:) = 0._DP
  zeta_store(:,:,:,:) = (0._DP,0._DP)
  !
  ALLOCATE(evc1(npwx,nbndval0x,nks))
  ALLOCATE(evc1_old(npwx,nbndval0x,nks))
  ALLOCATE(evc1_new(npwx,nbndval0x,nks))
  !
  SELECT CASE(wbse_calculation)
  CASE('l')
     !
     ! RESTART
     !
     CALL lanczos_restart_read(nipol_input,ipol_stopped,ilan_stopped)
     !
     ! 1) read ipol_stopped
     ! 2) read ilan_stopped
     ! 3) read store_beta, store_zeta
     !
     ipol_restart = ipol_stopped
     ilan_restart = ilan_stopped+1
     !
     ! 4) read d0psi evc1, evc_old, saved on files
     !
     CALL lanczos_d0psi_read()
     CALL lanczos_evcs_read(evc1,evc1_old)
     !
     l_from_scratch = .FALSE.
     !
  CASE('L')
     !
     ! FROM SCRATCH
     !
     CALL solve_e_psi()
     !
     CALL lanczos_d0psi_write()
     !
     ipol_restart = 1
     ilan_restart = 1
     !
     l_from_scratch = .TRUE.
     !
  CASE DEFAULT
     CALL errore('wbse_lanczos_diago','wrong wlzcos_calculation',1)
  END SELECT
  !
  WRITE(stdout,'(/,5x,"Lanczos linear-response absorption spectrum calculation")')
  WRITE(stdout,'(/,5x,"using Tamm-Dancoff Liouvillian operator")')
  WRITE(stdout,'(/,5x,"number of Lanczos iterations = ",i6)') n_lanczos
  !
  polarization_loop : DO ip = ipol_restart,nipol_input
     !
     pol_index = pol_index_input(ip)
     !
     IF(l_from_scratch) THEN
        CALL io_push_title('Starting new Lanczos loop at ipol: '//TRIM(pol_label_input(ip)))
        !
        evc1_old(:,:,:) = (0._DP,0._DP)
        evc1_new(:,:,:) = (0._DP,0._DP)
        evc1(:,:,:) = d0psi(:,:,:,pol_index)
     ELSE
        CALL io_push_title('Retarting Lanczos loop at ipol: '//TRIM(pol_label_input(ip)))
     ENDIF
     !
     CALL start_bar_type(barra,'lan_diago',n_lanczos-ilan_restart+1)
     !
     ! Loop on the Lanczos iterations
     !
     lancz_loop : DO iter = ilan_restart,n_lanczos
        !
        ! Application of the Liouvillian superoperator
        !
        CALL west_apply_liouvillian(evc1,evc1_new)
        !
        ! By construction <p|Lq>=0 should be 0, forcing this both conserves
        ! resources and increases stability.
        ! ( i.e., alpha = 0 )
        !
        ! Orthogonality requirement: <v|\bar{L}|v> = 1
        !
        CALL wbse_dot(evc1,evc1_new,nbndval0x,nks,wbse_dot_out)
        !
        beta(:) = REAL(wbse_dot_out,KIND=DP)
        !
        ! beta<0 is a serious error for the pseudo-Hermitian algorithm
        !
        DO is = 1,nspin
           IF(beta(is) < 0._DP) THEN
              CALL errore('wbse_lanczos_diago','negative beta',1)
           ELSE
              beta(is) = SQRT(beta(is))
              gamma(is) = beta(is)
           ENDIF
        ENDDO
        !
        beta_store(ip,iter,:) = beta
        !
        ! Renormalize q(i) and Lq(i)
        !
        DO iks = 1,nks
           current_spin = isk(iks)
           evc1(:,:,iks) = (1._DP/beta(current_spin))*evc1(:,:,iks)
           evc1_new(:,:,iks) = (1._DP/beta(current_spin))*evc1_new(:,:,iks)
        ENDDO
        !
        ! Calculation of zeta coefficients.
        ! See Eq.(35) in Malcioglu et al., Comput. Phys. Commun. 182, 1744 (2011).
        !
        IF(MOD(iter,2) == 0) THEN
           DO iip = 1,n_ipol
              CALL wbse_dot(d0psi(:,:,:,iip),evc1,nbndval0x,nks,wbse_dot_out)
              !
              zeta(:) = wbse_dot_out
              zeta_store(ip,iip,iter,:) = zeta
           ENDDO
        ELSE
           DO iip = 1,n_ipol
              zeta(:) = (0._DP,0._DP)
              zeta_store(ip,iip,iter,:) = zeta
           ENDDO
        ENDIF
        !
        DO iks = 1,nks
           current_spin = isk(iks)
           evc1_new(:,:,iks) = evc1_new(:,:,iks) - gamma(current_spin)*evc1_old(:,:,iks)
        ENDDO
        !
        ! Apply P_c|evc1_new>
        !
        DO iks = 1,nks
           !
           nbndval = nbnd_occ(iks)
           !
           ! ... Number of G vectors for PW expansion of wfs at k
           !
           npw = ngk(iks)
           !
           ! ... Read GS wavefunctions
           !
           IF(nks > 1) THEN
              IF(my_image_id == 0) CALL get_buffer(evc,lrwfc,iuwfc,iks)
              CALL mp_bcast(evc,0,inter_image_comm)
           ENDIF
           !
           CALL apply_alpha_pc_to_m_wfcs(nbndval,nbndval,evc1_new(1,1,iks),(1._DP,0._DP))
        ENDDO
        !
        ! Throw away q(i-1),and make q(i+1) to be the current vector,
        ! be ready for the next iteration. evc1_new will be free again after this step
        !
        evc1_old(:,:,:) = evc1
        evc1(:,:,:) = evc1_new
        !
        IF(n_steps_write_restart > 0 .AND. MOD(iter,n_steps_write_restart) == 0) THEN
           CALL lanczos_restart_write(nipol_input,ip,iter)
           CALL lanczos_evcs_write(evc1,evc1_old)
        ENDIF
        !
        CALL update_bar_type(barra,'lan_diago',1)
        !
     ENDDO lancz_loop
     !
     CALL lanczos_postpro_write(nipol_input,ip,pol_label_input(ip))
     !
     ilan_restart = 1
     l_from_scratch = .TRUE.
     !
     CALL stop_bar_type(barra,'lan_diago')
     !
  ENDDO polarization_loop
  !
  DEALLOCATE(pol_index_input)
  DEALLOCATE(pol_label_input)
  DEALLOCATE(d0psi)
  DEALLOCATE(evc1)
  DEALLOCATE(evc1_new)
  DEALLOCATE(evc1_old)
  DEALLOCATE(beta_store)
  DEALLOCATE(zeta_store)
  !
END SUBROUTINE
