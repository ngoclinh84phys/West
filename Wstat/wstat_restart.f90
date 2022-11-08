!
! Copyright (C) 2015-2022 M. Govoni
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file:
! Marco Govoni
!
!-----------------------------------------------------------------------
MODULE wstat_restart
  !----------------------------------------------------------------------------
  !
  USE kinds,       ONLY : DP
  USE json_module, ONLY : json_file
  !
  IMPLICIT NONE
  !
  INTERFACE wstat_restart_write
     MODULE PROCEDURE wstat_restart_write_real, wstat_restart_write_complex
  END INTERFACE
  !
  INTERFACE wstat_restart_read
     MODULE PROCEDURE wstat_restart_read_real, wstat_restart_read_complex
  END INTERFACE
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_write_real(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,my_pool_id,my_bgrp_id,me_bgrp,inter_image_comm,nimage
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : n_pdep_basis,ev,conv,dvg,dng,wstat_restart_dir
      USE mp,                   ONLY : mp_barrier,mp_get
      USE pdep_io,              ONLY : pdep_merge_and_write_G
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      REAL(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: local_j,global_j
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR
      !
      CALL my_mkdir(TRIM(wstat_restart_dir))
      !
      CALL start_clock('wstat_restart')
      time_spent(1) = get_clock('wstat_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add('dav_iter',dav_iter)
         CALL json%add('notcnv',notcnv)
         CALL json%add('nbase',nbase)
         CALL json%add('conv',conv(:))
         CALL json%add('ev',ev(:))
         CALL json%add('ew',ew(:))
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/'//TRIM('summary.json'))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! CREATE THE HR, VR FILE
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/hr_vr.dat',FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      ! CREATE THE EIGENVECTOR FILES
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         fname = TRIM(wstat_restart_dir)//'/V'//my_label//'.dat'
         CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
         fname = TRIM(wstat_restart_dir)//'/N'//my_label//'.dat'
         CALL pdep_merge_and_write_G(fname,dng(:,local_j))
         !
      ENDDO
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      time_spent(2) = get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(wstat_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_write_complex(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr,lastdone_iq)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : my_image_id,my_pool_id,my_bgrp_id,me_bgrp,inter_image_comm,nimage
      USE mp_world,             ONLY : mpime,root,world_comm
      USE io_global,            ONLY : stdout
      USE westcom,              ONLY : n_pdep_basis,ev,conv,dvg,dng,wstat_restart_dir
      USE mp,                   ONLY : mp_barrier,mp_get
      USE pdep_io,              ONLY : pdep_merge_and_write_G
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(IN) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(IN) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(IN) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(IN) :: vr_distr(n_pdep_basis,pert%nlocx)
      INTEGER,INTENT(IN),OPTIONAL :: lastdone_iq
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      CHARACTER(6) :: my_label
      INTEGER :: local_j,global_j
      INTEGER :: im
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      TYPE(json_file) :: json
      INTEGER :: iun
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! MKDIR
      !
      CALL my_mkdir(wstat_restart_dir)
      !
      CALL start_clock('wstat_restart')
      time_spent(1) = get_clock('wstat_restart')
      !
      ! CREATE THE SUMMARY FILE
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         !
         CALL json%add('dav_iter',dav_iter)
         CALL json%add('notcnv',notcnv)
         CALL json%add('nbase',nbase)
         CALL json%add('conv',conv(:))
         CALL json%add('ev',ev(:))
         CALL json%add('ew',ew(:))
         IF(PRESENT(lastdone_iq)) THEN
            CALL json%add('lastdone_iq',lastdone_iq)
         ENDIF
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/'//TRIM('summary.json'))
         CALL json%print(iun)
         CLOSE(iun)
         CALL json%destroy()
         !
      ENDIF
      !
      ! CREATE THE HR, VR FILE
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) THEN
         !
         OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/'//TRIM('hr_vr.dat'),FORM='unformatted')
         !
      ENDIF
      !
      DO im = 0,nimage-1
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,hr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(tmp_distr,vr_distr,my_image_id,0,im,im,inter_image_comm)
         IF(mpime == root) WRITE(iun) tmp_distr
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      ! CREATE THE EIGENVECTOR FILES
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         fname = TRIM(wstat_restart_dir)//'/V'//my_label//'.dat'
         IF(PRESENT(lastdone_iq)) THEN
            CALL pdep_merge_and_write_G(fname,dvg(:,local_j),lastdone_iq)
         ELSE
            CALL pdep_merge_and_write_G(fname,dvg(:,local_j))
         ENDIF
         fname = TRIM(wstat_restart_dir)//'/N'//my_label//'.dat'
         IF(PRESENT(lastdone_iq)) THEN
            CALL pdep_merge_and_write_G(fname,dng(:,local_j),lastdone_iq)
         ELSE
            CALL pdep_merge_and_write_G(fname,dng(:,local_j))
         ENDIF
         !
      ENDDO
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      time_spent(2) = get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART written in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location   : ",a)') TRIM(wstat_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_clear()
      !------------------------------------------------------------------------
      !
      USE mp_world,             ONLY : root,mpime,world_comm
      USE mp,                   ONLY : mp_barrier,mp_bcast
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE clib_wrappers,        ONLY : f_rmdir
      USE west_io,              ONLY : remove_if_present
      !
      IMPLICIT NONE
      !
      ! Workspace
      !
      CHARACTER(LEN=512) :: fname
      INTEGER :: ierr,ip
      CHARACTER(6) :: my_label
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      ! ... clear the main restart directory
      !
      IF(mpime == root) THEN
         CALL remove_if_present(TRIM(wstat_restart_dir)//'/summary.json')
         CALL remove_if_present(TRIM(wstat_restart_dir)//'/hr_vr.dat')
         DO ip = 1,n_pdep_basis
            WRITE(my_label,'(i6.6)') ip
            fname = 'V'//my_label//'.dat'
            CALL remove_if_present(TRIM(wstat_restart_dir)//'/'//TRIM(fname))
            fname = 'N'//my_label//'.dat'
            CALL remove_if_present(TRIM(wstat_restart_dir)//'/'//TRIM(fname))
         ENDDO
         ierr = f_rmdir(TRIM(wstat_restart_dir))
      ENDIF
      !
      CALL mp_bcast(ierr,root,world_comm)
      !
      CALL errore('wstat_restart','cannot clear restart',ierr)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_read_real(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : world_comm
      USE mp,                   ONLY : mp_barrier
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE io_global,            ONLY : stdout
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      REAL(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('wstat_restart')
      time_spent(1) = get_clock('wstat_restart')
      !
      CALL read_restart12_(dav_iter,notcnv,nbase,ew)
      !
      CALL read_restart3d_(hr_distr,vr_distr)
      !
      CALL read_restart4_(nbase)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      time_spent(2) = get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------")')
      WRITE(stdout,'(5x,"[I/O] RESTART read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
      WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wstat_restart_dir)
      WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------")')
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE wstat_restart_read_complex(dav_iter,notcnv,nbase,ew,hr_distr,vr_distr,lastdone_iq,iq)
      !------------------------------------------------------------------------
      !
      USE mp_global,            ONLY : world_comm
      USE mp,                   ONLY : mp_barrier
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE io_global,            ONLY : stdout
      USE distribution_center,  ONLY : pert
      USE types_bz_grid,        ONLY : q_grid
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      INTEGER,INTENT(OUT) :: lastdone_iq
      INTEGER,INTENT(IN) :: iq
      !
      ! Workspace
      !
      REAL(DP),EXTERNAL :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      INTEGER :: ipol
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('wstat_restart')
      time_spent(1) = get_clock('wstat_restart')
      !
      CALL read_restart12_(dav_iter,notcnv,nbase,ew,lastdone_iq)
      !
      CALL read_restart3z_(hr_distr,vr_distr)
      !
      CALL read_restart4_(nbase,lastdone_iq)
      !
      ! BARRIER
      !
      CALL mp_barrier(world_comm)
      !
      time_spent(2) = get_clock('wstat_restart')
      CALL stop_clock('wstat_restart')
      !
      IF(iq == lastdone_iq) THEN
         WRITE(stdout,'(1/,5x,"[I/O] -------------------------------------------------------------------")')
         WRITE(stdout,'(5x,"[I/O] Restarting from q(",i5,") = (",3f12.7,")")') &
              lastdone_iq,(q_grid%p_cryst(ipol,lastdone_iq),ipol=1,3)
         WRITE(stdout,'(5x,"[I/O] RESTART read in ",a20)') human_readable_time(time_spent(2)-time_spent(1))
         WRITE(stdout,'(5x,"[I/O] In location : ",a)') TRIM(wstat_restart_dir)
         WRITE(stdout,'(5x,"[I/O] -------------------------------------------------------------------")')
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart12_(dav_iter,notcnv,nbase,ew,iq)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : conv,n_pdep_eigen,n_pdep_basis,wstat_restart_dir,ev
      USE mp_world,             ONLY : world_comm,mpime,root
      USE mp,                   ONLY : mp_bcast
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER,INTENT(OUT) :: dav_iter,notcnv,nbase
      INTEGER,INTENT(OUT),OPTIONAL :: iq
      REAL(DP),INTENT(OUT) :: ew(n_pdep_basis)
      !
      ! Workspace
      !
      LOGICAL :: found
      TYPE(json_file) :: json
      REAL(DP),ALLOCATABLE :: rvals(:)
      LOGICAL,ALLOCATABLE :: lvals(:)
      INTEGER :: ival
      !
      IF(mpime == root) THEN
         !
         CALL json%initialize()
         CALL json%load(filename=TRIM(wstat_restart_dir)//'/'//TRIM('summary.json'))
         !
         CALL json%get('dav_iter',ival,found)
         IF(found) dav_iter = ival
         CALL json%get('notcnv',ival,found)
         IF(found) notcnv = ival
         CALL json%get('nbase',ival,found)
         IF(found) nbase = ival
         CALL json%get('conv',lvals,found)
         IF(found) conv(1:n_pdep_eigen) = lvals(1:n_pdep_eigen)
         CALL json%get('ev',rvals,found)
         IF(found) ev(:) = rvals(:)
         CALL json%get('ew',rvals,found)
         IF(found) ew(1:n_pdep_basis) = rvals(1:n_pdep_basis)
         IF(PRESENT(iq)) THEN
            CALL json%get('lastdone_iq',ival,found)
            IF(found) iq = ival
         ENDIF
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast(dav_iter,root,world_comm)
      CALL mp_bcast(notcnv,root,world_comm)
      CALL mp_bcast(nbase,root,world_comm)
      CALL mp_bcast(conv,root,world_comm)
      !
      CALL mp_bcast(ev,root,world_comm)
      CALL mp_bcast(ew,root,world_comm)
      IF(PRESENT(iq)) THEN
         CALL mp_bcast(iq,root,world_comm)
      ENDIF
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3d_(hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE mp_world,             ONLY : mpime,root
      USE mp,                   ONLY : mp_bcast,mp_get
      USE distribution_center,  ONLY : pert
      USE mp_global,            ONLY : nimage,my_pool_id,my_bgrp_id,me_bgrp,inter_image_comm,intra_image_comm,my_image_id
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      REAL(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      REAL(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: iun
      INTEGER :: im
      REAL(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/'//TRIM('hr_vr.dat'),FORM='unformatted')
      !
      DO im = 0,nimage-1
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      CALL mp_bcast(hr_distr,0,intra_image_comm)
      CALL mp_bcast(vr_distr,0,intra_image_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart3z_(hr_distr,vr_distr)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : n_pdep_basis,wstat_restart_dir
      USE mp_world,             ONLY : mpime,root
      USE mp,                   ONLY : mp_bcast,mp_get
      USE distribution_center,  ONLY : pert
      USE mp_global,            ONLY : nimage,my_pool_id,my_bgrp_id,me_bgrp,inter_image_comm,intra_image_comm,my_image_id
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      COMPLEX(DP),INTENT(OUT) :: hr_distr(n_pdep_basis,pert%nlocx)
      COMPLEX(DP),INTENT(OUT) :: vr_distr(n_pdep_basis,pert%nlocx)
      !
      ! Workspace
      !
      INTEGER :: iun
      INTEGER :: im
      COMPLEX(DP),ALLOCATABLE :: tmp_distr(:,:)
      !
      ALLOCATE(tmp_distr(n_pdep_basis,pert%nlocx))
      !
      IF(mpime == root) OPEN(NEWUNIT=iun,FILE=TRIM(wstat_restart_dir)//'/'//TRIM('hr_vr.dat'),FORM='unformatted')
      !
      DO im = 0,nimage-1
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(hr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
         IF(mpime == root) READ(iun) tmp_distr(:,:)
         IF(me_bgrp == 0 .AND. my_bgrp_id == 0 .AND. my_pool_id == 0) &
         & CALL mp_get(vr_distr,tmp_distr,my_image_id,im,0,im,inter_image_comm)
         !
      ENDDO
      !
      IF(mpime == root) CLOSE(iun)
      !
      DEALLOCATE(tmp_distr)
      !
      CALL mp_bcast(hr_distr,0,intra_image_comm)
      CALL mp_bcast(vr_distr,0,intra_image_comm)
      !
    END SUBROUTINE
    !
    !------------------------------------------------------------------------
    SUBROUTINE read_restart4_(nbase,iq)
      !------------------------------------------------------------------------
      !
      USE westcom,              ONLY : dvg,dng,npwqx,wstat_restart_dir
      USE pdep_io,              ONLY : pdep_read_G_and_distribute
      USE distribution_center,  ONLY : pert
      !
      IMPLICIT NONE
      !
      INTEGER,INTENT(IN) :: nbase
      INTEGER,INTENT(IN),OPTIONAL :: iq
      !
      INTEGER :: global_j,local_j
      CHARACTER(6) :: my_label
      CHARACTER(LEN=512) :: fname
      !
      IF(.NOT. ALLOCATED(dvg)) ALLOCATE(dvg(npwqx,pert%nlocx))
      IF(.NOT. ALLOCATED(dng)) ALLOCATE(dng(npwqx,pert%nlocx))
      dvg = 0._DP
      dng = 0._DP
      !
      DO local_j = 1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         WRITE(my_label,'(i6.6)') global_j
         IF(global_j > nbase) CYCLE
         !
         fname = TRIM(wstat_restart_dir)//'/V'//my_label//'.dat'
         IF(PRESENT(iq)) THEN
            CALL pdep_read_G_and_distribute(fname,dvg(:,local_j),iq)
         ELSE
            CALL pdep_read_G_and_distribute(fname,dvg(:,local_j))
         ENDIF
         fname = TRIM(wstat_restart_dir)//'/N'//my_label//'.dat'
         IF(PRESENT(iq)) THEN
            CALL pdep_read_G_and_distribute(fname,dng(:,local_j),iq)
         ELSE
            CALL pdep_read_G_and_distribute(fname,dng(:,local_j))
         ENDIF
         !
      ENDDO
      !
    END SUBROUTINE
    !
END MODULE
