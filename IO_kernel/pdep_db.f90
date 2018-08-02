!
! Copyright (C) 2015-2017 M. Govoni 
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
MODULE pdep_db
  !----------------------------------------------------------------------------
  !
  USE kinds,     ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  CONTAINS
    !
    !
    ! *****************************
    ! PDEP WRITE
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE pdep_db_write( iq, lprintinfo )
      !------------------------------------------------------------------------
      !
      USE mp,                   ONLY : mp_bcast,mp_barrier
      USE mp_world,             ONLY : mpime,root,world_comm
      USE mp_global,            ONLY : my_image_id
      USE io_global,            ONLY : stdout 
      USE westcom,              ONLY : wstat_calculation,n_pdep_times,n_pdep_eigen,n_pdep_maxiter,n_dfpt_maxiter, &
                                     & n_steps_write_restart,n_pdep_restart_from_itr,n_pdep_read_from_file,trev_pdep, &
                                     & tr2_dfpt,l_deflate,l_kinetic_only,ev,dvg,west_prefix,trev_pdep_rel, &
                                     & l_minimize_exx_if_active,l_use_ecutrho,wstat_save_dir,logfile 
      USE pdep_io,              ONLY : pdep_merge_and_write_G 
      USE io_push,              ONLY : io_push_bar
      USE distribution_center,  ONLY : pert
      USE types_bz_grid,        ONLY : q_grid 
      USE json_module,          ONLY : json_file, json_value, json_core
      USE cell_base,            ONLY : celldm,at,bg,tpiba
      USE gvect,                ONLY : ecutrho
      USE gvecw,                ONLY : ecutwfc
      !
      !
      IMPLICIT NONE
      !
      ! I/O 
      !
      INTEGER, INTENT(IN), OPTIONAL :: iq
      LOGICAL, INTENT(OUT), OPTIONAL :: lprintinfo
      !
      ! Workspace
      ! 
      ! optional 
      INTEGER, PARAMETER :: default_iq = 1 
      INTEGER            :: iq_
      LOGICAL, PARAMETER :: default_lprintinfo = .TRUE.
      LOGICAL            :: lprintinfo_
      ! labels  
      CHARACTER(LEN=9) :: label_j
      CHARACTER(LEN=9) :: label_q
      CHARACTER(LEN=9) :: label_i
      ! time 
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      ! scratch  
      INTEGER :: iunout,global_j,local_j
      INTEGER :: ierr
      ! json
      TYPE(json_core) :: jcor
      TYPE(json_file) :: json
      TYPE(json_value),POINTER :: jval
      INTEGER :: iunit, n_elements, ielement, myiq, write_element 
      LOGICAL :: found 
      ! files
      CHARACTER(LEN=:),ALLOCATABLE :: summary_file
      CHARACTER(LEN=:),ALLOCATABLE  :: eigenpot_filename(:)
      CHARACTER(LEN=:),ALLOCATABLE  :: fname
      LOGICAL :: lexists
      !
      ! Assign defaut to optional parameters 
      !
      IF(PRESENT(iq)) THEN 
         iq_ = iq 
      ELSE 
         iq_ = default_iq
      ENDIF
      IF(PRESENT(lprintinfo)) THEN 
         lprintinfo_ = lprintinfo
      ELSE 
         lprintinfo_ = default_lprintinfo
      ENDIF
      !
      ! MPI barrier
      !
      CALL mp_barrier(world_comm)
      !
      ! Start clock
      !
      CALL start_clock('pdep_db')
      time_spent(1)=get_clock('pdep_db')
      !
      ! Set filenames
      !
      IF(ALLOCATED(eigenpot_filename)) DEALLOCATE(eigenpot_filename)
      ALLOCATE( CHARACTER(LEN=25) :: eigenpot_filename(n_pdep_eigen) )
      DO global_j = 1, n_pdep_eigen
         WRITE(label_j,'(i9.9)') global_j
         WRITE(label_q,'(i9.9)') iq_
         eigenpot_filename(global_j) = "Q"//TRIM(ADJUSTL(label_q))//"E"//TRIM(ADJUSTL(label_j))//".json" 
      ENDDO
      IF(ALLOCATED(summary_file)) DEALLOCATE(summary_file)
      summary_file = TRIM(ADJUSTL(wstat_save_dir)) // "/summary.json"
      !
      ! Create summary file if it does not exist
      !
      IF ( mpime == root ) THEN 
        !
        INQUIRE(FILE=summary_file, EXIST=lexists)
        IF( (.NOT. lexists) ) THEN
          CALL json%initialize()
          CALL json%add('dielectric_matrix.domain.a1',celldm(1)*at(1:3,1))
          CALL json%add('dielectric_matrix.domain.a2',celldm(1)*at(1:3,2))
          CALL json%add('dielectric_matrix.domain.a3',celldm(1)*at(1:3,3))
          CALL json%add('dielectric_matrix.domain.b1',tpiba*bg(1:3,1))
          CALL json%add('dielectric_matrix.domain.b2',tpiba*bg(1:3,2))
          CALL json%add('dielectric_matrix.domain.b3',tpiba*bg(1:3,3))
          CALL jcor%create_array(jval,"pdep")
          CALL json%add("dielectric_matrix.pdep",jval)
          !
          OPEN(NEWUNIT= iunit, FILE= summary_file )
          CALL json%print_file( iunit )
          CLOSE( iunit )
          !
          CALL json%destroy()
        ENDIF
        !
      ENDIF 
      !
      ! Update summary file with current structure 
      !
      IF ( mpime == root ) THEN
         !
         CALL json%initialize()
         CALL json%load_file( filename = summary_file )
         !
         CALL json%info('dielectric_matrix.pdep',n_children=n_elements)
         write_element = n_elements + 1
         DO ielement = 1, n_elements
            WRITE(label_i,'(i9)') ielement
            CALL json%get('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').iq', myiq, found) 
            IF( found ) THEN  
               IF (myiq /= iq_ ) CYCLE 
               write_element = ielement
               EXIT
            ENDIF
         ENDDO
         WRITE(label_i,'(i9)') write_element
         CALL json%add('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').iq', iq_)
         CALL json%add('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').q', q_grid%p_cryst(1:3,iq_) )
         CALL json%add('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').eigenval' , ev(1:n_pdep_eigen) )
         CALL json%add('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').eigenvec' , eigenpot_filename(1:n_pdep_eigen))
         !
         OPEN( NEWUNIT=iunit, FILE=summary_file )
         CALL json%print_file( iunit )
         CLOSE( iunit )
         CALL json%destroy()
         !
      ENDIF
      !
      ! Dump eigenvectors
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         IF(global_j>n_pdep_eigen) CYCLE
         ! 
         fname = TRIM(ADJUSTL(wstat_save_dir)) // "/"// TRIM(ADJUSTL(eigenpot_filename(global_j))) 
         CALL pdep_merge_and_write_G(fname,dvg(:,local_j),iq_)
         !
      ENDDO
      !
      ! MPI barrier
      !
      CALL mp_barrier( world_comm )
      !
      ! timing
      !
      time_spent(2)=get_clock('pdep_db')
      CALL stop_clock('pdep_db')
      !
      IF (lprintinfo_) THEN
         WRITE(stdout,'(  5x," ")')
         CALL io_push_bar()
         WRITE(stdout, "(5x, 'SAVE written in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
         WRITE(stdout, "(5x, 'In location : ',a)") TRIM(ADJUSTL( wstat_save_dir ))  
         CALL io_push_bar()
      ENDIF 
      !
      IF(ALLOCATED(eigenpot_filename)) DEALLOCATE(eigenpot_filename)
      IF(ALLOCATED(summary_file)) DEALLOCATE(summary_file)
      IF(ALLOCATED(fname)) DEALLOCATE(fname)
      !
    END SUBROUTINE
    !
    !
    ! *****************************
    ! PDEP READ
    ! *****************************
    !
    !------------------------------------------------------------------------
    SUBROUTINE pdep_db_read( nglob_to_be_read, iq, lprintinfo )
      !------------------------------------------------------------------------
      !
      USE westcom,             ONLY : n_pdep_eigen,ev,dvg,west_prefix,npwqx,wstat_save_dir
      USE io_global,           ONLY : stdout 
      USE mp,                  ONLY : mp_bcast,mp_barrier
      USE mp_world,            ONLY : world_comm,mpime,root
      USE mp_global,           ONLY : my_image_id
      USE pdep_io,             ONLY : pdep_read_G_and_distribute
      USE io_push,             ONLY : io_push_bar
      USE distribution_center, ONLY : pert
      USE json_module,         ONLY : json_file, json_value, json_core
      !
      IMPLICIT NONE
      !
      ! I/O
      !
      INTEGER, INTENT(IN) :: nglob_to_be_read
      INTEGER, INTENT(IN), OPTIONAL :: iq
      LOGICAL, INTENT(IN), OPTIONAL :: lprintinfo
      ! 
      ! Workspace
      !
      ! optional 
      INTEGER, PARAMETER :: default_iq = 1 
      INTEGER            :: iq_  
      LOGICAL, PARAMETER :: default_lprintinfo = .TRUE.  
      LOGICAL            :: lprintinfo_
      ! labels
      CHARACTER(LEN=9)      :: label_j
      CHARACTER(LEN=9)      :: label_q
      CHARACTER(LEN=9)      :: label_i
      ! time 
      REAL(DP), EXTERNAL    :: GET_CLOCK
      REAL(DP) :: time_spent(2)
      CHARACTER(20),EXTERNAL :: human_readable_time
      ! scratch 
      INTEGER :: ierr, n_eigen_to_get
      INTEGER :: tmp_n_pdep_eigen
      INTEGER :: dime, iun, global_j, local_j
      REAL(DP),ALLOCATABLE :: tmp_ev(:)
      ! json managers
      TYPE(json_core) :: jcor
      TYPE(json_file) :: json
      TYPE(json_value),POINTER :: jval
      INTEGER :: iunit, n_elements, ielement, myiq 
      LOGICAL :: found
      INTEGER,ALLOCATABLE :: ilen(:)
      ! files
      CHARACTER(LEN=:),ALLOCATABLE :: eigenpot_filename(:) 
      CHARACTER(LEN=:),ALLOCATABLE :: fname
      !
      ! Assign defaut to optional parameters 
      !
      IF(PRESENT(iq)) THEN 
         iq_ = iq 
      ELSE 
         iq_ = default_iq
      ENDIF
      IF(PRESENT(lprintinfo)) THEN 
         lprintinfo_ = lprintinfo
      ELSE 
         lprintinfo_ = default_lprintinfo
      ENDIF
      !
      ! MPI barrier 
      !
      CALL mp_barrier(world_comm)
      !
      CALL start_clock('pdep_db')
      !
      ! TIMING
      !
      time_spent(1)=get_clock('pdep_db')
      !
      ! 1)  READ THE INPUT FILE
      !
      !
      IF ( mpime == root ) THEN
         !
         CALL json%initialize()
         CALL json%load_file( filename = TRIM(ADJUSTL(wstat_save_dir)) // "/summary.json" )
         IF( json%failed() ) THEN 
            CALL errore("", "Cannot open: " // TRIM(ADJUSTL(wstat_save_dir)) // "/summary.json", 1 )
         ENDIF 
         ! 
         !CALL json%get('dielectric_matrix.n_pdep_eigen', tmp_n_pdep_eigen, found)
         CALL json%info('dielectric_matrix.pdep',n_children=n_elements)
         !
         DO ielement = 1, n_elements
            WRITE(label_i,'(i9)') ielement
            CALL json%get('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').iq', myiq, found) 
            IF( found ) THEN
               IF( myiq /= iq_ ) CYCLE 
               CALL json%get('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').eigenval',tmp_ev)
               CALL json%get('dielectric_matrix.pdep('//TRIM(ADJUSTL(label_i))//').eigenvec',eigenpot_filename,ilen=ilen)
               tmp_n_pdep_eigen = SIZE(tmp_ev,1)
               EXIT
            ENDIF
         ENDDO
         !
         CALL json%destroy()
         !
      ENDIF
      !
      CALL mp_bcast( tmp_n_pdep_eigen, root, world_comm )
      !
      ! In case nglob_to_be_read is 0, overwrite it with the read value 
      !
      IF (nglob_to_be_read==0) THEN 
         n_eigen_to_get = tmp_n_pdep_eigen
         n_pdep_eigen=tmp_n_pdep_eigen
      ELSE
         n_eigen_to_get = MIN(tmp_n_pdep_eigen,nglob_to_be_read)
      ENDIF
      !
      ! 2)  READ THE EIGENVALUES FILE
      !
      IF(.NOT.ALLOCATED(ev)) ALLOCATE(ev(n_eigen_to_get))
      IF ( mpime==root ) ev(1:nglob_to_be_read) = tmp_ev(1:nglob_to_be_read)
      CALL mp_bcast( ev, root, world_comm )
      !
      IF(.NOT.ALLOCATED(eigenpot_filename)) ALLOCATE( CHARACTER(LEN=25) :: eigenpot_filename(n_eigen_to_get) )
      DO ielement = 1, n_eigen_to_get
         CALL mp_bcast(eigenpot_filename(ielement),root,world_comm)
      ENDDO
      !
      ! 3)  READ THE EIGENVECTOR FILES
      !
      IF(.NOT.ALLOCATED(dvg)) THEN
         ALLOCATE(dvg(npwqx,pert%nlocx))
         dvg = 0._DP
      ENDIF
      !
      DO local_j=1,pert%nloc
         !
         ! local -> global
         !
         global_j = pert%l2g(local_j)
         IF(global_j>n_eigen_to_get) CYCLE
         !
         fname = TRIM(ADJUSTL(wstat_save_dir)) // "/"// TRIM(ADJUSTL(eigenpot_filename(global_j)))
         CALL pdep_read_G_and_distribute(fname,dvg(:,local_j),iq_)
         !
      ENDDO
      !
      ! MPI BARRIER
      !
      CALL mp_barrier( world_comm )
      !
      ! TIMING
      !
      time_spent(2)=get_clock('pdep_db')
      CALL stop_clock('pdep_db')
      !
      IF (lprintinfo_) THEN
         WRITE(stdout,'(  5x," ")')
         CALL io_push_bar()
         WRITE(stdout, "(5x, 'SAVE read in ',a20)") human_readable_time(time_spent(2)-time_spent(1)) 
         WRITE(stdout, "(5x, 'In location : ',a)") TRIM(ADJUSTL( wstat_save_dir )) 
         WRITE(stdout, "(5x, 'Eigen. found : ',i12)") n_eigen_to_get
         CALL io_push_bar()
      ENDIF
      !
      IF(ALLOCATED(eigenpot_filename)) DEALLOCATE(eigenpot_filename)
      IF(ALLOCATED(fname)) DEALLOCATE(fname)
      !
    END SUBROUTINE
    !
END MODULE
