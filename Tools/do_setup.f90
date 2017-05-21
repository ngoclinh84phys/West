!
! Copyright (C) 2015-2016 M. Govoni 
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
SUBROUTINE do_setup
  !-----------------------------------------------------------------------
  !
  USE json_module,            ONLY : json_file
  USE pwcom,                  ONLY : npw,nbnd,nkstot,xk,wk,nspin,nelec,nelup,neldw,et,wg,&
                                   & lspinorb,domag,lsda,isk,nks,two_fermi_energies,ngk
  USE fixed_occ,              ONLY : tfixed_occ,f_inp
  USE kinds,                  ONLY : DP
  USE mp,                     ONLY : mp_sum
  USE mp_global,              ONLY : intra_bgrp_comm,npool,nbgrp
  USE mp_pools,               ONLY : intra_pool_comm, inter_pool_comm, &
                                     my_pool_id, nproc_pool, kunit
  USE io_global,              ONLY : stdout
  USE lsda_mod,               ONLY : current_spin,lsda
  USE constants,              ONLY : rytoev
  USE control_flags,          ONLY : gamma_only
  USE noncollin_module,       ONLY : noncolin,npol
  USE cell_base,              ONLY : omega,celldm,at
  USE fft_base,               ONLY : dfftp,dffts
  USE gvecs,                  ONLY : ngms_g, ngms
  USE gvect,                  ONLY : ngm_g, ngm
  USE gvecw,                  ONLY : ecutwfc
  USE io_push
  USE westcom,                ONLY : logfile  
  USE mp_world,               ONLY : mpime, root
  !
  IMPLICIT NONE
  !
  TYPE(json_file) :: json
  INTEGER :: iunit
  INTEGER :: auxi,ib
  INTEGER :: ipol,ik, npwx_g, nkbl, nkl, nkr, iks, ike, spin 
  INTEGER, ALLOCATABLE :: ngk_g(:)
  REAL(DP) :: xkg(3)
  REAL(DP) :: alat
  CHARACTER(LEN=6) :: cik
  !
  CALL start_clock('do_setup')
  !
  ! INIT PW
  !
  CALL init_pw_arrays(nbnd)
  CALL set_iks_l2g()
  !
  CALL set_dirs()
  !
  IF( mpime == root ) THEN 
     CALL json%initialize()
     CALL json%load_file(filename=TRIM(logfile))
  ENDIF
  !
  IF ( lsda ) THEN
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
     !IF ( .NOT. two_fermi_energies ) THEN
        DO iks = 1, nks
           spin = isk(iks)
           !
           SELECT CASE(spin)
           CASE(1)
              nelup = SUM( f_inp(:,1) )
           CASE(2)
              neldw = SUM( f_inp(:,2) )
           END SELECT
           !
        ENDDO
     ENDIF
     IF ( INT( nelup ) == 0 .AND. INT( neldw ) == 0 ) THEN
        CALL errore( 'do_setup','nelup = 0 and neldw = 0 ',1)
     ENDIF
  ENDIF
  !
  ! SYSTEM OVERVIEW
  !
  CALL io_push_title('System Overview')
  CALL io_push_value('gamma_only',gamma_only,20)
  IF( mpime == root ) CALL json%add('system.basis.gamma_only',gamma_only)
  CALL io_push_value('ecutwfc [Ry]',ecutwfc,20)
  IF( mpime == root ) CALL json%add('system.basis.ecutwfc:ry',ecutwfc)
  CALL io_push_es0('omega [au^3]',omega,20)
  IF( mpime == root ) CALL json%add('system.cell.omega:au',omega)
  IF ( gamma_only ) THEN
     auxi = npw
     CALL mp_sum(auxi,intra_bgrp_comm)
     CALL io_push_value('glob. #G',auxi,20)
     IF( mpime == root ) CALL json%add('system.basis.globg',auxi)
  ELSE
     ALLOCATE( ngk_g(nkstot) )
     !npool = nproc_image / nproc_pool
     nkbl = nkstot / kunit
     nkl = kunit * ( nkbl / npool )
     nkr = ( nkstot - nkl * npool ) / kunit
     IF ( my_pool_id < nkr ) nkl = nkl + kunit
     iks = nkl*my_pool_id + 1
     IF ( my_pool_id >= nkr ) iks = iks + nkr*kunit
     ike = iks + nkl - 1
     ngk_g = 0
     ngk_g(iks:ike) = ngk(1:nks)
     CALL mp_sum( ngk_g, inter_pool_comm )
     CALL mp_sum( ngk_g, intra_pool_comm )
     ngk_g = ngk_g / nbgrp
     npwx_g = MAXVAL( ngk_g(1:nkstot) )
     CALL io_push_value('glob. #PW',npwx_g,20)
     IF( mpime == root ) CALL json%add('system.basis.globpw',npwx_g)
     DEALLOCATE( ngk_g )
  ENDIF
  CALL io_push_value('nbnd',nbnd,20)
  IF( mpime == root ) CALL json%add('system.electron.nbnd',nbnd)
  CALL io_push_value('nkstot',nkstot,20)
  IF( mpime == root ) CALL json%add('system.electron.nkstot',nkstot)
  CALL io_push_value('nspin',nspin,20)
  IF( mpime == root ) CALL json%add('system.electron.nspin',nspin)
  CALL io_push_value('nelec',nelec,20)
  IF( mpime == root ) CALL json%add('system.electron.nelec',nelec)
  IF(nspin == 2) THEN
     CALL io_push_value('nelup',nelup,20)
     IF( mpime == root ) CALL json%add('system.electron.nelup',nelup)
     CALL io_push_value('neldw',neldw,20)
     IF( mpime == root ) CALL json%add('system.electron.neldw',neldw)
  ENDIF
  CALL io_push_value('npol',npol,20)
  IF( mpime == root ) CALL json%add('system.electron.npol',npol)
  CALL io_push_value('lsda',lsda,20)
  IF( mpime == root ) CALL json%add('system.electron.lsda',lsda)
  CALL io_push_value('noncolin',noncolin,20)
  IF( mpime == root ) CALL json%add('system.electron.noncolin',noncolin)
  CALL io_push_value('lspinorb',lspinorb,20)
  IF( mpime == root ) CALL json%add('system.electron.lspinorb',lspinorb)
  CALL io_push_value('domag',domag,20)
  IF( mpime == root ) CALL json%add('system.electron.domag',domag)
  CALL io_push_bar
  !
  alat = celldm(1)
  !
  WRITE( stdout, '(/5x,"sFFT G-space: ",i8," G-vectors", 5x, &
       &               "R-space: (",i4,",",i4,",",i4,")")') &
       &         ngms_g, dffts%nr1, dffts%nr2, dffts%nr3
  WRITE( stdout, '( 5x,"dFFT G-space: ",i8," G-vectors", 5x, &
       &               "R-space: (",i4,",",i4,",",i4,")")') &
       &         ngm_g, dfftp%nr1, dfftp%nr2, dfftp%nr3
  WRITE( stdout, '(/5x,"Cell [a.u.]          = ",3f14.6)') alat*at(1,1:3)
  WRITE( stdout, '( 5x,"                     = ",3f14.6)') alat*at(2,1:3)
  WRITE( stdout, '( 5x,"                     = ",3f14.6)') alat*at(3,1:3)
  WRITE( stdout, '( 5x," ")')
  IF( mpime == root ) THEN 
     CALL json%add('system.basis.sFFT.ngm',ngms_g)
     CALL json%add('system.basis.sFFT.nr1',dffts%nr1)
     CALL json%add('system.basis.sFFT.nr2',dffts%nr2)
     CALL json%add('system.basis.sFFT.nr3',dffts%nr3)
     CALL json%add('system.basis.dFFT.ngm',ngm_g)
     CALL json%add('system.basis.dFFT.nr1',dfftp%nr1)
     CALL json%add('system.basis.dFFT.nr2',dfftp%nr2)
     CALL json%add('system.basis.dFFT.nr3',dfftp%nr3)
     CALL json%add('system.cell.a1:au',alat*at(1:3,1))
     CALL json%add('system.cell.a2:au',alat*at(1:3,2))
     CALL json%add('system.cell.a3:au',alat*at(1:3,3))
     CALL json%add('system.cell.alat:au',alat)
  ENDIF
  !
  WRITE( stdout, '(5x,"number of ks points=",i6)') nkstot
  IF( mpime == root ) CALL json%add('system.kpt.nkstot',nkstot)
  WRITE( stdout, '(23x,"cart. coord. in units 2pi/alat")')
  DO ik = 1, nkstot
     WRITE( cik, '(i6)') ik 
     WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') ik, &
          (xk (ipol, ik) , ipol = 1, 3) , wk (ik)
     IF( mpime == root ) CALL json%add('system.kpt.k('//TRIM(ADJUSTL(cik))//').cartcoord:tpiba',xkg(1:3))
     IF( mpime == root ) CALL json%add('system.kpt.k('//TRIM(ADJUSTL(cik))//').weight',wk(ik))
  ENDDO
  WRITE( stdout, '(/23x,"cryst. coord.")')
  DO ik = 1, nkstot
     WRITE( cik, '(i6)') ik 
     DO ipol = 1, 3
        xkg(ipol) = at(1,ipol)*xk(1,ik) + at(2,ipol)*xk(2,ik) + &
                    at(3,ipol)*xk(3,ik)
        ! xkg are the component in the crystal RL basis
     ENDDO
     WRITE( stdout, '(8x,"k(",i5,") = (",3f12.7,"), wk =",f12.7)') &
          ik, (xkg (ipol) , ipol = 1, 3) , wk (ik)
     IF( mpime == root ) CALL json%add('system.kpt.k('//TRIM(ADJUSTL(cik))//').crystcoord',xkg(1:3))
  ENDDO
  WRITE( stdout, * )
  !
  IF( mpime == root ) THEN
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print_file( iunit )
     CLOSE( iunit )
     CALL json%destroy()
  ENDIF 
  !
  CALL stop_clock('do_setup')
  !
END SUBROUTINE 
