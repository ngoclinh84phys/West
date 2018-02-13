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
! Marco Govoni, Nicholas Brawand 
!
!----------------------------------------------------------------------------
SUBROUTINE do_sxx ( )
  !----------------------------------------------------------------------------
  !
  USE kinds,                 ONLY : DP
  USE uspp,                  ONLY : vkb,nkb
  USE io_global,             ONLY : stdout
  USE pwcom,                 ONLY : current_spin,wk,nks,nelup,neldw,isk,igk_k,xk,npw,npwx,lsda,nkstot,current_k,ngk,et
  USE io_push,               ONLY : io_push_title,io_push_bar
  USE westcom,               ONLY : iuwfc,lrwfc,westpp_range,westpp_save_dir,nbnd_occ,iks_l2g,westpp_epsinfty,dvg,ev,&
                                  & npwq,npwqx,fftdriver,logfile
  USE mp_global,             ONLY : inter_image_comm,my_image_id,intra_bgrp_comm
  USE mp,                    ONLY : mp_bcast,mp_sum
  USE fft_base,              ONLY : dfftp,dffts
  USE wvfct,                 ONLY : nbnd
  USE buffers,               ONLY : get_buffer
  USE wavefunctions_module,  ONLY : evc,psic,psic_nc
  USE bar,                   ONLY : bar_type,start_bar_type,update_bar_type,stop_bar_type
  USE fft_at_gamma,          ONLY : single_invfft_gamma,single_fwfft_gamma
  USE fft_at_k,              ONLY : single_invfft_k,single_fwfft_k
  USE distribution_center,   ONLY : pert
  USE control_flags,         ONLY : gamma_only 
  USE gvect,                 ONLY : g,nl,gstart,ngm_g,ngm
  USE cell_base,             ONLY : omega,at,alat
  USE noncollin_module,      ONLY : noncolin,npol 
  USE mp_world,              ONLY : mpime,root
  USE constants,             ONLY : rytoev
  USE json_module,           ONLY : json_file
  USE class_coulomb,         ONLY : coulomb
  USE types_coulomb,         ONLY : pot3D
  !
  IMPLICIT NONE
  !
  ! ... LOCAL variables
  !
  INTEGER :: ir, ip, ig, iks, ib, iv, ip_glob 
  COMPLEX(DP),ALLOCATABLE :: pertg(:),pertr(:),pertr_nc(:,:)
  TYPE(bar_type) :: barra
  REAL(DP),ALLOCATABLE :: sigma_exx( :, : ) 
  REAL(DP),ALLOCATABLE :: sigma_sxx( :, : ) 
  REAL(DP) :: peso
  REAL(DP), EXTERNAL :: DDOT
  CHARACTER(LEN=6) :: myglobk
  REAL(DP),ALLOCATABLE :: out_tab(:,:)
  COMPLEX(DP),ALLOCATABLE :: zproj(:,:)
  REAL(DP),ALLOCATABLE :: dproj(:,:)
  TYPE(json_file) :: json
  INTEGER :: iunit
  !
  CALL io_push_title("(S)creened eXact eXchange")
  !
  ALLOCATE( sigma_exx( westpp_range(1):westpp_range(2), nks) )
  ALLOCATE( sigma_sxx( westpp_range(1):westpp_range(2), nks) )
  !
  sigma_exx = 0._DP
  sigma_sxx = 0._DP
  !
  ALLOCATE( pertg(npwqx) )
  IF(noncolin) THEN 
     ALLOCATE( pertr_nc( dffts%nnr, npol ) )
  ELSE
     ALLOCATE( pertr( dffts%nnr ) )
  ENDIF
  !
  CALL pot3d%init('Wave','gb')
  !
  IF( gamma_only ) THEN 
     peso = 2._DP  
     ALLOCATE( dproj( 1, pert%nloc ) )
  ELSE
     peso = 1._DP
     ALLOCATE( zproj( 1, pert%nloc ) )
  ENDIF
  !
  CALL start_bar_type( barra, 'westpp', nks * (westpp_range(2)-westpp_range(1)+1)  ) 
  !
  DO iks = 1, nks  ! KPOINT-SPIN LOOP
     !
     ! ... Set k-point, spin, kinetic energy, needed by Hpsi
     !
     current_k = iks
     IF ( lsda ) current_spin = isk(iks)
     call g2_kin( iks )
     !
     ! ... More stuff needed by the hamiltonian: nonlocal projectors
     !
     !IF ( nkb > 0 ) CALL init_us_2( ngk(iks), igk_k(1,iks), xk(1,iks), vkb )
     npw = ngk(iks)
     !
     ! ... read in wavefunctions from the previous iteration
     !
     IF(nks>1) THEN
        !iuwfc = 20
        !lrwfc = nbnd * npwx * npol 
        !!CALL get_buffer( evc, nwordwfc, iunwfc, iks )
        IF(my_image_id==0) CALL get_buffer( evc, lrwfc, iuwfc, iks )
        !CALL mp_bcast(evc,0,inter_image_comm)
        !CALL davcio(evc,lrwfc,iuwfc,iks,-1)
        CALL mp_bcast(evc,0,inter_image_comm)
     ENDIF
     !
     !nbndval = nbnd_occ(iks)
     !
     DO ib = 1, nbnd
        !
        IF( ib < westpp_range(1) .OR. ib > westpp_range(2) ) CYCLE 
        !
        IF(gamma_only) THEN
           CALL single_invfft_gamma(dffts,npw,npwx,evc(1,ib),psic,'Wave') 
        ELSEIF(noncolin) THEN
           CALL single_invfft_k(dffts,npw,npwx,evc(1     ,ib),psic_nc(1,1),'Wave',igk_k(1,current_k))
           CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx,ib),psic_nc(1,2),'Wave',igk_k(1,current_k))
        ELSE
           CALL single_invfft_k(dffts,npw,npwx,evc(1,ib),psic,'Wave',igk_k(1,current_k))
        ENDIF
        !
        DO iv = 1, nbnd_occ(iks)
           !
           ! Bring it to R-space
           IF(gamma_only) THEN
              CALL single_invfft_gamma(dffts,npw,npwx,evc(1,iv),pertr,'Wave')
              DO ir=1,dffts%nnr 
                 pertr(ir)=psic(ir)*pertr(ir)
              ENDDO
              CALL single_fwfft_gamma(dffts,npwq,npwqx,pertr,pertg,TRIM(fftdriver))
           ELSEIF(noncolin) THEN
              CALL single_invfft_k(dffts,npw,npwx,evc(1     ,iv),pertr_nc(1,1),'Wave',igk_k(1,current_k))
              CALL single_invfft_k(dffts,npw,npwx,evc(1+npwx,iv),pertr_nc(1,2),'Wave',igk_k(1,current_k))
              DO ir=1,dffts%nnr 
                 pertr_nc(ir,1)=DCONJG(psic_nc(ir,1))*pertr_nc(ir,1)+DCONJG(psic_nc(ir,2))*pertr_nc(ir,2)
              ENDDO
              CALL single_fwfft_k(dffts,npwq,npwqx,pertr_nc(1,1),pertg,TRIM(fftdriver)) ! no igk
           ELSE
              CALL single_invfft_k(dffts,npw,npwx,evc(1,iv),pertr,'Wave',igk_k(1,current_k))
              DO ir=1,dffts%nnr 
                 pertr(ir)=DCONJG(psic(ir))*pertr(ir)
              ENDDO
              CALL single_fwfft_k(dffts,npwq,npwqx,pertr,pertg,TRIM(fftdriver)) ! no igk
           ENDIF 
           !
           DO ig = 1,npwq
              pertg(ig) = pertg(ig) * 3Dpot%sqvc(ig) 
           ENDDO
           sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - peso * DDOT( 2*npwq, pertg(1), 1, pertg(1), 1) / omega
           !IF(gstart==2) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) + REAL( pertg(1), KIND = DP )**2 / omega
           IF( ib == iv .AND. gstart == 2 ) sigma_exx( ib, iks ) = sigma_exx( ib, iks ) - pot3D%div
           !
           ! -- < SXX >
           IF( gamma_only ) THEN  
              CALL glbrak_gamma( pertg, dvg, dproj, npwq, npwqx, 1, pert%nloc, 1, npol)
              CALL mp_sum( dproj, intra_bgrp_comm )
              DO ip = 1, pert%nloc
                 ip_glob = pert%l2g(ip)
                 sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - dproj(1,ip)**2 * (ev(ip_glob)/(1._DP-ev(ip_glob))) / omega
              ENDDO
              IF( ib == iv ) sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - (1._DP/westpp_epsinfty-1._DP) * pot3D%div
           ELSE
              CALL glbrak_k( pertg, dvg, zproj, npwq, npwqx, 1, pert%nloc, 1, npol)
              CALL mp_sum( zproj, intra_bgrp_comm )
              DO ip = 1, pert%nloc
                 ip_glob = pert%l2g(ip)
                 sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - REAL(zproj(1,ip)*CONJG(zproj(1,ip)),KIND=DP) &
                                      & * (ev(ip_glob)/(1._DP-ev(ip_glob))) / omega
              ENDDO
              IF( ib == iv ) sigma_sxx( ib, iks ) = sigma_sxx( ib, iks ) - (1._DP/westpp_epsinfty-1._DP) * pot3D%div
           ENDIF 
           ! -- </ SXX > 
           !
        ENDDO
        !
        CALL update_bar_type( barra,'westpp', 1 )
        !
     ENDDO
     !
  ENDDO
  !
  CALL stop_bar_type( barra, 'westpp' )
  !
  CALL mp_sum( sigma_exx, intra_bgrp_comm )
  CALL mp_sum( sigma_sxx, inter_image_comm )
  !
  sigma_sxx = sigma_exx + sigma_sxx 
  !
  DEALLOCATE( pertg ) 
  IF( noncolin ) THEN 
    DEALLOCATE( pertr_nc ) 
  ELSE
    DEALLOCATE( pertr ) 
  ENDIF
  IF( gamma_only ) THEN
     DEALLOCATE( dproj )
  ELSE
     DEALLOCATE( zproj )
  ENDIF  
  !
  ! Output it per k-point
  !
  IF(mpime==root) THEN
     CALL json%initialize()
     CALL json%load_file(filename=TRIM(logfile))
  ENDIF
  !
  ! STDOUT
  !
  WRITE(stdout,"(5X)")  
  CALL io_push_bar()
  WRITE(stdout,"(5X,a,1X,a,1X,a,1X,a,1X,a,1X,a)") &
  & 'K     ', 'B     ', '      Eks [eV]', '       Sx [eV]', '      Sxx [eV]', '        Sxx/Sx'
  CALL io_push_bar()
  !
  ALLOCATE(out_tab(westpp_range(2)-westpp_range(1)+1,5))
  !
  DO iks=1,nks
     DO ib = westpp_range(1), westpp_range(2)
        out_tab( ib - westpp_range(1) + 1, 1) = REAL( ib, KIND=DP) 
        out_tab( ib - westpp_range(1) + 1, 2) = et(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 3) = sigma_exx(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 4) = sigma_sxx(ib,iks) * rytoev
        out_tab( ib - westpp_range(1) + 1, 5) = sigma_sxx(ib,iks) / sigma_exx(ib,iks)
        WRITE(stdout,"(5X,i6.6,1X,i6.6,1X,1f14.6,1X,1f14.6,1X,1f14.6,1X,1f14.6)") &
        & iks, ib, et(ib,iks)*rytoev, sigma_exx(ib,iks)*rytoev, sigma_sxx(ib,iks)*rytoev, sigma_sxx(ib,iks)/sigma_exx(ib,iks)
     ENDDO
     WRITE(myglobk,'(i6.6)') iks_l2g(iks)
     !
     !CALL serial_table_output(mpime==root,4000,'sxx_K'//myglobk,out_tab,&
     !& westpp_range(2)-westpp_range(1)+1,5,&
     !& (/'      band','    E0[eV]','    Sx[eV]','   Sxx[eV]','    Sxx/Sx'/))
     !
     IF(mpime==root) THEN 
        !
        CALL json%add('output.S.K'//TRIM(myglobk)//'.bandmap',out_tab(:,1))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.eks',out_tab(:,2))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.sx',out_tab(:,3))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.sxx',out_tab(:,4))
        CALL json%add('output.S.K'//TRIM(myglobk)//'.fraction',out_tab(:,5))
        !
     ENDIF
  ENDDO
  DEALLOCATE(out_tab)
  CALL io_push_bar()
  !
  IF( mpime == root ) THEN 
     !
     OPEN( NEWUNIT=iunit, FILE=TRIM(logfile) )
     CALL json%print_file( iunit )
     CLOSE( iunit )
     !
     CALL json%destroy()
     !
  ENDIF
  !
  DEALLOCATE( sigma_exx )
  DEALLOCATE( sigma_sxx ) 
  !
END SUBROUTINE
