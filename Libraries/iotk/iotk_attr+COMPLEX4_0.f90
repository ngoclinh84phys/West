! Input/Output Tool Kit (IOTK)
! Copyright (C) 2004-2006 Giovanni Bussi
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


#ifdef __IOTK_COMPLEX4
#if 0 <= __IOTK_MAXRANK

! This is needed as a workaround for bugged pack 
subroutine iotk_private_pack_COMPLEX4(out,in,n,l)
    use iotk_base
    implicit none
    integer,                                    intent(in)  :: n,l
    COMPLEX (kind=iotk_COMPLEX4), intent(out) :: out(n)
    COMPLEX (kind=iotk_COMPLEX4), intent(in)  :: in(n)
    out = in
end subroutine iotk_private_pack_COMPLEX4

subroutine iotk_write_COMPLEX4(val,string,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_fmt_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  COMPLEX(kind=iotk_COMPLEX4), intent(in) :: val(:)
#ifdef __IOTK_WORKAROUND6
  character(len=*)              :: string
#else
  character(len=*), intent(out) :: string
#endif
  character(len=*), intent(in)  :: fmt
  integer, intent(out) :: ierr
  character(len=100) :: tmpval
  integer :: index,iostat
  ierr = 0
  iostat = 0 
  string(1:1) = iotk_eos
  if(size(val)==0) return
  if(len(string)==0) then
    call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
    return
  end if
  do index=1,size(val)
    if (trim(fmt)=="!") then
       write(tmpval,trim(iotk_wfmt("COMPLEX",kind(val),1,-1," ")),iostat=iostat) val(index)
    else
       write(tmpval,fmt=fmt,iostat=iostat) val(index)
    end if
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
call iotk_error_msg(ierr,' ')
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    call iotk_strcat(string,trim(adjustl(tmpval))//" ",ierr)
    if(ierr/=0) then
      call iotk_error_issue(ierr,"iotk_write",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
      return
    end if
  end do
! the last blank is deleted
  string(iotk_strlen(string):iotk_strlen(string)) = iotk_eos
end subroutine iotk_write_COMPLEX4

subroutine iotk_read_COMPLEX4(val,string,index,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_xtox_interf
  use iotk_misc_interf
  implicit none
  COMPLEX(kind=iotk_COMPLEX4), intent(inout) :: val(:)
  character(len=*), intent(in) :: string
  integer, intent(inout) :: index
  integer, intent(out) :: ierr
  integer :: pos,pos1,iostat
  integer :: maxindex
  real(kind=iotk_COMPLEX4) :: tmpreal
  complex(kind=iotk_COMPLEX4) :: tmpcomplex
#ifdef __IOTK_WORKAROUND9
  character(len=100) :: tmpstr ! debug
#endif
  pos = 0
  pos1= 0
  ierr = 0
  iostat = 0
   maxindex = 2 * size(val)
! for the moment, commas are considered as blanks
  do
    pos = verify(string(pos1+1:)," ,")+pos1
    if(pos==pos1) exit
    pos = pos - 1
    pos1 = scan(string(pos+1:)," ,")+pos
    if(pos1==pos) pos1 = len(string) + 1
!READ string(pos+1:pos1-1)
    index = index+1
    if(index>maxindex) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
call iotk_error_msg(ierr,'Too many data')
    end if
#ifdef __IOTK_WORKAROUND9
    tmpstr = TRIM( string(pos+1:pos1-1) )
    read( tmpstr,"(G100.95)",iostat=iostat) tmpreal
#else
    read(string(pos+1:pos1-1),"(G100.95)",iostat=iostat) tmpreal
#endif
    if(modulo(index,2)==1) then
      tmpcomplex = cmplx(tmpreal,aimag((val((index+1)/2))),kind=iotk_COMPLEX4)
    else
      tmpcomplex = cmplx(real(val((index+1)/2)),tmpreal,kind=iotk_COMPLEX4)
    end if
    val((index+1)/2) = tmpcomplex
    if(iostat/=0) then
      call iotk_error_issue(ierr,"iotk_read",__FILE__,__LINE__)
call iotk_error_msg(ierr,"CVS Revision: 1.21 ")
call iotk_error_msg(ierr,'Error reading a COMPLEX number from string')
call iotk_error_write(ierr,"string",string(pos+1:pos1-1))
call iotk_error_write(ierr,"iostat",iostat)
      return
    end if
    if(pos1>=len(string)) exit
  end do
end subroutine iotk_read_COMPLEX4

subroutine iotk_write_attr_COMPLEX4_0(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX4), intent(in)  :: val 
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  character(*), optional, intent(in) :: fmt
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
  character(len=300) :: usefmt
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
!-<
  if (present(fmt)) then
    usefmt = fmt
  else
    usefmt = "!"
  end if
!->
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  delim = '"'
  call iotk_write((/val/),tmpval,usefmt,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_COMPLEX4_0

subroutine iotk_scan_attr_COMPLEX4_0(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=iotk_COMPLEX4)                        :: val 
#else
  COMPLEX(kind=iotk_COMPLEX4), intent(out)           :: val 
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX4), optional, intent(in)  :: default 
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  integer :: index
  COMPLEX(kind=iotk_COMPLEX4), allocatable :: tmpval (:)
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  allocate(tmpval(1))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  if(index/=2*1) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute size does not match')
call iotk_error_write(ierrl,"attr",valc)
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
  val = tmpval(1)
  deallocate(tmpval)
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    val = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_0

#endif
#endif

subroutine iotk_attr_dummy_COMPLEX4_0
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX4_0




!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


#ifdef __IOTK_COMPLEX4
#if 1 <= __IOTK_MAXRANK



subroutine iotk_write_attr_COMPLEX4_1(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX4), intent(in)  :: val (:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  character(*), optional, intent(in) :: fmt
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
  character(len=300) :: usefmt
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
!-<
  if (present(fmt)) then
    usefmt = fmt
  else
    usefmt = "!"
  end if
!->
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  delim = '"'
  call iotk_write(pack(val,mask=.true.),tmpval,usefmt,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_COMPLEX4_1

subroutine iotk_scan_attr_COMPLEX4_1(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=iotk_COMPLEX4)                        :: val (:)
#else
  COMPLEX(kind=iotk_COMPLEX4), intent(out)           :: val (:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX4), optional, intent(in)  :: default (:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  integer :: index
  COMPLEX(kind=iotk_COMPLEX4), allocatable :: tmpval (:)
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  if(index/=2*size(val)) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute size does not match')
call iotk_error_write(ierrl,"attr",valc)
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
  val = reshape (source=tmpval,shape=shape(val))
  deallocate(tmpval)
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    val = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_1

#endif
#endif

subroutine iotk_attr_dummy_COMPLEX4_1
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX4_1




!------------------------------------------------------------------------------!
! Inclusion of configuration file
#include "iotk_config.h"
!------------------------------------------------------------------------------!

#include "iotk_auxmacros.h"


#ifdef __IOTK_COMPLEX4
#if 2 <= __IOTK_MAXRANK



subroutine iotk_write_attr_COMPLEX4_2(attr,name,val,dummy,first,newline,fmt,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*), intent(inout) :: attr
  character(*), intent(in)    :: name
  COMPLEX(kind=iotk_COMPLEX4), intent(in)  :: val (:,:)
  type(iotk_dummytype), optional :: dummy
  logical, optional, intent(in)  :: first
  logical, optional, intent(in)  :: newline
  character(*), optional, intent(in) :: fmt
  integer, optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen
  integer :: vallen
  integer :: namlen
  character :: delim
  character(len=300) :: usefmt
  character(iotk_vallenx) :: tmpval
  logical :: nl
  if(present(newline)) then
    nl = newline
  else
    nl = .false.
  endif
!-<
  if (present(fmt)) then
    usefmt = fmt
  else
    usefmt = "!"
  end if
!->
  ierrl = 0
  if(present(first)) then
    if(first) attr(1:1) = iotk_eos
  end if
  attlen = iotk_strlen_trim(attr)
  namlen = iotk_strlen_trim(name)
  if(.not.iotk_check_name(name)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Wrong tag name')
call iotk_error_write(ierrl,"name",name(1:namlen))
    goto 1
  end if
  delim = '"'
  call iotk_write(pack(val,mask=.true.),tmpval,usefmt,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  vallen = iotk_strlen(tmpval)
  if(attlen+vallen+namlen+5>len(attr)) then
    call iotk_error_issue(ierrl,"iotk_write_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute dummy argument is too short')
    goto 1
  end if
  if(.not. nl) then
    attr(attlen+1:attlen+vallen+namlen+5) = " "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  else
    attr(attlen+1:attlen+vallen+namlen+len(iotk_newline)+5) &
       = iotk_newline//" "//name(1:namlen)//"="//delim//tmpval(1:vallen)//delim//iotk_eos
  endif
1 continue
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_write_attr_COMPLEX4_2

subroutine iotk_scan_attr_COMPLEX4_2(attr,name,val,dummy,found,default,eos,ierr)
  use iotk_base
  use iotk_error_interf
  use iotk_attr_interf, only : iotk_read
  use iotk_str_interf
  use iotk_misc_interf
  implicit none
  character(*),             intent(in)  :: attr
  character(*),             intent(in)  :: name
#ifdef __IOTK_WORKAROUND6
  COMPLEX(kind=iotk_COMPLEX4)                        :: val (:,:)
#else
  COMPLEX(kind=iotk_COMPLEX4), intent(out)           :: val (:,:)
#endif
  type(iotk_dummytype), optional :: dummy
  logical,        optional, intent(out) :: found
  COMPLEX(kind=iotk_COMPLEX4), optional, intent(in)  :: default (:,:)
  logical,        optional, intent(in)  :: eos
  integer,        optional, intent(out) :: ierr
  integer :: ierrl
  integer :: attlen,pos,equal,namlen
  character :: delim
  logical :: foundl
  character(iotk_vallenx) :: valc
  integer :: index
  COMPLEX(kind=iotk_COMPLEX4), allocatable :: tmpval (:)
  ierrl = 0
  attlen=iotk_strlen(attr)
  namlen=iotk_strlen_trim(name)
  foundl = .false.
  equal = 0
  do
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) exit
    equal = equal + pos
    pos = scan(attr(equal+1:attlen),"=")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'')
call iotk_error_write(ierrl,"attr",attr(equal+1:attlen))
      goto 1
    end if
    equal = equal + pos
    if(trim(attr(equal-pos:equal-1))==name(1:namlen)) foundl = .true.
    pos = verify(attr(equal+1:attlen)," ")
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    equal = equal + pos
    delim = attr(equal:equal)
    if(delim/="'" .and. delim/='"') then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    pos = scan(attr(equal+1:attlen),delim)
    if(pos<=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
    if(foundl) exit
    equal = equal + pos
  end do
  if(foundl) then
    call iotk_strcpy(valc,attr(equal+1:equal+pos-1),ierrl)
    if(ierrl/=0) then
      call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
      goto 1
    end if
  else
    goto 1
  end if
  allocate(tmpval(size(val)))
  index = 0
  call iotk_str_clean(valc)
  call iotk_read(tmpval,valc(1:iotk_strlen(valc)),index,ierrl)
  if(ierrl/=0) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
    goto 1
  end if
  if(index/=2*size(val)) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute size does not match')
call iotk_error_write(ierrl,"attr",valc)
call iotk_error_write(ierrl,"size",size(tmpval))
    goto 1
  end if
  val = reshape (source=tmpval,shape=shape(val))
  deallocate(tmpval)
1 continue
  if(ierrl/=0) foundl=.false.
  if(present(found)) found = foundl
  if(ierrl==0 .and. .not. present(found) .and. .not. present(default) .and. .not. foundl) then
    call iotk_error_issue(ierrl,"iotk_scan_attr",__FILE__,__LINE__)
call iotk_error_msg(ierrl,"CVS Revision: 1.21 ")
call iotk_error_msg(ierrl,'Attribute not found')
call iotk_error_write(ierrl,"name",name)
    ierrl = - ierrl
  end if
  if(present(default) .and. .not. foundl) then
    val = default
  end if
  if(present(ierr)) then
    ierr = ierrl
  else
    if(ierrl/=0) call iotk_error_handler(ierrl)
  end if
end subroutine iotk_scan_attr_COMPLEX4_2

#endif
#endif

subroutine iotk_attr_dummy_COMPLEX4_2
  write(0,*)
end subroutine iotk_attr_dummy_COMPLEX4_2


