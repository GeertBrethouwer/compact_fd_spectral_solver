! stand-in for hichan's mod_param (which needs MPI): only abort_run is used by compact.f90
module mod_param
  implicit none
contains
  subroutine abort_run(msg)
    character(len=*), intent(in) :: msg
    print*,'ABORT: ',msg; stop 1
  end subroutine
end module
