program host_model_simulator

!--------------------------------------------------------------------------------
! Program which simulates the host model for MICM
!--------------------------------------------------------------------------------
use chemSolve,          only: chemSolve_register,  chemSolve_init, chemSolve_run
use precision,          only: r8

implicit none

integer           :: nSpecies               ! Number of chemical species in run

real(r8),pointer  :: vmr_init(:)            ! Initial VMR
real(r8),pointer  :: vmr_curr(:)            ! Current VMR

real(r8)          :: timeStepSize = 2e+1_r8 ! seconds - this system is unstable, so 3e+1 fails

! convergence criteria will have to be set somewhere(Cafe) and passed to ode solver.
real(r8),pointer  :: relTol(:)              ! Relative tolerance
real(r8),pointer  :: absTol(:)              ! Absolute tolerance

integer           :: ierr                   ! Error code


! These will eventually be provided by the host model (are not currently used,  but  will be)

real(r8), parameter :: temperature = 273
real(r8), parameter :: pressure = 100000
real(r8), parameter :: mass_density = 1

!-----------------------------------------------
! Register the chemistry packages
!-----------------------------------------------

call chemSolve_register(nSpecies)

!-----------------------------------------------
! Allocate the local variables  (This will be done via CPF?)
!-----------------------------------------------
allocate (vmr_init(nSpecies))
allocate (vmr_curr(nSpecies))
allocate (absTol(nSpecies))
allocate (relTol(nSpecies))

!-----------------------------------------------
! Initialize the chemistry packages
!-----------------------------------------------

call chemSolve_init(absTol, relTol)

!-----------------------------------------------
! Explicitly specify the data which will come from the  cpf
!-----------------------------------------------

  vmr_init(:) =(/1._r8, 0._r8, 0._r8/)
  print *, 'initial value'
  print *, vmr_init


!-----------------------------------------------
! Simulate the XML file which CCPP will use to drive the model
! Only called at beginnning

  call  chemSolve_run(vmr_init, timeStepSize, absTol, relTol, vmr_curr, ierr)

  print *, 'final_value'
  print *, vmr_curr

!-----------------------------------
! some of these will be deallocated by CPF
!-----------------------------------

  deallocate (vmr_init)
  deallocate (vmr_curr)
  deallocate (absTol)
  deallocate (relTol)

end program host_model_simulator
