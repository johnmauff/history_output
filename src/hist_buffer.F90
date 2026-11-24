module hist_buffer
   ! Device residency (OpenACC builds only): when this file is compiled with
   ! OpenACC, every buffer's arrays (data, num_samples, var_buffer) are
   ! created on the device when first allocated and accumulated there. The
   ! input field may be on the device (accumulated in place) or on the host
   ! (copied in for the call). The host copies are brought up to date only
   ! by sync_host, which the hist_api read routines call before reading.
   ! Device kernels are plain-array routines with the host arithmetic and
   ! order, so results are bit-for-bit identical. Without OpenACC nothing
   ! changes.
   use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64
   use hist_hashable,                 only: hist_hashable_t
#ifdef _OPENACC
   use openacc,                       only: acc_is_present
#endif

   implicit none
   private

   ! Public interfaces
   public :: buffer_factory

   ! Accumulation types -- the integers and array positions below must match
   integer, parameter, public :: hist_accum_lst = 1 ! last sample
   integer, parameter, public :: hist_accum_min = 2 ! minimum sample
   integer, parameter, public :: hist_accum_max = 3 ! maximum sample
   integer, parameter, public :: hist_accum_avg = 4 ! sample average
   integer, parameter, public :: hist_accum_var = 5 ! sample standard deviation

   ! Number of accumulation types defined above
   integer, parameter, public :: hist_num_accum_types = 5

   integer, parameter         :: as_len = 36
   character(len=as_len), parameter, public ::                                &
        accum_strings(hist_num_accum_types) = [                               &
        'last sampled value                  ',                               &
        'minimum of sampled values           ',                               &
        'maximum of sampled values           ',                               &
        'average of sampled values           ',                               &
        'standard deviation of sampled values']

   character(len=3), parameter, public ::                                     &
        accum_abbrev(hist_num_accum_types) =                                  &
        ['lst', 'min', 'max', 'avg', 'var']

   type, abstract, public :: hist_buffer_t
      ! hist_buffer_t is an abstract base class for hist_outfld buffers
      class(hist_hashable_t), pointer              :: field_info => null()
      integer,                             private :: vol = -1 ! For host output
      integer,                             private :: horiz_axis_ind = 0
      integer,                             private :: rank = 0
      integer,                             private :: accum_type = 0
      integer,                allocatable, private :: field_shape(:)
      integer,                allocatable, private :: block_begs(:)
      integer,                allocatable, private :: block_ends(:)
      integer,                allocatable, private :: num_samples(:)
      character(len=:),       allocatable, private :: buff_type
      ! .true. once the buffer's arrays have device copies (OpenACC builds)
      logical,                             private :: on_device = .false.
      class(hist_buffer_t),   pointer              :: next => null()
   contains
      procedure                              :: field  => get_field_info
      procedure                              :: volume => get_volume
      procedure                              :: get_shape
      procedure                              :: horiz_axis_index
      procedure                              :: init_buffer
      procedure                              :: accum_string
      procedure                              :: buffer_type
      procedure                              :: check_status
      procedure                              :: has_blocks
      procedure                              :: get_num_samples
      procedure                              :: set_num_samples
      procedure                              :: is_on_device
      procedure(hist_buff_clear),   deferred :: clear
      procedure(hist_buff_init),    deferred :: initialize
      procedure(hist_buff_dev),     deferred :: sync_host
      procedure(hist_buff_dev),     deferred :: device_release
   end type hist_buffer_t

   type, public, extends(hist_buffer_t) :: hist_buff_1d_t
      real(REAL64), allocatable :: data(:)
      real(REAL64), allocatable :: var_buffer(:)
   contains
      procedure :: clear => buff_1d_clear
      procedure :: accumulate => buff_1d_accum
      procedure :: norm_value => buff_1d_value
      procedure :: initialize => init_buff_1d
      procedure :: sync_host => buff_1d_sync_host
      procedure :: device_release => buff_1d_device_release
      procedure :: get_var_buffer => get_var_buffer_1d
      procedure :: set_field_buffer => set_field_buffer_1d
      procedure :: set_var_buffer => set_var_buffer_1d
   end type hist_buff_1d_t

   type, public, extends(hist_buffer_t) :: hist_buff_2d_t
      real(REAL64), allocatable :: data(:,:)
      real(REAL64), allocatable :: var_buffer(:,:)
   contains
      procedure :: clear => buff_2d_clear
      procedure :: accumulate => buff_2d_accum
      procedure :: norm_value => buff_2d_value
      procedure :: initialize => init_buff_2d
      procedure :: sync_host => buff_2d_sync_host
      procedure :: device_release => buff_2d_device_release
      procedure :: check_fill_value => buff_2d_check_fill
      procedure :: get_var_buffer => get_var_buffer_2d
      procedure :: set_field_buffer => set_field_buffer_2d
      procedure :: set_var_buffer => set_var_buffer_2d
   end type hist_buff_2d_t

   type, public, extends(hist_buffer_t) :: hist_buff_3d_t
      real(REAL64), allocatable :: data(:,:,:)
      real(REAL64), allocatable :: var_buffer(:,:,:)
   contains
      procedure :: clear => buff_3d_clear
      procedure :: accumulate => buff_3d_accum
      procedure :: norm_value => buff_3d_value
      procedure :: initialize => init_buff_3d
      procedure :: sync_host => buff_3d_sync_host
      procedure :: device_release => buff_3d_device_release
      procedure :: check_fill_value => buff_3d_check_fill
      procedure :: get_var_buffer => get_var_buffer_3d
      procedure :: set_field_buffer => set_field_buffer_3d
      procedure :: set_var_buffer => set_var_buffer_3d
   end type hist_buff_3d_t

   ! Abstract interfaces for hist_buffer_t class
   abstract interface
      subroutine hist_buff_sub_log(this, logger)
         use hist_msg_handler, only: hist_log_messages
         import                   :: hist_buffer_t
         implicit none
         class(hist_buffer_t),              intent(inout) :: this
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_sub_log
   end interface

   abstract interface
      subroutine hist_buff_init(this, field_in, volume_in, horiz_axis_in,     &
           accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
         use hist_msg_handler,             only: hist_log_messages
         use, intrinsic :: ISO_FORTRAN_ENV, only: REAL64
         import                   :: hist_buffer_t
         import                   :: hist_hashable_t
         implicit none
         class(hist_buffer_t),              intent(inout) :: this
         class(hist_hashable_t),  pointer,  intent(in)    :: field_in
         integer,                           intent(in)    :: volume_in
         integer,                           intent(in)    :: horiz_axis_in
         integer,                           intent(in)    :: accum_type_in
         integer,                           intent(in)    :: shape_in(:)
         integer,                 optional, intent(in)    :: block_sizes_in(:)
         integer,                 optional, intent(in)    :: block_ind_in
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_init
   end interface

   abstract interface
      subroutine hist_buff_clear(this, logger)
         use hist_msg_handler, only: hist_log_messages
         import                  :: hist_buffer_t
         implicit none
         class(hist_buffer_t),              intent(inout) :: this
         type(hist_log_messages), optional, intent(inout) :: logger
      end subroutine hist_buff_clear
   end interface

   abstract interface
      subroutine hist_buff_dev(this)
         import                  :: hist_buffer_t
         implicit none
         class(hist_buffer_t), intent(inout) :: this
      end subroutine hist_buff_dev
   end interface

#ifdef _OPENACC
   ! Plain-array device helpers (device code never indexes the buffer
   ! types' allocatable components directly)
   interface dev_create
      module procedure dev_create_i1, dev_create_r1, dev_create_r2, dev_create_r3
   end interface dev_create
   interface dev_delete
      module procedure dev_delete_i1, dev_delete_r1, dev_delete_r2, dev_delete_r3
   end interface dev_delete
   interface dev_update_self
      module procedure dev_update_self_i1, dev_update_self_r1,                &
           dev_update_self_r2, dev_update_self_r3
   end interface dev_update_self
   interface dev_update_device
      module procedure dev_update_device_i1, dev_update_device_r1,            &
           dev_update_device_r2, dev_update_device_r3
   end interface dev_update_device
   interface dev_zero
      module procedure dev_zero_i1, dev_zero_r1, dev_zero_r2, dev_zero_r3
   end interface dev_zero
#endif

contains

   !#######################################################################

   function get_field_info(this) result(field_info)
      class(hist_buffer_t), intent(in) :: this
      class(hist_hashable_t), pointer  :: field_info

      field_info => this%field_info
   end function get_field_info

   !#######################################################################

   pure function get_volume(this) result(volume)
      class(hist_buffer_t), intent(in) :: this
      integer                          :: volume

      volume = this%vol
   end function get_volume

   !#######################################################################

   pure function get_shape(this) result(fshape)
      class(hist_buffer_t), intent(in) :: this
      integer, allocatable :: fshape(:)

      allocate(fshape(size(this%field_shape)))

      fshape = this%field_shape
   end function get_shape

   !#######################################################################

   pure function horiz_axis_index(this) result(axis_index)
      class(hist_buffer_t), intent(in) :: this
      integer                          :: axis_index

      axis_index = this%horiz_axis_ind
   end function horiz_axis_index
   !#######################################################################

   function check_status(this, logger, filename, line) result(is_valid)
      ! Check to see if this buffer is properly initialized
      use hist_msg_handler, only: hist_log_messages, hist_add_error, ERROR

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      character(len=*),        optional, intent(in)    :: filename
      integer,                 optional, intent(in)    :: line
      logical                                          :: is_valid
      ! Local variable
      character(len=*), parameter :: subname = 'check_status'

      is_valid = .true.
      if ((this%horiz_axis_index() < 1)      .or.                             &
           (.not. allocated(this%field_shape))) then
         is_valid = .false.
         call hist_add_error(subname,                                         &
              "buffer not properly initialized '", errors=logger)
         if (present(filename) .and. present(line) .and. present(logger)) then
            call logger%add_stack_frame(ERROR, filename, line)
         end if
      end if

   end function check_status

   !#######################################################################

   pure function has_blocks(this) result(blocked)
      ! Dummy argument
      class(hist_buffer_t),   intent(in)    :: this
      logical                               :: blocked

      blocked = allocated(this%block_begs) .and. allocated(this%block_ends)

   end function has_blocks

   !#######################################################################

   function get_num_samples(this, logger) result(samples)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      ! Dummy arguments
      class(hist_buffer_t),                 intent(in) :: this
      integer, allocatable                             :: samples(:)
      type(hist_log_messages), optional, intent(inout) :: logger

      ! Local variables
      integer :: ierr
      character(len=512) :: errmsg
      character(len=*), parameter :: subname = 'get_num_samples'

      allocate(samples(size(this%num_samples)), stat=ierr, errmsg=errmsg)
      if (ierr /= 0) then
         call hist_add_error(subname, 'failed to allocate num_samples; errmsg='//trim(errmsg), &
                 errors=logger)
      end if

      samples = this%num_samples

   end function get_num_samples

   !#######################################################################

   subroutine set_num_samples(this, samples)
      ! Dummy arguments
      class(hist_buffer_t), intent(inout) :: this
      integer,                 intent(in) :: samples(:)

      this%num_samples = samples
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%num_samples)
      end if
#endif

   end subroutine set_num_samples

   !#######################################################################

   pure function is_on_device(this) result(on_dev)
      ! .true. if the buffer's arrays are device-resident (OpenACC builds)
      class(hist_buffer_t), intent(in) :: this
      logical                          :: on_dev

      on_dev = this%on_device
   end function is_on_device

   !#######################################################################

   subroutine init_buffer(this, field_in, volume_in, horiz_axis_in,           &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error, error_length
      use hist_msg_handler, only: hist_add_error

      ! Dummy arguments
      class(hist_buffer_t),              intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: astat
      integer                     :: hsize
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'init_buffer'

      ! Sanity check
      this%rank = SIZE(shape_in, 1)
      if ((horiz_axis_in < 1) .or. (horiz_axis_in > this%rank)) then
         call hist_add_error(subname,                                         &
              "horiz_axis_in must be between 1 and '", errint1=this%rank,     &
              errors=logger)
      else
         this%field_info => field_in
         this%vol = volume_in
         this%horiz_axis_ind = horiz_axis_in
         this%accum_type = accum_type_in
         allocate(this%field_shape(size(shape_in, 1)), stat=astat, errmsg=errmsg)
         if (astat == 0) then
            this%field_shape(:) = shape_in(:)
         else
            call hist_add_alloc_error('field_shape, error: '//errmsg, __FILE__, __LINE__ - 4,  &
                 subname=subname, errors=logger)
         end if
      end if
   end subroutine init_buffer

   !#######################################################################

   function accum_string(this) result(ac_str)
      class(hist_buffer_t), intent(in)   :: this
      character(len=as_len), allocatable :: ac_str

      ac_str = accum_strings(this%accum_type)

   end function accum_string

   !#######################################################################

   function buffer_type(this) result(bt_str)
      class(hist_buffer_t), intent(in) :: this
      character(len=:), allocatable    :: bt_str

      if (allocated(this%buff_type)) then
         bt_str = this%buff_type
      else
         bt_str = 'UNKNOWN'
      end if

   end function buffer_type

   !#######################################################################

   subroutine buff_1d_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error, error_length

      ! Dummy arguments
      class(hist_buff_1d_t),       intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'buff_1d_clear'

      if (.not. allocated(this%data)) then
         allocate(this%data(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('data, error: '//errmsg, __FILE__, __LINE__ - 1,  &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%var_buffer) .and. this%accum_type == hist_accum_var) then
         allocate(this%var_buffer(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('var_buffer, error: '//errmsg, __FILE__, __LINE__ - 2, &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%num_samples)) then
         allocate(this%num_samples(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('num_samples, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL64
      this%num_samples = 0
      if (allocated(this%var_buffer)) then
         this%var_buffer = 0.0_REAL64
      end if
#ifdef _OPENACC
      ! First clear: create the device copies; every clear: zero them too
      if (.not. this%on_device) then
         call dev_create(this%data)
         call dev_create(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_create(this%var_buffer)
         end if
         this%on_device = .true.
      end if
      call dev_zero(this%data)
      call dev_zero(this%num_samples)
      if (allocated(this%var_buffer)) then
         call dev_zero(this%var_buffer)
      end if
#endif

   end subroutine buff_1d_clear

   !#######################################################################

   subroutine init_buff_1d(this, field_in, volume_in, horiz_axis_in, &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages

      class(hist_buff_1d_t),             intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in, accum_type_in, &
              shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_1d'

   end subroutine init_buff_1d

   !#######################################################################

   subroutine buff_1d_accum(this, field, cols_or_block, flag_xyfill, &
         fill_value, cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
      ! Dummy arguments
      class(hist_buff_1d_t),       intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: cols_or_block
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer            :: col_beg_use
      integer            :: col_end_use
      integer            :: ind1
      real(REAL64)       :: fld_val, tmp
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'buff_1d_accum'

      if (this%has_blocks()) then
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

#ifdef _OPENACC
      if (this%on_device) then
         call buff_1d_accum_device(this, field, col_beg_use, col_end_use, &
              flag_xyfill, fill_value, logger)
         return
      end if
#endif

      select case (this%accum_type)
      case (hist_accum_lst)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               if (fld_val == fill_value) then
                  ! Set num samples to 0 if this is a fill value
                  ! This will trigger logic in _value routine
                  this%data(ind1) = fill_value
                  this%num_samples(ind1) = 0
               else
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if
            else
               this%data(ind1) = fld_val
               this%num_samples(ind1) = 1
            end if
         end do
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               ! If the buffer (possibly) contains fill values,
               ! only check for minimum if not a fill value
               if (this%num_samples(ind1) == 0) then
                  if (fld_val /= fill_value) then
                     this%data(ind1) = fld_val
                     this%num_samples(ind1) = 1
                  else
                     ! Set to large positive number if this is a fill value
                     ! Also do not change num_samples - if num_samples
                     ! is zero at the end of the accumulation period,
                     ! the value will be overwritten to fill_value
                     this%data(ind1) = HUGE(REAL64)
                  end if
               else if (fld_val < this%data(ind1) .and. fld_val /= fill_value) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the minimum value for this col
            else
               if (this%num_samples(ind1) == 0 .or. fld_val < this%data(ind1)) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the minimum value for this col
            end if
         end do
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               ! If the buffer (possibly) contains fill values,
               ! only check for maximum if not a fill value
               if (this%num_samples(ind1) == 0) then
                  if (fld_val /= fill_value) then
                     this%data(ind1) = fld_val
                     this%num_samples(ind1) = 1
                  else
                     ! Set to large negative number if this is a fill value
                     ! Also do not change num_samples - if num_samples
                     ! is zero at the end of the accumulation period,
                     ! the value will be overwritten to fill_value
                     this%data(ind1) = -HUGE(REAL64)
                  end if
               else if (fld_val > this%data(ind1) .and. fld_val /= fill_value) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the maximum value for this col
            else
               if (this%num_samples(ind1) == 0 .or. fld_val > this%data(ind1)) then
                  this%data(ind1) = fld_val
                  this%num_samples(ind1) = 1
               end if ! No else, we already have the maximum value for this col
            end if
         end do
      case (hist_accum_avg)
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               ! Only include sample if not the fill value
               if (fld_val /= fill_value) then
                  this%data(ind1) = this%data(ind1) + fld_val
                  this%num_samples(ind1) = this%num_samples(ind1) + 1
               end if
            else
               this%data(ind1) = this%data(ind1) + fld_val
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
         end do
      case (hist_accum_var)
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = col_beg_use, col_end_use
            fld_val = field(ind1 - col_beg_use + 1)
            if (flag_xyfill) then
               if (fld_val /= fill_value) then
                  ! Only include the sample if it's not a fill value
                  if (this%num_samples(ind1) == 0) then
                     this%data(ind1) = fld_val
                     this%var_buffer(ind1) = 0.0_REAL64
                     this%num_samples(ind1) = this%num_samples(ind1) + 1
                  else
                     tmp = this%data(ind1)
                     this%num_samples(ind1) = this%num_samples(ind1) + 1
                     this%data(ind1) = this%data(ind1) + &
                             (fld_val - this%data(ind1)) / &
                             this%num_samples(ind1)
                     this%var_buffer(ind1) = this%var_buffer(ind1) + &
                             (fld_val - this%data(ind1)) * &
                             (fld_val - tmp)
                  end if
               end if
            else
               if (this%num_samples(ind1) == 0) then
                  this%data(ind1) = fld_val
                  this%var_buffer(ind1) = 0.0_REAL64
                  this%num_samples(ind1) = this%num_samples(ind1) + 1
               else
                  tmp = this%data(ind1)
                  this%num_samples(ind1) = this%num_samples(ind1) + 1
                  this%data(ind1) = this%data(ind1) + &
                          (fld_val - this%data(ind1)) / &
                          this%num_samples(ind1)
                  this%var_buffer(ind1) = this%var_buffer(ind1) + &
                          (fld_val - this%data(ind1)) * &
                          (fld_val - tmp)
               end if
            end if
         end do
      case default
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error(subname, errmsg, errors=logger)
      end select

   end subroutine buff_1d_accum

   !#######################################################################

   pure subroutine buff_1d_value(this, norm_val, flag_xyfill, fill_value)
      ! Dummy arguments
      class(hist_buff_1d_t),             intent(in)    :: this
      real(REAL64),                      intent(inout) :: norm_val(:)
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      ! Local variable
      integer :: ind1
      integer :: nacc
      real(REAL64) :: tmpfill, variance

      ! Initialize values that are not fill values:
      tmpfill = merge(fill_value, 0.0_REAL64, flag_xyfill)
      if (this%accum_type /= hist_accum_var) then
         do ind1 = 1, this%field_shape(1)
            nacc = this%num_samples(ind1)
            if (nacc > 0) then
               norm_val(ind1) = this%data(ind1) / nacc
            else
               norm_val(ind1) = tmpfill
            end if
         end do
      else
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = 1, this%field_shape(1)
            if (this%num_samples(ind1) > 0) then
               variance = this%var_buffer(ind1) / this%num_samples(ind1)
               norm_val(ind1) = sqrt(variance)
            else
               norm_val(ind1) = tmpfill
            end if
         end do
      end if

   end subroutine buff_1d_value

   !#######################################################################

   subroutine buff_2d_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error, error_length

      ! Dummy arguments
      class(hist_buff_2d_t),       intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'buff_2d_clear'

      if (.not. allocated(this%data)) then
         allocate(this%data(this%field_shape(1), this%field_shape(2)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('data, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%var_buffer) .and. this%accum_type == hist_accum_var) then
         allocate(this%var_buffer(this%field_shape(1), this%field_shape(2)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('var_buffer, error: '//errmsg, __FILE__, __LINE__ - 2,  &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%num_samples)) then
         allocate(this%num_samples(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('num_samples, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL64
      this%num_samples = 0
      if (allocated(this%var_buffer)) then
         this%var_buffer = 0.0_REAL64
      end if
#ifdef _OPENACC
      ! First clear: create the device copies; every clear: zero them too
      if (.not. this%on_device) then
         call dev_create(this%data)
         call dev_create(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_create(this%var_buffer)
         end if
         this%on_device = .true.
      end if
      call dev_zero(this%data)
      call dev_zero(this%num_samples)
      if (allocated(this%var_buffer)) then
         call dev_zero(this%var_buffer)
      end if
#endif

   end subroutine buff_2d_clear

   !#######################################################################

   subroutine init_buff_2d(this, field_in, volume_in, horiz_axis_in, &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages

      class(hist_buff_2d_t),             intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in, accum_type_in, &
              shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_2d'

   end subroutine init_buff_2d

   !#######################################################################

   subroutine buff_2d_accum(this, field, cols_or_block, flag_xyfill, &
         fill_value, cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
      ! Dummy arguments
      class(hist_buff_2d_t),       intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:)
      integer,                           intent(in)    :: cols_or_block
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer            :: col_beg_use
      integer            :: col_end_use
      integer            :: ind1, ind2
      real(REAL64)       :: fld_val, tmp
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'buff_2d_accum'

      if (this%has_blocks()) then
         ! For a blocked field, <cols_or_block> is a block index
         col_beg_use = this%block_begs(cols_or_block)
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

#ifdef _OPENACC
      if (this%on_device) then
         call buff_2d_accum_device(this, field, col_beg_use, col_end_use, &
              flag_xyfill, fill_value, logger)
         return
      end if
#endif

      select case (this%accum_type)
      case (hist_accum_lst)
         do ind1 = col_beg_use, col_end_use
            ! Only set samples for this column if not a fill value
            if (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value) then
               this%num_samples(ind1) = 0
            else
               this%num_samples(ind1) = 1
            end if
            do ind2 = 1, this%field_shape(2)
               fld_val = field(ind1 - col_beg_use + 1, ind2)
               if (flag_xyfill .and. fld_val == fill_value) then
                  this%data(ind1, ind2) = fill_value
               else
                  this%data(ind1, ind2) = fld_val
               end if
            end do
         end do
         if (flag_xyfill) then
            ! Check if assumption that columns are consistent wrt fill value was accurate
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            ! Has this column been sampled yet?
            if (this%num_samples(ind1) == 0) then
               ! First sample for this column: initialize all levels.
               ! Skip fill value columns (num_samples stays 0;
               ! buff_2d_value will write fill_value at output time)
               if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value)) then
                  do ind2 = 1, this%field_shape(2)
                     this%data(ind1, ind2) = field(ind1 - col_beg_use + 1, ind2)
                  end do
                  this%num_samples(ind1) = 1
               end if
            else
               ! Subsequent samples in column: update minimum per level
               do ind2 = 1, this%field_shape(2)
                  fld_val = field(ind1 - col_beg_use + 1, ind2)
                  if (flag_xyfill) then
                     if (fld_val < this%data(ind1, ind2) .and. &
                         fld_val /= fill_value) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  else
                     if (fld_val < this%data(ind1, ind2)) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  end if
               end do
            end if
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            ! Has this column been sampled yet?
            if (this%num_samples(ind1) == 0) then
               ! First sample for this column: initialize all levels.
               ! Skip fill value columns (num_samples stays 0;
               ! buff_2d_value will write fill_value at output time)
               if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value)) then
                  do ind2 = 1, this%field_shape(2)
                     this%data(ind1, ind2) = field(ind1 - col_beg_use + 1, ind2)
                  end do
                  this%num_samples(ind1) = 1
               end if
            else
               ! Subsequent samples: update maximum per level
               do ind2 = 1, this%field_shape(2)
                  fld_val = field(ind1 - col_beg_use + 1, ind2)
                  if (flag_xyfill) then
                     if (fld_val > this%data(ind1, ind2) .and. fld_val /= fill_value) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  else
                     if (fld_val > this%data(ind1, ind2)) then
                        this%data(ind1, ind2) = fld_val
                     end if
                  end if
               end do
            end if
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_avg)
         do ind1 = col_beg_use, col_end_use
            if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value)) then
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
            do ind2 = 1, this%field_shape(2)
               fld_val = field(ind1 - col_beg_use + 1, ind2)
               ! Compute running sum
               if (flag_xyfill) then
                  ! Only include sample if it is not the fill value
                  if (fld_val /= fill_value) then
                     this%data(ind1, ind2) = this%data(ind1, ind2) + fld_val
                  end if
               else
                  this%data(ind1, ind2) = this%data(ind1, ind2) + fld_val
               end if
            end do
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_var)
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = col_beg_use, col_end_use
            if (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1) == fill_value) then
               this%num_samples(ind1) = this%num_samples(ind1)
            else
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
            do ind2 = 1, this%field_shape(2)
               fld_val = field(ind1 - col_beg_use + 1, ind2)
               if (flag_xyfill) then
                  if (fld_val /= fill_value) then
                     ! Only include the sample if it's not a fill value
                     if (this%num_samples(ind1) == 1) then
                        this%data(ind1, ind2) = fld_val
                        this%var_buffer(ind1, ind2) = 0.0_REAL64
                     else
                        tmp = this%data(ind1, ind2)
                        this%data(ind1, ind2) = this%data(ind1, ind2) + &
                                (fld_val - this%data(ind1, ind2)) / &
                                this%num_samples(ind1)
                        this%var_buffer(ind1, ind2) = this%var_buffer(ind1, ind2) + &
                                (fld_val - this%data(ind1, ind2)) * &
                                (fld_val - tmp)
                     end if
                  end if
               else
                  if (this%num_samples(ind1) == 1) then
                     this%data(ind1, ind2) = fld_val
                     this%var_buffer(ind1, ind2) = 0.0_REAL64
                  else
                     tmp = this%data(ind1, ind2)
                     this%data(ind1, ind2) = this%data(ind1, ind2) + &
                             (fld_val - this%data(ind1, ind2)) / &
                             this%num_samples(ind1)
                     this%var_buffer(ind1, ind2) = this%var_buffer(ind1, ind2) + &
                             (fld_val - this%data(ind1, ind2)) * &
                             (fld_val - tmp)
                  end if
               end if
            end do
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case default
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error(subname, errmsg, errors=logger)
      end select

   end subroutine buff_2d_accum

   !#######################################################################

   subroutine buff_2d_check_fill(this, field, fill_value, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      ! Ensure that columns have fill value applied consistently across levels
      class(hist_buff_2d_t),             intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:)
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      character(len=512) :: errstr
      integer :: idx
      logical :: col_is_fill

      ! Every level in a column must match the fill status of its first level
      do idx = 1, size(field, 1)
         col_is_fill = (field(idx,1) == fill_value)
         if (any((field(idx,:) == fill_value) .neqv. col_is_fill)) then
            write(errstr, '(a,i0)') 'ERROR: fill value applied inconsistently for column ', idx
            call hist_add_error('buff_2d_check_fill', errstr, errors=logger)
            exit
         end if
      end do

   end subroutine buff_2d_check_fill

   !#######################################################################

   pure subroutine buff_2d_value(this, norm_val, flag_xyfill, fill_value)
      ! Dummy arguments
      class(hist_buff_2d_t),             intent(in)    :: this
      real(REAL64),                      intent(inout) :: norm_val(:,:)
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      ! Local variable
      integer :: ind1, ind2
      integer :: nacc
      real(REAL64) :: tmpfill, variance

      ! Initialize values that are not fill values:
      tmpfill = merge(fill_value, 0.0_REAL64, flag_xyfill)
      if (this%accum_type /= hist_accum_var) then
         do ind1 = 1, this%field_shape(1)
            do ind2 = 1, this%field_shape(2)
               nacc = this%num_samples(ind1)
               if (nacc > 0) then
                  norm_val(ind1, ind2) = this%data(ind1, ind2) / nacc
               else
                  norm_val(ind1, ind2) = tmpfill
               end if
            end do
         end do
      else
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = 1, this%field_shape(1)
            do ind2 = 1, this%field_shape(2)
               if (this%num_samples(ind1) > 0) then
                  variance = this%var_buffer(ind1, ind2) / this%num_samples(ind1)
                  norm_val(ind1, ind2) = sqrt(variance)
               else
                  norm_val(ind1, ind2) = tmpfill
               end if
            end do
         end do
      end if

   end subroutine buff_2d_value

   !#######################################################################

   subroutine buff_3d_clear(this, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_alloc_error, error_length

      ! Dummy arguments
      class(hist_buff_3d_t),       intent(inout) :: this
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer                     :: aerr
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'buff_3d_clear'

      if (.not. allocated(this%data)) then
         allocate(this%data(this%field_shape(1), this%field_shape(2), this%field_shape(3)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('data, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%var_buffer) .and. this%accum_type == hist_accum_var) then
         allocate(this%var_buffer(this%field_shape(1), this%field_shape(2), this%field_shape(3)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('var_buffer, error: '//errmsg, __FILE__, __LINE__ - 2,  &
                 subname=subname, errors=logger)
         end if
      end if
      if (.not. allocated(this%num_samples)) then
         allocate(this%num_samples(this%field_shape(1)), stat=aerr, errmsg=errmsg)
         if (aerr /= 0) then
            call hist_add_alloc_error('num_samples, error: '//errmsg, __FILE__, __LINE__ - 1, &
                 subname=subname, errors=logger)
         end if
      end if
      this%data = 0.0_REAL64
      this%num_samples = 0
      if (allocated(this%var_buffer)) then
         this%var_buffer = 0.0_REAL64
      end if
#ifdef _OPENACC
      ! First clear: create the device copies; every clear: zero them too
      if (.not. this%on_device) then
         call dev_create(this%data)
         call dev_create(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_create(this%var_buffer)
         end if
         this%on_device = .true.
      end if
      call dev_zero(this%data)
      call dev_zero(this%num_samples)
      if (allocated(this%var_buffer)) then
         call dev_zero(this%var_buffer)
      end if
#endif

   end subroutine buff_3d_clear

   !#######################################################################

   subroutine init_buff_3d(this, field_in, volume_in, horiz_axis_in, &
        accum_type_in, shape_in, block_sizes_in, block_ind_in, logger)
      use hist_msg_handler, only: hist_log_messages

      class(hist_buff_3d_t),             intent(inout) :: this
      class(hist_hashable_t),  pointer,  intent(in)    :: field_in
      integer,                           intent(in)    :: volume_in
      integer,                           intent(in)    :: horiz_axis_in
      integer,                           intent(in)    :: accum_type_in
      integer,                           intent(in)    :: shape_in(:)
      integer,                 optional, intent(in)    :: block_sizes_in(:)
      integer,                 optional, intent(in)    :: block_ind_in
      type(hist_log_messages), optional, intent(inout) :: logger

      call init_buffer(this, field_in, volume_in, horiz_axis_in, accum_type_in, &
              shape_in, block_sizes_in, block_ind_in, logger=logger)
      call this%clear(logger=logger)
      this%buff_type = 'buff_3d'

   end subroutine init_buff_3d

   !#######################################################################

   subroutine buff_3d_accum(this, field, cols_or_block, flag_xyfill, &
         fill_value, cole, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
      ! Dummy arguments
      class(hist_buff_3d_t),       intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:,:)
      integer,                           intent(in)    :: cols_or_block
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      integer,                 optional, intent(in)    :: cole
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      integer            :: col_beg_use
      integer            :: col_end_use
      integer            :: ind1, ind2, ind3
      real(REAL64)       :: fld_val, tmp
      character(len=error_length) :: errmsg
      character(len=*), parameter :: subname = 'buff_3d_accum'

      if (this%has_blocks()) then
         ! For a blocked field, <cols_or_block> is a block index
         col_beg_use = this%block_begs(cols_or_block)
         col_end_use = this%block_ends(cols_or_block)
      else
         ! Non blocked, <cols_or_block> is the first column index
         col_beg_use = cols_or_block
         if (present(cole)) then
            col_end_use = cole
         else
            col_end_use = col_beg_use +                                       &
                 this%field_shape(this%horiz_axis_ind) - 1
         end if
      end if

#ifdef _OPENACC
      if (this%on_device) then
         call buff_3d_accum_device(this, field, col_beg_use, col_end_use, &
              flag_xyfill, fill_value, logger)
         return
      end if
#endif

      select case (this%accum_type)
      case (hist_accum_lst)
         do ind1 = col_beg_use, col_end_use
            ! Only set samples for this column if not a fill value
            if (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1, 1) == fill_value) then
               this%num_samples(ind1) = 0
            else
               this%num_samples(ind1) = 1
            end if
            do ind2 = 1, this%field_shape(2)
               do ind3 = 1, this%field_shape(3)
                  fld_val = field(ind1 - col_beg_use + 1, ind2, ind3)
                  if (flag_xyfill .and. fld_val == fill_value) then
                     this%data(ind1, ind2, ind3) = fill_value
                  else
                     this%data(ind1, ind2, ind3) = fld_val
                  end if
               end do
            end do
         end do
         if (flag_xyfill) then
            ! Check if assumption that columns are consistent wrt fill value was accurate
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_min)
         do ind1 = col_beg_use, col_end_use
            ! Has this column been sampled yet?
            if (this%num_samples(ind1) == 0) then
               ! First sample for this column: initialize all levels.
               ! Skip fill value columns (num_samples stays 0;
               ! buff_3d_value will write fill_value at output time)
               if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1, 1) == fill_value)) then
                  do ind2 = 1, this%field_shape(2)
                     do ind3 = 1, this%field_shape(3)
                        this%data(ind1, ind2, ind3) = field(ind1 - col_beg_use + 1, ind2, ind3)
                     end do
                  end do
                  this%num_samples(ind1) = 1
               end if
            else
               ! Subsequent samples in column: update minimum per level
               do ind2 = 1, this%field_shape(2)
                  do ind3 = 1, this%field_shape(3)
                     fld_val = field(ind1 - col_beg_use + 1, ind2, ind3)
                     if (flag_xyfill) then
                        if (fld_val < this%data(ind1, ind2, ind3) .and. &
                            fld_val /= fill_value) then
                           this%data(ind1, ind2, ind3) = fld_val
                        end if
                     else
                        if (fld_val < this%data(ind1, ind2, ind3)) then
                           this%data(ind1, ind2, ind3) = fld_val
                        end if
                     end if
                  end do
               end do
            end if
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_max)
         do ind1 = col_beg_use, col_end_use
            ! Has this column been sampled yet?
            if (this%num_samples(ind1) == 0) then
               ! First sample for this column: initialize all levels.
               ! Skip fill value columns (num_samples stays 0;
               ! buff_3d_value will write fill_value at output time)
               if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1, 1) == fill_value)) then
                  do ind2 = 1, this%field_shape(2)
                     do ind3 = 1, this%field_shape(3)
                        this%data(ind1, ind2, ind3) = field(ind1 - col_beg_use + 1, ind2, ind3)
                     end do
                  end do
                  this%num_samples(ind1) = 1
               end if
            else
               ! Subsequent samples: update maximum per level
               do ind2 = 1, this%field_shape(2)
                  do ind3 = 1, this%field_shape(3)
                     fld_val = field(ind1 - col_beg_use + 1, ind2, ind3)
                     if (flag_xyfill) then
                        if (fld_val > this%data(ind1, ind2, ind3) .and. fld_val /= fill_value) then
                           this%data(ind1, ind2, ind3) = fld_val
                        end if
                     else
                        if (fld_val > this%data(ind1, ind2, ind3)) then
                           this%data(ind1, ind2, ind3) = fld_val
                        end if
                     end if
                  end do
               end do
            end if
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_avg)
         do ind1 = col_beg_use, col_end_use
            if (.not. (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1, 1) == fill_value)) then
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
            do ind2 = 1, this%field_shape(2)
               do ind3 = 1, this%field_shape(3)
                  fld_val = field(ind1 - col_beg_use + 1, ind2, ind3)
                  ! Compute running sum
                  if (flag_xyfill) then
                     ! Only include sample if it is not the fill value
                     if (fld_val /= fill_value) then
                        this%data(ind1, ind2, ind3) = this%data(ind1, ind2, ind3) + fld_val
                     end if
                  else
                     this%data(ind1, ind2, ind3) = this%data(ind1, ind2, ind3) + fld_val
                  end if
               end do
            end do
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case (hist_accum_var)
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = col_beg_use, col_end_use
            if (flag_xyfill .and. field(ind1 - col_beg_use + 1, 1, 1) == fill_value) then
               this%num_samples(ind1) = this%num_samples(ind1)
            else
               this%num_samples(ind1) = this%num_samples(ind1) + 1
            end if
            do ind2 = 1, this%field_shape(2)
               do ind3 = 1, this%field_shape(3)
                  fld_val = field(ind1 - col_beg_use + 1, ind2, ind3)
                  if (flag_xyfill) then
                     if (fld_val /= fill_value) then
                        ! Only include the sample if it's not a fill value
                        if (this%num_samples(ind1) == 1) then
                           this%data(ind1, ind2, ind3) = fld_val
                           this%var_buffer(ind1, ind2, ind3) = 0.0_REAL64
                        else
                           tmp = this%data(ind1, ind2, ind3)
                           this%data(ind1, ind2, ind3) = this%data(ind1, ind2, ind3) + &
                                   (fld_val - this%data(ind1, ind2, ind3)) / &
                                   this%num_samples(ind1)
                           this%var_buffer(ind1, ind2, ind3) = this%var_buffer(ind1, ind2, ind3) + &
                                   (fld_val - this%data(ind1, ind2, ind3)) * &
                                   (fld_val - tmp)
                        end if
                     end if
                  else
                     if (this%num_samples(ind1) == 1) then
                        this%data(ind1, ind2, ind3) = fld_val
                        this%var_buffer(ind1, ind2, ind3) = 0.0_REAL64
                     else
                        tmp = this%data(ind1, ind2, ind3)
                        this%data(ind1, ind2, ind3) = this%data(ind1, ind2, ind3) + &
                                (fld_val - this%data(ind1, ind2, ind3)) / &
                                this%num_samples(ind1)
                        this%var_buffer(ind1, ind2, ind3) = this%var_buffer(ind1, ind2, ind3) + &
                                (fld_val - this%data(ind1, ind2, ind3)) * &
                                (fld_val - tmp)
                     end if
                  end if
               end do
            end do
         end do
         if (flag_xyfill) then
            call this%check_fill_value(field, fill_value, logger)
         end if
      case default
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error(subname, errmsg, errors=logger)
      end select

   end subroutine buff_3d_accum

   !#######################################################################

   subroutine buff_3d_check_fill(this, field, fill_value, logger)
      use hist_msg_handler, only: hist_log_messages, hist_add_error
      ! Ensure that columns have fill value applied consistently across levels
      class(hist_buff_3d_t),             intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:,:)
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      character(len=512) :: errstr
      integer :: idx
      logical :: col_is_fill

      ! Every element in a column must match the fill status of its first element
      do idx = 1, size(field, 1)
         col_is_fill = (field(idx,1,1) == fill_value)
         if (any((field(idx,:,:) == fill_value) .neqv. col_is_fill)) then
            write(errstr, '(a,i0)') 'ERROR: fill value applied inconsistently for column ', idx
            call hist_add_error('buff_3d_check_fill', errstr, errors=logger)
            exit
         end if
      end do

   end subroutine buff_3d_check_fill

   !#######################################################################

   pure subroutine buff_3d_value(this, norm_val, flag_xyfill, fill_value)
      ! Dummy arguments
      class(hist_buff_3d_t),             intent(in)    :: this
      real(REAL64),                      intent(inout) :: norm_val(:,:,:)
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      ! Local variable
      integer :: ind1, ind2, ind3
      integer :: nacc
      real(REAL64) :: tmpfill, variance

      ! Initialize values that are not fill values:
      tmpfill = merge(fill_value, 0.0_REAL64, flag_xyfill)
      if (this%accum_type /= hist_accum_var) then
         do ind1 = 1, this%field_shape(1)
            do ind2 = 1, this%field_shape(2)
               do ind3 = 1, this%field_shape(3)
                  nacc = this%num_samples(ind1)
                  if (nacc > 0) then
                     norm_val(ind1, ind2, ind3) = this%data(ind1, ind2, ind3) / nacc
                  else
                     norm_val(ind1, ind2, ind3) = tmpfill
                  end if
               end do
            end do
         end do
      else
         ! Standard deviation using Welford's algorithm
         ! DOI: 10.1080/00401706.1962.10490022
         do ind1 = 1, this%field_shape(1)
            do ind2 = 1, this%field_shape(2)
               do ind3 = 1, this%field_shape(3)
                  if (this%num_samples(ind1) > 0) then
                     variance = this%var_buffer(ind1, ind2, ind3) / this%num_samples(ind1)
                     norm_val(ind1, ind2, ind3) = sqrt(variance)
                  else
                     norm_val(ind1, ind2, ind3) = tmpfill
                  end if
               end do
            end do
         end do
      end if

   end subroutine buff_3d_value

   !#######################################################################

   subroutine get_var_buffer_1d(this, var_buff)
      class(hist_buff_1d_t), intent(in) :: this
      real(REAL64), intent(out) :: var_buff(:)

      var_buff = this%var_buffer
   end subroutine get_var_buffer_1d

   !#######################################################################

   subroutine get_var_buffer_2d(this, var_buff)
      class(hist_buff_2d_t), intent(in) :: this
      real(REAL64), intent(out) :: var_buff(:,:)

      var_buff = this%var_buffer
   end subroutine get_var_buffer_2d

   !#######################################################################

   subroutine get_var_buffer_3d(this, var_buff)
      class(hist_buff_3d_t), intent(in) :: this
      real(REAL64), intent(out) :: var_buff(:,:,:)

      var_buff = this%var_buffer
   end subroutine get_var_buffer_3d

   !#######################################################################

   subroutine set_field_buffer_1d(this, field_buff)
      class(hist_buff_1d_t), intent(inout) :: this
      real(REAL64),             intent(in) :: field_buff(:)

      this%data = field_buff
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%data)
      end if
#endif
   end subroutine set_field_buffer_1d

   !#######################################################################

   subroutine set_field_buffer_2d(this, field_buff)
      class(hist_buff_2d_t), intent(inout) :: this
      real(REAL64),             intent(in) :: field_buff(:,:)

      this%data = field_buff
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%data)
      end if
#endif
   end subroutine set_field_buffer_2d

   !#######################################################################

   subroutine set_field_buffer_3d(this, field_buff)
      class(hist_buff_3d_t), intent(inout) :: this
      real(REAL64),             intent(in) :: field_buff(:,:,:)

      this%data = field_buff
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%data)
      end if
#endif
   end subroutine set_field_buffer_3d

   !#######################################################################

   subroutine set_var_buffer_1d(this, var_buff)
      class(hist_buff_1d_t), intent(inout) :: this
      real(REAL64),             intent(in) :: var_buff(:)

      this%var_buffer = var_buff
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%var_buffer)
      end if
#endif
   end subroutine set_var_buffer_1d

   !#######################################################################

   subroutine set_var_buffer_2d(this, var_buff)
      class(hist_buff_2d_t), intent(inout) :: this
      real(REAL64),             intent(in) :: var_buff(:,:)

      this%var_buffer = var_buff
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%var_buffer)
      end if
#endif
   end subroutine set_var_buffer_2d

   !#######################################################################

   subroutine set_var_buffer_3d(this, var_buff)
      class(hist_buff_3d_t), intent(inout) :: this
      real(REAL64),             intent(in) :: var_buff(:,:,:)

      this%var_buffer = var_buff
#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_device(this%var_buffer)
      end if
#endif
   end subroutine set_var_buffer_3d

   !#######################################################################

   subroutine buff_1d_sync_host(this)
      ! Bring the host copies of the buffer up to date (no-op unless the
      ! buffer is device-resident). Called before any host read.
      class(hist_buff_1d_t), intent(inout) :: this

#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_self(this%data)
         call dev_update_self(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_update_self(this%var_buffer)
         end if
      end if
#endif
   end subroutine buff_1d_sync_host

   !#######################################################################

   subroutine buff_1d_device_release(this)
      ! Delete the device copies (before the buffer is deallocated)
      class(hist_buff_1d_t), intent(inout) :: this

#ifdef _OPENACC
      if (this%on_device) then
         call dev_delete(this%data)
         call dev_delete(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_delete(this%var_buffer)
         end if
         this%on_device = .false.
      end if
#endif
   end subroutine buff_1d_device_release

   !#######################################################################

   subroutine buff_2d_sync_host(this)
      ! Bring the host copies of the buffer up to date (no-op unless the
      ! buffer is device-resident). Called before any host read.
      class(hist_buff_2d_t), intent(inout) :: this

#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_self(this%data)
         call dev_update_self(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_update_self(this%var_buffer)
         end if
      end if
#endif
   end subroutine buff_2d_sync_host

   !#######################################################################

   subroutine buff_2d_device_release(this)
      ! Delete the device copies (before the buffer is deallocated)
      class(hist_buff_2d_t), intent(inout) :: this

#ifdef _OPENACC
      if (this%on_device) then
         call dev_delete(this%data)
         call dev_delete(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_delete(this%var_buffer)
         end if
         this%on_device = .false.
      end if
#endif
   end subroutine buff_2d_device_release

   !#######################################################################

   subroutine buff_3d_sync_host(this)
      ! Bring the host copies of the buffer up to date (no-op unless the
      ! buffer is device-resident). Called before any host read.
      class(hist_buff_3d_t), intent(inout) :: this

#ifdef _OPENACC
      if (this%on_device) then
         call dev_update_self(this%data)
         call dev_update_self(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_update_self(this%var_buffer)
         end if
      end if
#endif
   end subroutine buff_3d_sync_host

   !#######################################################################

   subroutine buff_3d_device_release(this)
      ! Delete the device copies (before the buffer is deallocated)
      class(hist_buff_3d_t), intent(inout) :: this

#ifdef _OPENACC
      if (this%on_device) then
         call dev_delete(this%data)
         call dev_delete(this%num_samples)
         if (allocated(this%var_buffer)) then
            call dev_delete(this%var_buffer)
         end if
         this%on_device = .false.
      end if
#endif
   end subroutine buff_3d_device_release

#ifdef _OPENACC
   !#######################################################################

   subroutine buff_1d_accum_device(this, field, col_beg, col_end,           &
        flag_xyfill, fill_value, logger)
      ! Device version of buff_1d_accum: same arithmetic and order. A host
      ! <field> is copied to the device for the call.
      use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
      ! Dummy arguments
      class(hist_buff_1d_t),            intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:)
      integer,                           intent(in)    :: col_beg
      integer,                           intent(in)    :: col_end
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      logical                     :: field_on_device
      integer                     :: bad_col
      character(len=error_length) :: errmsg
      character(len=512)          :: errstr

      if ((this%accum_type < hist_accum_lst) .or.                             &
           (this%accum_type > hist_accum_var)) then
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error('buff_1d_accum', errmsg, errors=logger)
         return
      end if
      field_on_device = acc_is_present(field)
      if (.not. field_on_device) then
         !$acc enter data copyin(field)
      end if
      if (this%accum_type == hist_accum_var) then
         call accum_var_1d_dev(this%data, this%var_buffer, this%num_samples, &
              field, col_beg, col_end, flag_xyfill, fill_value)
      else
         call accum_1d_dev(this%accum_type, this%data, this%num_samples,     &
              field, col_beg, col_end, flag_xyfill, fill_value)
      end if
      if (.not. field_on_device) then
         !$acc exit data delete(field)
      end if

   end subroutine buff_1d_accum_device

   !#######################################################################

   subroutine buff_2d_accum_device(this, field, col_beg, col_end,           &
        flag_xyfill, fill_value, logger)
      ! Device version of buff_2d_accum: same arithmetic and order. A host
      ! <field> is copied to the device for the call.
      use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
      ! Dummy arguments
      class(hist_buff_2d_t),            intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:)
      integer,                           intent(in)    :: col_beg
      integer,                           intent(in)    :: col_end
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      logical                     :: field_on_device
      integer                     :: bad_col
      character(len=error_length) :: errmsg
      character(len=512)          :: errstr

      if ((this%accum_type < hist_accum_lst) .or.                             &
           (this%accum_type > hist_accum_var)) then
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error('buff_2d_accum', errmsg, errors=logger)
         return
      end if
      field_on_device = acc_is_present(field)
      if (.not. field_on_device) then
         !$acc enter data copyin(field)
      end if
      if (this%accum_type == hist_accum_var) then
         call accum_var_2d_dev(this%data, this%var_buffer, this%num_samples, &
              field, col_beg, col_end, this%field_shape(2), flag_xyfill, fill_value)
      else
         call accum_2d_dev(this%accum_type, this%data, this%num_samples,     &
              field, col_beg, col_end, this%field_shape(2), flag_xyfill, fill_value)
      end if
      if (flag_xyfill) then
         ! Same check as buff_2d_check_fill, as a device reduction
         call check_fill_2d_dev(field, fill_value, bad_col)
         if (bad_col > 0) then
            write(errstr, '(a,i0)') 'ERROR: fill value applied inconsistently for column ', bad_col
            call hist_add_error('buff_2d_check_fill', errstr, errors=logger)
         end if
      end if
      if (.not. field_on_device) then
         !$acc exit data delete(field)
      end if

   end subroutine buff_2d_accum_device

   !#######################################################################

   subroutine buff_3d_accum_device(this, field, col_beg, col_end,           &
        flag_xyfill, fill_value, logger)
      ! Device version of buff_3d_accum: same arithmetic and order. A host
      ! <field> is copied to the device for the call.
      use hist_msg_handler, only: hist_log_messages, hist_add_error, error_length
      ! Dummy arguments
      class(hist_buff_3d_t),            intent(inout) :: this
      real(REAL64),                      intent(in)    :: field(:,:,:)
      integer,                           intent(in)    :: col_beg
      integer,                           intent(in)    :: col_end
      logical,                           intent(in)    :: flag_xyfill
      real(REAL64),                      intent(in)    :: fill_value
      type(hist_log_messages), optional, intent(inout) :: logger
      ! Local variables
      logical                     :: field_on_device
      integer                     :: bad_col
      character(len=error_length) :: errmsg
      character(len=512)          :: errstr

      if ((this%accum_type < hist_accum_lst) .or.                             &
           (this%accum_type > hist_accum_var)) then
         ! This shouldn't trigger, so raise an error:
         write(errmsg, *) 'un-recognized accumulation type: ', this%accum_type
         call hist_add_error('buff_3d_accum', errmsg, errors=logger)
         return
      end if
      field_on_device = acc_is_present(field)
      if (.not. field_on_device) then
         !$acc enter data copyin(field)
      end if
      if (this%accum_type == hist_accum_var) then
         call accum_var_3d_dev(this%data, this%var_buffer, this%num_samples, &
              field, col_beg, col_end, this%field_shape(2), this%field_shape(3), flag_xyfill, fill_value)
      else
         call accum_3d_dev(this%accum_type, this%data, this%num_samples,     &
              field, col_beg, col_end, this%field_shape(2), this%field_shape(3), flag_xyfill, fill_value)
      end if
      if (flag_xyfill) then
         ! Same check as buff_3d_check_fill, as a device reduction
         call check_fill_3d_dev(field, fill_value, bad_col)
         if (bad_col > 0) then
            write(errstr, '(a,i0)') 'ERROR: fill value applied inconsistently for column ', bad_col
            call hist_add_error('buff_3d_check_fill', errstr, errors=logger)
         end if
      end if
      if (.not. field_on_device) then
         !$acc exit data delete(field)
      end if

   end subroutine buff_3d_accum_device

   !#######################################################################

   subroutine accum_1d_dev(accum_type, data, num_samples, field, col_beg,  &
        col_end, flag_xyfill, fill_value)
      ! lst/min/max/avg accumulation on the device (see buff_1d_accum)
      integer,      intent(in)    :: accum_type
      real(REAL64), intent(inout) :: data(:)
      integer,      intent(inout) :: num_samples(:)
      real(REAL64), intent(in)    :: field(:)
      integer,      intent(in)    :: col_beg
      integer,      intent(in)    :: col_end
      logical,      intent(in)    :: flag_xyfill
      real(REAL64), intent(in)    :: fill_value
      integer      :: i, ic
      real(REAL64) :: fld_val

      select case (accum_type)
      case (hist_accum_lst)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            fld_val = field(ic)
            if (flag_xyfill .and. fld_val == fill_value) then
               data(i) = fill_value
               num_samples(i) = 0
            else
               data(i) = fld_val
               num_samples(i) = 1
            end if
         end do
      case (hist_accum_min)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            fld_val = field(ic)
            if (flag_xyfill) then
               if (num_samples(i) == 0) then
                  if (fld_val /= fill_value) then
                     data(i) = fld_val
                     num_samples(i) = 1
                  else
                     data(i) = HUGE(REAL64)
                  end if
               else if (fld_val < data(i) .and. fld_val /= fill_value) then
                  data(i) = fld_val
                  num_samples(i) = 1
               end if
            else
               if (num_samples(i) == 0 .or. fld_val < data(i)) then
                  data(i) = fld_val
                  num_samples(i) = 1
               end if
            end if
         end do
      case (hist_accum_max)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            fld_val = field(ic)
            if (flag_xyfill) then
               if (num_samples(i) == 0) then
                  if (fld_val /= fill_value) then
                     data(i) = fld_val
                     num_samples(i) = 1
                  else
                     data(i) = -HUGE(REAL64)
                  end if
               else if (fld_val > data(i) .and. fld_val /= fill_value) then
                  data(i) = fld_val
                  num_samples(i) = 1
               end if
            else
               if (num_samples(i) == 0 .or. fld_val > data(i)) then
                  data(i) = fld_val
                  num_samples(i) = 1
               end if
            end if
         end do
      case (hist_accum_avg)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            fld_val = field(ic)
            if (flag_xyfill) then
               if (fld_val /= fill_value) then
                  data(i) = data(i) + fld_val
                  num_samples(i) = num_samples(i) + 1
               end if
            else
               data(i) = data(i) + fld_val
               num_samples(i) = num_samples(i) + 1
            end if
         end do
      end select

   end subroutine accum_1d_dev

   !#######################################################################

   subroutine accum_2d_dev(accum_type, data, num_samples, field, col_beg,  &
        col_end, n2, flag_xyfill, fill_value)
      ! lst/min/max/avg accumulation on the device (see buff_2d_accum)
      integer,      intent(in)    :: accum_type
      real(REAL64), intent(inout) :: data(:,:)
      integer,      intent(inout) :: num_samples(:)
      real(REAL64), intent(in)    :: field(:,:)
      integer,      intent(in)    :: col_beg
      integer,      intent(in)    :: col_end
      integer,      intent(in)    :: n2
      logical,      intent(in)    :: flag_xyfill
      real(REAL64), intent(in)    :: fill_value
      integer      :: i, ic, k
      real(REAL64) :: fld_val

      select case (accum_type)
      case (hist_accum_lst)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (flag_xyfill .and. field(ic, 1) == fill_value) then
               num_samples(i) = 0
            else
               num_samples(i) = 1
            end if
            !$acc loop seq
            do k = 1, n2
               fld_val = field(ic, k)
               if (flag_xyfill .and. fld_val == fill_value) then
                  data(i, k) = fill_value
               else
                  data(i, k) = fld_val
               end if
            end do
         end do
      case (hist_accum_min)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (num_samples(i) == 0) then
               if (.not. (flag_xyfill .and. field(ic, 1) == fill_value)) then
                  !$acc loop seq
                  do k = 1, n2
                     data(i, k) = field(ic, k)
                  end do
                  num_samples(i) = 1
               end if
            else
               !$acc loop seq
               do k = 1, n2
                  fld_val = field(ic, k)
                  if (flag_xyfill) then
                     if (fld_val < data(i, k) .and. fld_val /= fill_value) then
                        data(i, k) = fld_val
                     end if
                  else
                     if (fld_val < data(i, k)) then
                        data(i, k) = fld_val
                     end if
                  end if
               end do
            end if
         end do
      case (hist_accum_max)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (num_samples(i) == 0) then
               if (.not. (flag_xyfill .and. field(ic, 1) == fill_value)) then
                  !$acc loop seq
                  do k = 1, n2
                     data(i, k) = field(ic, k)
                  end do
                  num_samples(i) = 1
               end if
            else
               !$acc loop seq
               do k = 1, n2
                  fld_val = field(ic, k)
                  if (flag_xyfill) then
                     if (fld_val > data(i, k) .and. fld_val /= fill_value) then
                        data(i, k) = fld_val
                     end if
                  else
                     if (fld_val > data(i, k)) then
                        data(i, k) = fld_val
                     end if
                  end if
               end do
            end if
         end do
      case (hist_accum_avg)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (.not. (flag_xyfill .and. field(ic, 1) == fill_value)) then
               num_samples(i) = num_samples(i) + 1
            end if
            !$acc loop seq
            do k = 1, n2
               fld_val = field(ic, k)
               if (flag_xyfill) then
                  if (fld_val /= fill_value) then
                     data(i, k) = data(i, k) + fld_val
                  end if
               else
                  data(i, k) = data(i, k) + fld_val
               end if
            end do
         end do
      end select

   end subroutine accum_2d_dev

   !#######################################################################

   subroutine accum_3d_dev(accum_type, data, num_samples, field, col_beg,  &
        col_end, n2, n3, flag_xyfill, fill_value)
      ! lst/min/max/avg accumulation on the device (see buff_3d_accum)
      integer,      intent(in)    :: accum_type
      real(REAL64), intent(inout) :: data(:,:,:)
      integer,      intent(inout) :: num_samples(:)
      real(REAL64), intent(in)    :: field(:,:,:)
      integer,      intent(in)    :: col_beg
      integer,      intent(in)    :: col_end
      integer,      intent(in)    :: n2, n3
      logical,      intent(in)    :: flag_xyfill
      real(REAL64), intent(in)    :: fill_value
      integer      :: i, ic, k, m
      real(REAL64) :: fld_val

      select case (accum_type)
      case (hist_accum_lst)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (flag_xyfill .and. field(ic, 1, 1) == fill_value) then
               num_samples(i) = 0
            else
               num_samples(i) = 1
            end if
            !$acc loop seq
            do k = 1, n2
               !$acc loop seq
               do m = 1, n3
                  fld_val = field(ic, k, m)
                  if (flag_xyfill .and. fld_val == fill_value) then
                     data(i, k, m) = fill_value
                  else
                     data(i, k, m) = fld_val
                  end if
               end do
            end do
         end do
      case (hist_accum_min)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (num_samples(i) == 0) then
               if (.not. (flag_xyfill .and. field(ic, 1, 1) == fill_value)) then
                  !$acc loop seq
                  do k = 1, n2
                     !$acc loop seq
                     do m = 1, n3
                        data(i, k, m) = field(ic, k, m)
                     end do
                  end do
                  num_samples(i) = 1
               end if
            else
               !$acc loop seq
               do k = 1, n2
                  !$acc loop seq
                  do m = 1, n3
                     fld_val = field(ic, k, m)
                     if (flag_xyfill) then
                        if (fld_val < data(i, k, m) .and. fld_val /= fill_value) then
                           data(i, k, m) = fld_val
                        end if
                     else
                        if (fld_val < data(i, k, m)) then
                           data(i, k, m) = fld_val
                        end if
                     end if
                  end do
               end do
            end if
         end do
      case (hist_accum_max)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (num_samples(i) == 0) then
               if (.not. (flag_xyfill .and. field(ic, 1, 1) == fill_value)) then
                  !$acc loop seq
                  do k = 1, n2
                     !$acc loop seq
                     do m = 1, n3
                        data(i, k, m) = field(ic, k, m)
                     end do
                  end do
                  num_samples(i) = 1
               end if
            else
               !$acc loop seq
               do k = 1, n2
                  !$acc loop seq
                  do m = 1, n3
                     fld_val = field(ic, k, m)
                     if (flag_xyfill) then
                        if (fld_val > data(i, k, m) .and. fld_val /= fill_value) then
                           data(i, k, m) = fld_val
                        end if
                     else
                        if (fld_val > data(i, k, m)) then
                           data(i, k, m) = fld_val
                        end if
                     end if
                  end do
               end do
            end if
         end do
      case (hist_accum_avg)
         !$acc parallel loop gang vector private(ic, fld_val) &
         !$acc   present(data, num_samples, field)
         do i = col_beg, col_end
            ic = i - col_beg + 1
            if (.not. (flag_xyfill .and. field(ic, 1, 1) == fill_value)) then
               num_samples(i) = num_samples(i) + 1
            end if
            !$acc loop seq
            do k = 1, n2
               !$acc loop seq
               do m = 1, n3
                  fld_val = field(ic, k, m)
                  if (flag_xyfill) then
                     if (fld_val /= fill_value) then
                        data(i, k, m) = data(i, k, m) + fld_val
                     end if
                  else
                     data(i, k, m) = data(i, k, m) + fld_val
                  end if
               end do
            end do
         end do
      end select

   end subroutine accum_3d_dev

   !#######################################################################

   subroutine accum_var_1d_dev(data, var_buffer, num_samples, field,       &
        col_beg, col_end, flag_xyfill, fill_value)
      ! Standard deviation (Welford) on the device (see buff_1d_accum)
      real(REAL64), intent(inout) :: data(:)
      real(REAL64), intent(inout) :: var_buffer(:)
      integer,      intent(inout) :: num_samples(:)
      real(REAL64), intent(in)    :: field(:)
      integer,      intent(in)    :: col_beg
      integer,      intent(in)    :: col_end
      logical,      intent(in)    :: flag_xyfill
      real(REAL64), intent(in)    :: fill_value
      integer      :: i, ic
      real(REAL64) :: fld_val, tmp

      !$acc parallel loop gang vector private(ic, fld_val, tmp) &
      !$acc   present(data, var_buffer, num_samples, field)
      do i = col_beg, col_end
         ic = i - col_beg + 1
         fld_val = field(ic)
         if (flag_xyfill) then
            if (fld_val /= fill_value) then
               if (num_samples(i) == 0) then
                  data(i) = fld_val
                  var_buffer(i) = 0.0_REAL64
                  num_samples(i) = num_samples(i) + 1
               else
                  tmp = data(i)
                  num_samples(i) = num_samples(i) + 1
                  data(i) = data(i) + &
                          (fld_val - data(i)) / &
                          num_samples(i)
                  var_buffer(i) = var_buffer(i) + &
                          (fld_val - data(i)) * &
                          (fld_val - tmp)
               end if
            end if
         else
            if (num_samples(i) == 0) then
               data(i) = fld_val
               var_buffer(i) = 0.0_REAL64
               num_samples(i) = num_samples(i) + 1
            else
               tmp = data(i)
               num_samples(i) = num_samples(i) + 1
               data(i) = data(i) + &
                       (fld_val - data(i)) / &
                       num_samples(i)
               var_buffer(i) = var_buffer(i) + &
                       (fld_val - data(i)) * &
                       (fld_val - tmp)
            end if
         end if
      end do

   end subroutine accum_var_1d_dev

   !#######################################################################

   subroutine accum_var_2d_dev(data, var_buffer, num_samples, field,       &
        col_beg, col_end, n2, flag_xyfill, fill_value)
      ! Standard deviation (Welford) on the device (see buff_2d_accum)
      real(REAL64), intent(inout) :: data(:,:)
      real(REAL64), intent(inout) :: var_buffer(:,:)
      integer,      intent(inout) :: num_samples(:)
      real(REAL64), intent(in)    :: field(:,:)
      integer,      intent(in)    :: col_beg
      integer,      intent(in)    :: col_end
      integer,      intent(in)    :: n2
      logical,      intent(in)    :: flag_xyfill
      real(REAL64), intent(in)    :: fill_value
      integer      :: i, ic, k
      real(REAL64) :: fld_val, tmp

      !$acc parallel loop gang vector private(ic, fld_val, tmp) &
      !$acc   present(data, var_buffer, num_samples, field)
      do i = col_beg, col_end
         ic = i - col_beg + 1
         if (.not. (flag_xyfill .and. field(ic, 1) == fill_value)) then
            num_samples(i) = num_samples(i) + 1
         end if
         !$acc loop seq
         do k = 1, n2
            fld_val = field(ic, k)
            if (flag_xyfill) then
               if (fld_val /= fill_value) then
                  if (num_samples(i) == 1) then
                     data(i, k) = fld_val
                     var_buffer(i, k) = 0.0_REAL64
                  else
                     tmp = data(i, k)
                     data(i, k) = data(i, k) + &
                             (fld_val - data(i, k)) / &
                             num_samples(i)
                     var_buffer(i, k) = var_buffer(i, k) + &
                             (fld_val - data(i, k)) * &
                             (fld_val - tmp)
                  end if
               end if
            else
               if (num_samples(i) == 1) then
                  data(i, k) = fld_val
                  var_buffer(i, k) = 0.0_REAL64
               else
                  tmp = data(i, k)
                  data(i, k) = data(i, k) + &
                          (fld_val - data(i, k)) / &
                          num_samples(i)
                  var_buffer(i, k) = var_buffer(i, k) + &
                          (fld_val - data(i, k)) * &
                          (fld_val - tmp)
               end if
            end if
         end do
      end do

   end subroutine accum_var_2d_dev

   !#######################################################################

   subroutine accum_var_3d_dev(data, var_buffer, num_samples, field,       &
        col_beg, col_end, n2, n3, flag_xyfill, fill_value)
      ! Standard deviation (Welford) on the device (see buff_3d_accum)
      real(REAL64), intent(inout) :: data(:,:,:)
      real(REAL64), intent(inout) :: var_buffer(:,:,:)
      integer,      intent(inout) :: num_samples(:)
      real(REAL64), intent(in)    :: field(:,:,:)
      integer,      intent(in)    :: col_beg
      integer,      intent(in)    :: col_end
      integer,      intent(in)    :: n2, n3
      logical,      intent(in)    :: flag_xyfill
      real(REAL64), intent(in)    :: fill_value
      integer      :: i, ic, k, m
      real(REAL64) :: fld_val, tmp

      !$acc parallel loop gang vector private(ic, fld_val, tmp) &
      !$acc   present(data, var_buffer, num_samples, field)
      do i = col_beg, col_end
         ic = i - col_beg + 1
         if (.not. (flag_xyfill .and. field(ic, 1, 1) == fill_value)) then
            num_samples(i) = num_samples(i) + 1
         end if
         !$acc loop seq
         do k = 1, n2
            !$acc loop seq
            do m = 1, n3
               fld_val = field(ic, k, m)
               if (flag_xyfill) then
                  if (fld_val /= fill_value) then
                     if (num_samples(i) == 1) then
                        data(i, k, m) = fld_val
                        var_buffer(i, k, m) = 0.0_REAL64
                     else
                        tmp = data(i, k, m)
                        data(i, k, m) = data(i, k, m) + &
                                (fld_val - data(i, k, m)) / &
                                num_samples(i)
                        var_buffer(i, k, m) = var_buffer(i, k, m) + &
                                (fld_val - data(i, k, m)) * &
                                (fld_val - tmp)
                     end if
                  end if
               else
                  if (num_samples(i) == 1) then
                     data(i, k, m) = fld_val
                     var_buffer(i, k, m) = 0.0_REAL64
                  else
                     tmp = data(i, k, m)
                     data(i, k, m) = data(i, k, m) + &
                             (fld_val - data(i, k, m)) / &
                             num_samples(i)
                     var_buffer(i, k, m) = var_buffer(i, k, m) + &
                             (fld_val - data(i, k, m)) * &
                             (fld_val - tmp)
                  end if
               end if
            end do
         end do
      end do

   end subroutine accum_var_3d_dev

   !#######################################################################

   subroutine check_fill_2d_dev(field, fill_value, bad_col)
      ! First column (index into <field>) whose elements do not all share
      ! the fill status of its first element; 0 if none (as buff_2d_check_fill)
      real(REAL64), intent(in)  :: field(:,:)
      real(REAL64), intent(in)  :: fill_value
      integer,      intent(out) :: bad_col
      integer :: i, k, n1, n2, bad
      logical :: col_is_fill, inconsistent

      n1 = size(field, 1)
      n2 = size(field, 2)
      bad = HUGE(1)
      !$acc parallel loop gang vector private(col_is_fill, inconsistent) &
      !$acc   reduction(min:bad) present(field)
      do i = 1, n1
         col_is_fill = (field(i, 1) == fill_value)
         inconsistent = .false.
         !$acc loop seq
         do k = 1, n2
            if ((field(i, k) == fill_value) .neqv. col_is_fill) then
               inconsistent = .true.
            end if
         end do
         if (inconsistent) then
            bad = min(bad, i)
         end if
      end do
      if (bad == HUGE(1)) then
         bad_col = 0
      else
         bad_col = bad
      end if

   end subroutine check_fill_2d_dev

   !#######################################################################

   subroutine check_fill_3d_dev(field, fill_value, bad_col)
      ! First column (index into <field>) whose elements do not all share
      ! the fill status of its first element; 0 if none (as buff_3d_check_fill)
      real(REAL64), intent(in)  :: field(:,:,:)
      real(REAL64), intent(in)  :: fill_value
      integer,      intent(out) :: bad_col
      integer :: i, k, m, n1, n2, n3, bad
      logical :: col_is_fill, inconsistent

      n1 = size(field, 1)
      n2 = size(field, 2)
      n3 = size(field, 3)
      bad = HUGE(1)
      !$acc parallel loop gang vector private(col_is_fill, inconsistent) &
      !$acc   reduction(min:bad) present(field)
      do i = 1, n1
         col_is_fill = (field(i, 1, 1) == fill_value)
         inconsistent = .false.
         !$acc loop seq
         do k = 1, n2
            !$acc loop seq
            do m = 1, n3
               if ((field(i, k, m) == fill_value) .neqv. col_is_fill) then
                  inconsistent = .true.
               end if
            end do
         end do
         if (inconsistent) then
            bad = min(bad, i)
         end if
      end do
      if (bad == HUGE(1)) then
         bad_col = 0
      else
         bad_col = bad
      end if

   end subroutine check_fill_3d_dev

   !#######################################################################
   ! Plain-array device helpers

   subroutine dev_create_i1(a)
      integer, intent(inout) :: a(:)
      !$acc enter data create(a)
   end subroutine dev_create_i1

   subroutine dev_delete_i1(a)
      integer, intent(inout) :: a(:)
      !$acc exit data delete(a)
   end subroutine dev_delete_i1

   subroutine dev_update_self_i1(a)
      integer, intent(inout) :: a(:)
      !$acc update self(a)
   end subroutine dev_update_self_i1

   subroutine dev_update_device_i1(a)
      integer, intent(inout) :: a(:)
      !$acc update device(a)
   end subroutine dev_update_device_i1

   subroutine dev_zero_i1(a)
      integer, intent(inout) :: a(:)
      !$acc kernels present(a)
      a(:) = 0
      !$acc end kernels
   end subroutine dev_zero_i1

   subroutine dev_create_r1(a)
      real(REAL64), intent(inout) :: a(:)
      !$acc enter data create(a)
   end subroutine dev_create_r1

   subroutine dev_delete_r1(a)
      real(REAL64), intent(inout) :: a(:)
      !$acc exit data delete(a)
   end subroutine dev_delete_r1

   subroutine dev_update_self_r1(a)
      real(REAL64), intent(inout) :: a(:)
      !$acc update self(a)
   end subroutine dev_update_self_r1

   subroutine dev_update_device_r1(a)
      real(REAL64), intent(inout) :: a(:)
      !$acc update device(a)
   end subroutine dev_update_device_r1

   subroutine dev_zero_r1(a)
      real(REAL64), intent(inout) :: a(:)
      !$acc kernels present(a)
      a(:) = 0.0_REAL64
      !$acc end kernels
   end subroutine dev_zero_r1

   subroutine dev_create_r2(a)
      real(REAL64), intent(inout) :: a(:,:)
      !$acc enter data create(a)
   end subroutine dev_create_r2

   subroutine dev_delete_r2(a)
      real(REAL64), intent(inout) :: a(:,:)
      !$acc exit data delete(a)
   end subroutine dev_delete_r2

   subroutine dev_update_self_r2(a)
      real(REAL64), intent(inout) :: a(:,:)
      !$acc update self(a)
   end subroutine dev_update_self_r2

   subroutine dev_update_device_r2(a)
      real(REAL64), intent(inout) :: a(:,:)
      !$acc update device(a)
   end subroutine dev_update_device_r2

   subroutine dev_zero_r2(a)
      real(REAL64), intent(inout) :: a(:,:)
      !$acc kernels present(a)
      a(:,:) = 0.0_REAL64
      !$acc end kernels
   end subroutine dev_zero_r2

   subroutine dev_create_r3(a)
      real(REAL64), intent(inout) :: a(:,:,:)
      !$acc enter data create(a)
   end subroutine dev_create_r3

   subroutine dev_delete_r3(a)
      real(REAL64), intent(inout) :: a(:,:,:)
      !$acc exit data delete(a)
   end subroutine dev_delete_r3

   subroutine dev_update_self_r3(a)
      real(REAL64), intent(inout) :: a(:,:,:)
      !$acc update self(a)
   end subroutine dev_update_self_r3

   subroutine dev_update_device_r3(a)
      real(REAL64), intent(inout) :: a(:,:,:)
      !$acc update device(a)
   end subroutine dev_update_device_r3

   subroutine dev_zero_r3(a)
      real(REAL64), intent(inout) :: a(:,:,:)
      !$acc kernels present(a)
      a(:,:,:) = 0.0_REAL64
      !$acc end kernels
   end subroutine dev_zero_r3

#endif

   !#######################################################################

   function buffer_factory(buffer_type, logger) result(newbuf)
      ! Create a new buffer based on <buffer_type>.
      ! <buffer_type> has a format typekind_rank_accum
      ! Where:
      ! <typekind> is a lowercase string representation of a
      !    supported kind from the ISO_FORTRAN_ENV module.
      ! <rank> is the rank of the buffer (no leading zeros)
      ! <accum> is the accumulation type, one of:
      !    lst: Store the last value collected
      !    avg: Accumulate running average
      !    var: Accumulate standard deviation
      !    min: Accumulate smallest value
      !    max: Accumulate largest value

      use hist_msg_handler, only: hist_log_messages
      use hist_msg_handler, only: hist_add_error

      ! Arguments
      class(hist_buffer_t),    pointer                 :: newbuf
      character(len=*),                  intent(in)    :: buffer_type
      type(hist_log_messages), optional, intent(inout) :: logger

      ! Local variables
      character(len=*),   parameter :: subname = 'buffer_factory'
      character(len=512)            :: alloc_errmsg
      ! For buffer allocation
      integer                                    :: aerr
      type(hist_buff_1d_t), pointer :: real_1
      type(hist_buff_2d_t), pointer :: real_2
      type(hist_buff_3d_t), pointer :: real_3

      ! Clear relevant pointers
      nullify(newbuf)
      nullify(real_1)
      nullify(real_2)
      nullify(real_3)

      ! Create new buffer
      select case (trim(buffer_type))
      case ('real_1')
         allocate(real_1, stat=aerr, errmsg=alloc_errmsg)
         if (aerr == 0) then
            newbuf => real_1
         else
            call hist_add_error(subname, &
               'Failed to allocate real_1 buffer, errmsg = ', &
               errstr2=trim(alloc_errmsg), errors=logger)
         end if
      case ('real_2')
         allocate(real_2, stat=aerr, errmsg=alloc_errmsg)
         if (aerr == 0) then
            newbuf => real_2
         else
            call hist_add_error(subname, &
               'Failed to allocate real_2 buffer, errmsg = ', &
               errstr2=trim(alloc_errmsg), errors=logger)
         end if
      case ('real_3')
         allocate(real_3, stat=aerr, errmsg=alloc_errmsg)
         if (aerr == 0) then
            newbuf => real_3
         else
            call hist_add_error(subname, &
               'Failed to allocate real_3 buffer, errmsg = ', &
               errstr2=trim(alloc_errmsg), errors=logger)
         end if
      case default
         call hist_add_error(subname,                                         &
              "Invalid or unsupported buffer type, '",                        &
              errstr2=trim(buffer_type), errstr3="'", errors=logger)
      end select
   end function buffer_factory

   !#######################################################################

end module hist_buffer
