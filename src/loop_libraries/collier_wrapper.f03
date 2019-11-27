module Collier_wrapper
   use COLLIER
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env
   implicit none

contains

   function B0_dummy(p10, m02, m12, scl2) result(res) bind(C, name='B0_impl')

      ! inputs and output from this function
      ! we use the c++ equivalent type names from the iso_c_binding module
      complex(C_DOUBLE_COMPLEX), intent(in) :: p10
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      call Init_cll(2,2,'')
      call SetMuUV2_cll(scl2)

      call B0_cll(res, p10, m02, m12)
   end

   function B1_dummy(p10, m02, m12, scl2) result(res) bind(C, name='B1_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Bcoeff(:,:), Bcoeffuv(:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Bcoeff(0:rank/2, 0:rank))
      allocate(Bcoeffuv(0:rank/2, 0:rank))

      call B_cll(Bcoeff, Bcoeffuv, p10, m02, m12, rank)

      res = Bcoeff(0,1)

      deallocate(Bcoeff, Bcoeffuv)
   end

   function C0_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C0_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      call C0_cll(res, p10, p21, p20, m02, m12, m22)
   end

   function C1_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C1_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Ccoeff(0:rank/2, 0:rank, 0:rank))
      allocate(Ccoeffuv(0:rank/2, 0:rank, 0:rank))

      call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, rank)

      res = Ccoeff(0,1,0)

      deallocate(Ccoeff, Ccoeffuv)
   end

   function C2_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C2_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Ccoeff(0:rank/2, 0:rank, 0:rank))
      allocate(Ccoeffuv(0:rank/2, 0:rank, 0:rank))

      call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, rank)

      res = Ccoeff(0,0,1)

      deallocate(Ccoeff, Ccoeffuv)
   end

   function C00_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C00_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Ccoeff(0:rank/2, 0:rank, 0:rank))
      allocate(Ccoeffuv(0:rank/2, 0:rank, 0:rank))

      call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, rank)

      res = Ccoeff(1,0,0)

      deallocate(Ccoeff, Ccoeffuv)
   end

   function C11_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C11_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Ccoeff(0:rank/2, 0:rank, 0:rank))
      allocate(Ccoeffuv(0:rank/2, 0:rank, 0:rank))

      call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, rank)

      res = Ccoeff(0,2,0)

      deallocate(Ccoeff, Ccoeffuv)
   end

   function C12_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C12_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Ccoeff(0:rank/2, 0:rank, 0:rank))
      allocate(Ccoeffuv(0:rank/2, 0:rank, 0:rank))

      call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, rank)

      res = Ccoeff(0,1,1)

      deallocate(Ccoeff, Ccoeffuv)
   end

   function C22_dummy(p10, p21, p20, m02, m12, m22, scl2) result(res) bind(C, name='C22_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Ccoeff(0:rank/2, 0:rank, 0:rank))
      allocate(Ccoeffuv(0:rank/2, 0:rank, 0:rank))

      call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, rank)

      res = Ccoeff(0,0,2)

      deallocate(Ccoeff, Ccoeffuv)
   end

   function D0_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D0_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      call D0_cll(res,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32)
   end

   function D00_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D00_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(1,0,0,0)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D1_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D1_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,1,0,0)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D11_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D11_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,2,0,0)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D12_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D12_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,1,1,0)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D13_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D13_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,1,0,1)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D2_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D2_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,0,1,0)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D22_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D22_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,0,2,0)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D23_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D23_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,0,1,1)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D3_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D3_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,0,0,1)

      deallocate(Dcoeff, Dcoeffuv)
   end

   function D33_dummy(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,scl2) result(res) bind(C, name='D33_impl')

      complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31
      complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32
      real(C_DOUBLE), intent(in) :: scl2
      complex(C_DOUBLE_COMPLEX) :: res

      complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:)
      integer, parameter :: rank = 3

      call Init_cll(3,3,'')
      call SetMuUV2_cll(scl2)

      allocate(Dcoeff(0:rank/2, 0:rank, 0:rank, 0:rank))
      allocate(Dcoeffuv(0:rank/2, 0:rank, 0:rank, 0:rank))

      call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,rank)

      res = Dcoeff(0,0,0,2)

      deallocate(Dcoeff, Dcoeffuv)
   end

end module
