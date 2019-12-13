#define PASTE(a) a
#define CAT(a,b) PASTE(a)b
#define TO_STRING(ARG) _QUOTE_START_ ARG _QUOTE_END_

#define two_point(NAME,N1,N2) \
_NL_function CAT(NAME,_dummy)(p10, m02, m12) result(res) bind(C, name=TO_STRING(CAT(NAME,_impl)))\
_NL_   complex(C_DOUBLE_COMPLEX), intent(in) :: p10 \
_NL_   complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12 \
_NL_   complex(C_DOUBLE_COMPLEX) :: res \
_NL_   complex(REAL64), allocatable :: Bcoeff(:,:), Bcoeffuv(:,:) \
_NL_\
_NL_   allocate(Bcoeff(0:0, 0:1)) \
_NL_   allocate(Bcoeffuv(0:0, 0:1)) \
_NL_   call B_cll(Bcoeff, Bcoeffuv, p10, m02, m12, 1) \
_NL_\
_NL_   res = Bcoeff(N1,N2) \
_NL_\
_NL_   deallocate(Bcoeff, Bcoeffuv) \
_NL_end

#define three_point(NAME,N1,N2,N3) \
_NL_function CAT(NAME,_dummy)(p10, p21, p20, m02, m12, m22) result(res) bind(C, name=TO_STRING(CAT(NAME,_impl)))\
_NL_   complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p20 \
_NL_   complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22 \
_NL_   complex(C_DOUBLE_COMPLEX) :: res \
_NL_   complex(REAL64), allocatable :: Ccoeff(:,:,:), Ccoeffuv(:,:,:) \
_NL_\
_NL_   allocate(Ccoeff(0:1, 0:2, 0:2)) \
_NL_   allocate(Ccoeffuv(0:1, 0:2, 0:2)) \
_NL_   call C_cll(Ccoeff, Ccoeffuv, p10, p21, p20, m02, m12, m22, 2) \
_NL_\
_NL_   res = Ccoeff(N1,N2,N3) \
_NL_\
_NL_   deallocate(Ccoeff, Ccoeffuv) \
_NL_end

#define four_point(NAME,N1,N2,N3,N4) \
_NL_function CAT(NAME,_dummy)(p10,p21,p32,p30,p20,p31,m02,m12,m22,m32) result(res) bind(C, name=TO_STRING(CAT(NAME,_impl)))\
_NL_   complex(C_DOUBLE_COMPLEX), intent(in) :: p10, p21, p32, p30, p20, p31 \
_NL_   complex(C_DOUBLE_COMPLEX), intent(in) :: m02, m12, m22, m32 \
_NL_   complex(C_DOUBLE_COMPLEX) :: res \
_NL_   complex(REAL64), allocatable :: Dcoeff(:,:,:,:), Dcoeffuv(:,:,:,:) \
_NL_\
_NL_   allocate(Dcoeff(0:1, 0:2, 0:2, 0:2)) \
_NL_   allocate(Dcoeffuv(0:1, 0:2, 0:2, 0:2)) \
_NL_   call D_cll(Dcoeff,Dcoeffuv,p10,p21,p32,p30,p20,p31,m02,m12,m22,m32,2) \
_NL_\
_NL_   res = Dcoeff(N1,N2,N3,N4) \
_NL_\
_NL_   deallocate(Dcoeff, Dcoeffuv) \
_NL_end

module Collier_wrapper
   use COLLIER
   use, intrinsic :: iso_c_binding
   use, intrinsic :: iso_fortran_env
   implicit none

contains
   subroutine initialize_collier_dummy() bind(C, name='initialize_collier_impl')
      real(C_DOUBLE) :: mu2
      call Init_cll(4,2,'')
   end

   subroutine set_mu2_uv_dummy(scl2) bind(C, name='set_mu2_uv_impl')
      real(C_DOUBLE), intent(in) :: scl2
      call SetMuUV2_cll(scl2)
   end

   two_point(B0,0,0)
   two_point(B1,0,1)

   three_point(C0,0,0,0)
   three_point(C00,1,0,0)
   three_point(C1,0,1,0)
   three_point(C11,0,2,0)
   three_point(C12,0,1,1)
   three_point(C2,0,0,1)
   three_point(C22,0,0,2)

   four_point(D0,0,0,0,0)
   four_point(D00,1,0,0,0)
   four_point(D1,0,1,0,0)
   four_point(D11,0,2,0,0)
   four_point(D12,0,1,1,0)
   four_point(D13,0,1,0,1)
   four_point(D2,0,0,1,0)
   four_point(D22,0,0,2,0)
   four_point(D23,0,0,1,1)
   four_point(D3,0,0,0,1)
   four_point(D33,0,0,0,2)

end module
