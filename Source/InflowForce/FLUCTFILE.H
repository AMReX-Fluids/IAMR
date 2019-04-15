
!
! ::: Array dimension definitions which are always 3-d
!

#if __STDC__==1

#define FF_DIMS(a) a##_l1, a##_l2 , a##_l3, a##_h1, a##_h2, a##_h3
#define FF_DIMDEC(a) a##_l1, a##_l2 , a##_l3, a##_h1, a##_h2, a##_h3
#define FF_DIMV(a) a##_l1:a##_h1, a##_l2:a##_h2, a##_l3:a##_h3
#define FF_DIM1(a) a##_l1:a##_h1
#define FF_DIM2(a) a##_l2:a##_h2
#define FF_DIM3(a) a##_l3:a##_h3
#define FF_DIM12(a) a##_l1:a##_h1, a##_l2:a##_h2
#define FF_DIM23(a) a##_l2:a##_h2, a##_l3:a##_h3
#define FF_DIM13(a) a##_l1:a##_h1, a##_l3:a##_h3
#define FF_ARG_L1(a) a##_l1
#define FF_ARG_L2(a) a##_l2
#define FF_ARG_L3(a) a##_l3
#define FF_ARG_H1(a) a##_h1
#define FF_ARG_H2(a) a##_h2
#define FF_ARG_H3(a) a##_h3

#else

#define FF_DIMS(a) a/**/_l1, a/**/_l2 , a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define FF_DIMDEC(a) a/**/_l1, a/**/_l2 , a/**/_l3, a/**/_h1, a/**/_h2, a/**/_h3
#define FF_DIMV(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2, a/**/_l3:a/**/_h3
#define FF_DIM1(a) a/**/_l1:a/**/_h1
#define FF_DIM2(a) a/**/_l2:a/**/_h2
#define FF_DIM3(a) a/**/_l3:a/**/_h3
#define FF_DIM12(a) a/**/_l1:a/**/_h1, a/**/_l2:a/**/_h2
#define FF_DIM23(a) a/**/_l2:a/**/_h2, a/**/_l3:a/**/_h3
#define FF_DIM13(a) a/**/_l1:a/**/_h1, a/**/_l3:a/**/_h3
#define FF_ARG_L1(a) a/**/_l1
#define FF_ARG_L2(a) a/**/_l2
#define FF_ARG_L3(a) a/**/_l3
#define FF_ARG_H1(a) a/**/_h1
#define FF_ARG_H2(a) a/**/_h2
#define FF_ARG_H3(a) a/**/_h3

#endif
