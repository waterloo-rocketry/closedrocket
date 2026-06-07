#line 1 "GNC_codegen_SIL.c"
extern void profileEnd_GNC_codegen_SIL(unsigned); extern void profileStart_GNC_codegen_SIL(unsigned); 
#line 1 "C:\\Users\\trist\\Documents\\GitHub\\simulink-canards\\embedded-coder\\codegen\\lib\\GNC_codegen_SIL\\instrumented\\GNC_codegen_SIL.h"
#ifndef GNC_CODEGEN_SIL_H
#define GNC_CODEGEN_SIL_H
#line 1 "C:\\Users\\trist\\Documents\\GitHub\\simulink-canards\\embedded-coder\\codegen\\lib\\GNC_codegen_SIL\\instrumented\\GNC_codegen_SIL_types.h"
#ifndef GNC_CODEGEN_SIL_TYPES_H
#define GNC_CODEGEN_SIL_TYPES_H
#line 1 "C:\\Users\\trist\\Documents\\GitHub\\simulink-canards\\embedded-coder\\codegen\\lib\\GNC_codegen_SIL\\instrumented\\rtwtypes.h"
#ifndef RTWTYPES_H
#define RTWTYPES_H
#line 12 "C:\\Program Files\\MATLAB\\R2025b\\extern\\include\\tmwtypes.h"
#ifndef tmwtypes_h
#define tmwtypes_h

#ifndef __TMWTYPES__
#define __TMWTYPES__
#line 23 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\vcruntime.h"
#ifndef _VCRUNTIME_H
#define _VCRUNTIME_H
#line 36
#ifndef _UCRT
#define _UCRT
#endif /* _UCRT */
#line 148 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\sal.h"
#ifndef _SAL_VERSION
#define _SAL_VERSION 20
#endif /* _SAL_VERSION */

#ifndef __SAL_H_VERSION
#define __SAL_H_VERSION 180000000
#endif /* __SAL_H_VERSION */
#line 182
#ifndef _USE_DECLSPECS_FOR_SAL
#define _USE_DECLSPECS_FOR_SAL 0
#endif /* _USE_DECLSPECS_FOR_SAL */
#ifndef _USE_ATTRIBUTES_FOR_SAL
#define _USE_ATTRIBUTES_FOR_SAL 0
#endif /* _USE_ATTRIBUTES_FOR_SAL */
#line 224
#ifndef _SAL_L_Source_
#define _SAL_L_Source_(Name,args,annotes) _SA_annotes3(SAL_name, #Name, "", "2") _Group_(annotes _SAL_nop_impl_)

#endif /* _SAL_L_Source_ */
#line 708
#pragma region Input Buffer SAL 1 compatibility macros
#line 1472
#pragma endregion Input Buffer SAL 1 compatibility macros
#line 2363
#ifndef __nothrow
#define __nothrow
#endif /* __nothrow */
#line 16 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\concurrencysal.h"
#ifndef CONCURRENCYSAL_H
#define CONCURRENCYSAL_H
#line 293
#ifndef _Interlocked_operand_
#define _Interlocked_operand_
#endif /* _Interlocked_operand_ */
#line 394
#endif /* CONCURRENCYSAL_H */
#line 15 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\vadefs.h"
#pragma pack ( push, 8 )
#line 37
#ifndef _VCRUNTIME_EXTRA_DISABLED_WARNINGS
#define _VCRUNTIME_EXTRA_DISABLED_WARNINGS
#endif /* _VCRUNTIME_EXTRA_DISABLED_WARNINGS */



#ifndef _VCRUNTIME_DISABLED_WARNINGS
#define _VCRUNTIME_DISABLED_WARNINGS _VCRUNTIME_DISABLED_WARNING_4339 _VCRUNTIME_DISABLED_WARNING_4412 4514 4820 _VCRUNTIME_EXTRA_DISABLED_WARNINGS
#endif /* _VCRUNTIME_DISABLED_WARNINGS */

#pragma warning(push)
#pragma warning(disable: 4514 4820)
#line 58
#ifndef _UINTPTR_T_DEFINED
#define _UINTPTR_T_DEFINED

typedef unsigned __int64 uintptr_t; 



#endif /* _UINTPTR_T_DEFINED */

#ifndef _VA_LIST_DEFINED
#define _VA_LIST_DEFINED



typedef char *va_list; 

#endif /* _VA_LIST_DEFINED */
#line 155
void __cdecl __va_start(va_list *, ...); 
#line 207
#pragma warning(pop)
#pragma pack ( pop )
#line 60 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\vcruntime.h"
#pragma warning(push)
#pragma warning(disable: 4514 4820)
#line 96
__pragma( pack ( push, 8 )) 
#line 188
typedef unsigned __int64 size_t; 
typedef __int64 ptrdiff_t; 
typedef __int64 intptr_t; 
#line 204
typedef _Bool __vcrt_bool; 



#ifndef _SIZE_T_DEFINED
#define _SIZE_T_DEFINED
#endif /* _SIZE_T_DEFINED */

#ifndef _PTRDIFF_T_DEFINED
#define _PTRDIFF_T_DEFINED
#endif /* _PTRDIFF_T_DEFINED */

#ifndef _INTPTR_T_DEFINED
#define _INTPTR_T_DEFINED
#endif /* _INTPTR_T_DEFINED */


#ifndef _WCHAR_T_DEFINED
#define _WCHAR_T_DEFINED
typedef unsigned short wchar_t; 
#endif /* _WCHAR_T_DEFINED */
#line 378
void __cdecl __security_init_cookie(void); 
#line 387
void __cdecl __security_check_cookie(uintptr_t _StackCookie); 
__declspec(noreturn) void __cdecl __report_gsfailure(uintptr_t _StackCookie); 



extern uintptr_t __security_cookie; 
#line 400
__pragma( pack ( pop )) 

#pragma warning(pop)

#endif /* _VCRUNTIME_H */
#line 13 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\limits.h"
#pragma warning(push)
#pragma warning(disable: 4514 4820)

__pragma( pack ( push, 8 )) 
#line 70
#ifndef RSIZE_MAX
#define RSIZE_MAX (SIZE_MAX >> 1)
#endif /* RSIZE_MAX */

__pragma( pack ( pop )) 

#pragma warning(pop)
#line 8 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\stdbool.h"
#ifndef _STDBOOL
#define _STDBOOL
#line 21
#endif /* _STDBOOL */
#line 62 "C:\\Program Files\\MATLAB\\R2025b\\extern\\include\\tmwtypes.h"
#ifndef FLT_MANT_DIG
#define FLT_MANT_DIG 24
#endif /* FLT_MANT_DIG */
#ifndef DBL_MANT_DIG
#define DBL_MANT_DIG 53
#endif /* DBL_MANT_DIG */
#line 89
typedef unsigned char uchar_T; 
typedef unsigned short ushort_T; 
typedef unsigned long ulong_T; 
#line 97
typedef unsigned __int64 ulonglong_T; 
#line 222
typedef char int8_T; 
#line 239
typedef unsigned char uint8_T; 
#line 257
typedef short int16_T; 
#line 275
typedef unsigned short uint16_T; 
#line 293
typedef int int32_T; 
#line 311
typedef unsigned uint32_T; 
#line 372
typedef float real32_T; 
#line 386
typedef double real64_T; 
#line 436
#ifndef INT64_T
#define INT64_T __int64
#endif /* INT64_T */
#ifndef UINT64_T
#define UINT64_T unsigned __int64
#endif /* UINT64_T */
#ifndef FMT64
#define FMT64 "I64"
#endif /* FMT64 */
#line 465
typedef __int64 int64_T; 
#line 479
typedef unsigned __int64 uint64_T; 
#line 535
typedef real64_T real_T; 
#line 544
typedef real_T time_T; 
#line 556
typedef unsigned char boolean_T; 


#ifndef CHARACTER_T
#define CHARACTER_T char
#endif /* CHARACTER_T */
typedef char char_T; 


#ifndef INTEGER_T
#define INTEGER_T int
#endif /* INTEGER_T */
typedef int int_T; 


#ifndef UINTEGER_T
#define UINTEGER_T unsigned
#endif /* UINTEGER_T */
typedef unsigned uint_T; 


#ifndef BYTE_T
#define BYTE_T unsigned char
#endif /* BYTE_T */
typedef unsigned char byte_T; 
#line 592
typedef 
#line 590
struct { 
real32_T re, im; 
} creal32_T; 
#line 601
typedef 
#line 599
struct { 
real64_T re, im; 
} creal64_T; 
#line 610
typedef 
#line 608
struct { 
real_T re, im; 
} creal_T; 
#line 621
typedef 
#line 619
struct { 
int8_T re, im; 
} cint8_T; 
#line 630
typedef 
#line 628
struct { 
uint8_T re, im; 
} cuint8_T; 
#line 639
typedef 
#line 637
struct { 
int16_T re, im; 
} cint16_T; 
#line 648
typedef 
#line 646
struct { 
uint16_T re, im; 
} cuint16_T; 
#line 657
typedef 
#line 655
struct { 
int32_T re, im; 
} cint32_T; 
#line 666
typedef 
#line 664
struct { 
uint32_T re, im; 
} cuint32_T; 
#line 675
typedef 
#line 673
struct { 
int64_T re, im; 
} cint64_T; 
#line 684
typedef 
#line 682
struct { 
uint64_T re, im; 
} cuint64_T; 
#line 759
__forceinline uint64_T double_to_uint64_helper(double d) { 
union double_to_uint64_union_type { 
double dd; 
uint64_T i64; 
} di; 
di.dd = d; 
return ((((di.i64) & (0xfffffffffffffi64)) | (0x10000000000000i64)) << 11); 
} 
#line 9 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\stddef.h"
#ifndef _INC_STDDEF
#define _INC_STDDEF
#line 40 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt.h"
#ifndef _ARM_WINAPI_PARTITION_DESKTOP_SDK_AVAILABLE
#define _ARM_WINAPI_PARTITION_DESKTOP_SDK_AVAILABLE 1
#endif /* _ARM_WINAPI_PARTITION_DESKTOP_SDK_AVAILABLE */
#line 76
#ifndef _UCRT_EXTRA_DISABLED_WARNINGS
#define _UCRT_EXTRA_DISABLED_WARNINGS
#endif /* _UCRT_EXTRA_DISABLED_WARNINGS */
#line 92
#ifndef _UCRT_DISABLED_WARNINGS
#define _UCRT_DISABLED_WARNINGS 4324 _UCRT_DISABLED_WARNING_4412 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727 _UCRT_EXTRA_DISABLED_WARNINGS
#endif /* _UCRT_DISABLED_WARNINGS */
#line 121
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 
#line 144
#ifndef _ACRTIMP_ALT
#define _ACRTIMP_ALT _ACRTIMP
#endif /* _ACRTIMP_ALT */
#line 274
typedef _Bool __crt_bool; 
#line 328
#ifndef _CRT_UNUSED
#define _CRT_UNUSED(x) (void)x
#endif /* _CRT_UNUSED */
#line 371
void __cdecl _invalid_parameter_noinfo(void); 
__declspec(noreturn) void __cdecl _invalid_parameter_noinfo_noreturn(void); 

__declspec(noreturn) void __cdecl 
_invoke_watson(const wchar_t * _Expression, const wchar_t * _FunctionName, const wchar_t * _FileName, unsigned _LineNo, uintptr_t _Reserved); 
#line 498
#ifndef __STDC_WANT_SECURE_LIB__
#define __STDC_WANT_SECURE_LIB__ 1
#endif /* __STDC_WANT_SECURE_LIB__ */
#line 593
#ifndef _CRT_SECURE_CPP_NOTHROW
#define _CRT_SECURE_CPP_NOTHROW throw()
#endif /* _CRT_SECURE_CPP_NOTHROW */
#line 604
typedef int errno_t; 
typedef unsigned short wint_t; 
typedef unsigned short wctype_t; 
typedef long __time32_t; 
typedef __int64 __time64_t; 
#line 615
typedef 
#line 610
struct __crt_locale_data_public { 

const unsigned short *_locale_pctype; 
int _locale_mb_cur_max; 
unsigned _locale_lc_codepage; 
} __crt_locale_data_public; 
#line 621
typedef 
#line 617
struct __crt_locale_pointers { 

struct __crt_locale_data *locinfo; 
struct __crt_multibyte_data *mbcinfo; 
} __crt_locale_pointers; 

typedef __crt_locale_pointers *_locale_t; 
#line 629
typedef 
#line 625
struct _Mbstatet { 

unsigned long _Wchar; 
unsigned short _Byte, _State; 
} _Mbstatet; 

typedef _Mbstatet mbstate_t; 
#line 645
typedef __time64_t time_t; 




#ifndef _TIME_T_DEFINED
#define _TIME_T_DEFINED
#endif /* _TIME_T_DEFINED */


typedef size_t rsize_t; 
#line 2072
__pragma( pack ( pop )) 


#pragma warning(pop)
#line 14 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\stddef.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 
#line 35
int *__cdecl _errno(void); 


errno_t __cdecl _set_errno(int _Value); 
errno_t __cdecl _get_errno(int * _Value); 
#line 55
extern unsigned long __cdecl __threadid(void); 

extern uintptr_t __cdecl __threadhandle(void); 



__pragma( pack ( pop )) 

#pragma warning(pop)
#endif /* _INC_STDDEF */
#line 834 "C:\\Program Files\\MATLAB\\R2025b\\extern\\include\\tmwtypes.h"
typedef size_t mwSize; 
typedef size_t mwIndex; 
typedef ptrdiff_t mwSignedIndex; 





#ifndef SLSIZE_SLINDEX
#define SLSIZE_SLINDEX

typedef int64_T SLIndex; 
typedef int64_T SLSize; 




#endif /* SLSIZE_SLINDEX */
#line 880
typedef wchar_t CHAR16_T; 





#endif /* __TMWTYPES__ */

#endif /* tmwtypes_h */
#line 32 "C:\\Users\\trist\\Documents\\GitHub\\simulink-canards\\embedded-coder\\codegen\\lib\\GNC_codegen_SIL\\instrumented\\rtwtypes.h"
#endif /* RTWTYPES_H */
#line 6 "C:\\Users\\trist\\Documents\\GitHub\\simulink-canards\\embedded-coder\\codegen\\lib\\GNC_codegen_SIL\\instrumented\\GNC_codegen_SIL_types.h"
#ifndef typedef_struct1_T
#define typedef_struct1_T
#line 16
typedef 
#line 8
struct { 
double board_gyro[3]; 
double mti_gyro[3]; 
double ad_gyro[3]; 
double board_mag_earth[3]; 
double mti_mag_earth[3]; 
double board_baro; 
double mti_baro; 
} struct1_T; 
#endif /* typedef_struct1_T */

#ifndef typedef_struct2_T
#define typedef_struct2_T
#line 32
typedef 
#line 21
struct { 
double board_accel_f[3]; 
double board_gyro_f[3]; 
double mti_accel_f[3]; 
double mti_gyro_f[3]; 
double ad_accel_f[3]; 
double ad_gyro_f[3]; 
double board_baro_f; 
double board_mag_f[3]; 
double mti_baro_f; 
double mti_mag_f[3]; 
} struct2_T; 
#endif /* typedef_struct2_T */

#ifndef typedef_struct4_T
#define typedef_struct4_T



typedef 
#line 37
struct { 
double meas[3]; 
boolean_T status; 
} struct4_T; 
#endif /* typedef_struct4_T */

#ifndef typedef_struct5_T
#define typedef_struct5_T



typedef 
#line 45
struct { 
double meas; 
boolean_T status; 
} struct5_T; 
#endif /* typedef_struct5_T */

#ifndef typedef_struct3_T
#define typedef_struct3_T
#line 64
typedef 
#line 53
struct { 
struct4_T board_accel; 
struct4_T board_gyro; 
struct4_T mti_accel; 
struct4_T mti_gyro; 
struct4_T ad_accel; 
struct4_T ad_gyro; 
struct5_T board_baro; 
struct4_T board_mag; 
struct5_T mti_baro; 
struct4_T mti_mag; 
} struct3_T; 
#endif /* typedef_struct3_T */

#ifndef typedef_struct6_T
#define typedef_struct6_T
#line 76
typedef 
#line 69
struct { 
double pressure; 
double temperature; 
double density; 
double sonic_speed; 
double mach; 
double dynamic_pressure; 
} struct6_T; 
#endif /* typedef_struct6_T */

#ifndef typedef_struct0_T
#define typedef_struct0_T
#line 87
typedef 
#line 81
struct { 
double coeffs[2]; 
double w_old; 
double P_minus[4]; 
double d_old; 
double w_dot_old; 
} struct0_T; 
#endif /* typedef_struct0_T */

#ifndef typedef_struct_T
#define typedef_struct_T
#line 100
typedef 
#line 92
struct { 
double Cn_alpha; 
double J[9]; 
double Jinv[9]; 
double c_aero; 
double c_canard; 
double elevation; 
double g[3]; 
} struct_T; 
#endif /* typedef_struct_T */

#ifndef c_typedef_GNC_codegen_SILPersis
#define c_typedef_GNC_codegen_SILPersis





typedef 
#line 105
struct { 
struct_T param; 
struct_T b_param; 
struct_T c_param; 
struct_T d_param; 
} GNC_codegen_SILPersistentData; 
#endif /* c_typedef_GNC_codegen_SILPersis */

#ifndef c_typedef_GNC_codegen_SILStackD
#define c_typedef_GNC_codegen_SILStackD


typedef 
#line 115
struct { 
GNC_codegen_SILPersistentData *pd; 
} GNC_codegen_SILStackData; 
#endif /* c_typedef_GNC_codegen_SILStackD */

#endif /* GNC_CODEGEN_SIL_TYPES_H */
#line 9 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\stdlib.h"
#ifndef _INC_STDLIB
#define _INC_STDLIB
#line 13 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_malloc.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 
#line 58
__declspec(allocator) __declspec(restrict) void *__cdecl 
_calloc_base(size_t _Count, size_t _Size); 
#line 65
__declspec(allocator) __declspec(restrict) void *__cdecl 
calloc(size_t _Count, size_t _Size); 
#line 72
int __cdecl _callnewh(size_t _Size); 




__declspec(allocator) void *__cdecl 
_expand(void * _Block, size_t _Size); 
#line 84
void __cdecl _free_base(void * _Block); 




void __cdecl free(void * _Block); 




__declspec(allocator) __declspec(restrict) void *__cdecl 
_malloc_base(size_t _Size); 




__declspec(allocator) __declspec(restrict) void *__cdecl 
malloc(size_t _Size); 
#line 107
size_t __cdecl _msize_base(void * _Block); 
#line 113
size_t __cdecl _msize(void * _Block); 




__declspec(allocator) __declspec(restrict) void *__cdecl 
_realloc_base(void * _Block, size_t _Size); 
#line 125
__declspec(allocator) __declspec(restrict) void *__cdecl 
realloc(void * _Block, size_t _Size); 
#line 132
__declspec(allocator) __declspec(restrict) void *__cdecl 
_recalloc_base(void * _Block, size_t _Count, size_t _Size); 
#line 140
__declspec(allocator) __declspec(restrict) void *__cdecl 
_recalloc(void * _Block, size_t _Count, size_t _Size); 
#line 148
void __cdecl _aligned_free(void * _Block); 




__declspec(allocator) __declspec(restrict) void *__cdecl 
_aligned_malloc(size_t _Size, size_t _Alignment); 
#line 160
__declspec(allocator) __declspec(restrict) void *__cdecl 
_aligned_offset_malloc(size_t _Size, size_t _Alignment, size_t _Offset); 
#line 169
size_t __cdecl _aligned_msize(void * _Block, size_t _Alignment, size_t _Offset); 
#line 176
__declspec(allocator) __declspec(restrict) void *__cdecl 
_aligned_offset_realloc(void * _Block, size_t _Size, size_t _Alignment, size_t _Offset); 
#line 185
__declspec(allocator) __declspec(restrict) void *__cdecl 
_aligned_offset_recalloc(void * _Block, size_t _Count, size_t _Size, size_t _Alignment, size_t _Offset); 
#line 195
__declspec(allocator) __declspec(restrict) void *__cdecl 
_aligned_realloc(void * _Block, size_t _Size, size_t _Alignment); 
#line 203
__declspec(allocator) __declspec(restrict) void *__cdecl 
_aligned_recalloc(void * _Block, size_t _Count, size_t _Size, size_t _Alignment); 
#line 232
__pragma( pack ( pop )) 

#pragma warning(pop)
#line 16 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_search.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 


typedef int (__cdecl *_CoreCrtSecureSearchSortCompareFunction)(void *, const void *, const void *); 
typedef int (__cdecl *_CoreCrtNonSecureSearchSortCompareFunction)(const void *, const void *); 
#line 30
void *__cdecl bsearch_s(const void * _Key, const void * _Base, rsize_t _NumOfElements, rsize_t _SizeOfElements, _CoreCrtSecureSearchSortCompareFunction _CompareFunction, void * _Context); 
#line 39
void __cdecl qsort_s(void * _Base, rsize_t _NumOfElements, rsize_t _SizeOfElements, _CoreCrtSecureSearchSortCompareFunction _CompareFunction, void * _Context); 
#line 52
void *__cdecl bsearch(const void * _Key, const void * _Base, size_t _NumOfElements, size_t _SizeOfElements, _CoreCrtNonSecureSearchSortCompareFunction _CompareFunction); 
#line 60
void __cdecl qsort(void * _Base, size_t _NumOfElements, size_t _SizeOfElements, _CoreCrtNonSecureSearchSortCompareFunction _CompareFunction); 
#line 68
void *__cdecl _lfind_s(const void * _Key, const void * _Base, unsigned * _NumOfElements, size_t _SizeOfElements, _CoreCrtSecureSearchSortCompareFunction _CompareFunction, void * _Context); 
#line 78
void *__cdecl _lfind(const void * _Key, const void * _Base, unsigned * _NumOfElements, unsigned _SizeOfElements, _CoreCrtNonSecureSearchSortCompareFunction _CompareFunction); 
#line 87
void *__cdecl _lsearch_s(const void * _Key, void * _Base, unsigned * _NumOfElements, size_t _SizeOfElements, _CoreCrtSecureSearchSortCompareFunction _CompareFunction, void * _Context); 
#line 97
void *__cdecl _lsearch(const void * _Key, void * _Base, unsigned * _NumOfElements, unsigned _SizeOfElements, _CoreCrtNonSecureSearchSortCompareFunction _CompareFunction); 
#line 194
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _lfind. See online help for details." "")) void *__cdecl 
lfind(const void * _Key, const void * _Base, unsigned * _NumOfElements, unsigned _SizeOfElements, _CoreCrtNonSecureSearchSortCompareFunction _CompareFunction); 
#line 203
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _lsearch. See online help for detail" "s.")) void *__cdecl 
lsearch(const void * _Key, void * _Base, unsigned * _NumOfElements, unsigned _SizeOfElements, _CoreCrtNonSecureSearchSortCompareFunction _CompareFunction); 
#line 216
__pragma( pack ( pop )) 

#pragma warning(pop)
#line 13 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_wstdlib.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 
#line 54
errno_t __cdecl _itow_s(int _Value, wchar_t * _Buffer, size_t _BufferCount, int _Radix); 
#line 68
wchar_t *__cdecl _itow(int _Value, wchar_t * _Buffer, int _Radix); 
#line 77
errno_t __cdecl _ltow_s(long _Value, wchar_t * _Buffer, size_t _BufferCount, int _Radix); 
#line 91
wchar_t *__cdecl _ltow(long _Value, wchar_t * _Buffer, int _Radix); 
#line 99
errno_t __cdecl _ultow_s(unsigned long _Value, wchar_t * _Buffer, size_t _BufferCount, int _Radix); 
#line 113
wchar_t *__cdecl _ultow(unsigned long _Value, wchar_t * _Buffer, int _Radix); 
#line 121
double __cdecl wcstod(const wchar_t * _String, wchar_t ** _EndPtr); 
#line 127
double __cdecl _wcstod_l(const wchar_t * _String, wchar_t ** _EndPtr, _locale_t _Locale); 
#line 134
long __cdecl wcstol(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix); 
#line 141
long __cdecl _wcstol_l(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 149
__int64 __cdecl wcstoll(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix); 
#line 156
__int64 __cdecl _wcstoll_l(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 164
unsigned long __cdecl wcstoul(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix); 
#line 171
unsigned long __cdecl _wcstoul_l(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 179
unsigned __int64 __cdecl wcstoull(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix); 
#line 186
unsigned __int64 __cdecl _wcstoull_l(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 194
long double __cdecl wcstold(const wchar_t * _String, wchar_t ** _EndPtr); 
#line 200
long double __cdecl _wcstold_l(const wchar_t * _String, wchar_t ** _EndPtr, _locale_t _Locale); 
#line 207
float __cdecl wcstof(const wchar_t * _String, wchar_t ** _EndPtr); 
#line 213
float __cdecl _wcstof_l(const wchar_t * _String, wchar_t ** _EndPtr, _locale_t _Locale); 
#line 220
double __cdecl _wtof(const wchar_t * _String); 




double __cdecl _wtof_l(const wchar_t * _String, _locale_t _Locale); 
#line 231
int __cdecl _wtoi(const wchar_t * _String); 




int __cdecl _wtoi_l(const wchar_t * _String, _locale_t _Locale); 
#line 242
long __cdecl _wtol(const wchar_t * _String); 




long __cdecl _wtol_l(const wchar_t * _String, _locale_t _Locale); 
#line 253
__int64 __cdecl _wtoll(const wchar_t * _String); 




__int64 __cdecl _wtoll_l(const wchar_t * _String, _locale_t _Locale); 
#line 264
errno_t __cdecl _i64tow_s(__int64 _Value, wchar_t * _Buffer, size_t _BufferCount, int _Radix); 
#line 272
wchar_t *__cdecl _i64tow(__int64 _Value, wchar_t * _Buffer, int _Radix); 
#line 279
errno_t __cdecl _ui64tow_s(unsigned __int64 _Value, wchar_t * _Buffer, size_t _BufferCount, int _Radix); 
#line 287
wchar_t *__cdecl _ui64tow(unsigned __int64 _Value, wchar_t * _Buffer, int _Radix); 
#line 294
__int64 __cdecl _wtoi64(const wchar_t * _String); 




__int64 __cdecl _wtoi64_l(const wchar_t * _String, _locale_t _Locale); 
#line 305
__int64 __cdecl _wcstoi64(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix); 
#line 312
__int64 __cdecl _wcstoi64_l(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 320
unsigned __int64 __cdecl _wcstoui64(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix); 
#line 327
unsigned __int64 __cdecl _wcstoui64_l(const wchar_t * _String, wchar_t ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 339
__declspec(allocator) wchar_t *__cdecl _wfullpath(wchar_t * _Buffer, const wchar_t * _Path, size_t _BufferCount); 
#line 348
errno_t __cdecl _wmakepath_s(wchar_t * _Buffer, size_t _BufferCount, const wchar_t * _Drive, const wchar_t * _Dir, const wchar_t * _Filename, const wchar_t * _Ext); 
#line 366
void __cdecl _wmakepath(wchar_t * _Buffer, const wchar_t * _Drive, const wchar_t * _Dir, const wchar_t * _Filename, const wchar_t * _Ext); 
#line 375
void __cdecl _wperror(const wchar_t * _ErrorMessage); 




void __cdecl _wsplitpath(const wchar_t * _FullPath, wchar_t * _Drive, wchar_t * _Dir, wchar_t * _Filename, wchar_t * _Ext); 
#line 388
errno_t __cdecl _wsplitpath_s(const wchar_t * _FullPath, wchar_t * _Drive, size_t _DriveCount, wchar_t * _Dir, size_t _DirCount, wchar_t * _Filename, size_t _FilenameCount, wchar_t * _Ext, size_t _ExtCount); 
#line 409
errno_t __cdecl _wdupenv_s(wchar_t ** _Buffer, size_t * _BufferCount, const wchar_t * _VarName); 
#line 418
wchar_t *__cdecl _wgetenv(const wchar_t * _VarName); 
#line 424
errno_t __cdecl _wgetenv_s(size_t * _RequiredCount, wchar_t * _Buffer, size_t _BufferCount, const wchar_t * _VarName); 
#line 440
int __cdecl _wputenv(const wchar_t * _EnvString); 




errno_t __cdecl _wputenv_s(const wchar_t * _Name, const wchar_t * _Value); 




errno_t __cdecl _wsearchenv_s(const wchar_t * _Filename, const wchar_t * _VarName, wchar_t * _Buffer, size_t _BufferCount); 
#line 464
void __cdecl _wsearchenv(const wchar_t * _Filename, const wchar_t * _VarName, wchar_t * _ResultPath); 
#line 471
int __cdecl _wsystem(const wchar_t * _Command); 
#line 479
__pragma( pack ( pop )) 

#pragma warning(pop)
#line 18 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\stdlib.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 



#ifndef _countof
#define _countof __crt_countof
#endif /* _countof */
#line 38
void __cdecl _swab(char * _Buf1, char * _Buf2, int _SizeInBytes); 
#line 56
__declspec(noreturn) void __cdecl exit(int _Code); 
__declspec(noreturn) void __cdecl _exit(int _Code); 
__declspec(noreturn) void __cdecl _Exit(int _Code); 
__declspec(noreturn) void __cdecl quick_exit(int _Code); 
__declspec(noreturn) void __cdecl abort(void); 
#line 67
unsigned __cdecl _set_abort_behavior(unsigned _Flags, unsigned _Mask); 
#line 74
#ifndef _CRT_ONEXIT_T_DEFINED
#define _CRT_ONEXIT_T_DEFINED

typedef int (__cdecl *_onexit_t)(void); 



#endif /* _CRT_ONEXIT_T_DEFINED */
#line 144
int __cdecl atexit(void (__cdecl *)(void)); 
_onexit_t __cdecl _onexit(_onexit_t _Func); 


int __cdecl at_quick_exit(void (__cdecl *)(void)); 
#line 159
typedef void (__cdecl *_purecall_handler)(void); 


typedef void (__cdecl *_invalid_parameter_handler)(const wchar_t *, const wchar_t *, const wchar_t *, unsigned, uintptr_t); 
#line 171
_purecall_handler __cdecl _set_purecall_handler(_purecall_handler _Handler); 



_purecall_handler __cdecl _get_purecall_handler(void); 


_invalid_parameter_handler __cdecl _set_invalid_parameter_handler(_invalid_parameter_handler _Handler); 



_invalid_parameter_handler __cdecl _get_invalid_parameter_handler(void); 

_invalid_parameter_handler __cdecl _set_thread_local_invalid_parameter_handler(_invalid_parameter_handler _Handler); 



_invalid_parameter_handler __cdecl _get_thread_local_invalid_parameter_handler(void); 
#line 212
int __cdecl _set_error_mode(int _Mode); 




int *__cdecl _errno(void); 


errno_t __cdecl _set_errno(int _Value); 
errno_t __cdecl _get_errno(int * _Value); 

unsigned long *__cdecl __doserrno(void); 


errno_t __cdecl _set_doserrno(unsigned long _Value); 
errno_t __cdecl _get_doserrno(unsigned long * _Value); 


char **__cdecl __sys_errlist(void); 


int *__cdecl __sys_nerr(void); 


void __cdecl perror(const char * _ErrMsg); 
#line 242
char **__cdecl __p__pgmptr(void); 
wchar_t **__cdecl __p__wpgmptr(void); 
int *__cdecl __p__fmode(void); 
#line 259
errno_t __cdecl _get_pgmptr(char ** _Value); 


errno_t __cdecl _get_wpgmptr(wchar_t ** _Value); 

errno_t __cdecl _set_fmode(int _Mode); 

errno_t __cdecl _get_fmode(int * _PMode); 
#line 279
typedef 
#line 275
struct _div_t { 

int quot; 
int rem; 
} div_t; 
#line 285
typedef 
#line 281
struct _ldiv_t { 

long quot; 
long rem; 
} ldiv_t; 
#line 291
typedef 
#line 287
struct _lldiv_t { 

__int64 quot; 
__int64 rem; 
} lldiv_t; 

int __cdecl abs(int _Number); 
long __cdecl labs(long _Number); 
__int64 __cdecl llabs(__int64 _Number); 
__int64 __cdecl _abs64(__int64 _Number); 

unsigned short __cdecl _byteswap_ushort(unsigned short _Number); 
unsigned long __cdecl _byteswap_ulong(unsigned long _Number); 
unsigned __int64 __cdecl _byteswap_uint64(unsigned __int64 _Number); 

div_t __cdecl div(int _Numerator, int _Denominator); 
ldiv_t __cdecl ldiv(long _Numerator, long _Denominator); 
lldiv_t __cdecl lldiv(__int64 _Numerator, __int64 _Denominator); 



#pragma warning(push)
#pragma warning(disable: 6540)

unsigned __cdecl _rotl(unsigned _Value, int _Shift); 
#line 317
unsigned long __cdecl _lrotl(unsigned long _Value, int _Shift); 




unsigned __int64 __cdecl _rotl64(unsigned __int64 _Value, int _Shift); 




unsigned __cdecl _rotr(unsigned _Value, int _Shift); 
#line 333
unsigned long __cdecl _lrotr(unsigned long _Value, int _Shift); 




unsigned __int64 __cdecl _rotr64(unsigned __int64 _Value, int _Shift); 




#pragma warning(pop)
#line 350
void __cdecl srand(unsigned _Seed); 

int __cdecl rand(void); 
#line 394
#pragma pack ( push, 4 )



typedef 
#line 396
struct { 
unsigned char ld[10]; 
} _LDOUBLE; 
#pragma pack ( pop )
#line 418
typedef 
#line 416
struct { 
double x; 
} _CRT_DOUBLE; 




typedef 
#line 421
struct { 
float f; 
} _CRT_FLOAT; 
#line 432
typedef 
#line 430
struct { 
long double x; 
} _LONGDOUBLE; 



#pragma pack ( push, 4 )



typedef 
#line 438
struct { 
unsigned char ld12[12]; 
} _LDBL12; 
#pragma pack ( pop )
#line 450
double __cdecl atof(const char * _String); 
int __cdecl atoi(const char * _String); 
long __cdecl atol(const char * _String); 
__int64 __cdecl atoll(const char * _String); 
__int64 __cdecl _atoi64(const char * _String); 

double __cdecl _atof_l(const char * _String, _locale_t _Locale); 
int __cdecl _atoi_l(const char * _String, _locale_t _Locale); 
long __cdecl _atol_l(const char * _String, _locale_t _Locale); 
__int64 __cdecl _atoll_l(const char * _String, _locale_t _Locale); 
__int64 __cdecl _atoi64_l(const char * _String, _locale_t _Locale); 

int __cdecl _atoflt(_CRT_FLOAT * _Result, const char * _String); 
int __cdecl _atodbl(_CRT_DOUBLE * _Result, char * _String); 
int __cdecl _atoldbl(_LDOUBLE * _Result, char * _String); 


int __cdecl _atoflt_l(_CRT_FLOAT * _Result, const char * _String, _locale_t _Locale); 
#line 474
int __cdecl _atodbl_l(_CRT_DOUBLE * _Result, char * _String, _locale_t _Locale); 
#line 482
int __cdecl _atoldbl_l(_LDOUBLE * _Result, char * _String, _locale_t _Locale); 
#line 489
float __cdecl strtof(const char * _String, char ** _EndPtr); 
#line 495
float __cdecl _strtof_l(const char * _String, char ** _EndPtr, _locale_t _Locale); 
#line 502
double __cdecl strtod(const char * _String, char ** _EndPtr); 
#line 508
double __cdecl _strtod_l(const char * _String, char ** _EndPtr, _locale_t _Locale); 
#line 515
long double __cdecl strtold(const char * _String, char ** _EndPtr); 
#line 521
long double __cdecl _strtold_l(const char * _String, char ** _EndPtr, _locale_t _Locale); 
#line 528
long __cdecl strtol(const char * _String, char ** _EndPtr, int _Radix); 
#line 535
long __cdecl _strtol_l(const char * _String, char ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 543
__int64 __cdecl strtoll(const char * _String, char ** _EndPtr, int _Radix); 
#line 550
__int64 __cdecl _strtoll_l(const char * _String, char ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 558
unsigned long __cdecl strtoul(const char * _String, char ** _EndPtr, int _Radix); 
#line 565
unsigned long __cdecl _strtoul_l(const char * _String, char ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 573
unsigned __int64 __cdecl strtoull(const char * _String, char ** _EndPtr, int _Radix); 
#line 580
unsigned __int64 __cdecl _strtoull_l(const char * _String, char ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 588
__int64 __cdecl _strtoi64(const char * _String, char ** _EndPtr, int _Radix); 
#line 595
__int64 __cdecl _strtoi64_l(const char * _String, char ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 603
unsigned __int64 __cdecl _strtoui64(const char * _String, char ** _EndPtr, int _Radix); 
#line 610
unsigned __int64 __cdecl _strtoui64_l(const char * _String, char ** _EndPtr, int _Radix, _locale_t _Locale); 
#line 626
errno_t __cdecl _itoa_s(int _Value, char * _Buffer, size_t _BufferCount, int _Radix); 
#line 641
char *__cdecl _itoa(int _Value, char * _Buffer, int _Radix); 
#line 650
errno_t __cdecl _ltoa_s(long _Value, char * _Buffer, size_t _BufferCount, int _Radix); 
#line 664
char *__cdecl _ltoa(long _Value, char * _Buffer, int _Radix); 
#line 673
errno_t __cdecl _ultoa_s(unsigned long _Value, char * _Buffer, size_t _BufferCount, int _Radix); 
#line 687
char *__cdecl _ultoa(unsigned long _Value, char * _Buffer, int _Radix); 
#line 696
errno_t __cdecl _i64toa_s(__int64 _Value, char * _Buffer, size_t _BufferCount, int _Radix); 
#line 705
char *__cdecl _i64toa(__int64 _Value, char * _Buffer, int _Radix); 
#line 713
errno_t __cdecl _ui64toa_s(unsigned __int64 _Value, char * _Buffer, size_t _BufferCount, int _Radix); 
#line 721
char *__cdecl _ui64toa(unsigned __int64 _Value, char * _Buffer, int _Radix); 
#line 741
errno_t __cdecl _ecvt_s(char * _Buffer, size_t _BufferCount, double _Value, int _DigitCount, int * _PtDec, int * _PtSign); 
#line 760
char *__cdecl _ecvt(double _Value, int _DigitCount, int * _PtDec, int * _PtSign); 
#line 769
errno_t __cdecl _fcvt_s(char * _Buffer, size_t _BufferCount, double _Value, int _FractionalDigitCount, int * _PtDec, int * _PtSign); 
#line 790
char *__cdecl _fcvt(double _Value, int _FractionalDigitCount, int * _PtDec, int * _PtSign); 
#line 798
errno_t __cdecl _gcvt_s(char * _Buffer, size_t _BufferCount, double _Value, int _DigitCount); 
#line 814
char *__cdecl _gcvt(double _Value, int _DigitCount, char * _Buffer); 
#line 843
int __cdecl ___mb_cur_max_func(void); 


int __cdecl ___mb_cur_max_l_func(_locale_t _Locale); 
#line 852
int __cdecl mblen(const char * _Ch, size_t _MaxCount); 
#line 858
int __cdecl _mblen_l(const char * _Ch, size_t _MaxCount, _locale_t _Locale); 
#line 866
size_t __cdecl _mbstrlen(const char * _String); 
#line 872
size_t __cdecl _mbstrlen_l(const char * _String, _locale_t _Locale); 
#line 879
size_t __cdecl _mbstrnlen(const char * _String, size_t _MaxCount); 
#line 886
size_t __cdecl _mbstrnlen_l(const char * _String, size_t _MaxCount, _locale_t _Locale); 
#line 893
int __cdecl mbtowc(wchar_t * _DstCh, const char * _SrcCh, size_t _SrcSizeInBytes); 
#line 900
int __cdecl _mbtowc_l(wchar_t * _DstCh, const char * _SrcCh, size_t _SrcSizeInBytes, _locale_t _Locale); 
#line 908
errno_t __cdecl mbstowcs_s(size_t * _PtNumOfCharConverted, wchar_t * _DstBuf, size_t _SizeInWords, const char * _SrcBuf, size_t _MaxCount); 
#line 924
size_t __cdecl mbstowcs(wchar_t * _Dest, const char * _Source, size_t _MaxCount); 
#line 932
errno_t __cdecl _mbstowcs_s_l(size_t * _PtNumOfCharConverted, wchar_t * _DstBuf, size_t _SizeInWords, const char * _SrcBuf, size_t _MaxCount, _locale_t _Locale); 
#line 950
size_t __cdecl _mbstowcs_l(wchar_t * _Dest, const char * _Source, size_t _MaxCount, _locale_t _Locale); 
#line 963
int __cdecl wctomb(char * _MbCh, wchar_t _WCh); 
#line 969
int __cdecl _wctomb_l(char * _MbCh, wchar_t _WCh, _locale_t _Locale); 
#line 978
errno_t __cdecl wctomb_s(int * _SizeConverted, char * _MbCh, rsize_t _SizeInBytes, wchar_t _WCh); 
#line 988
errno_t __cdecl _wctomb_s_l(int * _SizeConverted, char * _MbCh, size_t _SizeInBytes, wchar_t _WCh, _locale_t _Locale); 
#line 996
errno_t __cdecl wcstombs_s(size_t * _PtNumOfCharConverted, char * _Dst, size_t _DstSizeInBytes, const wchar_t * _Src, size_t _MaxCountInBytes); 
#line 1012
size_t __cdecl wcstombs(char * _Dest, const wchar_t * _Source, size_t _MaxCount); 
#line 1020
errno_t __cdecl _wcstombs_s_l(size_t * _PtNumOfCharConverted, char * _Dst, size_t _DstSizeInBytes, const wchar_t * _Src, size_t _MaxCountInBytes, _locale_t _Locale); 
#line 1038
size_t __cdecl _wcstombs_l(char * _Dest, const wchar_t * _Source, size_t _MaxCount, _locale_t _Locale); 
#line 1068
__declspec(allocator) char *__cdecl _fullpath(char * _Buffer, const char * _Path, size_t _BufferCount); 
#line 1077
errno_t __cdecl _makepath_s(char * _Buffer, size_t _BufferCount, const char * _Drive, const char * _Dir, const char * _Filename, const char * _Ext); 
#line 1095
void __cdecl _makepath(char * _Buffer, const char * _Drive, const char * _Dir, const char * _Filename, const char * _Ext); 
#line 1105
void __cdecl _splitpath(const char * _FullPath, char * _Drive, char * _Dir, char * _Filename, char * _Ext); 
#line 1114
errno_t __cdecl _splitpath_s(const char * _FullPath, char * _Drive, size_t _DriveCount, char * _Dir, size_t _DirCount, char * _Filename, size_t _FilenameCount, char * _Ext, size_t _ExtCount); 
#line 1132
errno_t __cdecl getenv_s(size_t * _RequiredCount, char * _Buffer, rsize_t _BufferCount, const char * _VarName); 
#line 1144
int *__cdecl __p___argc(void); 
char ***__cdecl __p___argv(void); 
wchar_t ***__cdecl __p___wargv(void); 
#line 1158
char ***__cdecl __p__environ(void); 
wchar_t ***__cdecl __p__wenviron(void); 
#line 1184
char *__cdecl getenv(const char * _VarName); 
#line 1201
errno_t __cdecl _dupenv_s(char ** _Buffer, size_t * _BufferCount, const char * _VarName); 
#line 1211
int __cdecl system(const char * _Command); 





#pragma warning(push)
#pragma warning(disable: 6540)


int __cdecl _putenv(const char * _EnvString); 




errno_t __cdecl _putenv_s(const char * _Name, const char * _Value); 




#pragma warning(pop)

errno_t __cdecl _searchenv_s(const char * _Filename, const char * _VarName, char * _Buffer, size_t _BufferCount); 
#line 1247
void __cdecl _searchenv(const char * _Filename, const char * _VarName, char * _Buffer); 
#line 1255
__declspec(deprecated("This function or variable has been superceded by newer library or operating system functionality. Consider using SetErrorMode in" "stead. See online help for details.")) void __cdecl 
_seterrormode(int _Mode); 



__declspec(deprecated("This function or variable has been superceded by newer library or operating system functionality. Consider using Beep instead. S" "ee online help for details.")) void __cdecl 
_beep(unsigned _Frequency, unsigned _Duration); 




__declspec(deprecated("This function or variable has been superceded by newer library or operating system functionality. Consider using Sleep instead. " "See online help for details.")) void __cdecl 
_sleep(unsigned long _Duration); 
#line 1289
#pragma warning(push)
#pragma warning(disable: 4141)

__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _ecvt. See online help for details.")) char *__cdecl 
ecvt(double _Value, int _DigitCount, int * _PtDec, int * _PtSign); 
#line 1300
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _fcvt. See online help for details.")) char *__cdecl 
fcvt(double _Value, int _FractionalDigitCount, int * _PtDec, int * _PtSign); 
#line 1308
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _gcvt. See online help for details.")) char *__cdecl 
gcvt(double _Value, int _DigitCount, char * _DstBuf); 
#line 1315
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _itoa. See online help for details.")) char *__cdecl 
itoa(int _Value, char * _Buffer, int _Radix); 
#line 1322
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _ltoa. See online help for details.")) char *__cdecl 
ltoa(long _Value, char * _Buffer, int _Radix); 
#line 1330
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _swab. See online help for details.")) void __cdecl 
swab(char * _Buf1, char * _Buf2, int _SizeInBytes); 
#line 1337
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _ultoa. See online help for details." "")) char *__cdecl 
ultoa(unsigned long _Value, char * _Buffer, int _Radix); 
#line 1346
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _putenv. See online help for details" ".")) int __cdecl 
putenv(const char * _EnvString); 



#pragma warning(pop)

_onexit_t __cdecl onexit(_onexit_t _Func); 
#line 1359
__pragma( pack ( pop )) 

#pragma warning(pop)
#endif /* _INC_STDLIB */
#line 13 "C:\\Users\\trist\\Documents\\GitHub\\simulink-canards\\embedded-coder\\codegen\\lib\\GNC_codegen_SIL\\instrumented\\GNC_codegen_SIL.h"
extern void GNC_codegen_SIL_initialize(GNC_codegen_SILStackData * b_SD); 

extern void GNC_codegen_SIL_terminate(void); 

extern void controller_codegen_entry(GNC_codegen_SILStackData * b_SD, double b_time, double dt_ctrl, const double  xR[2], double pdyn, double delta, const struct0_T * ctrl_mem_in, double * u, double * r, struct0_T * ctrl_mem_out); 
#line 24
extern void navigation_codegen_entry(GNC_codegen_SILStackData * b_SD, double dt, boolean_T flight_phase, const double  x[11], const double  P[121], const struct1_T * bias, const struct2_T * sens_filt, const struct3_T * sens_input, double  x_ret[11], double  P_ret[121], struct1_T * bias_ret, struct2_T * sens_filt_ret, struct6_T * airdata, double  roll_state[2]); 
#line 35
#endif /* GNC_CODEGEN_SIL_H */
#line 9 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_math.h"
#ifndef _INC_MATH
#define _INC_MATH



#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 




struct _exception { 

int type; 
char *name; 
double arg1; 
double arg2; 
double retval; 
}; 



#ifndef _COMPLEX_DEFINED
#define _COMPLEX_DEFINED

struct _complex { 

double x, y; 
}; 





#endif /* _COMPLEX_DEFINED */
#line 59
typedef float float_t; 
typedef double double_t; 
#line 78
extern const double _HUGE; 





#ifndef _HUGE_ENUF
#define _HUGE_ENUF 1e+300
#endif /* _HUGE_ENUF */
#line 175
void __cdecl _fperrraise(int _Except); 

short __cdecl _dclass(double _X); 
short __cdecl _ldclass(long double _X); 
short __cdecl _fdclass(float _X); 

int __cdecl _dsign(double _X); 
int __cdecl _ldsign(long double _X); 
int __cdecl _fdsign(float _X); 

int __cdecl _dpcomp(double _X, double _Y); 
int __cdecl _ldpcomp(long double _X, long double _Y); 
int __cdecl _fdpcomp(float _X, float _Y); 

short __cdecl _dtest(double * _Px); 
short __cdecl _ldtest(long double * _Px); 
short __cdecl _fdtest(float * _Px); 

short __cdecl _d_int(double * _Px, short _Xexp); 
short __cdecl _ld_int(long double * _Px, short _Xexp); 
short __cdecl _fd_int(float * _Px, short _Xexp); 

short __cdecl _dscale(double * _Px, long _Lexp); 
short __cdecl _ldscale(long double * _Px, long _Lexp); 
short __cdecl _fdscale(float * _Px, long _Lexp); 

short __cdecl _dunscale(short * _Pex, double * _Px); 
short __cdecl _ldunscale(short * _Pex, long double * _Px); 
short __cdecl _fdunscale(short * _Pex, float * _Px); 

short __cdecl _dexp(double * _Px, double _Y, long _Eoff); 
short __cdecl _ldexp(long double * _Px, long double _Y, long _Eoff); 
short __cdecl _fdexp(float * _Px, float _Y, long _Eoff); 

short __cdecl _dnorm(unsigned short * _Ps); 
short __cdecl _fdnorm(unsigned short * _Ps); 

double __cdecl _dpoly(double _X, const double * _Tab, int _N); 
long double __cdecl _ldpoly(long double _X, const long double * _Tab, int _N); 
float __cdecl _fdpoly(float _X, const float * _Tab, int _N); 

double __cdecl _dlog(double _X, int _Baseflag); 
long double __cdecl _ldlog(long double _X, int _Baseflag); 
float __cdecl _fdlog(float _X, int _Baseflag); 

double __cdecl _dsin(double _X, unsigned _Qoff); 
long double __cdecl _ldsin(long double _X, unsigned _Qoff); 
float __cdecl _fdsin(float _X, unsigned _Qoff); 
#line 229
typedef 
#line 226
union { 
unsigned short _Sh[4]; 
double _Val; 
} _double_val; 
#line 236
typedef 
#line 233
union { 
unsigned short _Sh[2]; 
float _Val; 
} _float_val; 
#line 243
typedef 
#line 240
union { 
unsigned short _Sh[4]; 
long double _Val; 
} _ldouble_val; 
#line 251
typedef 
#line 246
union { 
unsigned short _Word[4]; 
float _Float; 
double _Double; 
long double _Long_double; 
} _float_const; 

extern const _float_const _Denorm_C, _Inf_C, _Nan_C, _Snan_C, _Hugeval_C; 
extern const _float_const _FDenorm_C, _FInf_C, _FNan_C, _FSnan_C; 
extern const _float_const _LDenorm_C, _LInf_C, _LNan_C, _LSnan_C; 

extern const _float_const _Eps_C, _Rteps_C; 
extern const _float_const _FEps_C, _FRteps_C; 
extern const _float_const _LEps_C, _LRteps_C; 

extern const double _Zero_C, _Xbig_C; 
extern const float _FZero_C, _FXbig_C; 
extern const long double _LZero_C, _LXbig_C; 
#line 470
int __cdecl abs(int _X); 
long __cdecl labs(long _X); 
__int64 __cdecl llabs(__int64 _X); 

double __cdecl acos(double _X); 
double __cdecl asin(double _X); 
double __cdecl atan(double _X); 
double __cdecl atan2(double _Y, double _X); 

double __cdecl cos(double _X); 
double __cdecl cosh(double _X); 
double __cdecl exp(double _X); 
double __cdecl fabs(double _X); 
double __cdecl fmod(double _X, double _Y); 
double __cdecl log(double _X); 
double __cdecl log10(double _X); 
double __cdecl pow(double _X, double _Y); 
double __cdecl sin(double _X); 
double __cdecl sinh(double _X); 
double __cdecl sqrt(double _X); 
double __cdecl tan(double _X); 
double __cdecl tanh(double _X); 

double __cdecl acosh(double _X); 
double __cdecl asinh(double _X); 
double __cdecl atanh(double _X); 
double __cdecl atof(const char * _String); 
double __cdecl _atof_l(const char * _String, _locale_t _Locale); 
double __cdecl _cabs(struct _complex _Complex_value); 
double __cdecl cbrt(double _X); 
double __cdecl ceil(double _X); 
double __cdecl _chgsign(double _X); 
double __cdecl copysign(double _Number, double _Sign); 
double __cdecl _copysign(double _Number, double _Sign); 
double __cdecl erf(double _X); 
double __cdecl erfc(double _X); 
double __cdecl exp2(double _X); 
double __cdecl expm1(double _X); 
double __cdecl fdim(double _X, double _Y); 
double __cdecl floor(double _X); 
double __cdecl fma(double _X, double _Y, double _Z); 
double __cdecl fmax(double _X, double _Y); 
double __cdecl fmin(double _X, double _Y); 
double __cdecl frexp(double _X, int * _Y); 
double __cdecl hypot(double _X, double _Y); 
double __cdecl _hypot(double _X, double _Y); 
int __cdecl ilogb(double _X); 
double __cdecl ldexp(double _X, int _Y); 
double __cdecl lgamma(double _X); 
__int64 __cdecl llrint(double _X); 
__int64 __cdecl llround(double _X); 
double __cdecl log1p(double _X); 
double __cdecl log2(double _X); 
double __cdecl logb(double _X); 
long __cdecl lrint(double _X); 
long __cdecl lround(double _X); 

int __cdecl _matherr(struct _exception * _Except); 

double __cdecl modf(double _X, double * _Y); 
double __cdecl nan(const char * _X); 
double __cdecl nearbyint(double _X); 
double __cdecl nextafter(double _X, double _Y); 
double __cdecl nexttoward(double _X, long double _Y); 
double __cdecl remainder(double _X, double _Y); 
double __cdecl remquo(double _X, double _Y, int * _Z); 
double __cdecl rint(double _X); 
double __cdecl round(double _X); 
double __cdecl scalbln(double _X, long _Y); 
double __cdecl scalbn(double _X, int _Y); 
double __cdecl tgamma(double _X); 
double __cdecl trunc(double _X); 
double __cdecl _j0(double _X); 
double __cdecl _j1(double _X); 
double __cdecl _jn(int _X, double _Y); 
double __cdecl _y0(double _X); 
double __cdecl _y1(double _X); 
double __cdecl _yn(int _X, double _Y); 

float __cdecl acoshf(float _X); 
float __cdecl asinhf(float _X); 
float __cdecl atanhf(float _X); 
float __cdecl cbrtf(float _X); 
float __cdecl _chgsignf(float _X); 
float __cdecl copysignf(float _Number, float _Sign); 
float __cdecl _copysignf(float _Number, float _Sign); 
float __cdecl erff(float _X); 
float __cdecl erfcf(float _X); 
float __cdecl expm1f(float _X); 
float __cdecl exp2f(float _X); 
float __cdecl fdimf(float _X, float _Y); 
float __cdecl fmaf(float _X, float _Y, float _Z); 
float __cdecl fmaxf(float _X, float _Y); 
float __cdecl fminf(float _X, float _Y); 
float __cdecl _hypotf(float _X, float _Y); 
int __cdecl ilogbf(float _X); 
float __cdecl lgammaf(float _X); 
__int64 __cdecl llrintf(float _X); 
__int64 __cdecl llroundf(float _X); 
float __cdecl log1pf(float _X); 
float __cdecl log2f(float _X); 
float __cdecl logbf(float _X); 
long __cdecl lrintf(float _X); 
long __cdecl lroundf(float _X); 
float __cdecl nanf(const char * _X); 
float __cdecl nearbyintf(float _X); 
float __cdecl nextafterf(float _X, float _Y); 
float __cdecl nexttowardf(float _X, long double _Y); 
float __cdecl remainderf(float _X, float _Y); 
float __cdecl remquof(float _X, float _Y, int * _Z); 
float __cdecl rintf(float _X); 
float __cdecl roundf(float _X); 
float __cdecl scalblnf(float _X, long _Y); 
float __cdecl scalbnf(float _X, int _Y); 
float __cdecl tgammaf(float _X); 
float __cdecl truncf(float _X); 
#line 595
float __cdecl _logbf(float _X); 
float __cdecl _nextafterf(float _X, float _Y); 
int __cdecl _finitef(float _X); 
int __cdecl _isnanf(float _X); 
int __cdecl _fpclassf(float _X); 

int __cdecl _set_FMA3_enable(int _Flag); 
int __cdecl _get_FMA3_enable(void); 
#line 615
float __cdecl acosf(float _X); 
float __cdecl asinf(float _X); 
float __cdecl atan2f(float _Y, float _X); 
float __cdecl atanf(float _X); 
float __cdecl ceilf(float _X); 
float __cdecl cosf(float _X); 
float __cdecl coshf(float _X); 
float __cdecl expf(float _X); 
#line 678
__inline float __cdecl fabsf(float _X) 
{ 
return (float)fabs(_X); 
} 
#line 687
float __cdecl floorf(float _X); 
float __cdecl fmodf(float _X, float _Y); 
#line 704
__inline float __cdecl frexpf(float _X, int *_Y) 
{ 
return (float)frexp(_X, _Y); 
} 

__inline float __cdecl hypotf(float _X, float _Y) 
{ 
return _hypotf(_X, _Y); 
} 

__inline float __cdecl ldexpf(float _X, int _Y) 
{ 
return (float)ldexp(_X, _Y); 
} 



float __cdecl log10f(float _X); 
float __cdecl logf(float _X); 
float __cdecl modff(float _X, float * _Y); 
float __cdecl powf(float _X, float _Y); 
float __cdecl sinf(float _X); 
float __cdecl sinhf(float _X); 
float __cdecl sqrtf(float _X); 
float __cdecl tanf(float _X); 
float __cdecl tanhf(float _X); 
#line 783
long double __cdecl acoshl(long double _X); 

__inline long double __cdecl acosl(long double _X) 
{ 
return acos((double)_X); 
} 

long double __cdecl asinhl(long double _X); 

__inline long double __cdecl asinl(long double _X) 
{ 
return asin((double)_X); 
} 

__inline long double __cdecl atan2l(long double _Y, long double _X) 
{ 
return atan2((double)_Y, (double)_X); 
} 

long double __cdecl atanhl(long double _X); 

__inline long double __cdecl atanl(long double _X) 
{ 
return atan((double)_X); 
} 

long double __cdecl cbrtl(long double _X); 

__inline long double __cdecl ceill(long double _X) 
{ 
return ceil((double)_X); 
} 

__inline long double __cdecl _chgsignl(long double _X) 
{ 
return _chgsign((double)_X); 
} 

long double __cdecl copysignl(long double _Number, long double _Sign); 

__inline long double __cdecl _copysignl(long double _Number, long double _Sign) 
{ 
return _copysign((double)_Number, (double)_Sign); 
} 

__inline long double __cdecl coshl(long double _X) 
{ 
return cosh((double)_X); 
} 

__inline long double __cdecl cosl(long double _X) 
{ 
return cos((double)_X); 
} 

long double __cdecl erfl(long double _X); 
long double __cdecl erfcl(long double _X); 

__inline long double __cdecl expl(long double _X) 
{ 
return exp((double)_X); 
} 

long double __cdecl exp2l(long double _X); 
long double __cdecl expm1l(long double _X); 

__inline long double __cdecl fabsl(long double _X) 
{ 
return fabs((double)_X); 
} 

long double __cdecl fdiml(long double _X, long double _Y); 

__inline long double __cdecl floorl(long double _X) 
{ 
return floor((double)_X); 
} 

long double __cdecl fmal(long double _X, long double _Y, long double _Z); 
long double __cdecl fmaxl(long double _X, long double _Y); 
long double __cdecl fminl(long double _X, long double _Y); 

__inline long double __cdecl fmodl(long double _X, long double _Y) 
{ 
return fmod((double)_X, (double)_Y); 
} 

__inline long double __cdecl frexpl(long double _X, int *_Y) 
{ 
return frexp((double)_X, _Y); 
} 

int __cdecl ilogbl(long double _X); 

__inline long double __cdecl _hypotl(long double _X, long double _Y) 
{ 
return _hypot((double)_X, (double)_Y); 
} 

__inline long double __cdecl hypotl(long double _X, long double _Y) 
{ 
return _hypot((double)_X, (double)_Y); 
} 

__inline long double __cdecl ldexpl(long double _X, int _Y) 
{ 
return ldexp((double)_X, _Y); 
} 

long double __cdecl lgammal(long double _X); 
__int64 __cdecl llrintl(long double _X); 
__int64 __cdecl llroundl(long double _X); 

__inline long double __cdecl logl(long double _X) 
{ 
return log((double)_X); 
} 

__inline long double __cdecl log10l(long double _X) 
{ 
return log10((double)_X); 
} 

long double __cdecl log1pl(long double _X); 
long double __cdecl log2l(long double _X); 
long double __cdecl logbl(long double _X); 
long __cdecl lrintl(long double _X); 
long __cdecl lroundl(long double _X); 

__inline long double __cdecl modfl(long double _X, long double *_Y) 
{ 
double _F, _I; 
_F = modf((double)_X, &_I); 
*_Y = _I; 
return _F; 
} 

long double __cdecl nanl(const char * _X); 
long double __cdecl nearbyintl(long double _X); 
long double __cdecl nextafterl(long double _X, long double _Y); 
long double __cdecl nexttowardl(long double _X, long double _Y); 

__inline long double __cdecl powl(long double _X, long double _Y) 
{ 
return pow((double)_X, (double)_Y); 
} 

long double __cdecl remainderl(long double _X, long double _Y); 
long double __cdecl remquol(long double _X, long double _Y, int * _Z); 
long double __cdecl rintl(long double _X); 
long double __cdecl roundl(long double _X); 
long double __cdecl scalblnl(long double _X, long _Y); 
long double __cdecl scalbnl(long double _X, int _Y); 

__inline long double __cdecl sinhl(long double _X) 
{ 
return sinh((double)_X); 
} 

__inline long double __cdecl sinl(long double _X) 
{ 
return sin((double)_X); 
} 

__inline long double __cdecl sqrtl(long double _X) 
{ 
return sqrt((double)_X); 
} 

__inline long double __cdecl tanhl(long double _X) 
{ 
return tanh((double)_X); 
} 

__inline long double __cdecl tanl(long double _X) 
{ 
return tan((double)_X); 
} 

long double __cdecl tgammal(long double _X); 
long double __cdecl truncl(long double _X); 
#line 984
extern double HUGE; 




__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _j0. See online help for details.")) double __cdecl j0(double _X); 
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _j1. See online help for details.")) double __cdecl j1(double _X); 
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _jn. See online help for details.")) double __cdecl jn(int _X, double _Y); 
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _y0. See online help for details.")) double __cdecl y0(double _X); 
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _y1. See online help for details.")) double __cdecl y1(double _X); 
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _yn. See online help for details.")) double __cdecl yn(int _X, double _Y); 




__pragma( pack ( pop )) 

#pragma warning(pop)
#endif /* _INC_MATH */
#line 9 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\string.h"
#ifndef _INC_STRING
#define _INC_STRING
#line 9 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\errno.h"
#ifndef _INC_ERRNO
#define _INC_ERRNO



#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 




int *__cdecl _errno(void); 


errno_t __cdecl _set_errno(int _Value); 
errno_t __cdecl _get_errno(int * _Value); 

unsigned long *__cdecl __doserrno(void); 


errno_t __cdecl _set_doserrno(unsigned long _Value); 
errno_t __cdecl _get_doserrno(unsigned long * _Value); 
#line 134
__pragma( pack ( pop )) 

#pragma warning(pop)
#endif /* _INC_ERRNO */
#line 12 "C:\\Program Files (x86)\\Microsoft Visual Studio\\2022\\BuildTools\\VC\\Tools\\MSVC\\14.41.34120\\include\\vcruntime_string.h"
#pragma warning(push)
#pragma warning(disable: 4514 4820)



__pragma( pack ( push, 8 )) 




void *__cdecl memchr(const void * _Buf, int _Val, size_t _MaxCount); 
#line 29
int __cdecl memcmp(const void * _Buf1, const void * _Buf2, size_t _Size); 
#line 43
void *__cdecl memcpy(void * _Dst, const void * _Src, size_t _Size); 
#line 50
void *__cdecl memmove(void * _Dst, const void * _Src, size_t _Size); 
#line 63
void *__cdecl memset(void * _Dst, int _Val, size_t _Size); 
#line 70
char *__cdecl strchr(const char * _Str, int _Val); 
#line 76
char *__cdecl strrchr(const char * _Str, int _Ch); 
#line 82
char *__cdecl strstr(const char * _Str, const char * _SubStr); 
#line 89
wchar_t *__cdecl wcschr(const wchar_t * _Str, wchar_t _Ch); 
#line 95
wchar_t *__cdecl wcsrchr(const wchar_t * _Str, wchar_t _Ch); 
#line 102
wchar_t *__cdecl wcsstr(const wchar_t * _Str, const wchar_t * _SubStr); 
#line 109
__pragma( pack ( pop )) 



#pragma warning(pop)
#line 14 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_memcpy_s.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 

#ifndef _CRT_MEMCPY_S_INLINE
#define _CRT_MEMCPY_S_INLINE static __inline
#endif /* _CRT_MEMCPY_S_INLINE */
#line 39
static __inline errno_t __cdecl memcpy_s(void *const 
_Destination, const rsize_t 
_DestinationSize, const void *const 
_Source, const rsize_t 
_SourceSize) 

{ { 
if (_SourceSize == (0)) 
{ 
return 0; 
}  } 

{ int _Expr_val = !(!(_Destination != (void *)0)); { if (!_Expr_val) { *_errno() = 22; _invalid_parameter_noinfo(); return 22; }  } } ; { 
if ((_Source == (void *)0) || (_DestinationSize < _SourceSize)) 
{ 
memset(_Destination, 0, _DestinationSize); 

{ int _Expr_val = !(!(_Source != (void *)0)); { if (!_Expr_val) { *_errno() = 22; _invalid_parameter_noinfo(); return 22; }  } } ; 
{ int _Expr_val = !(!(_DestinationSize >= _SourceSize)); { if (!_Expr_val) { *_errno() = 34; _invalid_parameter_noinfo(); return 34; }  } } ; 


return 22; 
}  } 
memcpy(_Destination, _Source, _SourceSize); 
return 0; 
} 


static __inline errno_t __cdecl memmove_s(void *const 
_Destination, const rsize_t 
_DestinationSize, const void *const 
_Source, const rsize_t 
_SourceSize) 

{ { 
if (_SourceSize == (0)) 
{ 
return 0; 
}  } 

{ int _Expr_val = !(!(_Destination != (void *)0)); { if (!_Expr_val) { *_errno() = 22; _invalid_parameter_noinfo(); return 22; }  } } ; 
{ int _Expr_val = !(!(_Source != (void *)0)); { if (!_Expr_val) { *_errno() = 22; _invalid_parameter_noinfo(); return 22; }  } } ; 
{ int _Expr_val = !(!(_DestinationSize >= _SourceSize)); { if (!_Expr_val) { *_errno() = 34; _invalid_parameter_noinfo(); return 34; }  } } ; 

memmove(_Destination, _Source, _SourceSize); 
return 0; 
} 
#line 92
#pragma warning(pop)
__pragma( pack ( pop )) 
#line 17 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_memory.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)




__pragma( pack ( push, 8 )) 




int __cdecl _memicmp(const void * _Buf1, const void * _Buf2, size_t _Size); 
#line 35
int __cdecl _memicmp_l(const void * _Buf1, const void * _Buf2, size_t _Size, _locale_t _Locale); 
#line 82
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _memccpy. See online help for detail" "s.")) void *__cdecl 
memccpy(void * _Dst, const void * _Src, int _Val, size_t _Size); 
#line 90
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _memicmp. See online help for detail" "s.")) int __cdecl 
memicmp(const void * _Buf1, const void * _Buf2, size_t _Size); 
#line 118
__pragma( pack ( pop )) 



#pragma warning(pop)
#line 14 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\corecrt_wstring.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)




__pragma( pack ( push, 8 )) 
#line 32
errno_t __cdecl wcscat_s(wchar_t * _Destination, rsize_t _SizeInWords, const wchar_t * _Source); 
#line 39
errno_t __cdecl wcscpy_s(wchar_t * _Destination, rsize_t _SizeInWords, const wchar_t * _Source); 
#line 46
errno_t __cdecl wcsncat_s(wchar_t * _Destination, rsize_t _SizeInWords, const wchar_t * _Source, rsize_t _MaxCount); 
#line 54
errno_t __cdecl wcsncpy_s(wchar_t * _Destination, rsize_t _SizeInWords, const wchar_t * _Source, rsize_t _MaxCount); 
#line 62
wchar_t *__cdecl wcstok_s(wchar_t * _String, const wchar_t * _Delimiter, wchar_t ** _Context); 
#line 83
__declspec(allocator) wchar_t *__cdecl _wcsdup(const wchar_t * _String); 
#line 100
wchar_t *__cdecl wcscat(wchar_t * _Destination, const wchar_t * _Source); 
#line 108
int __cdecl wcscmp(const wchar_t * _String1, const wchar_t * _String2); 
#line 119
wchar_t *__cdecl wcscpy(wchar_t * _Destination, const wchar_t * _Source); 
#line 126
size_t __cdecl wcscspn(const wchar_t * _String, const wchar_t * _Control); 
#line 132
size_t __cdecl wcslen(const wchar_t * _String); 
#line 145
size_t __cdecl wcsnlen(const wchar_t * _Source, size_t _MaxCount); 
#line 161
static __inline size_t __cdecl wcsnlen_s(const wchar_t *
_Source, size_t 
_MaxCount) 

{ 
return (_Source == (0)) ? (0) : wcsnlen(_Source, _MaxCount); 
} 
#line 178
wchar_t *__cdecl wcsncat(wchar_t * _Destination, const wchar_t * _Source, size_t _Count); 
#line 187
int __cdecl wcsncmp(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount); 
#line 200
wchar_t *__cdecl wcsncpy(wchar_t * _Destination, const wchar_t * _Source, size_t _Count); 
#line 209
wchar_t *__cdecl wcspbrk(const wchar_t * _String, const wchar_t * _Control); 
#line 215
size_t __cdecl wcsspn(const wchar_t * _String, const wchar_t * _Control); 
#line 221
wchar_t *__cdecl wcstok(wchar_t * _String, const wchar_t * _Delimiter, wchar_t ** _Context); 
#line 239
static __inline wchar_t *__cdecl _wcstok(wchar_t *const 
_String, const wchar_t *const 
_Delimiter) 

{ 
return wcstok(_String, _Delimiter, 0); 
} 
#line 268
wchar_t *__cdecl _wcserror(int _ErrorNumber); 




errno_t __cdecl _wcserror_s(wchar_t * _Buffer, size_t _SizeInWords, int _ErrorNumber); 
#line 288
wchar_t *__cdecl __wcserror(const wchar_t * _String); 



errno_t __cdecl __wcserror_s(wchar_t * _Buffer, size_t _SizeInWords, const wchar_t * _ErrorMessage); 
#line 304
int __cdecl _wcsicmp(const wchar_t * _String1, const wchar_t * _String2); 




int __cdecl _wcsicmp_l(const wchar_t * _String1, const wchar_t * _String2, _locale_t _Locale); 
#line 315
int __cdecl _wcsnicmp(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount); 
#line 321
int __cdecl _wcsnicmp_l(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount, _locale_t _Locale); 
#line 328
errno_t __cdecl _wcsnset_s(wchar_t * _Destination, size_t _SizeInWords, wchar_t _Value, size_t _MaxCount); 
#line 342
wchar_t *__cdecl _wcsnset(wchar_t * _String, wchar_t _Value, size_t _MaxCount); 
#line 350
wchar_t *__cdecl _wcsrev(wchar_t * _String); 



errno_t __cdecl _wcsset_s(wchar_t * _Destination, size_t _SizeInWords, wchar_t _Value); 
#line 366
wchar_t *__cdecl _wcsset(wchar_t * _String, wchar_t _Value); 
#line 373
errno_t __cdecl _wcslwr_s(wchar_t * _String, size_t _SizeInWords); 
#line 383
wchar_t *__cdecl _wcslwr(wchar_t * _String); 
#line 389
errno_t __cdecl _wcslwr_s_l(wchar_t * _String, size_t _SizeInWords, _locale_t _Locale); 
#line 401
wchar_t *__cdecl _wcslwr_l(wchar_t * _String, _locale_t _Locale); 
#line 409
errno_t __cdecl _wcsupr_s(wchar_t * _String, size_t _Size); 
#line 419
wchar_t *__cdecl _wcsupr(wchar_t * _String); 
#line 425
errno_t __cdecl _wcsupr_s_l(wchar_t * _String, size_t _Size, _locale_t _Locale); 
#line 437
wchar_t *__cdecl _wcsupr_l(wchar_t * _String, _locale_t _Locale); 
#line 446
size_t __cdecl wcsxfrm(wchar_t * _Destination, const wchar_t * _Source, size_t _MaxCount); 
#line 454
size_t __cdecl _wcsxfrm_l(wchar_t * _Destination, const wchar_t * _Source, size_t _MaxCount, _locale_t _Locale); 
#line 462
int __cdecl wcscoll(const wchar_t * _String1, const wchar_t * _String2); 
#line 468
int __cdecl _wcscoll_l(const wchar_t * _String1, const wchar_t * _String2, _locale_t _Locale); 
#line 475
int __cdecl _wcsicoll(const wchar_t * _String1, const wchar_t * _String2); 
#line 481
int __cdecl _wcsicoll_l(const wchar_t * _String1, const wchar_t * _String2, _locale_t _Locale); 
#line 488
int __cdecl _wcsncoll(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount); 
#line 495
int __cdecl _wcsncoll_l(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount, _locale_t _Locale); 
#line 503
int __cdecl _wcsnicoll(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount); 
#line 510
int __cdecl _wcsnicoll_l(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount, _locale_t _Locale); 
#line 569
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsdup. See online help for details" ".")) wchar_t *__cdecl 
wcsdup(const wchar_t * _String); 
#line 581
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsicmp. See online help for detail" "s.")) int __cdecl 
wcsicmp(const wchar_t * _String1, const wchar_t * _String2); 




__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsnicmp. See online help for detai" "ls.")) int __cdecl 
wcsnicmp(const wchar_t * _String1, const wchar_t * _String2, size_t _MaxCount); 
#line 594
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsnset. See online help for detail" "s.")) wchar_t *__cdecl 

wcsnset(wchar_t * _String, wchar_t _Value, size_t _MaxCount); 
#line 602
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsrev. See online help for details" ".")) wchar_t *__cdecl 

wcsrev(wchar_t * _String); 



__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsset. See online help for details" ".")) wchar_t *__cdecl 

wcsset(wchar_t * _String, wchar_t _Value); 




__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcslwr. See online help for details" ".")) wchar_t *__cdecl 

wcslwr(wchar_t * _String); 



__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsupr. See online help for details" ".")) wchar_t *__cdecl 

wcsupr(wchar_t * _String); 



__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _wcsicoll. See online help for detai" "ls.")) int __cdecl 
wcsicoll(const wchar_t * _String1, const wchar_t * _String2); 
#line 637
__pragma( pack ( pop )) 



#pragma warning(pop)
#line 19 "C:\\Program Files (x86)\\Windows Kits\\10\\include\\10.0.22621.0\\ucrt\\string.h"
#pragma warning(push)
#pragma warning(disable: 4324 4514 4574 4710 4793 4820 4995 4996 28719 28726 28727)


__pragma( pack ( push, 8 )) 
#line 32
errno_t __cdecl strcpy_s(char * _Destination, rsize_t _SizeInBytes, const char * _Source); 
#line 39
errno_t __cdecl strcat_s(char * _Destination, rsize_t _SizeInBytes, const char * _Source); 
#line 46
errno_t __cdecl strerror_s(char * _Buffer, size_t _SizeInBytes, int _ErrorNumber); 
#line 52
errno_t __cdecl strncat_s(char * _Destination, rsize_t _SizeInBytes, const char * _Source, rsize_t _MaxCount); 
#line 60
errno_t __cdecl strncpy_s(char * _Destination, rsize_t _SizeInBytes, const char * _Source, rsize_t _MaxCount); 
#line 68
char *__cdecl strtok_s(char * _String, const char * _Delimiter, char ** _Context); 
#line 76
void *__cdecl _memccpy(void * _Dst, const void * _Src, int _Val, size_t _MaxCount); 
#line 91
char *__cdecl strcat(char * _Destination, const char * _Source); 
#line 100
int __cdecl strcmp(const char * _Str1, const char * _Str2); 
#line 106
int __cdecl _strcmpi(const char * _String1, const char * _String2); 
#line 112
int __cdecl strcoll(const char * _String1, const char * _String2); 
#line 118
int __cdecl _strcoll_l(const char * _String1, const char * _String2, _locale_t _Locale); 
#line 130
char *__cdecl strcpy(char * _Destination, const char * _Source); 
#line 137
size_t __cdecl strcspn(const char * _Str, const char * _Control); 
#line 148
__declspec(allocator) char *__cdecl _strdup(const char * _Source); 
#line 159
char *__cdecl _strerror(const char * _ErrorMessage); 




errno_t __cdecl _strerror_s(char * _Buffer, size_t _SizeInBytes, const char * _ErrorMessage); 
#line 178
char *__cdecl strerror(int _ErrorMessage); 
#line 189
int __cdecl _stricmp(const char * _String1, const char * _String2); 
#line 195
int __cdecl _stricoll(const char * _String1, const char * _String2); 
#line 201
int __cdecl _stricoll_l(const char * _String1, const char * _String2, _locale_t _Locale); 
#line 208
int __cdecl _stricmp_l(const char * _String1, const char * _String2, _locale_t _Locale); 
#line 215
size_t __cdecl strlen(const char * _Str); 




errno_t __cdecl _strlwr_s(char * _String, size_t _Size); 
#line 230
char *__cdecl _strlwr(char * _String); 
#line 236
errno_t __cdecl _strlwr_s_l(char * _String, size_t _Size, _locale_t _Locale); 
#line 248
char *__cdecl _strlwr_l(char * _String, _locale_t _Locale); 
#line 262
char *__cdecl strncat(char * _Destination, const char * _Source, size_t _Count); 
#line 271
int __cdecl strncmp(const char * _Str1, const char * _Str2, size_t _MaxCount); 
#line 278
int __cdecl _strnicmp(const char * _String1, const char * _String2, size_t _MaxCount); 
#line 285
int __cdecl _strnicmp_l(const char * _String1, const char * _String2, size_t _MaxCount, _locale_t _Locale); 
#line 293
int __cdecl _strnicoll(const char * _String1, const char * _String2, size_t _MaxCount); 
#line 300
int __cdecl _strnicoll_l(const char * _String1, const char * _String2, size_t _MaxCount, _locale_t _Locale); 
#line 308
int __cdecl _strncoll(const char * _String1, const char * _String2, size_t _MaxCount); 
#line 315
int __cdecl _strncoll_l(const char * _String1, const char * _String2, size_t _MaxCount, _locale_t _Locale); 
#line 322
size_t __cdecl __strncnt(const char * _String, size_t _Count); 
#line 334
char *__cdecl strncpy(char * _Destination, const char * _Source, size_t _Count); 
#line 351
size_t __cdecl strnlen(const char * _String, size_t _MaxCount); 
#line 367
static __inline size_t __cdecl strnlen_s(const char *
_String, size_t 
_MaxCount) 

{ 
return (_String == (0)) ? (0) : strnlen(_String, _MaxCount); 
} 




errno_t __cdecl _strnset_s(char * _String, size_t _SizeInBytes, int _Value, size_t _MaxCount); 
#line 392
char *__cdecl _strnset(char * _Destination, int _Value, size_t _Count); 
#line 401
char *__cdecl strpbrk(const char * _Str, const char * _Control); 




char *__cdecl _strrev(char * _Str); 




errno_t __cdecl _strset_s(char * _Destination, size_t _DestinationSize, int _Value); 
#line 423
char *__cdecl _strset(char * _Destination, int _Value); 
#line 430
size_t __cdecl strspn(const char * _Str, const char * _Control); 
#line 436
char *__cdecl strtok(char * _String, const char * _Delimiter); 
#line 442
errno_t __cdecl _strupr_s(char * _String, size_t _Size); 
#line 452
char *__cdecl _strupr(char * _String); 
#line 458
errno_t __cdecl _strupr_s_l(char * _String, size_t _Size, _locale_t _Locale); 
#line 470
char *__cdecl _strupr_l(char * _String, _locale_t _Locale); 
#line 479
size_t __cdecl strxfrm(char * _Destination, const char * _Source, size_t _MaxCount); 
#line 487
size_t __cdecl _strxfrm_l(char * _Destination, const char * _Source, size_t _MaxCount, _locale_t _Locale); 
#line 531
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strdup. See online help for details" ".")) char *__cdecl 
strdup(const char * _String); 
#line 538
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strcmpi. See online help for detail" "s.")) int __cdecl 
strcmpi(const char * _String1, const char * _String2); 




__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _stricmp. See online help for detail" "s.")) int __cdecl 
stricmp(const char * _String1, const char * _String2); 




__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strlwr. See online help for details" ".")) char *__cdecl 
strlwr(char * _String); 



__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strnicmp. See online help for detai" "ls.")) int __cdecl 
strnicmp(const char * _String1, const char * _String2, size_t _MaxCount); 
#line 562
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strnset. See online help for detail" "s.")) char *__cdecl 
strnset(char * _String, int _Value, size_t _MaxCount); 
#line 569
__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strrev. See online help for details" ".")) char *__cdecl 
strrev(char * _String); 



__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strset. See online help for details" ".")) char *__cdecl 
strset(char * _String, int _Value); 



__declspec(deprecated("The POSIX name for this item is deprecated. Instead, use the ISO C and C++ conformant name: _strupr. See online help for details" ".")) char *__cdecl 
strupr(char * _String); 
#line 588
__pragma( pack ( pop )) 

#pragma warning(pop)

#endif /* _INC_STRING */
#line 6 "GNC_codegen_SIL.c"
static const double dv[9] = {(0.46), (0.0), (0.0), (0.0), (49.5), (0.0), (0.0), (0.0), (49.5)}; 

static const double dv1[9] = {(2.1739130434782608), (0.0), (0.0), (0.0), (0.020202020202020204), (0.0), (0.0), (0.0), (0.020202020202020204)}; 



static double airdata_atmos(double altitude, double * airdata_temperature, double * airdata_density, double * airdata_sonic_speed, double * airdata_mach, double * airdata_dynamic_pressure); 




static void b_ekf_correct(const double  x[11], const double  P[121], double y, double b, double  x_new[11], double  P_new[121]); 


static double b_norm(const double  x[3]); 

static void controller_codegen_entry_init(GNC_codegen_SILStackData * b_SD); 

static void dynamics_init(GNC_codegen_SILStackData * b_SD); 

static void dynamics_jacobian_init(GNC_codegen_SILStackData * b_SD); 

static void ekf_correct(const double  x[11], const double  P[121], const double  y[3], const double  b[3], const double  R[9], double  x_new[11], double  P_new[121]); 



static void inv(const double  x[9], double  y[9]); 

static void pad_filter_init(GNC_codegen_SILStackData * b_SD); 

static double airdata_atmos(double altitude, double *airdata_temperature, double *
airdata_density, double *
airdata_sonic_speed, double *airdata_mach, double *
airdata_dynamic_pressure) { 
double airdata_pressure; 
double layer_idx_1; 
double layer_idx_2; 
double layer_idx_3; 
double pressure; 
double temperature; 
int layer_idx_0; 
altitude = ((6.356766E+6) * altitude) / ((6.356766E+6) - altitude); 
layer_idx_0 = 0; 
layer_idx_1 = (101325.0); 
layer_idx_2 = (288.15); 
layer_idx_3 = (0.0065); { 
if (altitude > (11000.0)) { { 
if (altitude < (20000.0)) { 
layer_idx_0 = 11000; 
layer_idx_2 = (216.65); 
layer_idx_3 = (0.0); 
profileStart_GNC_codegen_SIL(1U); pressure = (22632.1) * exp((-(9.81) * (altitude - (11000.0))) / (62191.094035)); profileEnd_GNC_codegen_SIL(1U); 
} else { { 
if (altitude < (32000.0)) { 
layer_idx_0 = 20000; 
layer_idx_1 = (5474.9); 
layer_idx_2 = (216.65); 
layer_idx_3 = -(0.001); 
} else { 
layer_idx_0 = 32000; 
layer_idx_1 = (868.02); 
layer_idx_2 = (228.65); 
layer_idx_3 = -(0.0028); 
}  } 
profileStart_GNC_codegen_SIL(2U); pressure = layer_idx_1 * pow((1.0) - (layer_idx_3 / layer_idx_2) * (altitude - (double)layer_idx_0), (9.81) / ((287.0579) * layer_idx_3)); profileEnd_GNC_codegen_SIL(2U); 


}  } 
} else { 
profileStart_GNC_codegen_SIL(3U); pressure = layer_idx_1 * pow((1.0) - (layer_idx_3 / layer_idx_2) * (altitude - (double)layer_idx_0), (9.81) / ((287.0579) * layer_idx_3)); profileEnd_GNC_codegen_SIL(3U); 


}  } 
temperature = layer_idx_2 - layer_idx_3 * (altitude - (double)layer_idx_0); 
*airdata_density = pressure / ((287.0579) * temperature); 
profileStart_GNC_codegen_SIL(4U); *airdata_sonic_speed = sqrt((401.88106) * temperature); profileEnd_GNC_codegen_SIL(4U); 
airdata_pressure = pressure; 
*airdata_temperature = temperature; 
*airdata_mach = (0.0); 
*airdata_dynamic_pressure = (0.0); 
return airdata_pressure; 
} 

static void b_ekf_correct(const double x[11], const double P[121], double y, double 
b, double x_new[11], double P_new[121]) { 
double E[121]; 
double b_E[121]; 
double b_K[121]; 
double c_E[121]; 
double H[11]; 
double K[11]; 
double b_H[11]; 
double b_P[11]; 
double absxk; 
double airdata_altitude_pressure; 
double altitude; 
double altitude_ratio; 
double b_b; 
double b_expl_temp; 
double c_H; 
double c_b; 
double c_expl_temp; 
double d_expl_temp; 
double e_expl_temp; 
double expl_temp; 
double layer_idx_1; 
double layer_idx_2; 
double layer_idx_3; 
double q_mag; 
double scale; 
double t; 
double t0_pressure; 
int b_i; 
int i; 
int i1; 
int i10; 
int i11; 
int i12; 
int i13; 
int i2; 
int i3; 
int i4; 
int i5; 
int i6; 
int i7; 
int i8; 
int i9; 
int k; 
int layer_idx_0; 
char b_I[121]; 
profileStart_GNC_codegen_SIL(5U); t0_pressure = airdata_atmos(x[10], &expl_temp, &b_expl_temp, &c_expl_temp, &d_expl_temp, &e_expl_temp); profileEnd_GNC_codegen_SIL(5U); 

b_b = y - (t0_pressure + b); 
profileStart_GNC_codegen_SIL(6U); memset(&(H[0]), 0, (11U) * sizeof(double)); profileEnd_GNC_codegen_SIL(6U); 
altitude_ratio = (6.356766E+6) / ((6.356766E+6) - x[10]); 
altitude = altitude_ratio * x[10]; 
layer_idx_0 = 0; 
layer_idx_1 = (101325.0); 
layer_idx_2 = (288.15); 
layer_idx_3 = (0.0065); { 
if (altitude > (11000.0)) { { 
if (altitude < (20000.0)) { 
airdata_altitude_pressure = (-(3.5699790210323479) * (altitude_ratio * altitude_ratio)) * exp((-(9.81) * (altitude - (11000.0))) / (62191.094035)); 


} else { { 
if (altitude < (32000.0)) { 
layer_idx_0 = 20000; 
layer_idx_1 = (5474.9); 
layer_idx_2 = (216.65); 
layer_idx_3 = -(0.001); 
} else { 
layer_idx_0 = 32000; 
layer_idx_1 = (868.02); 
layer_idx_2 = (228.65); 
layer_idx_3 = -(0.0028); 
}  } 
airdata_altitude_pressure = (((-layer_idx_1 * (9.81)) / (layer_idx_2 * (287.0579))) * (altitude_ratio * altitude_ratio)) * pow((1.0) - (layer_idx_3 / layer_idx_2) * (altitude - (double)layer_idx_0), (9.81) / ((287.0579) * layer_idx_3) - (1.0)); 
#line 169
}  } 
} else { 
airdata_altitude_pressure = (((-layer_idx_1 * (9.81)) / (layer_idx_2 * (287.0579))) * (altitude_ratio * altitude_ratio)) * pow((1.0) - (layer_idx_3 / layer_idx_2) * (altitude - (double)layer_idx_0), (9.81) / ((287.0579) * layer_idx_3) - (1.0)); 




}  } 
H[10] = airdata_altitude_pressure; 
profileStart_GNC_codegen_SIL(7U); memset(&(b_H[0]), 0, (11U) * sizeof(double)); profileEnd_GNC_codegen_SIL(7U); 
c_H = (0.0); { 
for (i = 0; i < 11; i++) { 
double d; 
d = b_H[i]; { 
for (i1 = 0; i1 < 11; i1++) { 
d += H[i1] * P[i1 + 11 * i]; 
}  } 
b_H[i] = d; 
c_H += d * H[i]; 
}  } 
c_b = (1.0) / (c_H + (100.0)); 
profileStart_GNC_codegen_SIL(8U); memset(&(b_P[0]), 0, (11U) * sizeof(double)); profileEnd_GNC_codegen_SIL(8U); { 
for (i2 = 0; i2 < 11; i2++) { { 
for (i3 = 0; i3 < 11; i3++) { 
b_P[i3] += P[i3 + 11 * i2] * H[i2]; 
}  } 
}  } { 
for (i4 = 0; i4 < 11; i4++) { 
K[i4] = b_P[i4] * c_b; 
}  } 
profileStart_GNC_codegen_SIL(9U); memset(&(b_I[0]), 0, (121U) * sizeof(char)); profileEnd_GNC_codegen_SIL(9U); { 
for (k = 0; k < 11; k++) { 
b_I[k + 11 * k] = (1); 
}  } { 
for (i5 = 0; i5 < 11; i5++) { { 
for (i6 = 0; i6 < 11; i6++) { 
int E_tmp; 
E_tmp = i6 + 11 * i5; 
E[E_tmp] = (double)(b_I[E_tmp]) - K[i6] * H[i5]; 
}  } 
}  } 
profileStart_GNC_codegen_SIL(10U); memset(&(b_E[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(10U); { 
for (i7 = 0; i7 < 11; i7++) { { 
for (i8 = 0; i8 < 11; i8++) { 
double d1; 
d1 = P[i8 + 11 * i7]; { 
for (i10 = 0; i10 < 11; i10++) { 
int b_E_tmp; 
b_E_tmp = i10 + 11 * i7; 
b_E[b_E_tmp] += E[i10 + 11 * i8] * d1; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(11U); memset(&(c_E[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(11U); { 
for (i9 = 0; i9 < 11; i9++) { { 
for (i11 = 0; i11 < 11; i11++) { 
double d2; 
d2 = E[i9 + 11 * i11]; { 
for (i13 = 0; i13 < 11; i13++) { 
int c_E_tmp; 
c_E_tmp = i13 + 11 * i9; 
c_E[c_E_tmp] += b_E[i13 + 11 * i11] * d2; 
}  } 
b_K[i11 + 11 * i9] = (K[i11] * (100.0)) * K[i9]; 
}  } 
}  } { 
for (i12 = 0; i12 < 121; i12++) { 
P_new[i12] = c_E[i12] + b_K[i12]; 
}  } { 
for (b_i = 0; b_i < 11; b_i++) { 
x_new[b_i] = x[b_i] + K[b_i] * b_b; 
}  } 
scale = (3.3121686421112381E-170); 
profileStart_GNC_codegen_SIL(12U); absxk = fabs(x_new[0]); profileEnd_GNC_codegen_SIL(12U); { 
if (absxk > (3.3121686421112381E-170)) { 
q_mag = (1.0); 
scale = absxk; 
} else { 
t = absxk / (3.3121686421112381E-170); 
q_mag = t * t; 
}  } 
profileStart_GNC_codegen_SIL(13U); absxk = fabs(x_new[1]); profileEnd_GNC_codegen_SIL(13U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(14U); absxk = fabs(x_new[2]); profileEnd_GNC_codegen_SIL(14U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(15U); absxk = fabs(x_new[3]); profileEnd_GNC_codegen_SIL(15U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(16U); q_mag = scale * sqrt(q_mag); profileEnd_GNC_codegen_SIL(16U); 
x_new[0] /= q_mag; 
x_new[1] /= q_mag; 
x_new[2] /= q_mag; 
x_new[3] /= q_mag; 
} 

static double b_norm(const double x[3]) { 
double absxk; 
double scale; 
double t; 
double y; 
scale = (3.3121686421112381E-170); 
profileStart_GNC_codegen_SIL(17U); absxk = fabs(x[0]); profileEnd_GNC_codegen_SIL(17U); { 
if (absxk > (3.3121686421112381E-170)) { 
y = (1.0); 
scale = absxk; 
} else { 
t = absxk / (3.3121686421112381E-170); 
y = t * t; 
}  } 
profileStart_GNC_codegen_SIL(18U); absxk = fabs(x[1]); profileEnd_GNC_codegen_SIL(18U); { 
if (absxk > scale) { 
t = scale / absxk; 
y = (y * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
y += t * t; 
}  } 
profileStart_GNC_codegen_SIL(19U); absxk = fabs(x[2]); profileEnd_GNC_codegen_SIL(19U); { 
if (absxk > scale) { 
t = scale / absxk; 
y = (y * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
y += t * t; 
}  } 
return scale * sqrt(y); 
} 

static void controller_codegen_entry_init(GNC_codegen_SILStackData *b_SD) { 
int i; 
(b_SD->pd->param).Cn_alpha = (10.0); { 
for (i = 0; i < 9; i++) { 
(b_SD->pd->param).J[i] = dv[i]; 
(b_SD->pd->param).Jinv[i] = dv1[i]; 
}  } 
(b_SD->pd->param).c_aero = -(0.016182736457722724); 
(b_SD->pd->param).c_canard = (0.00061367415999999994); 
(b_SD->pd->param).elevation = (420.0); 
(b_SD->pd->param).g[0] = -(9.81); 
(b_SD->pd->param).g[1] = (0.0); 
(b_SD->pd->param).g[2] = (0.0); 
} 

static void dynamics_init(GNC_codegen_SILStackData *b_SD) { 
int i; 
(b_SD->pd->c_param).Cn_alpha = (10.0); { 
for (i = 0; i < 9; i++) { 
(b_SD->pd->c_param).J[i] = dv[i]; 
(b_SD->pd->c_param).Jinv[i] = dv1[i]; 
}  } 
(b_SD->pd->c_param).c_aero = -(0.016182736457722724); 
(b_SD->pd->c_param).c_canard = (0.00061367415999999994); 
(b_SD->pd->c_param).elevation = (420.0); 
(b_SD->pd->c_param).g[0] = -(9.81); 
(b_SD->pd->c_param).g[1] = (0.0); 
(b_SD->pd->c_param).g[2] = (0.0); 
} 

static void dynamics_jacobian_init(GNC_codegen_SILStackData *b_SD) { 
int i; 
(b_SD->pd->d_param).Cn_alpha = (10.0); { 
for (i = 0; i < 9; i++) { 
(b_SD->pd->d_param).J[i] = dv[i]; 
(b_SD->pd->d_param).Jinv[i] = dv1[i]; 
}  } 
(b_SD->pd->d_param).c_aero = -(0.016182736457722724); 
(b_SD->pd->d_param).c_canard = (0.00061367415999999994); 
(b_SD->pd->d_param).elevation = (420.0); 
(b_SD->pd->d_param).g[0] = -(9.81); 
(b_SD->pd->d_param).g[1] = (0.0); 
(b_SD->pd->d_param).g[2] = (0.0); 
} 

static void ekf_correct(const double x[11], const double P[121], const double 
y[3], const double b[3], const double R[9], double 
x_new[11], double P_new[121]) { 
static const char iv[9] = {(1), (0), (0), (0), (1), (0), (0), (0), (1)}; 
double E[121]; 
double b_E[121]; 
double c_E[121]; 
double c_K[121]; 
double H[33]; 
double K[33]; 
double b_H[33]; 
double b_K[33]; 
double b_P[33]; 
double y_tmp[33]; 
double b_dv[9]; 
double c_H[9]; 
double c_a[9]; 
double b_y[3]; 
double c_x[3]; 
double a; 
double absxk; 
double b_a; 
double b_x; 
double d12; 
double d13; 
double d14; 
double d15; 
double d16; 
double d17; 
double q_mag; 
double scale; 
double t; 
int i; 
int i1; 
int i10; 
int i11; 
int i12; 
int i13; 
int i14; 
int i15; 
int i16; 
int i17; 
int i18; 
int i19; 
int i2; 
int i20; 
int i21; 
int i22; 
int i23; 
int i24; 
int i25; 
int i26; 
int i27; 
int i28; 
int i29; 
int i3; 
int i30; 
int i4; 
int i5; 
int i6; 
int i7; 
int i8; 
int i9; 
int k; 
char b_I[121]; 
a = x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]); 
b_a = (2.0) * x[0]; 
profileStart_GNC_codegen_SIL(20U); memset(&(H[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(20U); 
b_x = (b[0] * x[1] + b[1] * x[2]) + b[2] * x[3]; 
c_x[0] = x[2] * b[2] - b[1] * x[3]; 
c_x[1] = b[0] * x[3] - x[1] * b[2]; 
c_x[2] = x[1] * b[1] - b[0] * x[2]; 
b_dv[0] = (0.0); 
b_dv[3] = x[0] * -b[2]; 
b_dv[6] = x[0] * b[1]; 
b_dv[1] = x[0] * b[2]; 
b_dv[4] = (0.0); 
b_dv[7] = x[0] * -b[0]; 
b_dv[2] = x[0] * -b[1]; 
b_dv[5] = x[0] * b[0]; 
b_dv[8] = (0.0); { 
for (i = 0; i < 3; i++) { 
double H_tmp; 
int b_H_tmp; 
int c_H_tmp; 
int d_H_tmp; 
H[i] = (2.0) * (x[0] * b[i] - c_x[i]); 
H_tmp = x[i + 1]; 
b_H_tmp = 3 * (i + 1); 
H[b_H_tmp] = (2.0) * (((b_x * (double)(iv[3 * i]) + x[1] * b[i]) - b[0] * H_tmp) + b_dv[3 * i]); 


c_H_tmp = 3 * i + 1; 
H[b_H_tmp + 1] = (2.0) * (((b_x * (double)(iv[c_H_tmp]) + x[2] * b[i]) - b[1] * H_tmp) + b_dv[c_H_tmp]); 


d_H_tmp = 3 * i + 2; 
H[b_H_tmp + 2] = (2.0) * (((b_x * (double)(iv[d_H_tmp]) + x[3] * b[i]) - b[2] * H_tmp) + b_dv[d_H_tmp]); 


}  } { 
for (i1 = 0; i1 < 3; i1++) { { 
for (i2 = 0; i2 < 11; i2++) { 
y_tmp[i2 + 11 * i1] = H[i1 + 3 * i2]; 
}  } 
}  } 
profileStart_GNC_codegen_SIL(21U); memset(&(b_H[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(21U); { 
for (i3 = 0; i3 < 11; i3++) { 
double d; 
int e_H_tmp; 
int f_H_tmp; 
d = b_H[3 * i3]; 
e_H_tmp = 3 * i3 + 1; 
f_H_tmp = 3 * i3 + 2; { 
for (i5 = 0; i5 < 11; i5++) { 
double d1; 
d1 = P[i5 + 11 * i3]; 
d += H[3 * i5] * d1; 
b_H[e_H_tmp] += H[3 * i5 + 1] * d1; 
b_H[f_H_tmp] += H[3 * i5 + 2] * d1; 
}  } 
b_H[3 * i3] = d; 
}  } 
profileStart_GNC_codegen_SIL(22U); memset(&(b_P[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(22U); { 
for (i4 = 0; i4 < 3; i4++) { { 
for (i6 = 0; i6 < 3; i6++) { 
double d2; 
d2 = (0.0); { 
for (i8 = 0; i8 < 11; i8++) { 
d2 += b_H[i4 + 3 * i8] * y_tmp[i8 + 11 * i6]; 
}  } 
int g_H_tmp; 
g_H_tmp = i4 + 3 * i6; 
c_H[g_H_tmp] = d2 + R[g_H_tmp]; 
}  } { 
for (i7 = 0; i7 < 11; i7++) { 
double d3; 
d3 = y_tmp[i7 + 11 * i4]; { 
for (i10 = 0; i10 < 11; i10++) { 
int P_tmp; 
P_tmp = i10 + 11 * i4; 
b_P[P_tmp] += P[i10 + 11 * i7] * d3; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(23U); inv(c_H, b_dv); profileEnd_GNC_codegen_SIL(23U); 
profileStart_GNC_codegen_SIL(24U); memset(&(K[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(24U); { 
for (i9 = 0; i9 < 3; i9++) { { 
for (i11 = 0; i11 < 3; i11++) { 
double d4; 
d4 = b_dv[i11 + 3 * i9]; { 
for (i12 = 0; i12 < 11; i12++) { 
int K_tmp; 
K_tmp = i12 + 11 * i9; 
K[K_tmp] += b_P[i12 + 11 * i11] * d4; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(25U); memset(&(b_I[0]), 0, (121U) * sizeof(char)); profileEnd_GNC_codegen_SIL(25U); { 
for (k = 0; k < 11; k++) { 
b_I[k + 11 * k] = (1); 
}  } { 
for (i13 = 0; i13 < 11; i13++) { 
double d5; 
double d6; 
double d7; 
d5 = K[i13]; 
d6 = K[i13 + 11]; 
d7 = K[i13 + 22]; { 
for (i15 = 0; i15 < 11; i15++) { 
int E_tmp; 
E_tmp = i13 + 11 * i15; 
E[E_tmp] = (double)(b_I[E_tmp]) - ((d5 * H[3 * i15] + d6 * H[3 * i15 + 1]) + d7 * H[3 * i15 + 2]); 

}  } 
}  } 
profileStart_GNC_codegen_SIL(26U); memset(&(b_E[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(26U); { 
for (i14 = 0; i14 < 11; i14++) { { 
for (i16 = 0; i16 < 11; i16++) { 
double d8; 
d8 = P[i16 + 11 * i14]; { 
for (i18 = 0; i18 < 11; i18++) { 
int b_E_tmp; 
b_E_tmp = i18 + 11 * i14; 
b_E[b_E_tmp] += E[i18 + 11 * i16] * d8; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(27U); memset(&(b_K[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(27U); { 
for (i17 = 0; i17 < 3; i17++) { { 
for (i19 = 0; i19 < 3; i19++) { 
double d9; 
d9 = R[i19 + 3 * i17]; { 
for (i20 = 0; i20 < 11; i20++) { 
int b_K_tmp; 
b_K_tmp = i20 + 11 * i17; 
b_K[b_K_tmp] += K[i20 + 11 * i19] * d9; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(28U); memset(&(c_E[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(28U); 
profileStart_GNC_codegen_SIL(29U); memset(&(c_K[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(29U); { 
for (i21 = 0; i21 < 11; i21++) { { 
for (i22 = 0; i22 < 11; i22++) { 
double d10; 
d10 = E[i21 + 11 * i22]; { 
for (i25 = 0; i25 < 11; i25++) { 
int c_E_tmp; 
c_E_tmp = i25 + 11 * i21; 
c_E[c_E_tmp] += b_E[i25 + 11 * i22] * d10; 
}  } 
}  } { 
for (i24 = 0; i24 < 3; i24++) { 
double d11; 
d11 = K[i21 + 11 * i24]; { 
for (i27 = 0; i27 < 11; i27++) { 
int c_K_tmp; 
c_K_tmp = i27 + 11 * i21; 
c_K[c_K_tmp] += b_K[i27 + 11 * i24] * d11; 
}  } 
}  } 
}  } { 
for (i23 = 0; i23 < 121; i23++) { 
P_new[i23] = c_E[i23] + c_K[i23]; 
}  } { 
for (i26 = 0; i26 < 3; i26++) { 
double a_tmp; 
int b_a_tmp; 
int c_a_tmp; 
a_tmp = x[i26 + 1]; 
c_a[3 * i26] = a * (double)(iv[3 * i26]) + ((2.0) * x[1]) * a_tmp; 
b_a_tmp = 3 * i26 + 1; 
c_a[b_a_tmp] = a * (double)(iv[b_a_tmp]) + ((2.0) * x[2]) * a_tmp; 
c_a_tmp = 3 * i26 + 2; 
c_a[c_a_tmp] = a * (double)(iv[c_a_tmp]) + ((2.0) * x[3]) * a_tmp; 
}  } 
b_dv[0] = (0.0); 
b_dv[3] = b_a * -x[3]; 
b_dv[6] = b_a * x[2]; 
b_dv[1] = b_a * x[3]; 
b_dv[4] = (0.0); 
b_dv[7] = b_a * -x[1]; 
b_dv[2] = b_a * -x[2]; 
b_dv[5] = b_a * x[1]; 
b_dv[8] = (0.0); { 
for (i28 = 0; i28 < 9; i28++) { 
c_a[i28] -= b_dv[i28]; 
}  } 
d12 = b[0]; 
d13 = b[1]; 
d14 = b[2]; { 
for (i29 = 0; i29 < 3; i29++) { 
b_y[i29] = y[i29] - ((c_a[i29] * d12 + c_a[i29 + 3] * d13) + c_a[i29 + 6] * d14); 

}  } 
d15 = b_y[0]; 
d16 = b_y[1]; 
d17 = b_y[2]; { 
for (i30 = 0; i30 < 11; i30++) { 
x_new[i30] = x[i30] + ((K[i30] * d15 + K[i30 + 11] * d16) + K[i30 + 22] * d17); 

}  } 
scale = (3.3121686421112381E-170); 
profileStart_GNC_codegen_SIL(30U); absxk = fabs(x_new[0]); profileEnd_GNC_codegen_SIL(30U); { 
if (absxk > (3.3121686421112381E-170)) { 
q_mag = (1.0); 
scale = absxk; 
} else { 
t = absxk / (3.3121686421112381E-170); 
q_mag = t * t; 
}  } 
profileStart_GNC_codegen_SIL(31U); absxk = fabs(x_new[1]); profileEnd_GNC_codegen_SIL(31U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(32U); absxk = fabs(x_new[2]); profileEnd_GNC_codegen_SIL(32U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(33U); absxk = fabs(x_new[3]); profileEnd_GNC_codegen_SIL(33U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(34U); q_mag = scale * sqrt(q_mag); profileEnd_GNC_codegen_SIL(34U); 
x_new[0] /= q_mag; 
x_new[1] /= q_mag; 
x_new[2] /= q_mag; 
x_new[3] /= q_mag; 
} 

static void inv(const double x[9], double y[9]) { 
double b_x[9]; 
double absx11; 
double absx21; 
double absx31; 
double t2; 
double t3; 
int p1; 
int p2; 
int p3; 
memcpy(&(b_x[0]), &(x[0]), (9U) * sizeof(double)); 
p1 = 0; 
p2 = 3; 
p3 = 6; 
profileStart_GNC_codegen_SIL(35U); absx11 = fabs(x[0]); profileEnd_GNC_codegen_SIL(35U); 
profileStart_GNC_codegen_SIL(36U); absx21 = fabs(x[1]); profileEnd_GNC_codegen_SIL(36U); 
profileStart_GNC_codegen_SIL(37U); absx31 = fabs(x[2]); profileEnd_GNC_codegen_SIL(37U); { 
if ((absx21 > absx11) && (absx21 > absx31)) { 
p1 = 3; 
p2 = 0; 
b_x[0] = x[1]; 
b_x[1] = x[0]; 
b_x[3] = x[4]; 
b_x[4] = x[3]; 
b_x[6] = x[7]; 
b_x[7] = x[6]; 
} else { { if (absx31 > absx11) { 
p1 = 6; 
p3 = 0; 
b_x[0] = x[2]; 
b_x[2] = x[0]; 
b_x[3] = x[5]; 
b_x[5] = x[3]; 
b_x[6] = x[8]; 
b_x[8] = x[6]; 
}  } }  } 
b_x[1] /= b_x[0]; 
b_x[2] /= b_x[0]; 
b_x[4] -= b_x[1] * b_x[3]; 
b_x[5] -= b_x[2] * b_x[3]; 
b_x[7] -= b_x[1] * b_x[6]; 
b_x[8] -= b_x[2] * b_x[6]; { 
if (fabs(b_x[5]) > fabs(b_x[4])) { 
double t1; 
int itmp; 
itmp = p2; 
p2 = p3; 
p3 = itmp; 
t1 = b_x[1]; 
b_x[1] = b_x[2]; 
b_x[2] = t1; 
t1 = b_x[4]; 
b_x[4] = b_x[5]; 
b_x[5] = t1; 
t1 = b_x[7]; 
b_x[7] = b_x[8]; 
b_x[8] = t1; 
}  } 
b_x[5] /= b_x[4]; 
b_x[8] -= b_x[5] * b_x[7]; 
t3 = (b_x[1] * b_x[5] - b_x[2]) / b_x[8]; 
t2 = -(b_x[1] + b_x[7] * t3) / b_x[4]; 
y[p1] = (((1.0) - b_x[3] * t2) - b_x[6] * t3) / b_x[0]; 
y[p1 + 1] = t2; 
y[p1 + 2] = t3; 
t3 = -b_x[5] / b_x[8]; 
t2 = ((1.0) - b_x[7] * t3) / b_x[4]; 
y[p2] = -(b_x[3] * t2 + b_x[6] * t3) / b_x[0]; 
y[p2 + 1] = t2; 
y[p2 + 2] = t3; 
t3 = (1.0) / b_x[8]; 
t2 = (-b_x[7] * t3) / b_x[4]; 
y[p3] = -(b_x[3] * t2 + b_x[6] * t3) / b_x[0]; 
y[p3 + 1] = t2; 
y[p3 + 2] = t3; 
} 

static void pad_filter_init(GNC_codegen_SILStackData *b_SD) { 
int i; 
(b_SD->pd->b_param).Cn_alpha = (10.0); { 
for (i = 0; i < 9; i++) { 
(b_SD->pd->b_param).J[i] = dv[i]; 
(b_SD->pd->b_param).Jinv[i] = dv1[i]; 
}  } 
(b_SD->pd->b_param).c_aero = -(0.016182736457722724); 
(b_SD->pd->b_param).c_canard = (0.00061367415999999994); 
(b_SD->pd->b_param).elevation = (420.0); 
(b_SD->pd->b_param).g[0] = -(9.81); 
(b_SD->pd->b_param).g[1] = (0.0); 
(b_SD->pd->b_param).g[2] = (0.0); 
} 

void GNC_codegen_SIL_initialize(GNC_codegen_SILStackData *b_SD) { 
controller_codegen_entry_init(b_SD); 
pad_filter_init(b_SD); 
dynamics_init(b_SD); 
dynamics_jacobian_init(b_SD); 
} 

void GNC_codegen_SIL_terminate(void) { } 

void controller_codegen_entry(GNC_codegen_SILStackData *b_SD, double b_time, double 
dt_ctrl, const double xR[2], double pdyn, double 
delta, const struct0_T *ctrl_mem_in, double *
u, double *r, struct0_T *ctrl_mem_out) { 
double P[4]; 
double b_K[4]; 
double b_dv[4]; 
double b_dv1[4]; 
double dv2[4]; 
double K[2]; 
double b_r[2]; 
double L_delta; 
double a; 
double b; 
double blend; 
double c_delta; 
double c_r; 
double d; 
double d1; 
double d10; 
double d11; 
double d12; 
double d2; 
double d3; 
double d4; 
double d5; 
double d8; 
double pdyn_params; 
double r_idx_0; 
double w_dot; 
double x; 
double x_tmp; 
int i; 
int i2; { 
if ((b_time >= (22.0)) && (b_time < (27.0))) { 
*r = (0.5); 
} else { { if ((b_time >= (27.0)) && (b_time < (32.0))) { 
*r = -(0.5); 
} else { { if ((b_time >= (32.0)) && (b_time < (39.0))) { 
*r = (0.5); 
} else { 
*r = (0.0); 
}  } }  } }  } 
pdyn_params = pdyn * (b_SD->pd->param).c_canard; 
c_delta = delta / (2.0); { 
if (fabs(c_delta) < (0.005)) { 
c_delta = (0.0); 
}  } 
c_delta = (0.75) * ctrl_mem_in->d_old + (0.25) * c_delta; 
w_dot = (0.75) * ctrl_mem_in->w_dot_old + ((0.25) * (xR[1] - ctrl_mem_in->w_old)) / dt_ctrl; 

r_idx_0 = pdyn_params * c_delta; 
P[0] = ctrl_mem_in->P_minus[0] + (1.0E-5); 
P[1] = ctrl_mem_in->P_minus[1]; 
P[2] = ctrl_mem_in->P_minus[2]; 
P[3] = ctrl_mem_in->P_minus[3] + (1.0E-9); 
profileStart_GNC_codegen_SIL(38U); memset(&(b_r[0]), 0, sizeof(double) << 1); profileEnd_GNC_codegen_SIL(38U); 
d = r_idx_0 * (ctrl_mem_in->P_minus[0] + (1.0E-5)); 
d1 = pdyn_params * (ctrl_mem_in->P_minus[3] + (1.0E-9)); 
c_r = ((b_r[0] + d) + pdyn_params * ctrl_mem_in->P_minus[1]) * r_idx_0 + ((b_r[1] + r_idx_0 * ctrl_mem_in->P_minus[2]) + d1) * pdyn_params; 

K[0] = (d + ctrl_mem_in->P_minus[2] * pdyn_params) / (c_r + (1.0)); 
K[1] = (ctrl_mem_in->P_minus[1] * r_idx_0 + d1) / (c_r + (1.0)); 
b = w_dot - (r_idx_0 * ctrl_mem_in->coeffs[0] + pdyn_params * ctrl_mem_in->coeffs[1]); 

ctrl_mem_out->coeffs[0] = ctrl_mem_in->coeffs[0] + K[0] * b; 
ctrl_mem_out->coeffs[1] = ctrl_mem_in->coeffs[1] + K[1] * b; 
b_dv[0] = (1.0) - K[0] * r_idx_0; 
b_dv[1] = (0.0) - K[1] * r_idx_0; 
b_dv[2] = (0.0) - K[0] * pdyn_params; 
b_dv[3] = (1.0) - K[1] * pdyn_params; 
profileStart_GNC_codegen_SIL(39U); memset(&(b_dv1[0]), 0, sizeof(double) << 2); profileEnd_GNC_codegen_SIL(39U); 
d2 = b_dv[0]; 
d3 = b_dv[1]; 
d4 = b_dv[2]; 
d5 = b_dv[3]; { 
for (i = 0; i < 2; i++) { 
double d6; 
double d7; 
double d9; 
int i1; 
i1 = i << 1; 
d6 = P[i1]; 
d7 = b_dv1[i1] + d2 * d6; 
d9 = b_dv1[i1 + 1] + d3 * d6; 
d6 = P[i1 + 1]; 
d7 += d4 * d6; 
b_dv1[i1] = d7; 
d9 += d5 * d6; 
b_dv1[i1 + 1] = d9; 
}  } 
profileStart_GNC_codegen_SIL(40U); memset(&(dv2[0]), 0, sizeof(double) << 2); profileEnd_GNC_codegen_SIL(40U); 
d8 = b_dv1[0]; 
d10 = b_dv1[1]; 
d11 = b_dv1[2]; 
d12 = b_dv1[3]; { 
for (i2 = 0; i2 < 2; i2++) { 
double d13; 
double d14; 
double d15; 
int i3; 
d13 = b_dv[i2]; 
i3 = i2 << 1; 
d14 = dv2[i3] + d8 * d13; 
d15 = dv2[i3 + 1] + d10 * d13; 
b_K[i3] = K[0] * K[i2]; 
d13 = b_dv[i2 + 2]; 
d14 += d11 * d13; 
dv2[i3] = d14; 
d15 += d12 * d13; 
dv2[i3 + 1] = d15; 
b_K[i3 + 1] = K[1] * K[i2]; 
}  } 
ctrl_mem_out->P_minus[0] = dv2[0] + b_K[0]; 
ctrl_mem_out->P_minus[1] = dv2[1] + b_K[1]; 
ctrl_mem_out->P_minus[2] = dv2[2] + b_K[2]; 
ctrl_mem_out->P_minus[3] = dv2[3] + b_K[3]; 
ctrl_mem_out->w_old = xR[1]; 
ctrl_mem_out->d_old = c_delta; 
ctrl_mem_out->w_dot_old = w_dot; 
L_delta = (ctrl_mem_out->coeffs[0] * pdyn_params) / (2.0); { 
if (fabs(L_delta) < (10.0)) { { 
if (L_delta >= (0.0)) { 
L_delta = (10.0); 
} else { 
L_delta = -(10.0); 
}  } 
}  } 
profileStart_GNC_codegen_SIL(41U); blend = fmax((0.0), fmin((1.0), (fabs(xR[1]) - (0.5)) / (0.5))); profileEnd_GNC_codegen_SIL(41U); 
a = -(1.0) / L_delta; 
x_tmp = ((1.0) - blend) * (5.0); 
profileStart_GNC_codegen_SIL(42U); x = sqrt(x_tmp); profileEnd_GNC_codegen_SIL(42U); 
K[0] = a * x; 
profileStart_GNC_codegen_SIL(43U); *u = fmin(fmax((K[0] * xR[0] + (a * sqrt((2.0) * x + (x_tmp + blend * (20.0)))) * xR[1]) + -K[0] * *r, -(0.3490658503988659)), (0.3490658503988659)); profileEnd_GNC_codegen_SIL(43U); { 




if (pdyn < (500.0)) { 
*u = (0.0); 
}  } 
} 

void navigation_codegen_entry(GNC_codegen_SILStackData *b_SD, double dt, boolean_T 
flight_phase, const double x[11], const double 
P[121], const struct1_T *bias, const struct2_T *
sens_filt, const struct3_T *
sens_input, double x_ret[11], double 
P_ret[121], struct1_T *bias_ret, struct2_T *
sens_filt_ret, struct6_T *airdata, double 
roll_state[2]) { 
static const double Q[121] = {(1.0E-10), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (1.0E-10), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (1.0E-10), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (1.0E-10), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.01), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.01), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.01), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0001), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0001), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0001), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.0), (0.001)}; 
#line 933
static const double R[9] = {(1.0E-9), (0.0), (0.0), (0.0), (1.0E-9), (0.0), (0.0), (0.0), (1.0E-9)}; 

static const double b_b[9] = {(1.0), (0.0), (0.0), (0.0), (1.0), (0.0), (0.0), (0.0), (1.0)}; 
double F[121]; 
double b_E[121]; 
double b_F[121]; 
double b_P_ret[121]; 
double c_E[121]; 
double c_K[121]; 
double c_P_ret[121]; 
double d_P_ret[121]; 
double e_P_ret[121]; 
double f_P_ret[121]; 
double g_P_ret[121]; 
double h_P_ret[121]; 
double i_P_ret[121]; 
double K[33]; 
double b_K[33]; 
double b_W_dt[16]; 
double c_b[16]; 
double d_b[16]; 
double b_x_ret[11]; 
double c_x_ret[11]; 
double d_x_ret[11]; 
double e_x_ret[11]; 
double f_x_ret[11]; 
double g_x_ret[11]; 
double h_x_ret[11]; 
double i_x_ret[11]; 
double b_dv[9]; 
double b_n_tilde[9]; 
double b_w_exp_tilde[9]; 
double c_skewed_exp_w_tmp[9]; 
double c_q[4]; 
double q[4]; 
double b_dt[3]; 
double c_dt[3]; 
double c_w_exp_tilde[3]; 
double dv3[3]; 
double airspeed; 
double b_expl_temp; 
double c_expl_temp; 
double d_expl_temp; 
double e_expl_temp; 
double expl_temp; 
double f_expl_temp; 
double g_expl_temp; 
double h_expl_temp; 
double i_expl_temp; 
double j_expl_temp; 
double k_expl_temp; 
double l_expl_temp; 
double t1_density; 
int b_k; 
int c_k; 
int i; 
int i1; 
int i10; 
int i11; 
int i12; 
int i15; 
int i17; 
int i18; 
int i19; 
int i2; 
int i20; 
int i21; 
int i22; 
int i23; 
int i24; 
int i25; 
int i28; 
int i29; 
int i3; 
int i30; 
int i31; 
int i32; 
int i33; 
int i34; 
int i36; 
int i37; 
int i38; 
int i39; 
int i4; 
int i40; 
int i41; 
int i42; 
int i43; 
int i44; 
int i45; 
int i46; 
int i47; 
int i48; 
int i49; 
int i5; 
int i50; 
int i51; 
int i52; 
int i53; 
int i54; 
int i55; 
int i56; 
int i57; 
int i58; 
int i59; 
int i6; 
int i60; 
int i61; 
int i62; 
int i7; 
int i8; 
int i9; 
int k; 
char c_I[121]; 
memcpy(&(P_ret[0]), &(P[0]), (121U) * sizeof(double)); 
*bias_ret = *bias; 
*sens_filt_ret = *sens_filt; { 
if (!flight_phase) { 
double ST[9]; 
double h_a[9]; 
double a[3]; 
double b_absxk; 
double b_scale; 
double b_t; 
double board_baro_f; 
double d11; 
double d12; 
double d5; 
double d6; 
double d9; 
double d_a; 
double f_a; 
double mti_baro_f; 
double qw; 
double qy; 
double qz; 
double y; 
board_baro_f = sens_filt->board_baro_f; 
mti_baro_f = sens_filt->mti_baro_f; { 
if ((sens_input->board_accel).status) { 
sens_filt_ret->board_accel_f[0] = (0.0005) * (sens_input->board_accel).meas[0] + (0.9995) * sens_filt->board_accel_f[0]; 


sens_filt_ret->board_accel_f[1] = (0.0005) * (sens_input->board_accel).meas[1] + (0.9995) * sens_filt->board_accel_f[1]; 


sens_filt_ret->board_accel_f[2] = (0.0005) * (sens_input->board_accel).meas[2] + (0.9995) * sens_filt->board_accel_f[2]; 


}  } { 
if ((sens_input->board_gyro).status) { 
sens_filt_ret->board_gyro_f[0] = (0.0005) * (sens_input->board_gyro).meas[0] + (0.9995) * sens_filt->board_gyro_f[0]; 

sens_filt_ret->board_gyro_f[1] = (0.0005) * (sens_input->board_gyro).meas[1] + (0.9995) * sens_filt->board_gyro_f[1]; 

sens_filt_ret->board_gyro_f[2] = (0.0005) * (sens_input->board_gyro).meas[2] + (0.9995) * sens_filt->board_gyro_f[2]; 

}  } { 
if ((sens_input->mti_accel).status) { 
sens_filt_ret->mti_accel_f[0] = (0.0005) * (sens_input->mti_accel).meas[0] + (0.9995) * sens_filt->mti_accel_f[0]; 

sens_filt_ret->mti_accel_f[1] = (0.0005) * (sens_input->mti_accel).meas[1] + (0.9995) * sens_filt->mti_accel_f[1]; 

sens_filt_ret->mti_accel_f[2] = (0.0005) * (sens_input->mti_accel).meas[2] + (0.9995) * sens_filt->mti_accel_f[2]; 

}  } { 
if ((sens_input->mti_gyro).status) { 
sens_filt_ret->mti_gyro_f[0] = (0.0005) * (sens_input->mti_gyro).meas[0] + (0.9995) * sens_filt->mti_gyro_f[0]; 

sens_filt_ret->mti_gyro_f[1] = (0.0005) * (sens_input->mti_gyro).meas[1] + (0.9995) * sens_filt->mti_gyro_f[1]; 

sens_filt_ret->mti_gyro_f[2] = (0.0005) * (sens_input->mti_gyro).meas[2] + (0.9995) * sens_filt->mti_gyro_f[2]; 

}  } { 
if ((sens_input->ad_accel).status) { 
sens_filt_ret->ad_accel_f[0] = (0.0005) * (sens_input->ad_accel).meas[0] + (0.9995) * sens_filt->ad_accel_f[0]; 

sens_filt_ret->ad_accel_f[1] = (0.0005) * (sens_input->ad_accel).meas[1] + (0.9995) * sens_filt->ad_accel_f[1]; 

sens_filt_ret->ad_accel_f[2] = (0.0005) * (sens_input->ad_accel).meas[2] + (0.9995) * sens_filt->ad_accel_f[2]; 

}  } { 
if ((sens_input->ad_gyro).status) { 
sens_filt_ret->ad_gyro_f[0] = (0.0005) * (sens_input->ad_gyro).meas[0] + (0.9995) * sens_filt->ad_gyro_f[0]; 

sens_filt_ret->ad_gyro_f[1] = (0.0005) * (sens_input->ad_gyro).meas[1] + (0.9995) * sens_filt->ad_gyro_f[1]; 

sens_filt_ret->ad_gyro_f[2] = (0.0005) * (sens_input->ad_gyro).meas[2] + (0.9995) * sens_filt->ad_gyro_f[2]; 

}  } { 
if ((sens_input->board_baro).status) { 
board_baro_f = (0.0005) * (sens_input->board_baro).meas + (0.9995) * sens_filt->board_baro_f; 

}  } { 
if ((sens_input->board_mag).status) { 
sens_filt_ret->board_mag_f[0] = (0.0005) * (sens_input->board_mag).meas[0] + (0.9995) * sens_filt->board_mag_f[0]; 

sens_filt_ret->board_mag_f[1] = (0.0005) * (sens_input->board_mag).meas[1] + (0.9995) * sens_filt->board_mag_f[1]; 

sens_filt_ret->board_mag_f[2] = (0.0005) * (sens_input->board_mag).meas[2] + (0.9995) * sens_filt->board_mag_f[2]; 

}  } { 
if ((sens_input->mti_baro).status) { 
mti_baro_f = (0.0005) * (sens_input->mti_baro).meas + (0.9995) * sens_filt->mti_baro_f; 

}  } { 
if ((sens_input->mti_mag).status) { 
sens_filt_ret->mti_mag_f[0] = (0.0005) * (sens_input->mti_mag).meas[0] + (0.9995) * sens_filt->mti_mag_f[0]; 

sens_filt_ret->mti_mag_f[1] = (0.0005) * (sens_input->mti_mag).meas[1] + (0.9995) * sens_filt->mti_mag_f[1]; 

sens_filt_ret->mti_mag_f[2] = (0.0005) * (sens_input->mti_mag).meas[2] + (0.9995) * sens_filt->mti_mag_f[2]; 

}  } 
sens_filt_ret->board_baro_f = board_baro_f; 
sens_filt_ret->mti_baro_f = mti_baro_f; 
a[0] = (0.0); 
a[1] = (0.0); 
a[2] = (0.0); { 
if ((sens_input->board_accel).status) { 
a[0] = sens_filt_ret->board_accel_f[0]; 
a[1] = sens_filt_ret->board_accel_f[1]; 
a[2] = sens_filt_ret->board_accel_f[2]; 
}  } { 
if ((sens_input->mti_accel).status) { 
a[0] += sens_filt_ret->mti_accel_f[0]; 
a[1] += sens_filt_ret->mti_accel_f[1]; 
a[2] += sens_filt_ret->mti_accel_f[2]; 
}  } { 
if ((sens_input->ad_accel).status) { 
a[0] += sens_filt_ret->ad_accel_f[0]; 
a[1] += sens_filt_ret->ad_accel_f[1]; 
a[2] += sens_filt_ret->ad_accel_f[2]; 
}  } 
profileStart_GNC_codegen_SIL(44U); d5 = b_norm(a) + (1.0E-6); profileEnd_GNC_codegen_SIL(44U); 
profileStart_GNC_codegen_SIL(45U); qw = sqrt((0.5) * (a[0] / d5) + (0.5)); profileEnd_GNC_codegen_SIL(45U); { 
if (qw == (0.0)) { 
qy = (1.0); 
qz = (0.0); 
} else { 
qy = ((0.5) * (a[2] / d5)) / qw; 
qz = (-(0.5) * (a[1] / d5)) / qw; 
}  } 
b_scale = (3.3121686421112381E-170); { 
if (qw > (3.3121686421112381E-170)) { 
y = (1.0); 
b_scale = qw; 
} else { 
b_t = qw / (3.3121686421112381E-170); 
y = b_t * b_t; 
}  } 
profileStart_GNC_codegen_SIL(46U); b_absxk = fabs(qy); profileEnd_GNC_codegen_SIL(46U); { 
if (b_absxk > b_scale) { 
b_t = b_scale / b_absxk; 
y = (y * b_t) * b_t + (1.0); 
b_scale = b_absxk; 
} else { 
b_t = b_absxk / b_scale; 
y += b_t * b_t; 
}  } 
profileStart_GNC_codegen_SIL(47U); b_absxk = fabs(qz); profileEnd_GNC_codegen_SIL(47U); { 
if (b_absxk > b_scale) { 
b_t = b_scale / b_absxk; 
y = (y * b_t) * b_t + (1.0); 
b_scale = b_absxk; 
} else { 
b_t = b_absxk / b_scale; 
y += b_t * b_t; 
}  } 
profileStart_GNC_codegen_SIL(48U); y = b_scale * sqrt(y); profileEnd_GNC_codegen_SIL(48U); 
d6 = qw / y; 
q[0] = d6; 
x_ret[0] = d6; 
d6 = (0.0) / y; 
q[1] = d6; 
x_ret[1] = d6; 
d6 = qy / y; 
q[2] = d6; 
x_ret[2] = d6; 
d6 = qz / y; 
q[3] = d6; 
x_ret[3] = d6; 
x_ret[10] = (b_SD->pd->b_param).elevation; 
x_ret[4] = (0.0); 
x_ret[7] = (0.0); 
bias_ret->board_gyro[0] = sens_filt_ret->board_gyro_f[0]; 
bias_ret->mti_gyro[0] = sens_filt_ret->mti_gyro_f[0]; 
bias_ret->ad_gyro[0] = sens_filt_ret->ad_gyro_f[0]; 
x_ret[5] = (0.0); 
x_ret[8] = (0.0); 
bias_ret->board_gyro[1] = sens_filt_ret->board_gyro_f[1]; 
bias_ret->mti_gyro[1] = sens_filt_ret->mti_gyro_f[1]; 
bias_ret->ad_gyro[1] = sens_filt_ret->ad_gyro_f[1]; 
x_ret[6] = (0.0); 
x_ret[9] = (0.0); 
bias_ret->board_gyro[2] = sens_filt_ret->board_gyro_f[2]; 
bias_ret->mti_gyro[2] = sens_filt_ret->mti_gyro_f[2]; 
bias_ret->ad_gyro[2] = sens_filt_ret->ad_gyro_f[2]; 
d_a = q[0] * q[0] - ((q[1] * q[1] + q[2] * q[2]) + d6 * d6); 
f_a = (2.0) * q[0]; { 
for (i3 = 0; i3 < 3; i3++) { 
double a_tmp; 
a_tmp = (2.0) * q[i3 + 1]; 
h_a[3 * i3] = d_a * b_b[i3] + a_tmp * q[1]; 
h_a[3 * i3 + 1] = d_a * b_b[i3 + 3] + a_tmp * q[2]; 
h_a[3 * i3 + 2] = d_a * b_b[i3 + 6] + a_tmp * d6; 
}  } 
b_dv[0] = (0.0); 
b_dv[1] = f_a * -d6; 
b_dv[2] = f_a * q[2]; 
b_dv[3] = f_a * d6; 
b_dv[4] = (0.0); 
b_dv[5] = f_a * -q[1]; 
b_dv[6] = f_a * -q[2]; 
b_dv[7] = f_a * q[1]; 
b_dv[8] = (0.0); { 
for (i5 = 0; i5 < 9; i5++) { 
ST[i5] = h_a[i5] - b_dv[i5]; 
}  } 
bias_ret->board_mag_earth[0] = (0.0); 
bias_ret->board_mag_earth[1] = (0.0); 
bias_ret->board_mag_earth[2] = (0.0); { 
for (i7 = 0; i7 < 3; i7++) { 
double d10; 
d10 = sens_filt_ret->board_mag_f[i7]; 
bias_ret->board_mag_earth[0] += ST[3 * i7] * d10; 
bias_ret->board_mag_earth[1] += ST[3 * i7 + 1] * d10; 
bias_ret->board_mag_earth[2] += ST[3 * i7 + 2] * d10; 
bias_ret->mti_mag_earth[i7] = (0.0); 
}  } 
d9 = bias_ret->mti_mag_earth[0]; 
d11 = bias_ret->mti_mag_earth[1]; 
d12 = bias_ret->mti_mag_earth[2]; { 
for (i8 = 0; i8 < 3; i8++) { 
double d13; 
d13 = sens_filt_ret->mti_mag_f[i8]; 
d9 += ST[3 * i8] * d13; 
d11 += ST[3 * i8 + 1] * d13; 
d12 += ST[3 * i8 + 2] * d13; 
}  } 
double t1_pressure; 
bias_ret->mti_mag_earth[2] = d12; 
bias_ret->mti_mag_earth[1] = d11; 
bias_ret->mti_mag_earth[0] = d9; 
t1_pressure = airdata_atmos((b_SD->pd->b_param).elevation, &e_expl_temp, &t1_density, &f_expl_temp, &g_expl_temp, &h_expl_temp); 


bias_ret->board_baro = board_baro_f - t1_pressure; 
bias_ret->mti_baro = mti_baro_f - t1_pressure; 
} else { 
double E[121]; 
double P_pred[121]; 
double W_dt[16]; 
double b_q[16]; 
double l_a[16]; 
double d_dt[12]; 
double x_pred[11]; 
double S[9]; 
double b_P_pred[9]; 
double b_skewed_exp_w_tmp[9]; 
double dv4[9]; 
double h_a[9]; 
double n_tilde[9]; 
double skewed_exp_w_tmp[9]; 
double w_exp_tilde[9]; 
double b_dv1[4]; 
double r_q_tmp[4]; 
double C_total_a[3]; 
double b_S[3]; 
double c_r_q_tmp[3]; 
double d_x[3]; 
double dn[3]; 
double dv2[3]; 
double C_ad_w_idx_0; 
double C_total_a_tmp; 
double C_total_a_tmp_tmp; 
double absxk; 
double b; 
double b_C_total_a_tmp_tmp; 
double b_a; 
double b_dphi_tmp; 
double b_q_mag; 
double b_r_q_tmp; 
double b_x; 
double c_C_total_a_tmp_tmp; 
double c_absxk; 
double c_scale; 
double c_t; 
double c_x; 
double d; 
double d1; 
double d16; 
double d17; 
double d18; 
double d19; 
double d2; 
double d20; 
double d21; 
double d22; 
double d23; 
double d24; 
double d25; 
double d26; 
double d27; 
double d3; 
double d30; 
double d31; 
double d33; 
double d34; 
double d4; 
double d72; 
double d73; 
double d74; 
double dphi; 
double dphi_tmp; 
double e_a; 
double g_a; 
double i_a; 
double j_a; 
double k_a; 
double m_a; 
double n_a; 
double n_idx_0; 
double n_idx_1; 
double n_idx_2; 
double o_a; 
double q_mag; 
double scale; 
double t; 
char b_I[16]; 
char w_exp_tilde_tmp[9]; 
d = (9.9999999999999981E+9) * (double)((sens_input->ad_gyro).status); 
C_total_a_tmp_tmp = (1.0000000000000002E+14) * (double)((sens_input->board_accel).status); 

b_C_total_a_tmp_tmp = (1.0000000000000002E+14) * (double)((sens_input->mti_accel).status); 

c_C_total_a_tmp_tmp = (1.0000000000000002E+14) * (double)((sens_input->ad_accel).status); 

C_total_a_tmp = (C_total_a_tmp_tmp + b_C_total_a_tmp_tmp) + c_C_total_a_tmp_tmp; 

C_total_a[0] = C_total_a_tmp; 
d1 = (9.9999999999999981E+9) * (double)((sens_input->board_gyro).status); 
d2 = (9.9999999999999981E+9) * (double)((sens_input->mti_gyro).status); 
d3 = d1 + d2; 
d4 = d3 + d; 
d /= d4; 
C_ad_w_idx_0 = d; 
C_total_a[1] = C_total_a_tmp; 
d = (0.0) / d3; 
C_total_a[2] = C_total_a_tmp; 
scale = (3.3121686421112381E-170); 
profileStart_GNC_codegen_SIL(49U); absxk = fabs(x[0]); profileEnd_GNC_codegen_SIL(49U); { 
if (absxk > (3.3121686421112381E-170)) { 
q_mag = (1.0); 
scale = absxk; 
} else { 
t = absxk / (3.3121686421112381E-170); 
q_mag = t * t; 
}  } 
profileStart_GNC_codegen_SIL(50U); absxk = fabs(x[1]); profileEnd_GNC_codegen_SIL(50U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(51U); absxk = fabs(x[2]); profileEnd_GNC_codegen_SIL(51U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(52U); absxk = fabs(x[3]); profileEnd_GNC_codegen_SIL(52U); { 
if (absxk > scale) { 
t = scale / absxk; 
q_mag = (q_mag * t) * t + (1.0); 
scale = absxk; 
} else { 
t = absxk / scale; 
q_mag += t * t; 
}  } 
profileStart_GNC_codegen_SIL(53U); q_mag = scale * sqrt(q_mag); profileEnd_GNC_codegen_SIL(53U); 
q[0] = x[0] / q_mag; 
q[1] = x[1] / q_mag; 
q[2] = x[2] / q_mag; 
q[3] = x[3] / q_mag; 
profileStart_GNC_codegen_SIL(54U); dphi_tmp = b_norm(&(x[4])); profileEnd_GNC_codegen_SIL(54U); 
b_dphi_tmp = dphi_tmp * dt; 
dphi = b_dphi_tmp / (2.0); { 
if (dphi_tmp == (0.0)) { 
dn[0] = (0.0); 
dn[1] = (0.0); 
dn[2] = (0.0); 
n_idx_0 = (0.0); 
n_idx_1 = (0.0); 
n_idx_2 = (0.0); 
} else { 
double b_dn_tmp; 
double c_dn_tmp; 
double dn_tmp; 
dn_tmp = x[4] / dphi_tmp; 
dn[0] = dn_tmp; 
b_dn_tmp = x[5] / dphi_tmp; 
dn[1] = b_dn_tmp; 
c_dn_tmp = x[6] / dphi_tmp; 
dn[2] = c_dn_tmp; 
n_idx_0 = dn_tmp; 
n_idx_1 = b_dn_tmp; 
n_idx_2 = c_dn_tmp; 
}  } 
profileStart_GNC_codegen_SIL(55U); b = sin(dphi); profileEnd_GNC_codegen_SIL(55U); 
n_tilde[0] = (0.0); 
n_tilde[3] = -n_idx_2; 
n_tilde[6] = n_idx_1; 
n_tilde[1] = n_idx_2; 
n_tilde[4] = (0.0); 
n_tilde[7] = -n_idx_0; 
n_tilde[2] = -n_idx_1; 
n_tilde[5] = n_idx_0; 
n_tilde[8] = (0.0); 
profileStart_GNC_codegen_SIL(56U); b_a = sin(b_dphi_tmp); profileEnd_GNC_codegen_SIL(56U); 
profileStart_GNC_codegen_SIL(57U); b_x = cos(b_dphi_tmp); profileEnd_GNC_codegen_SIL(57U); { 
for (i = 0; i < 9; i++) { 
w_exp_tilde_tmp[i] = (0); 
}  } 
profileStart_GNC_codegen_SIL(58U); memset(&(b_n_tilde[0]), 0, (9U) * sizeof(double)); profileEnd_GNC_codegen_SIL(58U); { 
for (k = 0; k < 3; k++) { 
double d7; 
int b_n_tilde_tmp; 
int n_tilde_tmp; 
w_exp_tilde_tmp[k + 3 * k] = (1); 
d7 = b_n_tilde[3 * k]; 
n_tilde_tmp = 3 * k + 1; 
b_n_tilde_tmp = 3 * k + 2; { 
for (i2 = 0; i2 < 3; i2++) { 
double d8; 
d8 = n_tilde[i2 + 3 * k]; 
d7 += n_tilde[3 * i2] * d8; 
b_n_tilde[n_tilde_tmp] += n_tilde[3 * i2 + 1] * d8; 
b_n_tilde[b_n_tilde_tmp] += n_tilde[3 * i2 + 2] * d8; 
}  } 
b_n_tilde[3 * k] = d7; 
}  } { 
for (i1 = 0; i1 < 9; i1++) { 
w_exp_tilde[i1] = ((double)(w_exp_tilde_tmp[i1]) - b_a * n_tilde[i1]) + ((1.0) - b_x) * b_n_tilde[i1]; 

}  } 
double c_a; 
profileStart_GNC_codegen_SIL(59U); c_a = b_norm(&(x[7])); profileEnd_GNC_codegen_SIL(59U); 
profileStart_GNC_codegen_SIL(60U); airdata_atmos(x[10], &expl_temp, &t1_density, &b_expl_temp, &c_expl_temp, &d_expl_temp); profileEnd_GNC_codegen_SIL(60U); 

e_a = ((0.5) * t1_density) * (c_a * c_a); 
g_a = (b_SD->pd->c_param).c_aero * (b_SD->pd->c_param).Cn_alpha; 
i_a = x[0] * x[0] - ((x[1] * x[1] + x[2] * x[2]) + x[3] * x[3]); 
j_a = (2.0) * x[0]; { 
for (i4 = 0; i4 < 3; i4++) { 
double b_a_tmp; 
int c_a_tmp; 
int d_a_tmp; 
b_a_tmp = x[i4 + 1]; 
h_a[3 * i4] = i_a * b_b[3 * i4] + ((2.0) * x[1]) * b_a_tmp; 
c_a_tmp = 3 * i4 + 1; 
h_a[c_a_tmp] = i_a * b_b[c_a_tmp] + ((2.0) * x[2]) * b_a_tmp; 
d_a_tmp = 3 * i4 + 2; 
h_a[d_a_tmp] = i_a * b_b[d_a_tmp] + ((2.0) * x[3]) * b_a_tmp; 
}  } 
b_dv[0] = (0.0); 
b_dv[3] = j_a * -x[3]; 
b_dv[6] = j_a * x[2]; 
b_dv[1] = j_a * x[3]; 
b_dv[4] = (0.0); 
b_dv[7] = j_a * -x[1]; 
b_dv[2] = j_a * -x[2]; 
b_dv[5] = j_a * x[1]; 
b_dv[8] = (0.0); { 
for (i6 = 0; i6 < 9; i6++) { 
S[i6] = h_a[i6] - b_dv[i6]; 
}  } 
b_q[0] = q[0]; 
b_q[4] = -q[1]; 
b_q[8] = -q[2]; 
b_q[12] = -q[3]; 
b_q[1] = q[1]; 
b_q[5] = q[0]; 
b_q[9] = -q[3]; 
b_q[13] = q[2]; 
b_q[2] = q[2]; 
b_q[6] = q[3]; 
b_q[10] = q[0]; 
b_q[14] = -q[1]; 
b_q[3] = q[3]; 
b_q[7] = -q[2]; 
b_q[11] = q[1]; 
b_q[15] = q[0]; 
profileStart_GNC_codegen_SIL(61U); b_dv1[0] = cos(dphi); profileEnd_GNC_codegen_SIL(61U); 
profileStart_GNC_codegen_SIL(62U); memset(&(b_w_exp_tilde[0]), 0, (9U) * sizeof(double)); profileEnd_GNC_codegen_SIL(62U); 
profileStart_GNC_codegen_SIL(63U); memset(&(c_w_exp_tilde[0]), 0, (3U) * sizeof(double)); profileEnd_GNC_codegen_SIL(63U); { 
for (i9 = 0; i9 < 3; i9++) { 
double d14; 
int b_w_exp_tilde_tmp; 
int c_w_exp_tilde_tmp; 
b_dv1[i9 + 1] = dn[i9] * b; 
d14 = b_w_exp_tilde[3 * i9]; 
b_w_exp_tilde_tmp = 3 * i9 + 1; 
c_w_exp_tilde_tmp = 3 * i9 + 2; { 
for (i10 = 0; i10 < 3; i10++) { 
double d15; 
d15 = (b_SD->pd->c_param).J[i10 + 3 * i9]; 
d14 += w_exp_tilde[3 * i10] * d15; 
b_w_exp_tilde[b_w_exp_tilde_tmp] += w_exp_tilde[3 * i10 + 1] * d15; 
b_w_exp_tilde[c_w_exp_tilde_tmp] += w_exp_tilde[3 * i10 + 2] * d15; 
}  } 
double d_w_exp_tilde_tmp; 
b_w_exp_tilde[3 * i9] = d14; 
d_w_exp_tilde_tmp = x[i9 + 4]; 
c_w_exp_tilde[0] += d14 * d_w_exp_tilde_tmp; 
c_w_exp_tilde[1] += b_w_exp_tilde[3 * i9 + 1] * d_w_exp_tilde_tmp; 
c_w_exp_tilde[2] += b_w_exp_tilde[3 * i9 + 2] * d_w_exp_tilde_tmp; 
}  } 
dv2[0] = (0.0); 
profileStart_GNC_codegen_SIL(64U); dv2[1] = e_a * (g_a * sin(atan2(x[9], x[7]))); profileEnd_GNC_codegen_SIL(64U); 
profileStart_GNC_codegen_SIL(65U); dv2[2] = e_a * (g_a * -sin(atan2(x[8], x[7]))); profileEnd_GNC_codegen_SIL(65U); 
profileStart_GNC_codegen_SIL(66U); memset(&(dv3[0]), 0, (3U) * sizeof(double)); profileEnd_GNC_codegen_SIL(66U); 
profileStart_GNC_codegen_SIL(67U); memset(&(b_dt[0]), 0, (3U) * sizeof(double)); profileEnd_GNC_codegen_SIL(67U); 
profileStart_GNC_codegen_SIL(68U); memset(&(c_dt[0]), 0, (3U) * sizeof(double)); profileEnd_GNC_codegen_SIL(68U); 
d16 = dv3[0]; 
d17 = dv3[1]; 
d18 = dv3[2]; 
d19 = b_dt[0]; 
d20 = b_dt[1]; 
d21 = b_dt[2]; 
d22 = c_dt[0]; 
d23 = c_dt[1]; 
d24 = c_dt[2]; 
d25 = x[7]; 
d26 = x[8]; 
d27 = x[9]; { 
for (i11 = 0; i11 < 3; i11++) { 
double d28; 
double d29; 
double d32; 
double d35; 
double d36; 
double d37; 
double d39; 
int i13; 
int i14; 
d28 = (b_SD->pd->c_param).Jinv[3 * i11]; 
d29 = c_w_exp_tilde[i11]; 
d16 += d28 * d29; 
d32 = dv2[i11]; 
d19 += (dt * d28) * d32; 
d35 = S[3 * i11]; 
d36 = (b_SD->pd->c_param).g[i11]; 
d22 += (dt * d35) * d36; 
d37 = d35 * d25; 
i13 = 3 * i11 + 1; 
d28 = (b_SD->pd->c_param).Jinv[i13]; 
d17 += d28 * d29; 
d20 += (dt * d28) * d32; 
d35 = S[i13]; 
d23 += (dt * d35) * d36; 
d37 += d35 * d26; 
i14 = 3 * i11 + 2; 
d28 = (b_SD->pd->c_param).Jinv[i14]; 
d18 += d28 * d29; 
d21 += (dt * d28) * d32; 
d35 = S[i14]; 
d24 += (dt * d35) * d36; 
d37 += d35 * d27; 
d39 = C_total_a[i11]; 
c_w_exp_tilde[i11] = ((w_exp_tilde[i11] * d25 + w_exp_tilde[i11 + 3] * d26) + w_exp_tilde[i11 + 6] * d27) + dt * (((C_total_a_tmp_tmp / d39) * (sens_input->board_accel).meas[i11] + (b_C_total_a_tmp_tmp / d39) * (sens_input->mti_accel).meas[i11]) + (c_C_total_a_tmp_tmp / d39) * (sens_input->ad_accel).meas[i11]); 
#line 1615
b_S[i11] = d37; 
}  } 
profileStart_GNC_codegen_SIL(69U); memset(&(c_q[0]), 0, sizeof(double) << 2); profileEnd_GNC_codegen_SIL(69U); 
d30 = c_q[0]; 
d31 = c_q[1]; 
d33 = c_q[2]; 
d34 = c_q[3]; { 
for (i12 = 0; i12 < 4; i12++) { 
double d38; 
int q_tmp; 
q_tmp = i12 << 2; 
d38 = b_dv1[i12]; 
d30 += b_q[q_tmp] * d38; 
d31 += b_q[q_tmp + 1] * d38; 
d33 += b_q[q_tmp + 2] * d38; 
d34 += b_q[q_tmp + 3] * d38; 
}  } 
double W_dt_tmp; 
double b_W_dt_tmp; 
double c_W_dt_tmp; 
double d_W_dt_tmp; 
double e_W_dt_tmp; 
double f_W_dt_tmp; 
x_pred[0] = d30; 
x_pred[1] = d31; 
x_pred[2] = d33; 
x_pred[3] = d34; 
x_pred[4] = d16 + d19; 
x_pred[7] = c_w_exp_tilde[0] + d22; 
x_pred[5] = d17 + d20; 
x_pred[8] = c_w_exp_tilde[1] + d23; 
x_pred[6] = d18 + d21; 
x_pred[9] = c_w_exp_tilde[2] + d24; 
x_pred[10] = x[10] + dt * b_S[0]; 
profileStart_GNC_codegen_SIL(70U); memset(&(F[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(70U); 
k_a = (0.5) * dt; 
W_dt[0] = (0.0); 
W_dt_tmp = k_a * -x[4]; 
W_dt[4] = W_dt_tmp; 
b_W_dt_tmp = k_a * -x[5]; 
W_dt[8] = b_W_dt_tmp; 
c_W_dt_tmp = k_a * -x[6]; 
W_dt[12] = c_W_dt_tmp; 
d_W_dt_tmp = k_a * x[4]; 
W_dt[1] = d_W_dt_tmp; 
W_dt[5] = (0.0); 
e_W_dt_tmp = k_a * x[6]; 
W_dt[9] = e_W_dt_tmp; 
W_dt[13] = b_W_dt_tmp; 
f_W_dt_tmp = k_a * x[5]; 
W_dt[2] = f_W_dt_tmp; 
W_dt[6] = c_W_dt_tmp; 
W_dt[10] = (0.0); 
W_dt[14] = d_W_dt_tmp; 
W_dt[3] = e_W_dt_tmp; 
W_dt[7] = f_W_dt_tmp; 
W_dt[11] = W_dt_tmp; 
W_dt[15] = (0.0); 
profileStart_GNC_codegen_SIL(71U); memset(&(c_b[0]), 0, sizeof(double) << 4); profileEnd_GNC_codegen_SIL(71U); { 
for (i15 = 0; i15 < 4; i15++) { 
int i16; 
i16 = i15 << 2; { 
for (i18 = 0; i18 < 4; i18++) { 
double d40; 
int b_tmp; 
d40 = W_dt[i18 + i16]; 
b_tmp = i18 << 2; 
c_b[i16] += W_dt[b_tmp] * d40; 
c_b[i16 + 1] += W_dt[b_tmp + 1] * d40; 
c_b[i16 + 2] += W_dt[b_tmp + 2] * d40; 
c_b[i16 + 3] += W_dt[b_tmp + 3] * d40; 
}  } 
}  } { 
for (i17 = 0; i17 < 16; i17++) { 
b_I[i17] = (0); 
}  } 
profileStart_GNC_codegen_SIL(72U); memset(&(b_W_dt[0]), 0, sizeof(double) << 4); profileEnd_GNC_codegen_SIL(72U); 
profileStart_GNC_codegen_SIL(73U); memset(&(d_b[0]), 0, sizeof(double) << 4); profileEnd_GNC_codegen_SIL(73U); { 
for (b_k = 0; b_k < 4; b_k++) { 
double d41; 
double d42; 
double d43; 
double d44; 
double d45; 
double d46; 
double d47; 
double d48; 
int I_tmp; 
I_tmp = b_k << 2; 
b_I[b_k + I_tmp] = (1); 
d41 = b_W_dt[I_tmp]; 
d42 = b_W_dt[I_tmp + 1]; 
d43 = b_W_dt[I_tmp + 2]; 
d44 = b_W_dt[I_tmp + 3]; 
d45 = d_b[I_tmp]; 
d46 = d_b[I_tmp + 1]; 
d47 = d_b[I_tmp + 2]; 
d48 = d_b[I_tmp + 3]; { 
for (i19 = 0; i19 < 4; i19++) { 
double d49; 
int g_W_dt_tmp; 
d49 = c_b[i19 + I_tmp]; 
g_W_dt_tmp = i19 << 2; 
d41 += W_dt[g_W_dt_tmp] * d49; 
d45 += c_b[g_W_dt_tmp] * d49; 
d42 += W_dt[g_W_dt_tmp + 1] * d49; 
d46 += c_b[g_W_dt_tmp + 1] * d49; 
d43 += W_dt[g_W_dt_tmp + 2] * d49; 
d47 += c_b[g_W_dt_tmp + 2] * d49; 
d44 += W_dt[g_W_dt_tmp + 3] * d49; 
d48 += c_b[g_W_dt_tmp + 3] * d49; 
}  } 
d_b[I_tmp + 3] = d48; 
d_b[I_tmp + 2] = d47; 
d_b[I_tmp + 1] = d46; 
d_b[I_tmp] = d45; 
b_W_dt[I_tmp + 3] = d44; 
b_W_dt[I_tmp + 2] = d43; 
b_W_dt[I_tmp + 1] = d42; 
b_W_dt[I_tmp] = d41; 
F[11 * b_k] = ((((double)(b_I[I_tmp]) + W_dt[I_tmp]) + (0.5) * c_b[I_tmp]) + (0.16666666666666666) * d41) + (0.041666666666666664) * d45; 


F[11 * b_k + 1] = ((((double)(b_I[I_tmp + 1]) + W_dt[I_tmp + 1]) + (0.5) * c_b[I_tmp + 1]) + (0.16666666666666666) * d42) + (0.041666666666666664) * d46; 



F[11 * b_k + 2] = ((((double)(b_I[I_tmp + 2]) + W_dt[I_tmp + 2]) + (0.5) * c_b[I_tmp + 2]) + (0.16666666666666666) * d43) + (0.041666666666666664) * d47; 



F[11 * b_k + 3] = ((((double)(b_I[I_tmp + 3]) + W_dt[I_tmp + 3]) + (0.5) * c_b[I_tmp + 3]) + (0.16666666666666666) * d44) + (0.041666666666666664) * d48; 



}  } 
double e_a_tmp; 
double f_a_tmp; 
double g_a_tmp; 
double h_a_tmp; 
double i_a_tmp; 
double j_a_tmp; 
double k_a_tmp; 
e_a_tmp = k_a * q[0]; 
l_a[0] = e_a_tmp; 
f_a_tmp = k_a * -q[1]; 
l_a[4] = f_a_tmp; 
g_a_tmp = k_a * -q[2]; 
l_a[8] = g_a_tmp; 
h_a_tmp = k_a * -q[3]; 
l_a[12] = h_a_tmp; 
i_a_tmp = k_a * q[1]; 
l_a[1] = i_a_tmp; 
l_a[5] = e_a_tmp; 
l_a[9] = h_a_tmp; 
j_a_tmp = k_a * q[2]; 
l_a[13] = j_a_tmp; 
l_a[2] = j_a_tmp; 
k_a_tmp = k_a * q[3]; 
l_a[6] = k_a_tmp; 
l_a[10] = e_a_tmp; 
l_a[14] = f_a_tmp; 
l_a[3] = k_a_tmp; 
l_a[7] = g_a_tmp; 
l_a[11] = i_a_tmp; 
l_a[15] = e_a_tmp; { 
for (i20 = 0; i20 < 3; i20++) { 
int F_tmp; 
int b_F_tmp; 
F_tmp = (i20 + 1) << 2; 
b_F_tmp = 11 * (i20 + 4); 
F[b_F_tmp] = l_a[F_tmp]; 
F[b_F_tmp + 1] = l_a[F_tmp + 1]; 
F[b_F_tmp + 2] = l_a[F_tmp + 2]; 
F[b_F_tmp + 3] = l_a[F_tmp + 3]; 
}  } 
m_a = ((0.5) * (b_SD->pd->d_param).c_aero) * (b_SD->pd->d_param).Cn_alpha; 
profileStart_GNC_codegen_SIL(74U); airdata_atmos(x[10], &i_expl_temp, &t1_density, &j_expl_temp, &k_expl_temp, &l_expl_temp); profileEnd_GNC_codegen_SIL(74U); { 

if (dphi_tmp == (0.0)) { 
n_idx_0 = (0.0); 
n_idx_1 = (0.0); 
n_idx_2 = (0.0); 
} else { 
n_idx_0 = x[4] / dphi_tmp; 
n_idx_1 = x[5] / dphi_tmp; 
n_idx_2 = x[6] / dphi_tmp; 
}  } 
n_tilde[0] = (0.0); 
n_tilde[3] = -n_idx_2; 
n_tilde[6] = n_idx_1; 
n_tilde[1] = n_idx_2; 
n_tilde[4] = (0.0); 
n_tilde[7] = -n_idx_0; 
n_tilde[2] = -n_idx_1; 
n_tilde[5] = n_idx_0; 
n_tilde[8] = (0.0); 
profileStart_GNC_codegen_SIL(75U); memset(&(b_n_tilde[0]), 0, (9U) * sizeof(double)); profileEnd_GNC_codegen_SIL(75U); { 
for (i21 = 0; i21 < 3; i21++) { 
double d50; 
int c_n_tilde_tmp; 
int d_n_tilde_tmp; 
d50 = b_n_tilde[3 * i21]; 
c_n_tilde_tmp = 3 * i21 + 1; 
d_n_tilde_tmp = 3 * i21 + 2; { 
for (i23 = 0; i23 < 3; i23++) { 
double d51; 
d51 = n_tilde[i23 + 3 * i21]; 
d50 += n_tilde[3 * i23] * d51; 
b_n_tilde[c_n_tilde_tmp] += n_tilde[3 * i23 + 1] * d51; 
b_n_tilde[d_n_tilde_tmp] += n_tilde[3 * i23 + 2] * d51; 
}  } 
b_n_tilde[3 * i21] = d50; 
}  } { 
for (i22 = 0; i22 < 9; i22++) { 
w_exp_tilde[i22] = ((double)(w_exp_tilde_tmp[i22]) - b_a * n_tilde[i22]) + ((1.0) - b_x) * b_n_tilde[i22]; 

}  } 
profileStart_GNC_codegen_SIL(76U); memset(&(b_dv[0]), 0, (9U) * sizeof(double)); profileEnd_GNC_codegen_SIL(76U); { 
for (i24 = 0; i24 < 3; i24++) { 
double d52; 
int i26; 
int i27; 
d52 = b_dv[3 * i24]; 
i26 = 3 * i24 + 1; 
i27 = 3 * i24 + 2; { 
for (i29 = 0; i29 < 3; i29++) { 
double d54; 
d54 = w_exp_tilde[i29 + 3 * i24]; 
d52 += (b_SD->pd->d_param).Jinv[3 * i29] * d54; 
b_dv[i26] += (b_SD->pd->d_param).Jinv[3 * i29 + 1] * d54; 
b_dv[i27] += (b_SD->pd->d_param).Jinv[3 * i29 + 2] * d54; 
F[(i29 + 11 * (i24 + 4)) + 4] = (0.0); 
}  } 
b_dv[3 * i24] = d52; 
}  } { 
for (i25 = 0; i25 < 3; i25++) { 
int F_tmp_tmp; 
F_tmp_tmp = 11 * (i25 + 4); { 
for (i28 = 0; i28 < 3; i28++) { 
double d53; 
d53 = (b_SD->pd->d_param).J[i28 + 3 * i25]; 
F[F_tmp_tmp + 4] += b_dv[3 * i28] * d53; 
F[F_tmp_tmp + 5] += b_dv[3 * i28 + 1] * d53; 
F[F_tmp_tmp + 6] += b_dv[3 * i28 + 2] * d53; 
}  } 
}  } 
b_dv[1] = t1_density * (m_a * x[9]); 
b_dv[4] = (0.0); 
b_dv[7] = t1_density * (m_a * x[7]); 
b_dv[2] = t1_density * (m_a * -x[8]); 
b_dv[5] = t1_density * (m_a * -x[7]); 
b_dv[8] = (0.0); 
c_x = (0.0); { 
for (i30 = 0; i30 < 3; i30++) { 
double d55; 
double d56; 
double d57; 
int c_F_tmp; 
b_dv[3 * i30] = (0.0); 
c_F_tmp = 11 * (i30 + 7); 
d55 = (0.0); 
d56 = (0.0); 
d57 = (0.0); { 
for (i31 = 0; i31 < 3; i31++) { 
double d58; 
d58 = b_dv[i31 + 3 * i30]; 
d55 += (dt * (b_SD->pd->d_param).Jinv[3 * i31]) * d58; 
d56 += (dt * (b_SD->pd->d_param).Jinv[3 * i31 + 1]) * d58; 
d57 += (dt * (b_SD->pd->d_param).Jinv[3 * i31 + 2]) * d58; 
}  } 
F[c_F_tmp + 6] = d57; 
F[c_F_tmp + 5] = d56; 
F[c_F_tmp + 4] = d55; 
c_x += x[i30 + 1] * (b_SD->pd->d_param).g[i30]; 
}  } 
d_x[0] = x[2] * (b_SD->pd->d_param).g[2] - (b_SD->pd->d_param).g[1] * x[3]; 
d_x[1] = (b_SD->pd->d_param).g[0] * x[3] - x[1] * (b_SD->pd->d_param).g[2]; 
d_x[2] = x[1] * (b_SD->pd->d_param).g[1] - (b_SD->pd->d_param).g[0] * x[2]; 
dv4[0] = (0.0); 
dv4[3] = x[0] * -(b_SD->pd->d_param).g[2]; 
dv4[6] = x[0] * (b_SD->pd->d_param).g[1]; 
dv4[1] = x[0] * (b_SD->pd->d_param).g[2]; 
dv4[4] = (0.0); 
dv4[7] = x[0] * -(b_SD->pd->d_param).g[0]; 
dv4[2] = x[0] * -(b_SD->pd->d_param).g[1]; 
dv4[5] = x[0] * (b_SD->pd->d_param).g[0]; 
dv4[8] = (0.0); 
skewed_exp_w_tmp[0] = (0.0); 
skewed_exp_w_tmp[3] = -x[9]; 
skewed_exp_w_tmp[6] = x[8]; 
skewed_exp_w_tmp[1] = x[9]; 
skewed_exp_w_tmp[4] = (0.0); 
skewed_exp_w_tmp[7] = -x[7]; 
skewed_exp_w_tmp[2] = -x[8]; 
skewed_exp_w_tmp[5] = x[7]; 
skewed_exp_w_tmp[8] = (0.0); 
b_skewed_exp_w_tmp[0] = (0.0); 
b_skewed_exp_w_tmp[3] = -x[6]; 
b_skewed_exp_w_tmp[6] = x[5]; 
b_skewed_exp_w_tmp[1] = x[6]; 
b_skewed_exp_w_tmp[4] = (0.0); 
b_skewed_exp_w_tmp[7] = -x[4]; 
b_skewed_exp_w_tmp[2] = -x[5]; 
b_skewed_exp_w_tmp[5] = x[4]; 
b_skewed_exp_w_tmp[8] = (0.0); 
n_a = (0.5) * (dt * dt); 
profileStart_GNC_codegen_SIL(77U); memset(&(c_skewed_exp_w_tmp[0]), 0, (9U) * sizeof(double)); profileEnd_GNC_codegen_SIL(77U); 
profileStart_GNC_codegen_SIL(78U); memset(&(b_dv[0]), 0, (9U) * sizeof(double)); profileEnd_GNC_codegen_SIL(78U); 
r_q_tmp[0] = x[0]; 
b_r_q_tmp = (0.0); { 
for (i32 = 0; i32 < 3; i32++) { 
double d59; 
double d60; 
double d_F_tmp; 
int e_F_tmp; 
int f_F_tmp; 
int g_F_tmp; 
F[i32 + 7] = dt * ((2.0) * (x[0] * (b_SD->pd->d_param).g[i32] - d_x[i32])); 
d_F_tmp = x[i32 + 1]; 
e_F_tmp = 11 * (i32 + 1); 
F[e_F_tmp + 7] = dt * ((2.0) * (((c_x * b_b[3 * i32] + x[1] * (b_SD->pd->d_param).g[i32]) - (b_SD->pd->d_param).g[0] * d_F_tmp) + dv4[3 * i32])); 



f_F_tmp = 3 * i32 + 1; 
F[e_F_tmp + 8] = dt * ((2.0) * (((c_x * b_b[f_F_tmp] + x[2] * (b_SD->pd->d_param).g[i32]) - (b_SD->pd->d_param).g[1] * d_F_tmp) + dv4[f_F_tmp])); 



g_F_tmp = 3 * i32 + 2; 
F[e_F_tmp + 9] = dt * ((2.0) * (((c_x * b_b[g_F_tmp] + x[3] * (b_SD->pd->d_param).g[i32]) - (b_SD->pd->d_param).g[2] * d_F_tmp) + dv4[g_F_tmp])); 



d59 = c_skewed_exp_w_tmp[3 * i32]; 
d60 = b_dv[3 * i32]; { 
for (i34 = 0; i34 < 3; i34++) { 
double d61; 
double d62; 
int b_skewed_exp_w_tmp_tmp; 
int i35; 
int skewed_exp_w_tmp_tmp; 
i35 = i34 + 3 * i32; 
d61 = b_skewed_exp_w_tmp[i35]; 
d62 = skewed_exp_w_tmp[i35]; 
d59 += skewed_exp_w_tmp[3 * i34] * d61; 
d60 += ((2.0) * b_skewed_exp_w_tmp[3 * i34]) * d62; 
skewed_exp_w_tmp_tmp = 3 * i34 + 1; 
c_skewed_exp_w_tmp[f_F_tmp] += skewed_exp_w_tmp[skewed_exp_w_tmp_tmp] * d61; 

b_dv[f_F_tmp] += ((2.0) * b_skewed_exp_w_tmp[skewed_exp_w_tmp_tmp]) * d62; 
b_skewed_exp_w_tmp_tmp = 3 * i34 + 2; 
c_skewed_exp_w_tmp[g_F_tmp] += skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp] * d61; 

b_dv[g_F_tmp] += ((2.0) * b_skewed_exp_w_tmp[b_skewed_exp_w_tmp_tmp]) * d62; 
}  } 
int h_F_tmp; 
int i_F_tmp; 
b_dv[3 * i32] = d60; 
c_skewed_exp_w_tmp[3 * i32] = d59; 
h_F_tmp = 11 * (i32 + 4); 
F[h_F_tmp + 7] = dt * skewed_exp_w_tmp[3 * i32] + n_a * (d59 - d60); 
i_F_tmp = 11 * (i32 + 7); 
F[i_F_tmp + 7] = w_exp_tilde[3 * i32]; 
F[h_F_tmp + 8] = dt * skewed_exp_w_tmp[f_F_tmp] + n_a * (c_skewed_exp_w_tmp[f_F_tmp] - b_dv[f_F_tmp]); 

F[i_F_tmp + 8] = w_exp_tilde[f_F_tmp]; 
F[h_F_tmp + 9] = dt * skewed_exp_w_tmp[g_F_tmp] + n_a * (c_skewed_exp_w_tmp[g_F_tmp] - b_dv[g_F_tmp]); 

F[i_F_tmp + 9] = w_exp_tilde[g_F_tmp]; 
r_q_tmp[i32 + 1] = -d_F_tmp; 
b_r_q_tmp += -d_F_tmp * x[i32 + 7]; 
}  } 
c_r_q_tmp[0] = r_q_tmp[2] * x[9] - r_q_tmp[3] * x[8]; 
c_r_q_tmp[1] = r_q_tmp[3] * x[7] - r_q_tmp[1] * x[9]; 
c_r_q_tmp[2] = r_q_tmp[1] * x[8] - r_q_tmp[2] * x[7]; { 
for (i33 = 0; i33 < 3; i33++) { 
double b_dt_tmp; 
double dt_tmp; 
int c_dt_tmp; 
int d_dt_tmp; 
int e_dt_tmp; 
dt_tmp = x[i33 + 7]; 
d_dt[i33] = dt * ((2.0) * (r_q_tmp[0] * dt_tmp - c_r_q_tmp[i33])); 
b_dt_tmp = r_q_tmp[i33 + 1]; 
c_dt_tmp = 3 * (i33 + 1); 
d_dt[c_dt_tmp] = dt * ((2.0) * (((b_r_q_tmp * b_b[3 * i33] + r_q_tmp[1] * dt_tmp) - x[7] * b_dt_tmp) + r_q_tmp[0] * skewed_exp_w_tmp[3 * i33])); 



d_dt_tmp = 3 * i33 + 1; 
d_dt[c_dt_tmp + 1] = dt * ((2.0) * (((b_r_q_tmp * b_b[d_dt_tmp] + r_q_tmp[2] * dt_tmp) - x[8] * b_dt_tmp) + r_q_tmp[0] * skewed_exp_w_tmp[d_dt_tmp])); 



e_dt_tmp = 3 * i33 + 2; 
d_dt[c_dt_tmp + 2] = dt * ((2.0) * (((b_r_q_tmp * b_b[e_dt_tmp] + r_q_tmp[3] * dt_tmp) - x[9] * b_dt_tmp) + r_q_tmp[0] * skewed_exp_w_tmp[e_dt_tmp])); 



}  } 
double p_a; 
F[10] = d_dt[0]; 
F[21] = d_dt[3]; 
F[32] = d_dt[6]; 
F[43] = d_dt[9]; 
o_a = r_q_tmp[0] * r_q_tmp[0] - ((r_q_tmp[1] * r_q_tmp[1] + r_q_tmp[2] * r_q_tmp[2]) + r_q_tmp[3] * r_q_tmp[3]); 


p_a = (2.0) * r_q_tmp[0]; 
b_dv[0] = (0.0); 
b_dv[3] = p_a * -r_q_tmp[3]; 
b_dv[6] = p_a * r_q_tmp[2]; 
b_dv[1] = p_a * r_q_tmp[3]; 
b_dv[4] = (0.0); 
b_dv[7] = p_a * -r_q_tmp[1]; 
b_dv[2] = p_a * -r_q_tmp[2]; 
b_dv[5] = p_a * r_q_tmp[1]; 
b_dv[8] = (0.0); { 
for (i36 = 0; i36 < 3; i36++) { 
F[11 * (i36 + 7) + 10] = dt * ((o_a * b_b[3 * i36] + ((2.0) * r_q_tmp[1]) * r_q_tmp[i36 + 1]) - b_dv[3 * i36]); 


}  } 
F[120] = (1.0); 
profileStart_GNC_codegen_SIL(79U); memset(&(b_F[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(79U); { 
for (i37 = 0; i37 < 11; i37++) { { 
for (i38 = 0; i38 < 11; i38++) { 
double d63; 
d63 = P[i38 + 11 * i37]; { 
for (i41 = 0; i41 < 11; i41++) { 
int j_F_tmp; 
j_F_tmp = i41 + 11 * i37; 
b_F[j_F_tmp] += F[i41 + 11 * i38] * d63; 
}  } 
}  } 
}  } { 
for (i39 = 0; i39 < 11; i39++) { { 
for (i40 = 0; i40 < 11; i40++) { 
double d64; 
d64 = (0.0); { 
for (i43 = 0; i43 < 11; i43++) { 
d64 += b_F[i39 + 11 * i43] * F[i40 + 11 * i43]; 
}  } 
int c_P_pred_tmp; 
c_P_pred_tmp = i39 + 11 * i40; 
P_pred[c_P_pred_tmp] = d64 + Q[c_P_pred_tmp]; 
}  } 
}  } { 
for (i42 = 0; i42 < 3; i42++) { 
int P_pred_tmp; 
int b_P_pred_tmp; 
int d_P_pred_tmp; 
P_pred_tmp = 11 * (i42 + 4); 
b_P_pred[3 * i42] = P_pred[P_pred_tmp + 4] + R[3 * i42]; 
b_P_pred_tmp = 3 * i42 + 1; 
b_P_pred[b_P_pred_tmp] = P_pred[P_pred_tmp + 5] + R[b_P_pred_tmp]; 
d_P_pred_tmp = 3 * i42 + 2; 
b_P_pred[d_P_pred_tmp] = P_pred[P_pred_tmp + 6] + R[d_P_pred_tmp]; 
}  } 
profileStart_GNC_codegen_SIL(80U); inv(b_P_pred, b_dv); profileEnd_GNC_codegen_SIL(80U); 
profileStart_GNC_codegen_SIL(81U); memset(&(K[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(81U); { 
for (i44 = 0; i44 < 3; i44++) { { 
for (i45 = 0; i45 < 3; i45++) { 
double d65; 
d65 = b_dv[i45 + 3 * i44]; { 
for (i46 = 0; i46 < 11; i46++) { 
int K_tmp; 
K_tmp = i46 + 11 * i44; 
K[K_tmp] += P_pred[i46 + 11 * (i45 + 4)] * d65; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(82U); memset(&(c_I[0]), 0, (121U) * sizeof(char)); profileEnd_GNC_codegen_SIL(82U); { 
for (c_k = 0; c_k < 11; c_k++) { 
c_I[c_k + 11 * c_k] = (1); 
}  } { 
for (i47 = 0; i47 < 44; i47++) { 
E[i47] = c_I[i47]; 
}  } { 
for (i48 = 0; i48 < 33; i48++) { 
E[i48 + 44] = (double)(c_I[i48 + 44]) - K[i48]; 
}  } { 
for (i49 = 0; i49 < 44; i49++) { 
E[i49 + 77] = c_I[i49 + 77]; 
}  } 
profileStart_GNC_codegen_SIL(83U); memset(&(b_E[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(83U); { 
for (i50 = 0; i50 < 11; i50++) { { 
for (i51 = 0; i51 < 11; i51++) { 
double d66; 
d66 = P_pred[i51 + 11 * i50]; { 
for (i53 = 0; i53 < 11; i53++) { 
int E_tmp; 
E_tmp = i53 + 11 * i50; 
b_E[E_tmp] += E[i53 + 11 * i51] * d66; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(84U); memset(&(b_K[0]), 0, (33U) * sizeof(double)); profileEnd_GNC_codegen_SIL(84U); { 
for (i52 = 0; i52 < 3; i52++) { { 
for (i54 = 0; i54 < 3; i54++) { 
double d67; 
d67 = R[i54 + 3 * i52]; { 
for (i55 = 0; i55 < 11; i55++) { 
int b_K_tmp; 
b_K_tmp = i55 + 11 * i52; 
b_K[b_K_tmp] += K[i55 + 11 * i54] * d67; 
}  } 
}  } 
}  } 
profileStart_GNC_codegen_SIL(85U); memset(&(c_E[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(85U); 
profileStart_GNC_codegen_SIL(86U); memset(&(c_K[0]), 0, (121U) * sizeof(double)); profileEnd_GNC_codegen_SIL(86U); { 
for (i56 = 0; i56 < 11; i56++) { { 
for (i57 = 0; i57 < 11; i57++) { 
double d68; 
d68 = E[i56 + 11 * i57]; { 
for (i60 = 0; i60 < 11; i60++) { 
int b_E_tmp; 
b_E_tmp = i60 + 11 * i56; 
c_E[b_E_tmp] += b_E[i60 + 11 * i57] * d68; 
}  } 
}  } { 
for (i59 = 0; i59 < 3; i59++) { 
double d70; 
d70 = K[i56 + 11 * i59]; { 
for (i61 = 0; i61 < 11; i61++) { 
int c_K_tmp; 
c_K_tmp = i61 + 11 * i56; 
c_K[c_K_tmp] += b_K[i61 + 11 * i59] * d70; 
}  } 
}  } 
}  } { 
for (i58 = 0; i58 < 121; i58++) { 
P_ret[i58] = c_E[i58] + c_K[i58]; 
}  } 
double d69; 
double d71; 
d69 = d1 / d3; 
d71 = d2 / d3; 
d72 = (((d1 / d4) * ((sens_input->board_gyro).meas[0] - bias->board_gyro[0]) + (d2 / d4) * ((sens_input->mti_gyro).meas[0] - bias->mti_gyro[0])) + C_ad_w_idx_0 * ((sens_input->ad_gyro).meas[0] - bias->ad_gyro[0])) - x_pred[4]; 



d73 = ((d69 * ((sens_input->board_gyro).meas[1] - bias->board_gyro[1]) + d71 * ((sens_input->mti_gyro).meas[1] - bias->mti_gyro[1])) + d * ((sens_input->ad_gyro).meas[1] - bias->ad_gyro[1])) - x_pred[5]; 



d74 = ((d69 * ((sens_input->board_gyro).meas[2] - bias->board_gyro[2]) + d71 * ((sens_input->mti_gyro).meas[2] - bias->mti_gyro[2])) + d * ((sens_input->ad_gyro).meas[2] - bias->ad_gyro[2])) - x_pred[6]; { 



for (i62 = 0; i62 < 11; i62++) { 
x_ret[i62] = x_pred[i62] + ((K[i62] * d72 + K[i62 + 11] * d73) + K[i62 + 22] * d74); 

}  } 
c_scale = (3.3121686421112381E-170); 
profileStart_GNC_codegen_SIL(87U); c_absxk = fabs(x_ret[0]); profileEnd_GNC_codegen_SIL(87U); { 
if (c_absxk > (3.3121686421112381E-170)) { 
b_q_mag = (1.0); 
c_scale = c_absxk; 
} else { 
c_t = c_absxk / (3.3121686421112381E-170); 
b_q_mag = c_t * c_t; 
}  } 
profileStart_GNC_codegen_SIL(88U); c_absxk = fabs(x_ret[1]); profileEnd_GNC_codegen_SIL(88U); { 
if (c_absxk > c_scale) { 
c_t = c_scale / c_absxk; 
b_q_mag = (b_q_mag * c_t) * c_t + (1.0); 
c_scale = c_absxk; 
} else { 
c_t = c_absxk / c_scale; 
b_q_mag += c_t * c_t; 
}  } 
profileStart_GNC_codegen_SIL(89U); c_absxk = fabs(x_ret[2]); profileEnd_GNC_codegen_SIL(89U); { 
if (c_absxk > c_scale) { 
c_t = c_scale / c_absxk; 
b_q_mag = (b_q_mag * c_t) * c_t + (1.0); 
c_scale = c_absxk; 
} else { 
c_t = c_absxk / c_scale; 
b_q_mag += c_t * c_t; 
}  } 
profileStart_GNC_codegen_SIL(90U); c_absxk = fabs(x_ret[3]); profileEnd_GNC_codegen_SIL(90U); { 
if (c_absxk > c_scale) { 
c_t = c_scale / c_absxk; 
b_q_mag = (b_q_mag * c_t) * c_t + (1.0); 
c_scale = c_absxk; 
} else { 
c_t = c_absxk / c_scale; 
b_q_mag += c_t * c_t; 
}  } 
profileStart_GNC_codegen_SIL(91U); b_q_mag = c_scale * sqrt(b_q_mag); profileEnd_GNC_codegen_SIL(91U); 
x_ret[0] /= b_q_mag; 
x_ret[1] /= b_q_mag; 
x_ret[2] /= b_q_mag; 
x_ret[3] /= b_q_mag; { 
if ((sens_input->board_baro).status) { 
memcpy(&(b_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(b_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
memcpy(&(e_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(d_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
b_ekf_correct(e_x_ret, d_P_ret, (sens_input->board_baro).meas, bias->board_baro, x_ret, P_ret); 

}  } { 
if ((sens_input->mti_baro).status) { 
memcpy(&(c_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(c_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
memcpy(&(g_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(f_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
b_ekf_correct(g_x_ret, f_P_ret, (sens_input->mti_baro).meas, bias->mti_baro, x_ret, P_ret); 

}  } { 
if ((sens_input->board_mag).status) { 
memcpy(&(d_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(e_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
memcpy(&(h_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(h_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
profileStart_GNC_codegen_SIL(92U); ekf_correct(h_x_ret, h_P_ret, (sens_input->board_mag).meas, bias->board_mag_earth, b_b, x_ret, P_ret); profileEnd_GNC_codegen_SIL(92U); 

}  } { 
if ((sens_input->mti_mag).status) { 
memcpy(&(f_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(g_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
memcpy(&(i_x_ret[0]), &(x_ret[0]), (11U) * sizeof(double)); 
memcpy(&(i_P_ret[0]), &(P_ret[0]), (121U) * sizeof(double)); 
profileStart_GNC_codegen_SIL(93U); ekf_correct(i_x_ret, i_P_ret, (sens_input->mti_mag).meas, bias->mti_mag_earth, b_b, x_ret, P_ret); profileEnd_GNC_codegen_SIL(93U); 

}  } 
}  } 
profileStart_GNC_codegen_SIL(94U); airdata->pressure = airdata_atmos(x_ret[10], &(airdata->temperature), &(airdata->density), &(airdata->sonic_speed), &(airdata->mach), &(airdata->dynamic_pressure)); profileEnd_GNC_codegen_SIL(94U); 


profileStart_GNC_codegen_SIL(95U); airspeed = b_norm(&(x_ret[7])); profileEnd_GNC_codegen_SIL(95U); 
airdata->mach = airspeed / airdata->sonic_speed; 
airdata->dynamic_pressure = ((0.5) * airdata->density) * (airspeed * airspeed); 
profileStart_GNC_codegen_SIL(96U); roll_state[0] = atan2((2.0) * (x_ret[2] * x_ret[3] + x_ret[0] * x_ret[1]), ((x_ret[0] * x_ret[0] - x_ret[1] * x_ret[1]) - x_ret[2] * x_ret[2]) + x_ret[3] * x_ret[3]); profileEnd_GNC_codegen_SIL(96U); 



roll_state[1] = x_ret[4]; 
} 
