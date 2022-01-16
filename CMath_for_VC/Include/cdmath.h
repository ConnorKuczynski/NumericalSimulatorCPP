/*	CDMATH.H

	Complex library for the languages C and C++.

	This header file contains all definitions for
	double-precision complex numbers (complex double).

	Copyright (c) 1996-2020 by OptiCode - Dr. Martin Sander Software Dev.
	Address of the author:
		OptiCode - Dr. Martin Sander Software Dev.
		Brahmsstr. 6
		D-32756 Detmold
		Germany
		optivec@gmx.de
		http://www.optivec.com
*/


#ifndef __CDMATH_H
#define __CDMATH_H

#if !defined( _CMATH_DEFS )
	#if (defined __BORLANDC__)
		#if (!defined __clang__)
			#pragma option -a-
		#else  /* bcc32c or 64-bit */
			#pragma pack( push,8 )
		#endif
	#else /* Visual C++, GCC, Watcom */
		#pragma pack( push,1 )
	#endif /* avoid insertion of dummy bytes  */
	typedef struct {float	Re,  Im;}  fComplex;
	typedef struct {double	Re,  Im;}  dComplex;
	typedef struct {float	Mag, Arg;} fPolar;
	typedef struct {double	Mag, Arg;} dPolar;
    #if ((defined __BORLANDC__) && !(defined _WIN64)) || (defined __GNUC__) || (defined __linux__)
        #define __vLDBL_SUPPORT   /* GCC, 32-bit BC, Linux compilers support IEEE 80-bit real numbers */
		typedef long double	extended;
		#if (!defined __clang__)
			typedef struct {extended Re,  Im;}  eComplex;
			typedef struct {extended Mag, Arg;} ePolar;
		#else  /* automatic struct-member alignment of bcc32c does not work properly */
			typedef struct {extended Re;  short pad1; extended Im;  short pad2;}  eComplex;
			typedef struct {extended Mag; short pad1; extended Arg; short pad2;}  ePolar;
		#endif
		#if ((defined __GNUC__) && !(defined __clang__)) && ((defined __linux__) || (defined _WIN64))
			typedef __float128 great;
			typedef short      half;   // place holder for half-float
			typedef struct {great  Re,  Im;}  gComplex;
			typedef struct {great  Mag, Arg;} gPolar;
			typedef struct {half   Re,  Im;}  hComplex;
			typedef struct {half   Mag, Arg;} hPolar;
			typedef gComplex gcomplex;
			typedef hComplex hcomplex;
			typedef gPolar   gpolar;
			typedef hPolar   hpolar;
		#endif
	#else /* Visual C++, BC 64-bit */
		typedef  double extended; /* no 80-bit IEEE numbers. So make
								     extended equal to double	*/
		typedef dComplex  eComplex;
		typedef dPolar	ePolar;
	#endif	
	#if (defined __BORLANDC__) && (!defined __clang__)
		#pragma option -a.
	#else
		#pragma pack( pop )
	#endif /* restore default data packing */
	typedef fComplex fcomplex;
	typedef dComplex dcomplex;
	typedef eComplex ecomplex;
	typedef fPolar fpolar;
	typedef dPolar dpolar;
	typedef ePolar epolar;  // tolerate all-lower case
	#define _CMATH_DEFS
#endif
#ifdef __BORLANDC__
	#include <_defs.h>
	#if (__BORLANDC__ >= 0x450)
		 #define __cmf _RTLENTRY _EXPFUNC
	#else
		 #define __cmf  _Cdecl _FARFUNC
	#endif
	#if __BORLANDC__ < 0x500
		#define VBOOL int
	#else
		#define VBOOL bool
	#endif
#else  /* Visual C++, Watcom, linux */
	#if defined __GNUC__  /* Standard for gcc is cdecl, other attributes ignored anyway */
		#define __cmf
	#else   /*  Visual C++, Watcom, Clang:  */
		#define __cmf __cdecl
	#endif
	#define VBOOL int
#endif

/*  first the constructors:  */
#ifdef __cplusplus
	/* since dComplex and dPolar are declared as a structs instead of classes,
	the constructors cannot get the names "dComplex" and "dPolar" here.	*/
  #ifndef _DCPLX_DEFINED
  inline dComplex __cmf dcplx( double __ReVal )
  {	dComplex Result;
	Result.Re = __ReVal;
	Result.Im = 0.0;
	return Result;
  }

  inline dPolar __cmf dpolr( double __MagVal )
  {	dPolar Result;
	Result.Mag = __MagVal;
	Result.Arg = 0.0;
	return Result;
  }

	// up-conversions from single precision
  inline dComplex __cmf dcplx( fComplex const & __zf )
  {	dComplex Result;
	Result.Re = __zf.Re;
	Result.Im = __zf.Im;
	return Result;
  }

  inline dPolar __cmf dpolr( fPolar const & __pf )
  {	dPolar Result;
	Result.Mag = __pf.Mag;
	Result.Arg = __pf.Arg;
	return Result;
  }
	// down-conversion from extended precision
	// (with OVERFLOW error handling):
  #if defined __vLDBL_SUPPORT  /* 80-bit IEEE numbers supported */
	dComplex __cmf dcplx( eComplex __ze );
	dPolar	__cmf dpolr( ePolar	__pe );
	//  identity: 32-bit BC++ does not like referenced definition
	inline dComplex __cmf dcplx( dComplex __zd )
	{	return __zd;  }
	inline dPolar	__cmf dpolr( dPolar __pd )
	{	return __pd;  }
  
  #else  // MSVC and BCC64: eComplex == dComplex
	inline dComplex __cmf dcplx( eComplex const & __ze )
	{	return (dComplex) __ze;  }
	inline dPolar __cmf dpolr( ePolar const & __pe )
	{	return (dPolar) __pe;  }

		//  identity: referenced definition possible
	inline dComplex & __cmf dcplx( dComplex & __zd )
	{	return __zd;  }
	inline dPolar	& __cmf dpolr( dPolar & __pd )
	{	return __pd;  }
  #endif  /* 80-bit real supported or not */

	// interconversions of Cartesian and polar
  dComplex __cmf dcplx( fPolar	__pf );
  dComplex __cmf dcplx( dPolar	__pd );
  dPolar	__cmf dpolr( fComplex __zf );
  dPolar	__cmf dpolr( dComplex __zd );
  #if defined __vLDBL_SUPPORT  /* 80-bit IEEE numbers supported */
	dComplex __cmf dcplx( ePolar	__pe );
	dPolar	__cmf dpolr( eComplex __ze );
  #else
	#ifdef _CMATH_CLASSDEFS
		dComplex __cmf dcplx( polar<long double> __pe );
		dPolar __cmf dpolr( complex<long double> __ze );
	#endif
  #endif

  #define _DCPLX_DEFINED
  #endif  // _DCPLX_DEFINED
#endif  /* __cplusplus */


#if defined __cplusplus && !defined _CMATH_CLASSDEFS
extern "C" {  // the following functions cannot be "extern C",
#endif		// if dComplex is a class

#if !defined _CMATH_CLASSDEFS /* already declared for C++ in <newcplx.h>? */
	/* basic form of constructor for C and C++ : */
	dComplex __cmf dcplx( double __ReVal,  double __ImVal);
	dPolar	__cmf dpolr( double __MagVal, double __ArgVal );

	/* plain-C version of conversions: 
		 (with OVERFLOW handling, where appropriate)	*/

	dComplex __cmf cftocd( fComplex __zf );
	dPolar	__cmf pftopd( fPolar	__pf );
	dComplex __cmf pftocd( fPolar	__pf );
	dComplex __cmf pdtocd( dPolar	__pd );
	dPolar	__cmf cftopd( fComplex __zf );
	dPolar	__cmf cdtopd( dComplex __zd );
    #if defined __vLDBL_SUPPORT  /* 80-bit IEEE numbers supported */
		dComplex __cmf cetocd( eComplex __ze );
		dPolar	__cmf petopd( ePolar	__pe );
		dComplex __cmf petocd( ePolar	__pe );
		dPolar	__cmf cetopd( eComplex __ze );
	#else	/* extended = double */
		#define cetocd(z) (z)
		#define petopd(z) (z)
		#define petocd pdtocd
		#define cetopd cdtopd
	#endif
#endif

	/* Basic complex operations. They are defined both
	for C and C++. However, for C++ you may as well use the
	overloaded operators and functions defined further below. */
#define		 cd_real( z )  (z).Re
#define		 cd_imag( z )  (z).Im
#define		 pd_abs(  p )  (p).Mag
#define		 pd_arg(  p )  (p).Arg
dComplex __cmf  cd_neg(  dComplex __z );
dComplex __cmf  cd_conj( dComplex __z );
dPolar	__cmf  pd_neg(  dPolar	__p );
dPolar	__cmf  pd_conj( dPolar	__p );
dPolar	__cmf  pd_principal( dPolar __p );
#if defined __cplusplus && !defined _CMATH_CLASSDEFS
}
#endif
#ifdef __cplusplus  // even if _CMATH_CLASSDEFS
	extern "C" {
#endif
double __cmf  cd_norm( dComplex __z );
double __cmf  cd_arg(  dComplex __z );
double __cmf  pd_norm( dPolar __p );
double __cmf  pd_real( dPolar __p );
double __cmf  pd_imag( dPolar __p );
#ifdef __cplusplus
} // end of extern "C"
#endif
#if defined __cplusplus && !defined _CMATH_CLASSDEFS
extern "C" {  // the following functions cannot be "extern C",
#endif		// if dComplex identical to complex<double>
dComplex __cmf  cd_polar( double __mag, double __angle );
dComplex __cmf  cd_magargtoc( double __mag, double __angle ); /* same as cd_polar */
dPolar	__cmf  pd_complex( double __re,  double __im );
dPolar	__cmf  pd_reimtop( double __re,  double __im ); /* same as pd_complex */

dComplex __cmf  cd_add(	dComplex __x, dComplex __y );
dComplex __cmf  cd_addRe( dComplex __x, double __yRe );
dComplex __cmf  cd_sub(	dComplex __x, dComplex __y );
dComplex __cmf  cd_subRe( dComplex __x, double __yRe );  /*  x - yRe  */
dComplex __cmf  cd_subrRe( dComplex __x, double __yRe ); /*  yRe - x  */
dComplex __cmf  cd_mul(	dComplex __x, dComplex __y );
dComplex __cmf  cd_mulconj(	dComplex __x, dComplex __y );
dComplex __cmf  cd_mulRe( dComplex __x, double __yRe );
dComplex __cmf  cd_div(	dComplex __x, dComplex __y );
dComplex __cmf  cd_divRe( dComplex __x, double __yRe );	/*  x / yRe  */
dComplex __cmf  cd_divrRe( dComplex __x, double __yRe );  /*  yRe / x  */

dPolar	__cmf  pd_mul(	dPolar	__x, dPolar	__y );
dPolar	__cmf  pd_mulconj(	dPolar	__x, dPolar	__y );
dPolar	__cmf  pd_mulRe(  dPolar	__x, double __yRe );
dPolar	__cmf  pd_div(	dPolar	__x, dPolar	__y );
dPolar	__cmf  pd_divRe(  dPolar	__x, double __yRe );	/*  x / yRe  */
dPolar	__cmf  pd_divrRe( dPolar	__x, double __yRe );  /*  yRe / x  */

/*  mathematical functions with error handling through _matherr: */
#ifdef __cplusplus
	extern "C" double	__cmf  cd_abs(  dComplex __z );
#else
	double	__cmf  cd_abs(  dComplex __z );
#endif
dComplex __cmf  cd_acos( dComplex __z );
dComplex __cmf  cd_asin( dComplex __z );
dComplex __cmf  cd_atan( dComplex __z );
dComplex __cmf  cd_cos(  dComplex __z );
dComplex __cmf  cd_cosh( dComplex __z );
dComplex __cmf  cd_cubic( dComplex __z );  /* raise to the third power */
dComplex __cmf  cd_exp(  dComplex __z );
dPolar	__cmf  cd_exptop(  dComplex __z ); /* exp with result as dPolar */
dComplex __cmf  cd_inv(  dComplex __z );	/*	1.0 / z	*/
dComplex __cmf  cd_ipow( dComplex __z, int __exponent );
											 /* raise z to integer power */
dComplex __cmf  cd_ln(	dComplex __z );
dComplex __cmf  cd_log(	dComplex __z );  /* sane as cd_ln */
dComplex __cmf  cd_log2(  dComplex __z );
dComplex __cmf  cd_log10( dComplex __z );
dComplex __cmf  cd_pow( dComplex __base, dComplex __exponent );
dComplex __cmf  cd_powReBase( double __base, dComplex __exponent ); /* power of real base */
dComplex __cmf  cd_powReExpo( dComplex __base, double __exponent ); /* raise z to real power */
						 /* for integer exponents, use cd_ipow ! */
dComplex __cmf  cd_quartic( dComplex __z );  /* raise to the fourth power */
dComplex __cmf  cd_sin(  dComplex __z );
dComplex __cmf  cd_sinh( dComplex __z );
dComplex __cmf  cd_square( dComplex __z );
dComplex __cmf  cd_sqrt( dComplex __z );
dComplex __cmf  cd_tan(  dComplex __z );
dComplex __cmf  cd_tanh( dComplex __z );

		/* some of the mathematical functions also for polar: */
dPolar	__cmf  pd_cubic( dPolar	__p );  /* raise to the third power */
dPolar	__cmf  pd_inv(  dPolar	__p );	/*	1.0 / z	*/
dPolar	__cmf  pd_ipow( dPolar	__p, int __exponent );
											 /* raise z to integer power */
dComplex __cmf  pd_lntoc(	dPolar	__p );
dComplex __cmf  pd_logtoc(	dPolar	__p ); /* same as pd_lntocd */
dComplex __cmf  pd_log2toc(  dPolar	__p );
dComplex __cmf  pd_log10toc( dPolar	__p );
dPolar	__cmf  pd_powReExpo( dPolar	__base, double __exponent ); /* raise z to real power */
						 /* for integer exponents, use pd_ipow ! */
dPolar	__cmf  pd_quartic( dPolar	__p );  /* raise to the fourth power */
dPolar	__cmf  pd_square( dPolar	__p );
dPolar	__cmf  pd_sqrt( dPolar	__p );


#if defined __cplusplus && !defined _CMATH_CLASSDEFS
}	//  end of the extern "C" statement
dPolar	__cmf  pd_principal( double __mag, double __angle );
		// this overload cannot be extern C
#endif

#if defined __cplusplus && !defined __STD_COMPLEX && !defined __COMPLEX_H && !defined __NEWCPLX_H
	/* in addition to the basic operations defined above for C,
	here is the same complete set of overloaded operators and
	functions as offered by <newcplx.h> for the complex classes.  */

	inline double real( dComplex const & __z )
	{
		return __z.Re;
	}

	inline double imag( dComplex const & __z )
	{
		return __z.Im;
	}

	inline dComplex neg( dComplex const & __z1 )
	{	dComplex Result;
		Result.Re = -__z1.Re;
		Result.Im = -__z1.Im;
		return Result;
	}

	inline dComplex conj( dComplex const & __z)
	{	dComplex Result;
		Result.Re =  __z.Re;
		Result.Im = -__z.Im;
		return Result;
	}

	double	__cmf  norm( dComplex __z );
	double	__cmf  arg(  dComplex __z );
	dComplex __cmf  magargtoc( double Mag, double Angle );


  //  unary operators:

  inline dComplex & operator +( dComplex & __z1 )
  {
	return __z1;
  }

  inline dComplex operator -( dComplex const & __z1 )
  {	dComplex Result;
	Result.Re = -__z1.Re;
	Result.Im = -__z1.Im;
	return Result;
  }


  //  binary operators:

  inline dComplex operator +( dComplex const & __z1, dComplex const & __z2 )
  {	dComplex Result;
	Result.Re = __z1.Re + __z2.Re;
	Result.Im = __z1.Im + __z2.Im;
	return Result;
  }

  inline dComplex operator +( dComplex const & __z1, double __z2Re )
  {	dComplex Result;
	Result.Re = __z1.Re + __z2Re;
	Result.Im = __z1.Im;
	return Result;
  }

  inline dComplex operator +( double __z1Re, dComplex const & __z2 )
  {	dComplex Result;
	Result.Re = __z1Re + __z2.Re;
	Result.Im = __z2.Im;
	return Result;
  }

  inline dComplex operator -( dComplex const & __z1, dComplex const & __z2 )
  {	dComplex Result;
	Result.Re = __z1.Re - __z2.Re;
	Result.Im = __z1.Im - __z2.Im;
	return Result;
  }

  inline dComplex operator -( dComplex const & __z1, double __z2Re )
  {	dComplex Result;
	Result.Re = __z1.Re - __z2Re;
	Result.Im = __z1.Im;
	return Result;
  }

  inline dComplex operator -( double __z1Re, dComplex const & __z2 )
  {	dComplex Result;
	Result.Re = __z1Re - __z2.Re;
	Result.Im = -__z2.Im;
	return Result;
  }

  inline dComplex operator *( dComplex const & __z1, dComplex const & __z2 )
  {	dComplex Result;
	Result.Re  = __z1.Re * __z2.Re - __z1.Im * __z2.Im;
	Result.Im  = __z1.Re * __z2.Im + __z1.Im * __z2.Re;
	return Result;
  }

  inline dComplex operator *( dComplex const & __z1, double __z2Re )
  {	dComplex Result;
	Result.Re = __z1.Re * __z2Re;
	Result.Im = __z1.Im * __z2Re;
	return Result;
  }

  inline dComplex operator *( double __z1Re, dComplex const & __z2 )
  {	dComplex Result;
	Result.Re  = __z1Re * __z2.Re;
	Result.Im  = __z1Re * __z2.Im;
	return Result;
  }

  inline dComplex operator /( dComplex const & __z1, double __z2Re )
  {	dComplex Result;
	Result.Re = __z1.Re / __z2Re;
	Result.Im = __z1.Im / __z2Re;
	return Result;
  }


  dComplex __cmf operator /( dComplex const & __z1, dComplex const & __z2 );
  dComplex __cmf operator /( double __z1Re, dComplex const & __z2 );

	/* compound-assignment operators:  */
  inline dComplex & operator +=( dComplex & __z1, dComplex const & __z2 )
  {
	__z1.Re += __z2.Re;
	__z1.Im += __z2.Im;
	return __z1;
  }

  inline dComplex & operator +=( dComplex & __z1, double __z2Re )
  {
	__z1.Re += __z2Re;
	return __z1;
  }

  inline dComplex & operator -=( dComplex & __z1, dComplex const & __z2 )
  {
	__z1.Re -= __z2.Re;
	__z1.Im -= __z2.Im;
	return __z1;
  }

  inline dComplex & operator -=( dComplex & __z1, double __z2Re )
  {
	__z1.Re -= __z2Re;
	return __z1;
  }

  inline dComplex & operator *=( dComplex & __z1, dComplex const & __z2 )
  {
	double tmpRe;
	tmpRe	= __z1.Re * __z2.Re - __z1.Im * __z2.Im;
	__z1.Im = __z1.Re * __z2.Im + __z1.Im * __z2.Re;
	__z1.Re = tmpRe;
	return __z1;
  }

  inline dComplex & operator *=( dComplex & __z1, double __z2Re )
  {
	__z1.Re *= __z2Re;
	__z1.Im *= __z2Re;
	return __z1;
  }

  dComplex & __cmf operator /=( dComplex & __z1, dComplex const & __z2 );

  inline dComplex & operator /=( dComplex & __z1, double __z2Re )
  {
	__z1.Re /= __z2Re;
	__z1.Im /= __z2Re;
	return __z1;
  }


  inline VBOOL operator ==( dComplex const & __z1, dComplex const & __z2 )
  {
	return (__z1.Re == __z2.Re) && (__z1.Im == __z2.Im );
  }

  inline VBOOL operator ==( dComplex const & __z1, double __z2Re )
  {
	return (__z1.Re == __z2Re) && (__z1.Im == 0.0 );
  }

  inline VBOOL operator !=( dComplex const & __z1, dComplex const & __z2 )
  {
	return (__z1.Re != __z2.Re) || (__z1.Im != __z2.Im );
  }

  inline VBOOL operator !=( dComplex const & __z1, double __z2Re )
  {
	return (__z1.Im != 0.0 ) || (__z1.Re != __z2Re);
  }


  /* now the operators for polar: */

	double  __cmf real(  dPolar __p );
	double  __cmf imag(  dPolar __p );
	dPolar __cmf reimtop( double __Re, double __Im );
	inline double norm( dPolar const & __p )
	{
		return __p.Mag * __p.Mag;
	}
	inline double arg( dPolar const & __p )
	{
		return __p.Arg;
	}

	dPolar  __cmf neg(  dPolar __p );

	inline dPolar conj( dPolar const & __p)
	{	dPolar Result;
		Result.Mag =  __p.Mag;
		Result.Arg = -__p.Arg;
		return Result;
	}


  //  unary operators:

  inline dPolar & operator +( dPolar & __p1 )
  {
	return __p1;
  }

  inline dPolar operator -( dPolar const & __p1 )
  {
	return neg( __p1 );
  }

  //  binary operators:

  inline dPolar operator *( dPolar const & __p1, dPolar const & __p2 )
  {	dPolar Result;
	Result.Mag  = __p1.Mag * __p2.Mag;
	Result.Arg  = __p1.Arg + __p2.Arg;
	return Result;
  }

  inline dPolar operator *( dPolar const & __p1, double __z2Re )
  {	dPolar Result;
	Result.Mag  = __p1.Mag * __z2Re;
	Result.Arg  = __p1.Arg;
	return Result;
  }

  inline dPolar operator *( double __z1Re, dPolar const & __p2 )
  {	dPolar Result;
	Result.Mag  = __p2.Mag * __z1Re;
	Result.Arg  = __p2.Arg;
	return Result;
  }

  inline dPolar operator /( dPolar const & __p1, dPolar const & __p2 )
  {	dPolar Result;
	Result.Mag  = __p1.Mag / __p2.Mag;
	Result.Arg  = __p1.Arg - __p2.Arg;
	return Result;
  }

  inline dPolar operator /( dPolar const & __p1, double __z2Re )
  {	dPolar Result;
	Result.Mag  = __p1.Mag / __z2Re;
	Result.Arg  = __p1.Arg;
	return Result;
  }

  dPolar __cmf operator /( double __p1Re, dPolar const & __p2 );

	/*  compound assinment operators:  */
  inline dPolar & operator *=( dPolar & __p1, dPolar const & __p2 )
  {
	__p1.Mag *= __p2.Mag;
	__p1.Arg += __p2.Arg;
	return __p1;
  }

  inline dPolar & __cmf operator *=( dPolar & __p1, double __z2Re )
  {
	__p1.Mag *= __z2Re;
	return __p1;
  }

  inline dPolar & operator /=( dPolar & __p1, dPolar const & __p2 )
  {
	__p1.Mag /= __p2.Mag;
	__p1.Arg -= __p2.Arg;
	return __p1;
  }

  inline dPolar & __cmf operator /=( dPolar & __p1, double __z2Re )
  {
	__p1.Mag /= __z2Re;
	return __p1;
  }

  inline VBOOL operator ==( dPolar const & __p1, dPolar const & __p2 )
  {
	return (__p1.Mag == __p2.Mag) &&
			 (pd_arg(__p1) == pd_arg(__p2) );
  }

  inline VBOOL operator ==( dPolar const & __p1, double __p2Re )
  {
	return (__p1.Mag == __p2Re) && (pd_arg(__p1) == 0.0 );
  }

  inline VBOOL operator !=( dPolar const & __p1, dPolar const & __p2 )
  {
	return (__p1.Mag != __p2.Mag) || (pd_arg(__p1) != pd_arg(__p2) );
  }

  inline VBOOL operator !=( dPolar const & __p1, double __p2Re )
  {
	return (__p1.Mag != __p2Re) || (pd_arg(__p1) != 0.0 );
  }

  /*  C++ version of the mathematical functions defined above.
	They use the same code as the C versions. In case of an error,
	you get a message in which the name of the C version is
	stated.
	Note that these functions require complex arguments to be
	passed by value, not by reference, as it is done in the member
	functions of the class complex. In terms of efficiency, this
	is about the same. (The math functions of the class complex
	store complex results at intermediate addresses and copy them
	to the desired address afterwards. This final copy is not
	necessary here.)											 */

  double	__cmf  abs(  dComplex __z );
  dComplex __cmf  acos( dComplex __z );
  dComplex __cmf  asin( dComplex __z );
  dComplex __cmf  atan( dComplex __z );
  dComplex __cmf  cos(  dComplex __z );
  dComplex __cmf  cosh( dComplex __z );
  dComplex __cmf  cubic( dComplex __z );  /* raise to the third power */
  dComplex __cmf  exp(  dComplex __z );
  dPolar	__cmf  exptop(  dComplex __z );
  dComplex __cmf  inv(  dComplex __z );	/*	1.0 / z	*/
  dComplex __cmf  ipow( dComplex __z, int __exponent );
											/* raise z to integer power */
  dComplex __cmf  ln(  dComplex __z );
  dComplex __cmf  log(  dComplex __z );  /* sane as ln */
  dComplex __cmf  log2( dComplex __z );
  dComplex __cmf  log10( dComplex __z );
  dComplex __cmf  pow( dComplex __z, dComplex __exponent );
  dComplex __cmf  pow( dComplex __z,  double __exponent ); // identical to powReExpo
  dComplex __cmf  pow( double __base,  dComplex __exponent ); // identical to powReBase
  dComplex __cmf  powReBase( double __base, dComplex __exponent ); // power of real base
  dComplex __cmf  powReExpo( dComplex __z, double __exponent );	// raise z to real power
							// for integer exponents, use ipow !
  dComplex __cmf  quartic( dComplex __z );  // raise to the fourth power
  dComplex __cmf  sin(  dComplex __z );
  dComplex __cmf  sinh( dComplex __z );
  dComplex __cmf  square( dComplex __z );
  dComplex __cmf  sqrt( dComplex __z );
  dComplex __cmf  tan(  dComplex __z );
  dComplex __cmf  tanh( dComplex __z );

		/*  polar math functions:  */
  inline double abs( dPolar const & __p )
  {
	return __p.Mag;
  }
  dPolar	__cmf  cubic( dPolar	__p );  /* raise to the third power */
  dPolar	__cmf  inv(  dPolar	__p );	/*	1.0 / z	*/
  dPolar	__cmf  ipow( dPolar	__p, int __exponent );
											 /* raise z to integer power */
  dComplex __cmf  lntoc(	dPolar	__p );
  dComplex __cmf  logtoc(	dPolar	__p ); /* same as lntoc */
  dComplex __cmf  log2toc(  dPolar	__p );
  dComplex __cmf  log10toc( dPolar	__p );
  dPolar	__cmf  pow( dPolar	__base, double __exponent ); /* raise z to real power */
  dPolar	__cmf  powReExpo( dPolar	__base, double __exponent ); /* raise z to real power */
						 /* for integer exponents, use ipow ! */
  dPolar	__cmf  quartic( dPolar	__p );  /* raise to the fourth power */
  dPolar	__cmf  square( dPolar	__p );
  dPolar	__cmf  sqrt( dPolar	__p );
#endif //  __cplusplus, not __STD_COMPLEX , __COMPLEX_H, __NEWCPLX_H


/********** Floating-point math error handling ***************

    The following constants and functions are relevant only for
	the OptiVec / CMATH debug libraries. All others treat errors 
	silently and continue execution after setting default results.
*/
#ifdef __cplusplus
	extern "C" {
#endif

#if !defined __VFPHAND_DEFINES
typedef enum
{
        fperrIgnore          = 0,
        fperrNoteDOMAIN      = 0x0001,
        fperrNoteSING        = 0x0002,
        fperrNoteOVERFLOW    = 0x0004,
        fperrNoteTLOSS       = 0x0010,
        fperrAbortDOMAIN     = 0x0101,
        fperrAbortSING       = 0x0202,
        fperrAbortOVERFLOW   = 0x0404,
        fperrAbortTLOSS      = 0x1010,
        fperrDefaultHandling = 0x0107   /* default = fperrAbortDOMAIN + fperrNoteSING + fperrNoteOVERFLOW */
}   V_fphand;
#define __VFPHAND_DEFINES
#endif

/* error handling functions, borrowed from VectorLib:  */
void  __cmf  V_setFPErrorHandling( int handlingMode );
void  __cmf  V_noteError( const char *fname, unsigned why );
void  __cmf  V_printErrorMsg( const char *ErrMsg );
void  __cmf  V_setErrorEventFile( const char *filename,  unsigned ScreenAndFile );
void  __cmf  V_closeErrorEventFile( void );

        /* library identification */
const char * __cmf V_getLibVersion( void );

#ifdef __cplusplus
}	// end of extern "C"
#endif
#endif	/*  __CDMATH_H  */
