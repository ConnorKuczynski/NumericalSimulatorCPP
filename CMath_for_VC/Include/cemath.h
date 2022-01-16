/*	CEMATH.H

	Complex library for the languages C and C++.

	This header file contains all definitions for
	extended-precision complex numbers (complex long double).

	Copyright (c) 1996-2020 by OptiCode - Dr. Martin Sander Software Dev.
	Address of the author:
			OptiCode - Dr. Martin Sander Software Dev.
			Brahmsstr. 6
			D-32756 Detmold
			Germany
			optivec@gmx.de
			http://www.optivec.com

*/


#ifndef __CEMATH_H
#define __CEMATH_H

#if !defined( _CMATH_DEFS )
	#if (defined __BORLANDC__)
		#if (!defined __clang__)
			#pragma option -a-
		#else  /* bcc32c or 64-bit */
			#pragma pack( push,1 )
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
#else  /* Visual C++, Watcom, Linux */
	#if defined __GNUC__ || (defined __linux__) /* Standard for gcc is cdecl, other attributes ignored anyway */
		#define __cmf
	#else   /*  Visual C++, Watcom, Clang:  */
		#define __cmf __cdecl
	#endif
	#define VBOOL int
	#ifndef __CDMATH_H
		#include <cdmath.h>
	#endif
#endif

#ifdef __vLDBL_SUPPORT  /* 80-bit IEEE numbers supported */
	/* 80-bit IEEE numbers supported: the following 340 lines apply
	only to GCC, BC++ 32-bit, and to Linux (64 bit) */
#ifdef __cplusplus
	/*  first the constructors:
	    since eComplex and ePolar are declared as a structs instead of classes,
	    the constructors cannot get the names "eComplex" and "ePolar" here.	*/
  #ifndef _ECPLX_DEFINED
  inline eComplex __cmf ecplx( extended __ReVal )
  {	eComplex Result;
	Result.Re = __ReVal;
	Result.Im = 0.0;
	return Result;
  }

  inline ePolar __cmf epolr( extended __MagVal )
  {	ePolar Result;
	Result.Mag = __MagVal;
	Result.Arg = 0.0f;
	return Result;
  }

	//  identity:
  inline eComplex __cmf ecplx( eComplex __ze )
  {	return __ze;  } 
  inline ePolar __cmf epolr( ePolar __pe )
  {	return __pe;  }

	// up-conversions from single and double precision:
  inline eComplex __cmf ecplx( fComplex const & __zf )
  {	eComplex Result;
	Result.Re = __zf.Re;
	Result.Im = __zf.Im;
	return Result;
  }

  inline eComplex __cmf ecplx( dComplex const & __zd )
  {	eComplex Result;
	Result.Re = __zd.Re;
	Result.Im = __zd.Im;
	return Result;
  }

  inline ePolar __cmf epolr( fPolar const & __pf )
  {	ePolar Result;
	Result.Mag = __pf.Mag;
	Result.Arg = __pf.Arg;
	return Result;
  }

  inline ePolar __cmf epolr( dPolar const & __pd )
  {	ePolar Result;
	Result.Mag = __pd.Mag;
	Result.Arg = __pd.Arg;
	return Result;
  }
  #ifdef __COMPLEX_H
		// conversion from class complex
	 inline eComplex __cmf ecplx( complex const & __zc )
	 {	eComplex Result;
		 Result.Re = real(__zc);
		 Result.Im = imag(__zc);
		 return Result;
	 }
  #endif  // __COMPLEX_H

	// interconversions of Cartesian and polar
  eComplex __cmf ecplx( fPolar	__pf );
  eComplex __cmf ecplx( dPolar	__pd );
  eComplex __cmf ecplx( ePolar	__pe );
  ePolar	__cmf epolr( fComplex __zf );
  ePolar	__cmf epolr( dComplex __zd );
  ePolar	__cmf epolr( eComplex __ze );
  #define _ECPLX_DEFINED
  #endif  // _ECPLX_DEFINED
#endif	/* __cplusplus */

#if defined __cplusplus && !defined _CMATH_CLASSDEFS
extern "C" {  // the following functions cannot be "extern C",
#endif		// if eComplex is a class

#if !defined _CMATH_CLASSDEFS /* already declared for C++ in <nexcplx.h>? */
		/* basic form of constructor for C and C++ : */
	eComplex __cmf ecplx( extended __ReVal,  extended __ImVal);
	ePolar	__cmf epolr( extended __MagVal, extended __ArgVal );

	 /* plain-C versions of conversions */
	eComplex __cmf cftoce( fComplex __zf );
	eComplex __cmf cdtoce( dComplex __zd );
	ePolar	__cmf pftope( fPolar	__zf );
	ePolar	__cmf pdtope( dPolar	__zd );
	eComplex __cmf pftoce( fPolar	__pf );
	eComplex __cmf pdtoce( dPolar	__pd );
	eComplex __cmf petoce( ePolar	__pe );
	ePolar	__cmf cftope( fComplex __zf );
	ePolar	__cmf cdtope( dComplex __zd );
	ePolar	__cmf cetope( eComplex __ze );
#endif

	/* Basic complex operations. They are defined both
	for C and C++. However, for C++ you may as well use the
	overloaded operators and functions defined further below. */
#define		 ce_real( z )  (z).Re
#define		 pe_abs(  p )  (p).Mag
#if defined __clang__
	/* work-around for pointer bug in CLang (both Win and Linux) and bcc32c: */ 
	extended __cmf  ce_imag( eComplex __z );
	extended __cmf  pe_arg( ePolar __p );
#else		
	#define		 ce_imag( z )  (z).Im
	#define		 pe_arg(  p )  (p).Arg
#endif
eComplex __cmf  ce_neg(  eComplex __z );
eComplex __cmf  ce_conj( eComplex __z );
ePolar	__cmf  pe_neg(  ePolar	__p );
ePolar	__cmf  pe_conj( ePolar	__p );
ePolar	__cmf  pe_principal( ePolar __p );
#if defined __cplusplus && !defined _CMATH_CLASSDEFS
}
#endif
#ifdef __cplusplus  // even if _CMATH_CLASSDEFS
	extern "C" {
#endif
extended __cmf  ce_norm( eComplex __z );
extended __cmf  ce_arg(  eComplex __z );
extended __cmf  pe_norm( ePolar __p );
extended __cmf  pe_real( ePolar __p );
extended __cmf  pe_imag( ePolar __p );
#ifdef __cplusplus
	}
#endif
#if defined __cplusplus && !defined _CMATH_CLASSDEFS
extern "C" {  // the following functions cannot be "extern C",
#endif		// if eComplex identical with complex<extended>
eComplex __cmf  ce_polar( extended __mag, extended __angle );
eComplex __cmf  ce_magargtoc( extended __mag, extended __angle ); /* same as ce_polar */
ePolar	__cmf  pe_complex( extended __re,  extended __im );
ePolar	__cmf  pe_reimtop( extended __re,  extended __im ); /* same as pe_complex */

eComplex __cmf  ce_add(	eComplex __x, eComplex __y );
eComplex __cmf  ce_addRe( eComplex __x, extended __yRe );
eComplex __cmf  ce_sub(	eComplex __x, eComplex __y );
eComplex __cmf  ce_subRe( eComplex __x, extended __yRe );  /* x - yRe */
eComplex __cmf  ce_subrRe( eComplex __x, extended __yRe ); /* yRe - x */
eComplex __cmf  ce_mul(	eComplex __x, eComplex __y );
eComplex __cmf  ce_mulconj(	eComplex __x, eComplex __y );
eComplex __cmf  ce_mulRe( eComplex __x, extended __yRe );
eComplex __cmf  ce_div(	eComplex __x, eComplex __y );
eComplex __cmf  ce_divRe( eComplex __x, extended __yRe );  /*  x / yRe  */
eComplex __cmf  ce_divrRe( eComplex __x, extended __yRe ); /* yRe / x	*/

ePolar	__cmf  pe_mul(	ePolar	__x, ePolar	__y );
ePolar	__cmf  pe_mulconj(	ePolar	__x, ePolar	__y );
ePolar	__cmf  pe_mulRe(  ePolar	__x, extended __yRe );
ePolar	__cmf  pe_div(	ePolar	__x, ePolar	__y );
ePolar	__cmf  pe_divRe(  ePolar	__x, extended __yRe );	/*  x / yRe  */
ePolar	__cmf  pe_divrRe( ePolar	__x, extended __yRe );  /*  yRe / x  */

/*  mathematical functions with error handling through _matherr: */
#ifdef __cplusplus  // even if _CMATH_CLASSDEFS
	extern "C" extended __cmf  ce_abs(  eComplex __z );
#else
	extended __cmf  ce_abs(  eComplex __z );
#endif
eComplex __cmf  ce_acos( eComplex __z );
eComplex __cmf  ce_asin( eComplex __z );
eComplex __cmf  ce_atan( eComplex __z );
eComplex __cmf  ce_cos(  eComplex __z );
eComplex __cmf  ce_cosh( eComplex __z );
eComplex __cmf  ce_cubic( eComplex __z );  /* raise to the third power */
eComplex __cmf  ce_exp(  eComplex __z );
ePolar	__cmf  ce_exptop(  eComplex __z ); /* exp with result as ePolar */
eComplex __cmf  ce_inv(  eComplex __z );	/*	1.0 / z	*/
eComplex __cmf  ce_ipow( eComplex __z, int __exponent );
											 /* raise z to integer power */
eComplex __cmf  ce_ln(	eComplex __z );
eComplex __cmf  ce_log(	eComplex __z ); /* same as ce_ln */
eComplex __cmf  ce_log2(  eComplex __z );
eComplex __cmf  ce_log10( eComplex __z );
eComplex __cmf  ce_pow( eComplex __base, eComplex __exponent );
eComplex __cmf  ce_powReBase( extended __base, eComplex __exponent ); /* power of real base */
eComplex __cmf  ce_powReExpo( eComplex __base, extended __exponent ); /* raise z to real power */
						 /* for integer exponents, use ce_ipow ! */
eComplex __cmf  ce_quartic( eComplex __z );  /* raise to the fourth power */
eComplex __cmf  ce_sin(  eComplex __z );
eComplex __cmf  ce_sinh( eComplex __z );
eComplex __cmf  ce_square( eComplex __z );
eComplex __cmf  ce_sqrt( eComplex __z );
eComplex __cmf  ce_tan(  eComplex __z );
eComplex __cmf  ce_tanh( eComplex __z );

		/* some of the mathematical functions also for polar: */
ePolar	__cmf  pe_cubic( ePolar	__p );  /* raise to the third power */
ePolar	__cmf  pe_inv(  ePolar	__p );	/*	1.0 / z	*/
ePolar	__cmf  pe_ipow( ePolar	__p, int __exponent );
											 /* raise z to integer power */
eComplex __cmf  pe_lntoc(	ePolar	__p );
eComplex __cmf  pe_logtoc(	ePolar	__p ); /* same as pe_lntocf */
eComplex __cmf  pe_log2toc(  ePolar	__p );
eComplex __cmf  pe_log10toc( ePolar	__p );
ePolar	__cmf  pe_powReExpo( ePolar	__base, extended __exponent ); /* raise z to real power */
						 /* for integer exponents, use pe_ipow ! */
ePolar	__cmf  pe_quartic( ePolar	__p );  /* raise to the fourth power */
ePolar	__cmf  pe_square( ePolar	__p );
ePolar	__cmf  pe_sqrt( ePolar	__p );

#if defined __cplusplus && !defined _CMATH_CLASSDEFS
}	//  end of the extern "C" statement
ePolar	__cmf  pe_principal( extended __mag, extended __angle );
		// this overload cannot be extern C
#endif

#if defined __cplusplus && !defined __STD_COMPLEX && !defined __NEWCPLX_H
	/* in addition to the basic operations defined above for C,
	here is the same complete set of overloaded operators and
	functions as offered by <newcplx.h> for the complex classes.  */

	inline extended real( eComplex const & __z )
	{
		return __z.Re;
	}

	inline extended imag( eComplex const & __z )
	{
		return __z.Im;
	}

	inline eComplex neg( eComplex const & __z1 )
	{	eComplex Result;
		Result.Re = -__z1.Re;
		Result.Im = -__z1.Im;
		return Result;
	}

	inline eComplex conj( eComplex const & __z)
	{	eComplex Result;
		Result.Re =  __z.Re;
		Result.Im = -__z.Im;
		return Result;
	}

	extended __cmf  norm( eComplex __z );
	extended __cmf  arg(  eComplex __z );
	eComplex __cmf  magargtoc( extended Mag, extended Angle );  // identical
	//  unary operators:

	inline eComplex & operator +( eComplex & __z1 )
	{
		return __z1;
	}

	inline eComplex operator -( eComplex const & __z1 )
	{	eComplex Result;
		Result.Re = -__z1.Re;
		Result.Im = -__z1.Im;
		return Result;
	}

	//  binary operators:

	inline eComplex operator +( eComplex const & __z1, eComplex const & __z2 )
	{	eComplex Result;
		Result.Re = __z1.Re + __z2.Re;
		Result.Im = __z1.Im + __z2.Im;
		return Result;
	}

	inline eComplex operator +( eComplex const & __z1, extended __z2Re )
	{	eComplex Result;
		Result.Re = __z1.Re + __z2Re;
		Result.Im = __z1.Im;
		return Result;
	}

	inline eComplex operator +( extended __z1Re, eComplex const & __z2 )
	{	eComplex Result;
		Result.Re = __z1Re + __z2.Re;
		Result.Im = __z2.Im;
		return Result;
	}

	inline eComplex operator -( eComplex const & __z1, eComplex const & __z2 )
	{	eComplex Result;
		Result.Re = __z1.Re - __z2.Re;
		Result.Im = __z1.Im - __z2.Im;
		return Result;
	}

	inline eComplex operator -( eComplex const & __z1, extended __z2Re )
	{	eComplex Result;
		Result.Re = __z1.Re - __z2Re;
		Result.Im = __z1.Im;
		return Result;
	}

	inline eComplex operator -( extended __z1Re, eComplex const & __z2 )
	{	eComplex Result;
		Result.Re = __z1Re - __z2.Re;
		Result.Im = -__z2.Im;
		return Result;
	}

	inline eComplex operator *( eComplex const & __z1, eComplex const & __z2 )
	{	eComplex Result;
		Result.Re  = __z1.Re * __z2.Re - __z1.Im * __z2.Im;
		Result.Im  = __z1.Re * __z2.Im + __z1.Im * __z2.Re;
		return Result;
	}

	inline eComplex operator *( eComplex const & __z1, extended __z2Re )
	{	eComplex Result;
		Result.Re = __z1.Re * __z2Re;
		Result.Im = __z1.Im * __z2Re;
		return Result;
	}

	inline eComplex operator *( extended __z1Re, eComplex const & __z2 )
	{	eComplex Result;
		Result.Re  = __z1Re * __z2.Re;
		Result.Im  = __z1Re * __z2.Im;
		return Result;
	}

	eComplex __cmf operator /( eComplex const & __z1, eComplex const & __z2 );
	eComplex __cmf operator /( extended __z1Re, eComplex const & __z2 );
			// cannot be safely inlined for extended precision

	inline eComplex operator /( eComplex const & __z1, extended __z2Re )
	{	eComplex Result;
		Result.Re = __z1.Re / __z2Re;
		Result.Im = __z1.Im / __z2Re;
		return Result;
	}

		/*  compound-assignment operators:  */
	inline eComplex & operator +=( eComplex & __z1, eComplex const & __z2 )
	{
		__z1.Re += __z2.Re;
		__z1.Im += __z2.Im;
		return __z1;
	}

	inline eComplex & operator +=( eComplex & __z1, extended __z2Re )
	{
		__z1.Re += __z2Re;
		return __z1;
	}

	inline eComplex & operator -=( eComplex & __z1, eComplex const & __z2 )
	{
		__z1.Re -= __z2.Re;
		__z1.Im -= __z2.Im;
		return __z1;
	}

	inline eComplex & operator -=( eComplex & __z1, extended __z2Re )
	{
		__z1.Re -= __z2Re;
		return __z1;
	}

	inline eComplex & operator *=( eComplex & __z1, eComplex const & __z2 )
	{
		extended tmpRe;
		tmpRe	= __z1.Re * __z2.Re - __z1.Im * __z2.Im;
		__z1.Im = __z1.Re * __z2.Im + __z1.Im * __z2.Re;
		__z1.Re = tmpRe;
		return __z1;
	}

	inline eComplex & operator *=( eComplex & __z1, extended __z2Re )
	{
		__z1.Re *= __z2Re;
		__z1.Im *= __z2Re;
		return __z1;
	}

	eComplex & __cmf operator /=( eComplex & __z1, eComplex const & __z2 );

	inline eComplex & operator /=( eComplex & __z1, extended __z2Re )
	{
		__z1.Re /= __z2Re;
		__z1.Im /= __z2Re;
		return __z1;
	}


	inline VBOOL operator ==( eComplex const & __z1, eComplex const & __z2 )
	{
		return (__z1.Re == __z2.Re) && (__z1.Im == __z2.Im );
	}

	inline VBOOL operator ==( eComplex const & __z1, extended __z2Re )
	{
		return (__z1.Re == __z2Re) && (__z1.Im == 0.0 );
	}

	inline VBOOL operator !=( eComplex const & __z1, eComplex const & __z2 )
	{
		return (__z1.Re != __z2.Re) || (__z1.Im != __z2.Im );
	}

	inline VBOOL operator !=( eComplex const & __z1, extended __z2Re )
	{
		return (__z1.Im != 0.0 ) || (__z1.Re != __z2Re);
	}

  /* now the operators for polar: */

	extended  __cmf real(  ePolar __p );
	extended  __cmf imag(  ePolar __p );
	ePolar	__cmf reimtop( extended __Re, extended __Im );
	inline extended norm( ePolar const & __p )
	{
		return __p.Mag * __p.Mag;
	}
	inline extended arg( ePolar const & __p )
	{
		return __p.Arg;
	}
	ePolar	__cmf neg( ePolar __p );

	inline ePolar conj( ePolar const & __p)
	{	ePolar Result;
		Result.Mag =  __p.Mag;
		Result.Arg = -__p.Arg;
		return Result;
	}


  //  unary operators:

  inline ePolar & operator +( ePolar & __p1 )
  {
	return __p1;
  }

  inline ePolar operator -( ePolar const & __p1 )
  {
	return neg( __p1 );
  }

  //  binary operators:

  inline ePolar operator *( ePolar const & __p1, ePolar const & __p2 )
  {	ePolar Result;
	Result.Mag  = __p1.Mag * __p2.Mag;
	Result.Arg  = __p1.Arg + __p2.Arg;
	return Result;
  }

  inline ePolar operator *( ePolar const & __p1, extended __z2Re )
  {	ePolar Result;
	Result.Mag  = __p1.Mag * __z2Re;
	Result.Arg  = __p1.Arg;
	return Result;
  }

  inline ePolar operator *( extended __z1Re, ePolar const & __p2 )
  {	ePolar Result;
	Result.Mag  = __p2.Mag * __z1Re;
	Result.Arg  = __p2.Arg;
	return Result;
  }

  inline ePolar operator /( ePolar const & __p1, ePolar const & __p2 )
  {	ePolar Result;
	Result.Mag  = __p1.Mag / __p2.Mag;
	Result.Arg  = __p1.Arg - __p2.Arg;
	return Result;
  }

  inline ePolar operator /( ePolar const & __p1, extended __z2Re )
  {	ePolar Result;
	Result.Mag  = __p1.Mag / __z2Re;
	Result.Arg  = __p1.Arg;
	return Result;
  }

  ePolar __cmf operator /( extended __p1Re, ePolar const & __p2 );

	 /*  compound assignment operators:  */
  inline ePolar & operator *=( ePolar & __p1, ePolar const & __p2 )
  {
	__p1.Mag *= __p2.Mag;
	__p1.Arg += __p2.Arg;
	return __p1;
  }

  inline ePolar & __cmf operator *=( ePolar & __p1, extended __p2Re )
  {
	__p1.Mag *= __p2Re;
	return __p1;
  }

  inline ePolar & operator /=( ePolar & __p1, ePolar const & __p2 )
  {
	__p1.Mag /= __p2.Mag;
	__p1.Arg -= __p2.Arg;
	return __p1;
  }

  inline ePolar & __cmf operator /=( ePolar & __p1, extended __p2Re )
  {
	__p1.Mag /= __p2Re;
	return __p1;
  }

  inline VBOOL operator ==( ePolar const & __p1, ePolar const & __p2 )
  {
	return (__p1.Mag == __p2.Mag) &&
			 (__p1.Arg == __p2.Arg );
  }

  inline VBOOL operator ==( ePolar const & __p1, extended __p2Re )
  {
	return (__p1.Mag == __p2Re) && (pe_arg(__p1) == 0.0 );
  }

  inline VBOOL operator !=( ePolar const & __p1, ePolar const & __p2 )
  {
	return (__p1.Mag != __p2.Mag) || (pe_arg(__p1) != pe_arg(__p2) );
  }

  inline VBOOL operator !=( ePolar const & __p1, extended __p2Re )
  {
	return (__p1.Mag != __p2Re) || (pe_arg(__p1) != 0.0 );
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
	necessary here.)

  */

  extended __cmf  abs(  eComplex __z );
  eComplex __cmf  acos( eComplex __z );
  eComplex __cmf  asin( eComplex __z );
  eComplex __cmf  atan( eComplex __z );
  eComplex __cmf  cos(  eComplex __z );
  eComplex __cmf  cosh( eComplex __z );
  eComplex __cmf  cubic( eComplex __z );  /* raise to the third power */
  eComplex __cmf  exp(  eComplex __z );
  ePolar	__cmf  exptop(  eComplex __z );
  eComplex __cmf  inv(  eComplex __z );	/*	1.0 / z	*/
  eComplex __cmf  ipow( eComplex __z, int __exponent );
											/* raise z to integer power */
  eComplex __cmf  ln(	eComplex __z );
  eComplex __cmf  log(  eComplex __z );  /* same as ln */
  eComplex __cmf  log2( eComplex __z );
  eComplex __cmf  log10( eComplex __z );
  eComplex __cmf  pow( eComplex __z,  eComplex __exponent );
  eComplex __cmf  pow( eComplex __z,  extended __exponent ); // identical to powReExpo
  eComplex __cmf  pow( extended __base,  eComplex __exponent ); // identical to powReBase
  eComplex __cmf  powReBase( extended __base, eComplex __exponent ); // power of real base
  eComplex __cmf  powReExpo( eComplex __z, extended __exponent );	// raise z to real power
							// for integer exponents, use ipow !
  eComplex __cmf  quartic( eComplex __z );  // raise to the fourth power
  eComplex __cmf  sin(  eComplex __z );
  eComplex __cmf  sinh( eComplex __z );
  eComplex __cmf  square( eComplex __z );
  eComplex __cmf  sqrt( eComplex __z );
  eComplex __cmf  tan(  eComplex __z );
  eComplex __cmf  tanh( eComplex __z );

		/*  polar math functions:  */
  inline extended abs( ePolar const & __p )
  {
	return __p.Mag;
  }
  ePolar	__cmf  cubic( ePolar	__p );  /* raise to the third power */
  ePolar	__cmf  inv(  ePolar	__p );	/*	1.0 / p	*/
  ePolar	__cmf  ipow( ePolar	__p, int __exponent );
											 /* raise p to integer power */
  eComplex __cmf  lntoc(	ePolar	__p );
  eComplex __cmf  logtoc(	ePolar	__p ); /* same as lntoc */
  eComplex __cmf  log2toc(  ePolar	__p );
  eComplex __cmf  log10toc( ePolar	__p );
  ePolar	__cmf  pow( ePolar	__base, extended __exponent ); /* raise p to real power */
  ePolar	__cmf  powReExpo( ePolar	__base, extended __exponent ); /* raise p to real power */
						 /* for integer exponents, use ipow ! */
  ePolar	__cmf  quartic( ePolar	__p );  /* raise to the fourth power */
  ePolar	__cmf  square( ePolar	__p );
  ePolar	__cmf  sqrt( ePolar	__p );
#endif //  __cplusplus, not __STD_COMPLEX , __NEWCPLX_H

#else /* no 80-bit IEEE number support:
		The following 50 lines apply only to
		Visual C++, BC++ 64-bit, and Watcom:		*/

#include <cdmath.h>
#define ecplx	 dcplx
#if !defined _CMATH_CLASSDEFS /* declared as friend functions in <nexcplx.h> */
	#define cftoce	cftocd
	#define cdtoce(z) (z)
	#define cftope	cftopd
	#define cdtope	cdtopd
	#define cetope	cdtopd
	#define pftoce	pftocd
	#define pdtoce	pdtocd
	#define petoce	pdtocd
#endif
#define ce_real( z )  (z).Re
#define ce_imag( z )  (z).Im
#define ce_neg	cd_neg
#define ce_conj	 cd_conj
#define ce_norm	 cd_norm
#define ce_arg	cd_arg
#define ce_polar	cd_polar
#define ce_add	cd_add
#define ce_addRe	cd_addRe
#define ce_sub	cd_sub
#define ce_subRe	cd_subRe
#define ce_subrRe	cd_subrRe
#define ce_mul	cd_mul
#define ce_mulconj	cd_mulconj
#define ce_mulRe	cd_mulRe
#define ce_div	cd_div
#define ce_divRe	cd_divRe
#define ce_divrRe	cd_divrRe

#define ce_abs	cd_abs
#define ce_acos	 cd_acos
#define ce_asin	 cd_asin
#define ce_atan	 cd_atan
#define ce_cos	cd_cos
#define ce_cosh	 cd_cosh
#define ce_cubic	cd_cubic
#define ce_exp	cd_exp
#define ce_exptop	cd_exptop
#define ce_inv	cd_inv
#define ce_ipow	 cd_ipow
#define ce_ln		cd_ln
#define ce_log	cd_log
#define ce_log2	 cd_log2
#define ce_log10	cd_log10
#define ce_pow	cd_pow
#define ce_powReBase	cd_powReBase
#define ce_powReExpo	cd_powReExpo
#define ce_quartic  cd_quartic
#define ce_sin	cd_sin
#define ce_sinh	 cd_sinh
#define ce_square	cd_square
#define ce_sqrt	 cd_sqrt
#define ce_tan	cd_tan
#define ce_tanh	 cd_tanh

#define epolr		 dpolr
#if !defined _CMATH_CLASSDEFS /* declared as friend functions in <nexcplx.h> */
	#define pftope		pftopd
	#define pdtope(p)	 (p)
#endif
#define pe_abs(p)	 (p).Mag
#define pe_arg(p)	 (p).Arg
#define pe_complex	pd_complex
#define pe_neg		pd_neg
#define pe_conj		pd_conj
#define pe_norm		pd_norm
#define pe_real		pd_real
#define pe_imag		pd_imag
#define pe_principal  pd_principal
#define pe_complex	pd_complex
#define pe_mul		pd_mul
#define pe_mulconj	pd_mulconj
#define pe_mulRe	pd_mulRe
#define pe_div		pd_div
#define pe_divRe	pd_divRe
#define pe_divrRe	 pd_divrRe
#define pe_cubic	pd_cubic
#define pe_inv		pd_inv
#define pe_ipow		pd_ipow
#define pe_lntoc	pd_lntoc
#define pe_logtoc	 pd_logtoc
#define pe_log2toc	pd_log2toc
#define pe_log10toc	pd_log10toc
#define pe_pow		pd_pow
#define pe_powReExpo  pd_powReExpo
#define pe_quartic	pd_quartic
#define pe_square	 pd_square
#define pe_sqrt		pd_sqrt

#endif  /* Borland C++, Visual C++, or Watcom */


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
#endif /*  __CEMATH_H  */
