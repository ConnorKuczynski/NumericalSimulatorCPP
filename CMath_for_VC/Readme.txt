      *****************************************************************
      *                                                               *
      *                   C M A T H  for  Visual C++                  *
      *                                                               *
      *                      Shareware Version 8                      *
      *****************************************************************

Contents
    1. Introduction
    2. System requirements
    3. Installation
    4. Running the example
    5. Documentation
    6. What's New?
    7. Copyright
    8. Registered Version

1. Introduction
---------------
CMATH is a comprehensive library for complex-number arithmetis and
mathematics, both in cartesian and in polar coordinates.
All functions may alternatively be called from classic C with
type-specific function names (like cf_sin, cd_exp, pe_sqrt),
or from C++ with overloaded function names and operators.

Superior speed, accuracy and safety are achieved through the implementation
in Assembly language (as opposed to the C++ code of other complex class
libraries).

Each of the three floating-point accuracies (float, double, and extended)
is given its own, optimized version of each function.


2. System requirements
----------------------
This version of CMATH is designed for PC systems,
equipped with a single- or multi-core CPU.

This package includes the libraries to be used  with the 
Microsoft Visual C++ compiler series (all versions of Visual Studio, 
down-compatible until VS 2005), 32-bit and/or 64-bit.

If you need versions for other compilers, please visit 
http://www.optivec.com/download 

At present, CMATH is available for the following target compilers:
- C++ Builder  (Embarcadero RAD Studio; formerly Borland C++)
- Microsoft Visual C++
- GCC, the GNU Compiler Collection
- LLVM Clang
- Delphi (Embarcadero RAD Studio)
- Lazarus / FreePascal


3. Installation
---------------

a) To install the Shareware version of CMATH, you need to execute INSTALL.EXE,
   contained in the ZIP file you downloaded. Only then, programs using the
   CMATH library can run. It is not possible to just copy the .lib and .h
   files, e.g., from another computer.

b) In case you wish to install CMATH for several target compilers
   (e.g., Visual Studio and RAD Studio), these versions may be installed 
   one after the other. You are free to choose the same directory or
   different directories, just as you like.

c) Windows has a habit of sometimes warning users after the installation
   of third-party products: "Setup possibly incomplete". If INSTALL.EXE
   did finish and displayed this README.TXT file, you know that the
   installation is complete.

d) After you completed the installation, you must set the search paths for
   include-files and for libraries according to your CMATH directory choice:
   Say, your CMATH directory is C:\CMATH. Then, these additional search 
   paths are:
   C:\CMATH\INCLUDE      for the include-files.
   C:\CMATH\LIB          for the libraries.

   Add these paths to the standard settings in the menues
   "Project / Settings / Configuration Settings / C/C++ / General / Additional include directories" and 
   "Project / Settings / Configuration Settings / Linker / General / Additional directories" 
   in MS Visual C++; or 
   "Extras / Options / Directories" 
   in older MS Visual C++ versions

e) You always have to include two CMATH libraries in your project. 
   One of these libraries is an interface or "base" library and specific for the 
   individual version of the target COMPILER. The other contains the bulk of the 
   high-performance functions and is specific for the target PROCESSOR. Details are
   described below for the example files. 
   For MS Visual C++, you also have to include the Windows API import library:
   In the menu   Project / (Configuration) Settings / Linker /
   Object and Library modules,
   you must add    user32.lib   (if it is not yet there).
   Otherwise you would get the linker error
   LNK2001: Unresolved external symbol __imp__MessageBoxA@16??

f) Each new release Visual Studio comes with its own specific version of the
   runtime DLL (with the exception of VS2017 and VS 2019, which do have their
   own runtimes, but also still allow to use the same runtime as VS2015). 
   Consequently, in the configurations with DYNAMIC runtime(code generation 
   options /MDd or /MD), any third-party library has to link exactly to that
   version of the runtime. This is why you find one set of CMATH base libraries
   for each version of Visual Studio.
   Only if you use the STATIC runtime (code generation options /MTd or /MT), 
   one and the same CMATH base library for a given configuration can be used
   with all versions of Visual Studio.
   There is a certain inconsistency in the description of the configurations
   in Visual Studio: The default configurations "Debug" and "Release" use the
   runtime library and MFC as DLL. This means that they are actually the 
   "Debug DLL" and "Release DLL" configurations.
   With these explanations, here are the libraries you have to choose:
   1. Configurations with STATIC VS runtime ("DebugStatic" configuration 
      in the examples):
      One CMATH base library covers all supported VS versions.
      For 64-bit, this is CMVCx64MTD.LIB, for 32-bit, it is CMVCMTD.LIB.
   2. Configurations with VS runtime as DLL:
                 64-bit Debug         32-bit Debug
      VS 2019    CMVC16x64MDD.LIB     CMVC16MDD.LIB
      VS 2017    CMVC15x64MDD.LIB     CMVC15MDD.LIB
      VS 2015    CMVC14x64MDD.LIB     CMVC14MDD.LIB
      VS 2013    CMVC12x64MDD.LIB     CMVC12MDD.LIB
      VS 2012    CMVC11x64MDD.LIB     CMVC11MDD.LIB
      VS 2010    *                    CMVC10MDD.LIB
      VS 2008    *                    CMVC9MDD.LIB
      VS 2005    CMVC8x64MDD.LIB      CMVC8MDD.LIB

      * The 64-bit libraries to be linked with the dynamic Runtime are 
        not available for VS 2008 and VS 2010.


4. Running the example
----------------------
Check your installation by compiling and running the small demo files.

 -  open the project map  CDEMO_VS20??.sln  (32-bit)  
    or CDEMO64_VS20??.sln (64-bit)
    (e.g., VS 2005 32-bit:   CDEMO_VS2005.sln, 
           VS 2015 64-bit:   CDEMO64_VS2015.sln, etc.)
    Depending on the exact version, Visual Studio may have to convert 
    this project map and its parts into a newer format. When prompted, 
    answer "Yes to all" to accept this automatic conversion.

 -  The project map contains two projects:
    CDEMO  and MANDEL.

    CDEMO   shows how to use CMATH functions for complex-number applications,
    MANDEL  shows a simple way how a Mandelbrodt plot can be coded with CMATH
           (actually, this is not the original Mandelbrodt, but a modified formula
            with several free parameters to play with.)

    Compile and run these projects separately.


5. Documentation
----------------
The full CMATH documentation is to be found in the file CMATH.HTM
to be read with a browser like Firefox, IE, Edge, Chrome, Opera, etc.


6. What's New?
--------------
Version 8.0
- Bug fixes in cd_arg, cd_abs

Version 7.3.6
- Bug fix in pd_mulRe (64-bit)

Version 7.3.4
-  Bug fix in the class version of cf_mul and cf_mulconj

Version 7.3.3
-  Bug fix in cf_neg

Version 7.3.2
-  Compatibility with the latest compiler release, 
   Visual Studio 2019

Version 7.3
- New: V_setFPErrorHandling
     This function allows to determine which types of math errors
     lead to a notification (popup window or console message,
     depending on setting in V_setErrorEventFile),
     and which types of math errors lead to the program being aborted.
     This function acts only in the debug libraries, as the production
     libraries always treat math errors silently and never actively
     abort program execution.
     The introduction of this function became necessary by the deprecation
     of the _matherr mechanism in all modern compilers. Should you
     still define _matherr in order to influence the behaviour of CMATH
     functions, you will need to switch over to V_setFPErrorHandling.
     
- New: Master License available, covering all supported target compilers
       (bundle, comprising all inividual "CMATH for xxx" products)

- bug fixes in
  c?_exptop (both 32-bit and 64-bit),
  c?_atan (64-bit)

Version 7.1
-  Compatibility with the latest compiler releases,
   Visual Studio 2017 and C++ Builder 10.2 Tokyo 
   
-  Bug fixes in cd_mulRe amd pf_powReExpo (64-bit only)

Version 7.0
-  Compatibility with latest compiler versions: 
   Visual Studio 2015, C++ Builder 10 Seattle and Berlin

-  Bug fix in pd_logtoC, pd_log2toC, pd_log10toC and their pe_ counterparts (32-bit only:
	would overwrite ebx on exit)

Version 6.5.8
-  Bug fixes in the 64-bit version for RAD Studio, affecting:
	pow( complex<float>, complex<float> )
	ipow( complex<float>, int )
	log( complex<float> ), log2, log10
	acos( complex<float> )
	inv( polar<float> )
   (The non-class variants, i.e. cf_pow, cf_ipow, cf_log etc., are
   not affected; the 32-bit and Visual C++ versions were ok, too.)
	
Version 6:
-  64-bit version (Visual Studio: at least 2005;  
                   C++ Builder / RAD Studio: at least XE3)

-  Compatibility with latest versions of target compilers

Version 5:
-  Improved thread safety,

-  Compatibility with latest versions of target compilers

Version 4:
-  In the functions for data-type down-conversion (cdtocf etc.),
   OVERFLOW errors do no longer lead to an error message, but are
   silently treated by setting the result to the maximum value possible
   (with the correct sign).

-   a few bug fixes and minor improvements of accuracy

Version 3:
The classes "polar" of float, double, and extended precision are introduced
along with a whole range of functions and operators in polar coordinates.
As a consequence, the member function "polar" of the cartesian complex
classes had to be replaced by "magargtoc".



7. Copyright
------------
The copyright owner of this product as a whole and of all its constituent
parts is
         OptiCode
         Dr. Martin Sander Software Development
         Brahmsstr. 6
         D-32756 Detmold
         Germany
         e-mail: optivec@gmx.de
         http://www.optivec.com

This Shareware version of CMATH is freely distributable in unchanged form.
For the distribution of applications created using CMATH, you need the
registered version. The detailed licence conditions are described in
chapter 1.2 of the file CMATH.HTM.


8. Registered Version
---------------------
If you like CMATH and decide to use it, consider registering:
The registered version of CMATH

-  has individually optimized libraries for each degree of processor
   backward-compatibility:
      Core2xxx, Core i3, i5 etc., AMD x64  (64-bit version only)
      486DX/Pentium+

-  entitles you to two years of free updates (by download from our web site).

-  costs USD 60 or EUR 59 for the commercial edition,
         USD 39 or EUR 39 for the educational edition,
   and can be ordered most conveniently through ShareIt:

   CMATH for MS Visual C++ / Visual Studio:
      https://order.shareit.com/product?productid=103422

   alternatively:
   CMATH Master License for all supported compilers:
      https://order.shareit.com/product?productid=300860338


See chapter 1.3 of the file CMATH.HTM for more ordering options and
volume discounts.


    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

Copyright (C) OptiCode - Dr. Martin Sander Software Dev. 1996-2020
