/********************   CDEMO.CPP   **************************************
*                                                                        *
*                Simple Demo Program for                                 *
*                       C M A T H                                        *
*        with C++ Builder, Visual C++, GCC, or CLang                     *
*                                                                        *
*   Copyright 1996-2020 by OptiCode - Dr. Martin Sander Software Dev.    *
*                                                                        *
*       This sample program provides a very basic test                   *
*       for the correct installation of CMATH                            *
*       and for correct error handling.                                  *
*                                                                        *
**************************************************************************/
/*
Visual C++:
    Open the project map  VDEMO_VS20??.sln / CDEMO_VS20??.sln (32-bit)  
    or VDEMO64_VS20??.sln / CDEMO_VS20??.sln (64-bit)
	Select the project CDemo
	Compile and run

C++ Builder:
	Open 
    CDEMO.cbproj  (RAD Studio 2009, 2010, the XE series, and 10.x)
    CDEMO.bdsproj (BDS 2006 and 2007)
    CDEMOB6.BPR   (BC++ Builder 6+)
    On the command-line, type
	    BCC32 -Iinclude cdemo.cpp lib\vcf4d.lib lib\vcfs.lib
	or
	    BCC32 -Iinclude cdemo.cpp lib\cmath4d.lib lib\cmathfs.lib

GCC, CLang:
    Use GNU make to build the target CDemo in the accompanying makefile.

*/

#include <newcplx.h>
#include <stdio.h>

int main( void )
{
    fComplex x, y;
    char     DataText[121];

    // V_setErrorEventFile( "NULL", 2 );
	     // as this is a console program, we might activate this line in order 
	     // to direct output to the screen instead of to popup windows

    V_printErrorMsg( "Welcome to the CMATH installation check!\n\n"
                     "Do not wonder to see a piece of error handling now\n"
					 "(for the debug library only)!\n");

	V_setFPErrorHandling( fperrNoteSING + fperrNoteOVERFLOW );
	x = sin( fComplex( 3.f, 5.f));
    y = inv( fComplex(0,0));  // will produce an error message in debug mode
    #if defined _MSC_VER && _MSC_VER >= 1400 
	    /* safe function of Visual Studio 2005 + and CLang */
        sprintf_s( DataText, 120,
                   "The sine of  { 3.0, 5.0 } is {%f, %f}\n"
                   "Complex 1/0 situations yield {%g, %g}\n",
                    x.Re, x.Im, y.Re, y.Im );
    #else
        sprintf( DataText, "The sine of  { 3.0, 5.0 } is {%f, %f}\n"
                           "Complex 1/0 situations yield {%g, %g}\n",
                           x.Re, x.Im, y.Re, y.Im );
    #endif
    V_printErrorMsg( DataText );
    V_printErrorMsg( "Thank you!\n"
                     "Enjoy using CMATH in your programs!" );
    return 0;
}
