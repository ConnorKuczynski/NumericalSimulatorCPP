/***********************  MANDELW.CPP  ***********************************
*                                                                        *
*               Simple Mandelbrodt program                               *
*                         using                                          *
*                       C M A T H                                        *
*        with C++ Builder, Visual C++, GCC, or CLang                     *
*                                                                        *
*   Copyright 1996-2020 by OptiCode - Dr. Martin Sander Software Dev.    *
*                                                                        *
*                                                                        *
*       This sample program is meant to demonstrate how to use           *
*       CMATH within your Windows programs. No efforts are made          *
*       to provide a comfortable user interface etc.                     *
*       Appearance could be considerably improved by refinement          *
*       of color definition.                                             *
*                                                                        *
**************************************************************************/
/*

Microsoft Visual C++:
    Open the project map  CDEMO_VS20??.sln  (32-bit)  
    or CDEMO64_VS20??.sln (64-bit)
	Select the project Mandel
	Compile and run

C++ Builder:
	Open 
    MANDEL.cbproj  (RAD Studio 2009, 2010, the XE series, and 10.x)
    MANDEL.bdsproj (BDS 2006 and 2007)
    MANDELB6.BPR   (BC++ Builder 6+)
    MANDELB.BPR    (BC++ Builder 4)
    On the command-line, type 
	    BCC32 -W -Iinclude mandelw.cpp lib\cmathf4d.lib lib\cmathfs.lib

GCC, CLang:
    Use GNU make to build the target Mandel in the accompanying makefile.
*/

#include <windows.h>                    /* Compiler's include files */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <winbase.h>

#include <newcplx.h>                      /* CMath include file */
HWND     hWndMain;

   /* Mandelbrodt parameters: */
fComplex  Offset(-0.5f,0), Asymmetry(-0.2f, +0.2f);
float     zoom = 1.0f;
unsigned  maxit = 120;
    // play around with these values!

#ifdef _WIN64
LRESULT CALLBACK
#else // preserve compatibility to older MSVC versions
LONG FAR PASCAL
#endif
MainMessageHandler (HWND, UINT, WPARAM, LPARAM);

const char AppName[] = "MandelDemo";     /* Name for the window */

int PASCAL WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance,
                   LPSTR /* lpCmdLine */, int /* nCmdShow */)
{
    MSG       msg;                      /* MSG structure to pass to windows proc */
    WNDCLASS  wc;

    if(!hPrevInstance)
    {
        wc.style      = CS_HREDRAW | CS_VREDRAW;
        wc.lpfnWndProc= MainMessageHandler;
        wc.cbClsExtra = 0;
        wc.cbWndExtra = 0;
        wc.hInstance  = hInstance;
        wc.hIcon      = LoadIcon (hInstance, AppName);
        wc.hCursor    = LoadCursor (NULL, IDC_ARROW);
        wc.hbrBackground  = (HBRUSH) GetStockObject (WHITE_BRUSH);
        wc.lpszMenuName   = AppName;
        wc.lpszClassName  = AppName;
        RegisterClass (&wc);
    }
    
                             /* create application's Main window:  */
    hWndMain = CreateWindow (AppName,
                             "CMATH Demo: Madelbrodt",
                             WS_OVERLAPPEDWINDOW,
                             CW_USEDEFAULT,     /* Use default X, Y, and width  */
                             CW_USEDEFAULT,
                             CW_USEDEFAULT,
                             CW_USEDEFAULT,
                             NULL,              /* Parent window's handle      */
                             NULL,              /* Default to Class Menu       */
                             hInstance,         /* Instance of window          */
                             NULL);             /* Create struct for WM_CREATE */


    if (hWndMain == NULL)
    {
        MessageBox(NULL, "Could not create window in WinMain", NULL, MB_ICONEXCLAMATION);
        return (1);
    }

    ShowWindow(hWndMain, SW_SHOWMAXIMIZED);     /* Display main window      */
    UpdateWindow(hWndMain);

    while(GetMessage(&msg, NULL, 0, 0)) /* Main message loop */
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

    UnregisterClass (AppName, hInstance);
    return (int)(msg.wParam);
}

#ifdef _WIN64
LRESULT CALLBACK
#else // preserve compatibility to older MSVC versions
LONG FAR PASCAL
#endif
MainMessageHandler(HWND hWnd, UINT Message, WPARAM wParam, LPARAM lParam)
{
    HDC         vDC;
    RECT        mrect;
    PAINTSTRUCT ps;
    HPEN        Pens[16];

        // set the Mandelbrodt parameters at the top of this file
        // (below the #include directives)! 
        // A more advanced user interface would have to allow for live
        // input of these parameters. For our purposes, setting them in
        // in the source code has to be enough.
    unsigned  i, j, kk;
    float     scal; 
    int       ii, jj;
    unsigned  pixx, pixy, width, width_2, height, height_2;
    fComplex  CM, ZM;

    switch (Message)                    /* Windows Message Loop           */
    {
        case WM_CREATE:
            break;

        case WM_PAINT:
            vDC = BeginPaint(hWndMain, &ps);
            GetClientRect( hWnd, &mrect );
            Pens[0]   = CreatePen( PS_SOLID, 1, 0x00000000 );  // Black
            Pens[1]   = CreatePen( PS_SOLID, 1, 0x00800000 );  // Dark Blue
            Pens[2]   = CreatePen( PS_SOLID, 1, 0x00008000 );  // Dark Green
            Pens[3]   = CreatePen( PS_SOLID, 1, 0x00000080 );  // Brown
            Pens[4]   = CreatePen( PS_SOLID, 1, 0x000000A0 );  // Dark Red
            Pens[5]   = CreatePen( PS_SOLID, 1, 0x00808000 );  // Cyan
            Pens[6]   = CreatePen( PS_SOLID, 1, 0x00800080 );  // Magenta
            Pens[7]   = CreatePen( PS_SOLID, 1, 0x00808080 );  // Dark Grey
            Pens[8]   = CreatePen( PS_SOLID, 1, 0x00B0B0B0 );  // Light Grey
            Pens[9]   = CreatePen( PS_SOLID, 1, 0x00FF0000 );  // Light Blue
            Pens[10]  = CreatePen( PS_SOLID, 1, 0x0000FF00 );  // Light Green
            Pens[11]  = CreatePen( PS_SOLID, 1, 0x000000FF );  // Light Red
            Pens[12]  = CreatePen( PS_SOLID, 1, 0x00FFFF00 );  // Light Cyan
            Pens[13]  = CreatePen( PS_SOLID, 1, 0x00FF00FF );  // Light Magenta
            Pens[14]  = CreatePen( PS_SOLID, 1, 0x0000FFFF );  // Yellow
            Pens[15]  = CreatePen( PS_SOLID, 1, 0x00FFFFFF );  // White
            width_2  = (width  = mrect.right - mrect.left) / 2;
            height_2 = (height = mrect.bottom - mrect.top) / 2;
            scal     = 1.5f / zoom / height_2;
            for (i = 0; i < width ; i++)
            {
                ii = i - width_2;
                for (j = 0; j < height; j++)
                {
                    jj = j - height_2;
                    ZM = fComplex((float)ii, (float)jj) * scal + Offset;
                    CM = ZM + Asymmetry;
                    for (kk = 0; kk < maxit; kk++)
                    {     // iterate until |Z| > 2, or a maximum of maxit times
                          // iterative Mandelbrodt formula
                        ZM = (ZM*ZM) + CM;
                        if( abs(ZM) > 2. ) break;
                    }  // kk now can be mapped to the plotting color
                       // for more details, use some nonlinear mapping function, like:
                    SelectObject( vDC, Pens[(int)(16.*fabs(sin( (double)kk )))] );
                    pixx = ii + width_2;  pixy = jj + height_2;
                    MoveToEx( vDC, pixx, pixy, NULL ); LineTo( vDC, pixx+1, pixy+1 );
                }
            }
            for( i=0; i<16; i++ ) DeleteObject( Pens[i] );
            EndPaint(hWndMain, &ps);
            break;

        case WM_CLOSE:
            DestroyWindow(hWnd);
            if (hWnd == hWndMain)
                PostQuitMessage(0);
            break;

        default:
            return (DefWindowProc(hWnd, Message, wParam, lParam));
    }  // end of switch( message )
    return (0);
}

