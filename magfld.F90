! ---------------------------------------------------------------
! magcof_mod: replaces COMMON /MAGCOF/ used in COFRM, DYPOL, FELDG.
! ---------------------------------------------------------------
    module magcof_mod
        implicit none
        integer, save :: NMAX = 0
        real,    save :: GB(255) = 0.0
        real,    save :: GV(225) = 0.0
        integer, save :: ICHG = -99999
    end module magcof_mod

    SUBROUTINE COFRM (DATE)
!          Define the International Geomagnetic Reference Field (IGRF) as a
!          scalar potential field using a truncated series expansion with
!          Schmidt semi-normalized associated Legendre functions of degree n and
!          order m.  The polynomial coefficients are a function of time and are
!          interpolated between five year epochs or extrapolated at a constant
!          rate after the last epoch.

!          INPUTS:
!            DATE = yyyy.fraction (UT)
!          OUTPUTS (in magcof_mod):
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed in COFRM

!          It is fatal to supply a DATE before the first epoch.  A warning is
!          issued to Fortran unit 0 (stderr) if DATE is later than the
!          recommended limit, five years after the last epoch.

!          HISTORY (blame):
!          Apr 1983:  Written by Vincent B. Wickwar (Utah State Univ.) including
!          secular variation acceleration rate set to zero in case the IGRF
!          definition includes such second time derivitives.  The maximum degree
!          (n) defined was 10.

!          Jun 1986:  Updated coefficients adding Definitive Geomagnetic Reference
!          Field (DGRF) for 1980 and IGRF for 1985 (EOS Volume 7 Number 24).  The
!          designation DGRF means coefficients will not change in the future
!          whereas IGRF coefficients are interim pending incorporation of new
!          magnetometer data.  Common block MAG was replaced by MAGCOF, thus
!          removing variables not used in subroutine FELDG.  (Roy Barnes)

!          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as given in EOS
!          Volume 73 Number 16 April 21 1992.  Other changes were made so future
!          updates should:  (1) Increment NDGY; (2) Append to EPOCH the next IGRF
!          year; (3) Append the next DGRF coefficients to G1DIM and H1DIM; and (4)
!          replace the IGRF initial values (G0, GT) and rates of change indices
!          (H0, HT).

!          Apr 1994 (Art Richmond): Computation of GV added, for finding magnetic
!          potential.

!          Aug 1995 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
!          obtained by anonymous ftp to geomag.gsfc.nasa.gov (cd pub, mget table*)
!          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov) with
!          problems reported to baldwin@geomag.gsfc.nasa.gov.

!          Oct 1995 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see email
!          in folder).  Also found bug whereby coefficients were not being updated
!          in FELDG when IENTY did not change so ICHG was added to flag date
!          changes.  Also, a vestigial switch (IS) was removed from COFRM; it was
!          always zero and involved 3 branch if statements in the main polynomial
!          construction loop now numbered 200.

!          Feb 1999 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid the
!          possibility of compiler or loader options initializing memory to
!          something else (e.g., indefinite).  Also simplify the algebra in COFRM
!          with no effect on results.

!          Mar 1999 (Barnes):  Removed three branch if's from FELDG and changed
!          statement labels to ascending order.

!          Jun 1999 (Barnes):  Corrected RTOD definition in GD2CART.

!          May 2000 (Barnes):  Replace IGRF 1995, add IGRF 2000, and extend the
!          earlier DGRF's back to 1900.  The coefficients came from an NGDC web
!          page.  Related documentation is in $APXROOT/docs/igrf.2000.*  where
!          $APXROOT, defined by 'source envapex', is traditionally ~bozo/apex).

!          Mar 2004 (Barnes):  Replace 1995 and 2000 coefficients; now both are
!          DGRF.  Coefficients for 2000 are degree 13 with precision increased to
!          tenths nT and accommodating this instigated changes:  (1) degree (NMAX)
!          is now a function of epoch (NMXE) to curtail irrelevant looping over
!          unused high order terms (n > 10 in epochs before 2000) when calculating
!          GB; (2) expand coefficients data statement layout for G1D and H1D,
!          formerly G1DIM and H1DIM; (3) omit secular variation acceleration terms
!          which were always zero; (4) increase array dimensions in common block
!          MAGCOF and associated arrays G and H in FELDG; (5) change earth's shape
!          in CONVRT from the IAU-1966 to the WGS-1984 spheroid; (6) eliminate
!          reference to 'definitive' in variables in COFRM which were not always
!          definitive; (7) change G to GB in COFRM s.t. arrays GB and GV in common
!          block MAGCOF are consistently named in all subroutines; (8) remove
!          unused constants in all five subroutines.  See EOS Volume 84 Number 46
!          November 18 2003, www.ngdc.noaa.gov/IAGA/vmod/igrf.html or local files
!          $APXROOT/docs/igrf.2004.*

!          Sept. 2005 (Maute): update with IGRF10 from
!          http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html use script
!          ~maute/apex.d/apex_update/igrf2f Note that the downloaded file the start
!          column of the year in the first line has to be before the start of each
!          number in the same column

!          Jan. 2010 (Maute) update with IGRF11 (same instructions as Sep. 2005
!          comment)

!          Jan. 2015 (Maute) update with IGRF12 (same instructions as Sep. 2005
!          comment)

    use magcof_mod
    implicit none

    ! Argument
    real, intent(in) :: DATE

    ! Local variables (previously implicit)
    integer          :: I, IY, IY1, N, M, NN, MM, I1
    real             :: TIME, T, TO5, F, RNN
    double precision :: F0

    ! ICHG initialised in magcof_mod (= -99999); DATEL saved here.
    real, save :: DATEL = -999.
     
!          NEPO = Number of epochs
!          NGH  = Single dimensioned array size of 2D version (GYR or HYR)
!          NGHT = Single dimensioned array size of 2D version (GT  or HT)
    integer,PARAMETER :: NEPO = 24, NGH = 225*NEPO, NGHT = 225
    real :: GYR(15,15,NEPO), HYR(15,15,NEPO), EPOCH(NEPO), &
    GT (15,15),      HT (15,15),       NMXE(NEPO), &
    GY1D(NGH),       HY1D(NGH), &
    GT1D(NGHT),      HT1D(NGHT)
    EQUIVALENCE (GYR(1,1,1),GY1D(1)), (HYR(1,1,1),HY1D(1)), &
    (GT (1,1),  GT1D(1)), (HT (1,1),  HT1D(1))

    SAVE EPOCH, NMXE, GYR, HYR, GT, HT, GY1D, HY1D, GT1D, HT1D
    DATA &
    EPOCH / 1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940, &
    &              1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, &
    &              1990, 1995, 2000, 2005, 2010, 2015/, &
    NMXE  /   10,   10,   10,   10,   10,   10,   10,   10,   10, &
    &                10,   10,   10,   10,   10,   10,   10,   10,   10, &
    &                10,   10,   13,   13,   13,   13/

!          g(n,m) for 1900
!          Fields across a line are (degree) n=1,13; lines are (order) m=0,13 as indicated
!          in column 6; e.g., for 1965 g(n=3,m=0) = 1297 or g(n=6,m=6) = -111

!           1       2       3      4      5      6      7      8      9
!                                        10     11     12     13          (n)
    DATA (GY1D(I),I=1,145) /0, &
    -31543,   -677,  1022,   876,  -184,    63,    70,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2298,   2905, -1469,   628,   328,    61,   -55,     8,    10, &
    -4,     0,     0,     0,   3*0, &
    &              924,  1256,   660,   264,   -11,     0,    -4,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     572,  -361,     5,  -217,    34,    -9,   -11, &
    -5,     0,     0,     0,   5*0, &
    &                            134,   -86,   -58,   -41,     1,    12, &
    -2,     0,     0,     0,   6*0, &
    -16,    59,   -21,     2,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -90,    18,    -9,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     5,     2, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,    -1, &
    &                                     2,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=146,225) / &
    &                                     2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1905
    DATA (GY1D(I),I=226,370) /0, &
    -31464,   -728,  1037,   880,  -192,    62,    70,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2298,   2928, -1494,   643,   328,    60,   -54,     8,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1041,  1239,   653,   259,   -11,     0,    -4,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     635,  -380,    -1,  -221,    33,    -9,   -11, &
    -5,     0,     0,     0,   5*0, &
    &                            146,   -93,   -57,   -41,     1,    12, &
    -2,     0,     0,     0,   6*0, &
    -26,    57,   -20,     2,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -92,    18,    -8,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     5,     2, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,     0, &
    &                                     2,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=371,450) / &
    &                                     2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1910
    DATA (GY1D(I),I=451,595) /0, &
    -31354,   -769,  1058,   884,  -201,    62,    71,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2297,   2948, -1524,   660,   327,    58,   -54,     8,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1176,  1223,   644,   253,   -11,     1,    -4,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     705,  -400,    -9,  -224,    32,    -9,   -11, &
    -5,     0,     0,     0,   5*0, &
    &                            160,  -102,   -54,   -40,     1,    12, &
    -2,     0,     0,     0,   6*0, &
    -38,    54,   -19,     2,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -95,    18,    -8,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     5,     2, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,     0, &
    &                                     2,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=596,675) / &
    &                                     2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1915
    DATA (GY1D(I),I=676,820) /0, &
    -31212,   -802,  1084,   887,  -211,    61,    72,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2306,   2956, -1559,   678,   327,    57,   -54,     8,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1309,  1212,   631,   245,   -10,     2,    -4,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     778,  -416,   -16,  -228,    31,    -9,   -11, &
    -5,     0,     0,     0,   5*0, &
    &                            178,  -111,   -51,   -38,     2,    12, &
    -2,     0,     0,     0,   6*0, &
    -51,    49,   -18,     3,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -98,    19,    -8,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     6,     2, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,     0, &
    &                                     1,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=821,900) / &
    &                                     2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1920
    DATA (GY1D(I),I=901,1045) /0, &
    -31060,   -839,  1111,   889,  -221,    61,    73,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2317,   2959, -1600,   695,   326,    55,   -54,     7,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1407,  1205,   616,   236,   -10,     2,    -3,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     839,  -424,   -23,  -233,    29,    -9,   -11, &
    -5,     0,     0,     0,   5*0, &
    &                            199,  -119,   -46,   -37,     2,    12, &
    -2,     0,     0,     0,   6*0, &
    -62,    44,   -16,     4,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -101,    19,    -7,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     6,     2, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,     0, &
    &                                     1,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=1046,1125) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1925
    DATA (GY1D(I),I=1126,1270) /0, &
    -30926,   -893,  1140,   891,  -230,    61,    73,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2318,   2969, -1645,   711,   326,    54,   -54,     7,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1471,  1202,   601,   226,    -9,     3,    -3,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     881,  -426,   -28,  -238,    27,    -9,   -11, &
    -5,     0,     0,     0,   5*0, &
    &                            217,  -125,   -40,   -35,     2,    12, &
    -2,     0,     0,     0,   6*0, &
    -69,    39,   -14,     4,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -103,    19,    -7,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     7,     2, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,     0, &
    &                                     1,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=1271,1350) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1930
    DATA (GY1D(I),I=1351,1495) /0, &
    -30805,   -951,  1172,   896,  -237,    60,    74,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2316,   2980, -1692,   727,   327,    53,   -54,     7,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1517,  1205,   584,   218,    -9,     4,    -3,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     907,  -422,   -32,  -242,    25,    -9,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            234,  -131,   -32,   -34,     2,    12, &
    -2,     0,     0,     0,   6*0, &
    -74,    32,   -12,     5,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -104,    18,    -6,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     8,     3, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          8,     0, &
    &                                     1,     0,     0,     0,  10*0, &
    -2/
    DATA (GY1D(I),I=1496,1575) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1935
    DATA (GY1D(I),I=1576,1720) /0, &
    -30715,  -1018,  1206,   903,  -241,    59,    74,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2306,   2984, -1740,   744,   329,    53,   -53,     7,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1550,  1215,   565,   211,    -8,     4,    -3,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     918,  -415,   -33,  -246,    23,    -9,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            249,  -136,   -25,   -33,     1,    11, &
    -2,     0,     0,     0,   6*0, &
    -76,    25,   -11,     6,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -106,    18,    -6,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   6,     8,     3, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          7,     0, &
    &                                     2,     0,     0,     0,  10*0, &
    -2/
    DATA (GY1D(I),I=1721,1800) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1940
    DATA (GY1D(I),I=1801,1945) /0, &
    -30654,  -1106,  1240,   914,  -241,    57,    74,    11,     8, &
    -3,     0,     0,     0,   2*0, &
    -2292,   2981, -1790,   762,   334,    54,   -53,     7,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1566,  1232,   550,   208,    -7,     4,    -3,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     916,  -405,   -33,  -249,    20,   -10,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            265,  -141,   -18,   -31,     1,    11, &
    -2,     0,     0,     0,   6*0, &
    -76,    18,    -9,     6,     1, &
    &                                     6,     0,     0,     0,   7*0, &
    -107,    17,    -5,    -2, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   5,     9,     3, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          7,     1, &
    &                                     2,     0,     0,     0,  10*0, &
    -2/
    DATA (GY1D(I),I=1946,2025) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1945
    DATA (GY1D(I),I=2026,2170) /0, &
    -30594,  -1244,  1282,   944,  -253,    59,    70,    13,     5, &
    -3,     0,     0,     0,   2*0, &
    -2285,   2990, -1834,   776,   346,    57,   -40,     7,   -21, &
    &                                    11,     0,     0,     0,   3*0, &
    &             1578,  1255,   544,   194,     6,     0,    -8,     1, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     913,  -421,   -20,  -246,     0,    -5,   -11, &
    &                                     2,     0,     0,     0,   5*0, &
    &                            304,  -142,   -25,   -29,     9,     3, &
    -5,     0,     0,     0,   6*0, &
    -82,    21,   -10,     7,    16, &
    -1,     0,     0,     0,   7*0, &
    -104,    15,   -10,    -3, &
    &                                     8,     0,     0,     0,   8*0, &
    &                                                  29,     7,    -4, &
    -1,     0,     0,     0,   9*0, &
    &                                                          2,    -3, &
    -3,     0,     0,     0,  10*0, &
    -4/
    DATA (GY1D(I),I=2171,2250) / &
    &                                     5,     0,     0,     0,  11*0, &
    -2,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1950
    DATA (GY1D(I),I=2251,2395) /0, &
    -30554,  -1341,  1297,   954,  -240,    54,    65,    22,     3, &
    -8,     0,     0,     0,   2*0, &
    -2250,   2998, -1889,   792,   349,    57,   -55,    15,    -7, &
    &                                     4,     0,     0,     0,   3*0, &
    &             1576,  1274,   528,   211,     4,     2,    -4,    -1, &
    -1,     0,     0,     0,   4*0, &
    &                     896,  -408,   -20,  -247,     1,    -1,   -25, &
    &                                    13,     0,     0,     0,   5*0, &
    &                            303,  -147,   -16,   -40,    11,    10, &
    -4,     0,     0,     0,   6*0, &
    -76,    12,    -7,    15,     5, &
    &                                     4,     0,     0,     0,   7*0, &
    -105,     5,   -13,    -5, &
    &                                    12,     0,     0,     0,   8*0, &
    &                                                  19,     5,    -2, &
    &                                     3,     0,     0,     0,   9*0, &
    -1,     3, &
    &                                     2,     0,     0,     0,  10*0, &
    &                                                                 8/
    DATA (GY1D(I),I=2396,2475) / &
    &                                    10,     0,     0,     0,  11*0, &
    &                                     3,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1955
    DATA (GY1D(I),I=2476,2620) /0, &
    -30500,  -1440,  1302,   958,  -229,    47,    65,    11,     4, &
    -3,     0,     0,     0,   2*0, &
    -2215,   3003, -1944,   796,   360,    57,   -56,     9,     9, &
    -5,     0,     0,     0,   3*0, &
    &             1581,  1288,   510,   230,     3,     2,    -6,    -4, &
    -1,     0,     0,     0,   4*0, &
    &                     882,  -397,   -23,  -247,    10,   -14,    -5, &
    &                                     2,     0,     0,     0,   5*0, &
    &                            290,  -152,    -8,   -32,     6,     2, &
    -3,     0,     0,     0,   6*0, &
    -69,     7,   -11,    10,     4, &
    &                                     7,     0,     0,     0,   7*0, &
    -107,     9,    -7,     1, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                  18,     6,     2, &
    -2,     0,     0,     0,   9*0, &
    &                                                          9,     2, &
    &                                     6,     0,     0,     0,  10*0, &
    &                                                                 5/
    DATA (GY1D(I),I=2621,2700) / &
    -2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1960
    DATA (GY1D(I),I=2701,2845) /0, &
    -30421,  -1555,  1302,   957,  -222,    46,    67,    15,     4, &
    &                                     1,     0,     0,     0,   2*0, &
    -2169,   3002, -1992,   800,   362,    58,   -56,     6,     6, &
    -3,     0,     0,     0,   3*0, &
    &             1590,  1289,   504,   242,     1,     5,    -4,     0, &
    &                                     4,     0,     0,     0,   4*0, &
    &                     878,  -394,   -26,  -237,    15,   -11,    -9, &
    &                                     0,     0,     0,     0,   5*0, &
    &                            269,  -156,    -1,   -32,     2,     1, &
    -1,     0,     0,     0,   6*0, &
    -63,    -2,    -7,    10,     4, &
    &                                     4,     0,     0,     0,   7*0, &
    -113,    17,    -5,    -1, &
    &                                     6,     0,     0,     0,   8*0, &
    &                                                   8,    10,    -2, &
    &                                     1,     0,     0,     0,   9*0, &
    &                                                          8,     3, &
    -1,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=2846,2925) / &
    &                                     2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1965
    DATA (GY1D(I),I=2926,3070) /0, &
    -30334,  -1662,  1297,   957,  -219,    45,    75,    13,     8, &
    -2,     0,     0,     0,   2*0, &
    -2119,   2997, -2038,   804,   358,    61,   -57,     5,    10, &
    -3,     0,     0,     0,   3*0, &
    &             1594,  1292,   479,   254,     8,     4,    -4,     2, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     856,  -390,   -31,  -228,    13,   -14,   -13, &
    -5,     0,     0,     0,   5*0, &
    &                            252,  -157,     4,   -26,     0,    10, &
    -2,     0,     0,     0,   6*0, &
    -62,     1,    -6,     8,    -1, &
    &                                     4,     0,     0,     0,   7*0, &
    -111,    13,    -1,    -1, &
    &                                     4,     0,     0,     0,   8*0, &
    &                                                   1,    11,     5, &
    &                                     0,     0,     0,     0,   9*0, &
    &                                                          4,     1, &
    &                                     2,     0,     0,     0,  10*0, &
    -2/
    DATA (GY1D(I),I=3071,3150) / &
    &                                     2,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1970
    DATA (GY1D(I),I=3151,3295) /0, &
    -30220,  -1781,  1287,   952,  -216,    43,    72,    14,     8, &
    -3,     0,     0,     0,   2*0, &
    -2068,   3000, -2091,   800,   359,    64,   -57,     6,    10, &
    -3,     0,     0,     0,   3*0, &
    &             1611,  1278,   461,   262,    15,     1,    -2,     2, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     838,  -395,   -42,  -212,    14,   -13,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            234,  -160,     2,   -22,    -3,    10, &
    -1,     0,     0,     0,   6*0, &
    -56,     3,    -2,     5,    -1, &
    &                                     6,     0,     0,     0,   7*0, &
    -112,    13,     0,     0, &
    &                                     4,     0,     0,     0,   8*0, &
    -2,    11,     3, &
    &                                     1,     0,     0,     0,   9*0, &
    &                                                          3,     1, &
    &                                     0,     0,     0,     0,  10*0, &
    -1/
    DATA (GY1D(I),I=3296,3375) / &
    &                                     3,     0,     0,     0,  11*0, &
    -1,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1975
    DATA (GY1D(I),I=3376,3520) /0, &
    -30100,  -1902,  1276,   946,  -218,    45,    71,    14,     7, &
    -3,     0,     0,     0,   2*0, &
    -2013,   3010, -2144,   791,   356,    66,   -56,     6,    10, &
    -3,     0,     0,     0,   3*0, &
    &             1632,  1260,   438,   264,    28,     1,    -1,     2, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     830,  -405,   -59,  -198,    16,   -12,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            216,  -159,     1,   -14,    -8,    10, &
    -2,     0,     0,     0,   6*0, &
    -49,     6,     0,     4,    -1, &
    &                                     5,     0,     0,     0,   7*0, &
    -111,    12,     0,    -1, &
    &                                     4,     0,     0,     0,   8*0, &
    -5,    10,     4, &
    &                                     1,     0,     0,     0,   9*0, &
    &                                                          1,     1, &
    &                                     0,     0,     0,     0,  10*0, &
    -2/
    DATA (GY1D(I),I=3521,3600) / &
    &                                     3,     0,     0,     0,  11*0, &
    -1,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1980
    DATA (GY1D(I),I=3601,3745) /0, &
    -29992,  -1997,  1281,   938,  -218,    48,    72,    18,     5, &
    -4,     0,     0,     0,   2*0, &
    -1956,   3027, -2180,   782,   357,    66,   -59,     6,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1663,  1251,   398,   261,    42,     2,     0,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     833,  -419,   -74,  -192,    21,   -11,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            199,  -162,     4,   -12,    -7,     9, &
    -2,     0,     0,     0,   6*0, &
    -48,    14,     1,     4,    -3, &
    &                                     5,     0,     0,     0,   7*0, &
    -108,    11,     3,    -1, &
    &                                     3,     0,     0,     0,   8*0, &
    -2,     6,     7, &
    &                                     1,     0,     0,     0,   9*0, &
    -1,     2, &
    &                                     2,     0,     0,     0,  10*0, &
    -5/
    DATA (GY1D(I),I=3746,3825) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1985
    DATA (GY1D(I),I=3826,3970) /0, &
    -29873,  -2072,  1296,   936,  -214,    53,    74,    21,     5, &
    -4,     0,     0,     0,   2*0, &
    -1905,   3044, -2208,   780,   355,    65,   -62,     6,    10, &
    -4,     0,     0,     0,   3*0, &
    &             1687,  1247,   361,   253,    51,     3,     0,     1, &
    &                                     3,     0,     0,     0,   4*0, &
    &                     829,  -424,   -93,  -185,    24,   -11,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            170,  -164,     4,    -6,    -9,     9, &
    -2,     0,     0,     0,   6*0, &
    -46,    16,     4,     4,    -3, &
    &                                     5,     0,     0,     0,   7*0, &
    -102,    10,     4,    -1, &
    &                                     3,     0,     0,     0,   8*0, &
    &                                                   0,     4,     7, &
    &                                     1,     0,     0,     0,   9*0, &
    -4,     1, &
    &                                     2,     0,     0,     0,  10*0, &
    -5/
    DATA (GY1D(I),I=3971,4050) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1990
    DATA (GY1D(I),I=4051,4195) /0, &
    -29775,  -2131,  1314,   939,  -214,    61,    77,    23,     4, &
    -3,     0,     0,     0,   2*0, &
    -1848,   3059, -2239,   780,   353,    65,   -64,     5,     9, &
    -4,     0,     0,     0,   3*0, &
    &             1686,  1248,   325,   245,    59,     2,    -1,     1, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     802,  -423,  -109,  -178,    26,   -10,   -12, &
    -5,     0,     0,     0,   5*0, &
    &                            141,  -165,     3,    -1,   -12,     9, &
    -2,     0,     0,     0,   6*0, &
    -36,    18,     5,     3,    -4, &
    &                                     4,     0,     0,     0,   7*0, &
    -96,     9,     4,    -2, &
    &                                     3,     0,     0,     0,   8*0, &
    &                                                   0,     2,     7, &
    &                                     1,     0,     0,     0,   9*0, &
    -6,     1, &
    &                                     3,     0,     0,     0,  10*0, &
    -6/
    DATA (GY1D(I),I=4196,4275) / &
    &                                     3,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 1995
    DATA (GY1D(I),I=4276,4420) /0, &
    -29692,  -2200,  1335,   940,  -214,    68,    77,    25,     4, &
    -3,     0,     0,     0,   2*0, &
    -1784,   3070, -2267,   780,   352,    67,   -72,     6,     9, &
    -6,     0,     0,     0,   3*0, &
    &             1681,  1249,   290,   235,    68,     1,    -6,     3, &
    &                                     2,     0,     0,     0,   4*0, &
    &                     759,  -418,  -118,  -170,    28,    -9,   -10, &
    -4,     0,     0,     0,   5*0, &
    &                            122,  -166,    -1,     5,   -14,     8, &
    -1,     0,     0,     0,   6*0, &
    -17,    19,     4,     9,    -8, &
    &                                     4,     0,     0,     0,   7*0, &
    -93,     8,     6,    -1, &
    &                                     2,     0,     0,     0,   8*0, &
    -2,    -5,    10, &
    &                                     2,     0,     0,     0,   9*0, &
    -7,    -2, &
    &                                     5,     0,     0,     0,  10*0, &
    -8/
    DATA (GY1D(I),I=4421,4500) / &
    &                                     1,     0,     0,     0,  11*0, &
    &                                     0,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          g(n,m) for 2000
    DATA (GY1D(I),I=4501,4645) /0, &
    -29619.4,-2267.7,1339.6, 932.3,-218.8,  72.3,  79.0,  24.4,   5.0, &
    -2.6,   2.7,  -2.2,  -0.2,   2*0, &
    -1728.2, 3068.4,-2288.0, 786.8, 351.4,  68.2, -74.0,   6.6,  9.4, &
    -6.0,  -1.7,  -0.3,  -0.9,   3*0, &
    &           1670.9,1252.1, 250.0, 222.3,  74.2,   0.0,  -9.2,   3.0, &
    &                                   1.7,  -1.9,   0.2,   0.3,   4*0, &
    &                   714.5,-403.0,-130.4,-160.9,  33.3,  -7.9,  -8.4, &
    -3.1,   1.5,   0.9,   0.1,   5*0, &
    &                          111.3,-168.6,  -5.9,   9.1, -16.6,   6.3, &
    -0.5,  -0.1,  -0.2,  -0.4,   6*0, &
    -12.9,  16.9,   6.9,   9.1,  -8.9, &
    &                                   3.7,   0.1,   0.9,   1.3,   7*0, &
    -90.4,   7.3,   7.0,  -1.5, &
    &                                   1.0,  -0.7,  -0.5,  -0.4,   8*0, &
    -1.2,  -7.9,   9.3, &
    &                                   2.0,   0.7,   0.3,   0.7,   9*0, &
    -7.0,  -4.3, &
    &                                   4.2,   1.7,  -0.3,  -0.4,  10*0, &
    -8.2/
    DATA (GY1D(I),I=4646,4725) / &
    &                                   0.3,   0.1,  -0.4,   0.3,  11*0, &
    -1.1,   1.2,  -0.1,  -0.1,  12*0, &
    &                                          4.0,  -0.2,   0.4,  13*0, &
    -0.4,   0.0,  14*0, &
    &                                                        0.1,  16*0/
!          g(n,m) for 2005
    DATA (GY1D(I),I=4726,4870) /0, &
    -29554.63,-2337.24,1336.30,920.55,-227.00, 73.60,79.88,24.80,5.58, &
    -2.17,  2.95, -2.15, -0.16,   2*0, &
    -1669.05,3047.69,-2305.83,797.96,354.41, 69.56,-74.46,  7.62,9.76, &
    -6.12, -1.60, -0.29, -0.88,   3*0, &
    &          1657.76,1246.39,210.65,208.95, 76.74, -1.65,-11.73, 3.58, &
    &                                  1.42, -1.88,  0.21,  0.30,   4*0, &
    &                672.51,-379.86,-136.54,-151.34, 38.73, -6.88,-6.94, &
    -2.35,  1.44,  0.89,  0.28,   5*0, &
    &                         100.00,-168.05,-14.58, 12.30,-18.11, 5.01, &
    -0.15, -0.31, -0.38, -0.43,   6*0, &
    -13.55, 14.58,  9.37, 10.17,-10.76, &
    &                                  3.06,  0.29,  0.96,  1.18,   7*0, &
    -86.36,  5.42,  9.36, -1.25, &
    &                                  0.29, -0.79, -0.30, -0.37,   8*0, &
    &                                                1.94,-11.25,  8.76, &
    &                                  2.06,  0.53,  0.46,  0.75,   9*0, &
    -4.87, -6.66, &
    &                                  3.77,  1.80, -0.35, -0.26,  10*0, &
    -9.22/
    DATA (GY1D(I),I=4871,4950) / &
    -0.21,  0.16, -0.36,  0.35,  11*0, &
    -2.09,  0.96,  0.08, -0.05,  12*0, &
    &                                         3.99, -0.49,  0.41,  13*0, &
    -0.08, -0.10,  14*0, &
    -0.18,  16*0/
!          g(n,m) for 2010
    DATA (GY1D(I),I=4951,5095) /0, &
    -29496.57,-2396.06,1339.85,912.66,-230.87, 72.78,80.44,24.41,5.50, &
    -1.94,  3.05, -2.12, -0.09,   2*0, &
    -1586.42,3026.34,-2326.54,808.97,357.29, 68.69,-75.00,  8.21,9.45, &
    -6.24, -1.48, -0.21, -0.89,   3*0, &
    &         1668.17,1232.10,166.58,200.26, 75.92, -4.55,-14.50,  3.45, &
    &                                  0.89, -2.03,  0.30,  0.31,   4*0, &
    &               633.73,-356.83,-141.05,-141.40, 45.24, -5.59, -5.27, &
    -1.07,  1.65,  1.04,  0.42,   5*0, &
    &                         89.40,-163.17,-22.83, 14.00,-19.34,  3.13, &
    -0.16, -0.51, -0.63, -0.45,   6*0, &
    -8.03, 13.10, 10.46, 11.61,-12.38, &
    &                                  2.45,  0.54,  0.95,  1.08,   7*0, &
    -78.09,  1.64, 10.85, -0.76, &
    -0.33, -0.79, -0.11, -0.31,   8*0, &
    &                                                4.92,-14.05,  8.43, &
    &                                  2.13,  0.37,  0.52,  0.78,   9*0, &
    -3.54, -8.42, &
    &                                  3.09,  1.79, -0.39, -0.18,  10*0, &
    -10.08/
    DATA (GY1D(I),I=5096,5175) / &
    -1.03,  0.12, -0.37,  0.38,  11*0, &
    -2.80,  0.75,  0.21,  0.02,  12*0, &
    &                                         3.75, -0.77,  0.42,  13*0, &
    &                                                0.04, -0.26,  14*0, &
    -0.26,  16*0/
!          g(n,m) for 2015
    DATA (GY1D(I),I=5176,5320) /0, &
    -29442.0,-2445.1,1350.7, 907.6,-232.6,  70.0,  81.6,  24.2,   5.4, &
    -1.9,   3.1,  -1.9,   0.0,   2*0, &
    -1501.0, 3012.9,-2352.3, 813.7, 360.1,  67.7, -76.1,   8.8,  8.8, &
    -6.3,  -1.5,  -0.2,  -0.9,   3*0, &
    &           1676.7,1225.6, 120.4, 192.4,  72.7,  -6.8, -16.9,   3.1, &
    &                                   0.1,  -2.3,   0.4,   0.4,   4*0, &
    &                   582.0,-334.9,-140.9,-129.9,  51.8,  -3.2,  -3.3, &
    &                                   0.5,   2.0,   1.2,   0.5,   5*0, &
    &                           70.4,-157.5, -28.9,  15.0, -20.6,   0.7, &
    -0.5,  -0.8,  -0.8,  -0.5,   6*0, &
    &                                   4.1,  13.2,   9.4,  13.4, -13.3, &
    &                                   1.8,   0.6,   0.9,   1.0,   7*0, &
    -70.9,  -2.8,  11.7,  -0.1, &
    -0.7,  -0.7,   0.1,  -0.2,   8*0, &
    &                                                 6.8, -15.9,   8.7, &
    &                                   2.1,   0.2,   0.5,   0.8,   9*0, &
    -2.0,  -9.1, &
    &                                   2.4,   1.7,  -0.3,  -0.1,  10*0, &
    -10.5/
    DATA (GY1D(I),I=5321,5400) / &
    -1.8,  -0.2,  -0.4,   0.3,  11*0, &
    -3.6,   0.4,   0.2,   0.1,  12*0, &
    &                                          3.5,  -0.9,   0.5,  13*0, &
    &                                                 0.0,  -0.4,  14*0, &
    -0.3,  16*0/
!          h(n,m) for 1900
    DATA (HY1D(I),I=1,145) /16*0, &
    &     5922,  -1061,  -330,   195,  -210,    -9,   -45,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &             1121,     3,   -69,    53,    83,   -13,   -14,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     523,  -210,   -33,     2,   -10,     7,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -75,  -124,   -35,    -1,   -13,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                     3,    36,    28,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -69,   -12,    16,     8, &
    &                                     0,     0,     0,     0,   8*0, &
    -22,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -18,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=146,225) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1905
    DATA (HY1D(I),I=226,370) /16*0, &
    &     5909,  -1086,  -357,   203,  -193,    -7,   -46,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &             1065,    34,   -77,    56,    86,   -14,   -15,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     480,  -201,   -32,     4,   -11,     7,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -65,  -125,   -32,     0,   -13,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    11,    32,    28,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -67,   -12,    16,     8, &
    &                                     0,     0,     0,     0,   8*0, &
    -22,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -18,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=371,450) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1910
    DATA (HY1D(I),I=451,595) /16*0, &
    &     5898,  -1128,  -389,   211,  -172,    -5,   -47,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &             1000,    62,   -90,    57,    89,   -14,   -15,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     425,  -189,   -33,     5,   -12,     6,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -55,  -126,   -29,     1,   -13,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    21,    28,    28,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -65,   -13,    16,     8, &
    &                                     0,     0,     0,     0,   8*0, &
    -22,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -18,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=596,675) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1915
    DATA (HY1D(I),I=676,820) /16*0, &
    &     5875,  -1191,  -421,   218,  -148,    -2,   -48,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &              917,    84,  -109,    58,    93,   -14,   -15,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     360,  -173,   -34,     8,   -12,     6,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -51,  -126,   -26,     2,   -13,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    32,    23,    28,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -62,   -15,    16,     8, &
    &                                     0,     0,     0,     0,   8*0, &
    -22,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -18,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=821,900) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1920
    DATA (HY1D(I),I=901,1045) /16*0, &
    &     5845,  -1259,  -445,   220,  -122,     0,   -49,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &              823,   103,  -134,    58,    96,   -14,   -15,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     293,  -153,   -38,    11,   -13,     6,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -57,  -125,   -22,     4,   -14,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    43,    18,    28,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -57,   -16,    17,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -22,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -19,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=1046,1125) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1925
    DATA (HY1D(I),I=1126,1270) /16*0, &
    &     5817,  -1334,  -462,   216,   -96,     3,   -50,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &              728,   119,  -163,    58,    99,   -14,   -15,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     229,  -130,   -44,    14,   -14,     6,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -70,  -122,   -18,     5,   -14,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    51,    13,    29,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -52,   -17,    17,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -21,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -19,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=1271,1350) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1930
    DATA (HY1D(I),I=1351,1495) /16*0, &
    &     5808,  -1424,  -480,   205,   -72,     4,   -51,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &              644,   133,  -195,    60,   102,   -15,   -15,    14, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     166,  -109,   -53,    19,   -14,     5,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -90,  -118,   -16,     6,   -14,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    58,     8,    29,     5,    -2, &
    -4,     0,     0,     0,   7*0, &
    -46,   -18,    18,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -20,    -5,    10, &
    -2,     0,     0,     0,   9*0, &
    -19,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=1496,1575) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1935
    DATA (HY1D(I),I=1576,1720) /16*0, &
    &     5812,  -1520,  -494,   188,   -51,     4,   -52,     8,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    &              586,   146,  -226,    64,   104,   -17,   -15,    15, &
    &                                     1,     0,     0,     0,   4*0, &
    &                     101,   -90,   -64,    25,   -14,     5,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -114,  -115,   -15,     7,   -15,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    64,     4,    29,     5,    -3, &
    -4,     0,     0,     0,   7*0, &
    -40,   -19,    18,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -19,    -5,    11, &
    -1,     0,     0,     0,   9*0, &
    -19,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=1721,1800) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1940
    DATA (HY1D(I),I=1801,1945) /16*0, &
    &     5821,  -1614,  -499,   169,   -33,     4,   -52,     8,   -21, &
    &                                     2,     0,     0,     0,   3*0, &
    &              528,   163,  -252,    71,   105,   -18,   -14,    15, &
    &                                     1,     0,     0,     0,   4*0, &
    &                      43,   -72,   -75,    33,   -14,     5,     5, &
    &                                     2,     0,     0,     0,   5*0, &
    -141,  -113,   -15,     7,   -15,    -3, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    69,     0,    29,     5,    -3, &
    -4,     0,     0,     0,   7*0, &
    -33,   -20,    19,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -19,    -5,    11, &
    -1,     0,     0,     0,   9*0, &
    -19,    -2, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=1946,2025) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1945
    DATA (HY1D(I),I=2026,2170) /16*0, &
    &     5810,  -1702,  -499,   144,   -12,     6,   -45,    12,   -27, &
    &                                     5,     0,     0,     0,   3*0, &
    &              477,   186,  -276,    95,   100,   -18,   -21,    17, &
    &                                     1,     0,     0,     0,   4*0, &
    -11,   -55,   -67,    16,     2,   -12,    29, &
    -20,     0,     0,     0,   5*0, &
    -178,  -119,    -9,     6,    -7,    -9, &
    -1,     0,     0,     0,   6*0, &
    &                                    82,   -16,    28,     2,     4, &
    -6,     0,     0,     0,   7*0, &
    -39,   -17,    18,     9, &
    &                                     6,     0,     0,     0,   8*0, &
    -22,     3,     6, &
    -4,     0,     0,     0,   9*0, &
    -11,     1, &
    -2,     0,     0,     0,  10*0, &
    &                                                                 8/
    DATA (HY1D(I),I=2171,2250) / &
    &                                     0,     0,     0,     0,  11*0, &
    -2,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1950
    DATA (HY1D(I),I=2251,2395) /16*0, &
    &     5815,  -1810,  -476,   136,     3,    -1,   -35,     5,   -24, &
    &                                    13,     0,     0,     0,   3*0, &
    &              381,   206,  -278,   103,    99,   -17,   -22,    19, &
    -2,     0,     0,     0,   4*0, &
    -46,   -37,   -87,    33,     0,     0,    12, &
    -10,     0,     0,     0,   5*0, &
    -210,  -122,   -12,    10,   -21,     2, &
    &                                     2,     0,     0,     0,   6*0, &
    &                                    80,   -12,    36,    -8,     2, &
    -3,     0,     0,     0,   7*0, &
    -30,   -18,    17,     8, &
    &                                     6,     0,     0,     0,   8*0, &
    -16,    -4,     8, &
    -3,     0,     0,     0,   9*0, &
    -17,   -11, &
    &                                     6,     0,     0,     0,  10*0, &
    -7/
    DATA (HY1D(I),I=2396,2475) / &
    &                                    11,     0,     0,     0,  11*0, &
    &                                     8,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1955
    DATA (HY1D(I),I=2476,2620) /16*0, &
    &     5820,  -1898,  -462,   133,    15,    -9,   -50,    10,   -11, &
    -4,     0,     0,     0,   3*0, &
    &              291,   216,  -274,   110,    96,   -24,   -15,    12, &
    &                                     0,     0,     0,     0,   4*0, &
    -83,   -23,   -98,    48,    -4,     5,     7, &
    -8,     0,     0,     0,   5*0, &
    -230,  -121,   -16,     8,   -23,     6, &
    -2,     0,     0,     0,   6*0, &
    &                                    78,   -12,    28,     3,    -2, &
    -4,     0,     0,     0,   7*0, &
    -24,   -20,    23,    10, &
    &                                     1,     0,     0,     0,   8*0, &
    -18,    -4,     7, &
    -3,     0,     0,     0,   9*0, &
    -13,    -6, &
    &                                     7,     0,     0,     0,  10*0, &
    &                                                                 5/
    DATA (HY1D(I),I=2621,2700) / &
    -1,     0,     0,     0,  11*0, &
    -3,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1960
    DATA (HY1D(I),I=2701,2845) /16*0, &
    &     5791,  -1967,  -414,   135,    16,   -10,   -55,    11,   -18, &
    &                                     4,     0,     0,     0,   3*0, &
    &              206,   224,  -278,   125,    99,   -28,   -14,    12, &
    &                                     1,     0,     0,     0,   4*0, &
    -130,     3,  -117,    60,    -6,     7,     2, &
    &                                     0,     0,     0,     0,   5*0, &
    -255,  -114,   -20,     7,   -18,     0, &
    &                                     2,     0,     0,     0,   6*0, &
    &                                    81,   -11,    23,     4,    -3, &
    -5,     0,     0,     0,   7*0, &
    -17,   -18,    23,     9, &
    &                                     1,     0,     0,     0,   8*0, &
    -17,     1,     8, &
    -1,     0,     0,     0,   9*0, &
    -20,     0, &
    &                                     6,     0,     0,     0,  10*0, &
    &                                                                 5/
    DATA (HY1D(I),I=2846,2925) / &
    &                                     0,     0,     0,     0,  11*0, &
    -7,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1965
    DATA (HY1D(I),I=2926,3070) /16*0, &
    &     5776,  -2016,  -404,   148,    19,   -11,   -61,     7,   -22, &
    &                                     2,     0,     0,     0,   3*0, &
    &              114,   240,  -269,   128,   100,   -27,   -12,    15, &
    &                                     1,     0,     0,     0,   4*0, &
    -165,    13,  -126,    68,    -2,     9,     7, &
    &                                     2,     0,     0,     0,   5*0, &
    -269,   -97,   -32,     6,   -16,    -4, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    81,    -8,    26,     4,    -5, &
    -4,     0,     0,     0,   7*0, &
    -7,   -23,    24,    10, &
    &                                     0,     0,     0,     0,   8*0, &
    -12,    -3,    10, &
    -2,     0,     0,     0,   9*0, &
    -17,    -4, &
    &                                     3,     0,     0,     0,  10*0, &
    &                                                                 1/
    DATA (HY1D(I),I=3071,3150) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1970
    DATA (HY1D(I),I=3151,3295) /16*0, &
    &     5737,  -2047,  -366,   167,    26,   -12,   -70,     7,   -21, &
    &                                     1,     0,     0,     0,   3*0, &
    &               25,   251,  -266,   139,   100,   -27,   -15,    16, &
    &                                     1,     0,     0,     0,   4*0, &
    -196,    26,  -139,    72,    -4,     6,     6, &
    &                                     3,     0,     0,     0,   5*0, &
    -279,   -91,   -37,     8,   -17,    -4, &
    &                                     4,     0,     0,     0,   6*0, &
    &                                    83,    -6,    23,     6,    -5, &
    -4,     0,     0,     0,   7*0, &
    &                                            1,   -23,    21,    10, &
    &                                     0,     0,     0,     0,   8*0, &
    -11,    -6,    11, &
    -1,     0,     0,     0,   9*0, &
    -16,    -2, &
    &                                     3,     0,     0,     0,  10*0, &
    &                                                                 1/
    DATA (HY1D(I),I=3296,3375) / &
    &                                     1,     0,     0,     0,  11*0, &
    -4,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1975
    DATA (HY1D(I),I=3376,3520) /16*0, &
    &     5675,  -2067,  -333,   191,    31,   -13,   -77,     6,   -21, &
    &                                     1,     0,     0,     0,   3*0, &
    -68,   262,  -265,   148,    99,   -26,   -16,    16, &
    &                                     1,     0,     0,     0,   4*0, &
    -223,    39,  -152,    75,    -5,     4,     7, &
    &                                     3,     0,     0,     0,   5*0, &
    -288,   -83,   -41,    10,   -19,    -4, &
    &                                     4,     0,     0,     0,   6*0, &
    &                                    88,    -4,    22,     6,    -5, &
    -4,     0,     0,     0,   7*0, &
    &                                           11,   -23,    18,    10, &
    -1,     0,     0,     0,   8*0, &
    -12,   -10,    11, &
    -1,     0,     0,     0,   9*0, &
    -17,    -3, &
    &                                     3,     0,     0,     0,  10*0, &
    &                                                                 1/
    DATA (HY1D(I),I=3521,3600) / &
    &                                     1,     0,     0,     0,  11*0, &
    -5,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1980
    DATA (HY1D(I),I=3601,3745) /16*0, &
    &     5604,  -2129,  -336,   212,    46,   -15,   -82,     7,   -21, &
    &                                     1,     0,     0,     0,   3*0, &
    -200,   271,  -257,   150,    93,   -27,   -18,    16, &
    &                                     0,     0,     0,     0,   4*0, &
    -252,    53,  -151,    71,    -5,     4,     9, &
    &                                     3,     0,     0,     0,   5*0, &
    -297,   -78,   -43,    16,   -22,    -5, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    92,    -2,    18,     9,    -6, &
    -4,     0,     0,     0,   7*0, &
    &                                           17,   -23,    16,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -10,   -13,    10, &
    -1,     0,     0,     0,   9*0, &
    -15,    -6, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=3746,3825) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1985
    DATA (HY1D(I),I=3826,3970) /16*0, &
    &     5500,  -2197,  -310,   232,    47,   -16,   -83,     8,   -21, &
    &                                     1,     0,     0,     0,   3*0, &
    -306,   284,  -249,   150,    88,   -27,   -19,    15, &
    &                                     0,     0,     0,     0,   4*0, &
    -297,    69,  -154,    69,    -2,     5,     9, &
    &                                     3,     0,     0,     0,   5*0, &
    -297,   -75,   -48,    20,   -23,    -6, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    95,    -1,    17,    11,    -6, &
    -4,     0,     0,     0,   7*0, &
    &                                           21,   -23,    14,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -7,   -15,     9, &
    -1,     0,     0,     0,   9*0, &
    -11,    -7, &
    &                                     4,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=3971,4050) / &
    &                                     0,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1990
    DATA (HY1D(I),I=4051,4195) /16*0, &
    &     5406,  -2279,  -284,   247,    46,   -16,   -80,    10,   -20, &
    &                                     2,     0,     0,     0,   3*0, &
    -373,   293,  -240,   154,    82,   -26,   -19,    15, &
    &                                     1,     0,     0,     0,   4*0, &
    -352,    84,  -153,    69,     0,     6,    11, &
    &                                     3,     0,     0,     0,   5*0, &
    -299,   -69,   -52,    21,   -22,    -7, &
    &                                     6,     0,     0,     0,   6*0, &
    &                                    97,     1,    17,    12,    -7, &
    -4,     0,     0,     0,   7*0, &
    &                                           24,   -23,    12,     9, &
    &                                     0,     0,     0,     0,   8*0, &
    -4,   -16,     8, &
    -2,     0,     0,     0,   9*0, &
    -10,    -7, &
    &                                     3,     0,     0,     0,  10*0, &
    &                                                                 2/
    DATA (HY1D(I),I=4196,4275) / &
    -1,     0,     0,     0,  11*0, &
    -6,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 1995
    DATA (HY1D(I),I=4276,4420) /16*0, &
    &     5306,  -2366,  -262,   262,    46,   -17,   -69,    11,   -20, &
    &                                     1,     0,     0,     0,   3*0, &
    -413,   302,  -236,   165,    72,   -25,   -21,    15, &
    &                                     0,     0,     0,     0,   4*0, &
    -427,    97,  -143,    67,     4,     8,    12, &
    &                                     4,     0,     0,     0,   5*0, &
    -306,   -55,   -58,    24,   -23,    -6, &
    &                                     5,     0,     0,     0,   6*0, &
    &                                   107,     1,    17,    15,    -8, &
    -5,     0,     0,     0,   7*0, &
    &                                           36,   -24,    11,     8, &
    -1,     0,     0,     0,   8*0, &
    -6,   -16,     5, &
    -2,     0,     0,     0,   9*0, &
    -4,    -8, &
    &                                     1,     0,     0,     0,  10*0, &
    &                                                                 3/
    DATA (HY1D(I),I=4421,4500) / &
    -2,     0,     0,     0,  11*0, &
    -7,     0,     0,     0,  12*0, &
    &                                            0,     0,     0,  13*0, &
    &                                                   0,     0,  14*0, &
    &                                                          0,  16*0/
!          h(n,m) for 2000
    DATA (HY1D(I),I=4501,4645) /16*0, &
    &   5186.1,-2481.6,-227.6, 272.6,  43.8, -17.4, -64.6,  11.9, -19.7, &
    &                                   1.7,   0.1,  -0.4,  -0.9,   3*0, &
    -458.0, 293.4,-231.9, 171.9,  63.7, -24.2, -21.5,  13.4, &
    &                                   0.0,   1.3,   0.3,   0.2,   4*0, &
    -491.1, 119.8,-133.1,  65.1,   6.2,   8.5,  12.5, &
    &                                   4.0,  -0.9,   2.5,   1.8,   5*0, &
    -303.8, -39.3, -61.2,  24.0, -21.5,  -6.2, &
    &                                   4.9,  -2.6,  -2.6,  -0.4,   6*0, &
    &                                 106.3,   0.7,  14.8,  15.5,  -8.4, &
    -5.9,   0.9,   0.7,  -1.0,   7*0, &
    &                                         43.8, -25.4,   8.9,   8.4, &
    -1.2,  -0.7,   0.3,  -0.1,   8*0, &
    -5.8, -14.9,   3.8, &
    -2.9,  -2.8,   0.0,   0.7,   9*0, &
    -2.1,  -8.2, &
    &                                   0.2,  -0.9,   0.0,   0.3,  10*0, &
    &                                                               4.8/
    DATA (HY1D(I),I=4646,4725) / &
    -2.2,  -1.2,   0.3,   0.6,  11*0, &
    -7.4,  -1.9,  -0.9,   0.3,  12*0, &
    -0.9,  -0.4,  -0.2,  13*0, &
    &                                                 0.8,  -0.5,  14*0, &
    -0.9,  16*0/
!          h(n,m) for 2005
    DATA (HY1D(I),I=4726,4870) /16*0, &
    & 5077.99,-2594.50,-198.86,282.07, 42.72,-20.33,-61.14,11.20,-20.11, &
    &                                  2.19,  0.26, -0.55, -0.76,   3*0, &
    -515.43,269.72,-225.23,180.25, 54.75,-22.57,-20.88, 12.69, &
    &                                  0.10,  1.44,  0.23,  0.33,   4*0, &
    -524.72,145.15,-123.45, 63.63,  6.82,  9.83, 12.67, &
    &                                  4.46, -0.77,  2.38,  1.72,   5*0, &
    -305.36,-19.57,-63.53, 25.35,-19.71, -6.72, &
    &                                  4.76, -2.27, -2.63, -0.54,   6*0, &
    &                                103.85,  0.24, 10.93, 16.22, -8.16, &
    -6.58,  0.90,  0.61, -1.07,   7*0, &
    &                                        50.94,-26.32,  7.61,  8.10, &
    -1.01, -0.58,  0.40, -0.04,   8*0, &
    -4.64,-12.76,  2.92, &
    -3.47, -2.69,  0.01,  0.63,   9*0, &
    -0.06, -7.73, &
    -0.86, -1.08,  0.02,  0.21,  10*0, &
    &                                                              6.01/
    DATA (HY1D(I),I=4871,4950) / &
    -2.31, -1.58,  0.28,  0.53,  11*0, &
    -7.93, -1.90, -0.87,  0.38,  12*0, &
    -1.39, -0.34, -0.22,  13*0, &
    &                                                0.88, -0.57,  14*0, &
    -0.82,  16*0/
!          h(n,m) for 2010
    DATA (HY1D(I),I=4951,5095) /16*0, &
    & 4944.26,-2708.54,-160.40,286.48,44.58,-20.90,-57.80, 10.84,-20.54, &
    &                                  2.73,  0.13, -0.87, -0.87,   3*0, &
    -575.73,251.75,-211.03,189.01, 44.18,-21.20,-20.03, 11.51, &
    -0.10,  1.67,  0.27,  0.30,   4*0, &
    -537.03,164.46,-118.06, 61.54,  6.54, 11.83, 12.75, &
    &                                  4.71, -0.66,  2.13,  1.66,   5*0, &
    -309.72, -0.01,-66.26, 24.96,-17.41, -7.14, &
    &                                  4.44, -1.76, -2.49, -0.59,   6*0, &
    &                                101.04,  3.02,  7.03, 16.71, -7.42, &
    -7.22,  0.85,  0.49, -1.14,   7*0, &
    &                                        55.40,-27.61,  6.96,  7.97, &
    -0.96, -0.39,  0.59, -0.07,   8*0, &
    -3.28,-10.74,  2.14, &
    -3.95, -2.51,  0.00,  0.54,   9*0, &
    &                                                       1.64, -6.08, &
    -1.99, -1.27,  0.13,  0.10,  10*0, &
    &                                                              7.01/
    DATA (HY1D(I),I=5096,5175) / &
    -1.97, -2.11,  0.27,  0.49,  11*0, &
    -8.31, -1.94, -0.86,  0.44,  12*0, &
    -1.86, -0.23, -0.25,  13*0, &
    &                                                0.87, -0.53,  14*0, &
    -0.79,  16*0/
!          h(n,m) for 2015
    DATA (HY1D(I),I=5176,5320) /16*0, &
    &   4797.1,-2845.6,-115.3, 283.3,  47.3, -20.8, -54.1,  10.1, -21.6, &
    &                                   3.2,  -0.1,  -1.1,  -0.9,   3*0, &
    -641.9, 244.9,-188.7, 197.0,  33.2, -19.5, -18.3,  10.8, &
    -0.4,   2.0,   0.4,   0.4,   4*0, &
    -538.4, 180.9,-119.3,  58.9,   5.7,  13.3,  11.8, &
    &                                   4.6,  -0.7,   1.9,   1.6,   5*0, &
    -329.5,  16.0, -66.7,  24.4, -14.6,  -6.8, &
    &                                   4.4,  -1.1,  -2.2,  -0.5,   6*0, &
    &                                 100.2,   7.3,   3.4,  16.2,  -6.9, &
    -7.9,   0.8,   0.3,  -1.2,   7*0, &
    &                                         62.6, -27.4,   5.7,   7.8, &
    -0.6,  -0.2,   0.7,  -0.1,   8*0, &
    -2.2,  -9.1,   1.0, &
    -4.2,  -2.2,  -0.1,   0.4,   9*0, &
    &                                                        2.1,  -4.0, &
    -2.8,  -1.4,   0.3,  -0.1,  10*0, &
    &                                                               8.4/
    DATA (HY1D(I),I=5321,5400) / &
    -1.2,  -2.5,   0.2,   0.4,  11*0, &
    -8.7,  -2.0,  -0.9,   0.5,  12*0, &
    -2.4,  -0.1,  -0.3,  13*0, &
    &                                                 0.7,  -0.4,  14*0, &
    -0.8,  16*0/
!          Secular variation rates are nominally okay through 2020
    DATA (GT1D(I),I=1,145) /0, &
    &     10.3,   -8.7,   3.4,  -0.7,  -0.2,  -0.3,   0.3,   0.2,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   2*0, &
    &     18.1,   -3.3,  -5.5,   0.2,   0.5,  -0.1,  -0.2,   0.0,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   3*0, &
    &              2.1,  -0.7,  -9.1,  -1.3,  -0.7,  -0.5,  -0.6,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   4*0, &
    -10.1,   4.1,  -0.1,   2.1,   1.3,   0.5,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   5*0, &
    -4.3,   1.4,  -1.2,   0.1,  -0.2,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   6*0, &
    &                                   3.9,   0.3,  -0.6,   0.4,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   7*0, &
    &                                          1.6,  -0.8,   0.1,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   8*0, &
    &                                                 0.2,  -0.4,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   9*0, &
    &                                                        0.3,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,  10*0, &
    &                                                               0.0/
    DATA (GT1D(I),I=146,225) / &
    &                                   0.0,   0.0,   0.0,   0.0,  11*0, &
    &                                   0.0,   0.0,   0.0,   0.0,  12*0, &
    &                                          0.0,   0.0,   0.0,  13*0, &
    &                                                 0.0,   0.0,  14*0, &
    &                                                        0.0,  16*0/
    DATA (HT1D(I),I=1,145) /16*0, &
    -26.6,  -27.4,   8.2,  -1.3,   0.6,   0.0,   0.8,  -0.3,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   3*0, &
    -14.1,  -0.4,   5.3,   1.7,  -2.1,   0.4,   0.3,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   4*0, &
    &                     1.8,   2.9,  -1.2,  -0.7,  -0.2,   0.1,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   5*0, &
    -5.2,   3.4,   0.2,  -0.3,   0.5,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   6*0, &
    &                                   0.0,   0.9,  -0.6,  -0.2,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   7*0, &
    &                                          1.0,   0.1,  -0.3,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   8*0, &
    -0.2,   0.3,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,   9*0, &
    &                                                        0.0,   0.0, &
    &                                   0.0,   0.0,   0.0,   0.0,  10*0, &
    &                                                               0.0/
    DATA (HT1D(I),I=146,225) / &
    &                                   0.0,   0.0,   0.0,   0.0,  11*0, &
    &                                   0.0,   0.0,   0.0,   0.0,  12*0, &
    &                                          0.0,   0.0,   0.0,  13*0, &
    &                                                 0.0,   0.0,  14*0, &
    &                                                        0.0,  16*0/

!          Do not need to load new coefficients if date has not changed
    ICHG = 0
    IF (DATE == DATEL) GO TO 300
    DATEL = DATE
    ICHG = 1

!          Trap out of range date:
    IF (DATE < EPOCH(1)) GO TO 9100
    IF (DATE > EPOCH(NEPO)+5.) WRITE(0,9200) DATE, EPOCH(NEPO) + 5.
     
    DO 100 I=1,NEPO
        IF (DATE < EPOCH(I)) GO TO 110
        IY = I
    100 END DO
    110 CONTINUE
     
    NMAX  = NMXE(IY)
    TIME  = DATE
    T     = TIME-EPOCH(IY)
    TO5   = T/5.
    IY1   = IY + 1
    GB(1) = 0.0
    GV(1) = 0.0
    I  = 2
    F0 = -1.0D-5
    DO 200 N=1,NMAX
        F0 = F0 * REAL(N)/2.
        F  = F0 / SQRT(2.0)
        NN = N+1
        MM = 1
        IF (IY < NEPO) GB(I) = (GYR(NN,MM,IY) +                              & ! interpolate (m=0 terms)
        (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F0
        IF (IY == NEPO) GB(I) = (GYR(NN,MM,IY) + GT(NN,MM)    *T  ) * F0     ! extrapolate (m=0 terms)
        GV(I) = GB(I) / REAL(NN)
        I = I+1
        DO 200 M=1,N
            F  = F / SQRT( REAL(N-M+1) / REAL(N+M) )
            NN = N+1
            MM = M+1
            I1 = I+1
            IF (IY < NEPO) THEN                                                ! interpolate (m>0 terms)
                GB(I)  = (GYR(NN,MM,IY) + &
                (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F
                GB(I1) = (HYR(NN,MM,IY) + &
                (HYR(NN,MM,IY1)-HYR(NN,MM,IY))*TO5) * F
            ELSE                                                                  ! extrapolate (m>0 terms)
                GB(I)  = (GYR(NN,MM,IY) +GT (NN,MM)    *T  ) * F
                GB(I1) = (HYR(NN,MM,IY) +HT (NN,MM)    *T  ) * F
            ENDIF
            RNN = REAL(NN)
            GV(I)  = GB(I)  / RNN
            GV(I1) = GB(I1) / RNN
            I = I+2
    200 END DO
     
    300 CONTINUE

    RETURN
     
!          Error trap diagnostics:
    9100 WRITE (0,'(''COFRM:  DATE'',F9.3,'' precedes earliest available ('' &
    ,F6.1,'')'')') DATE, EPOCH(1)
    ERROR STOP 1
    9200 FORMAT('COFRM:  DATE',F9.3,' is after the last recommended for ext &
    rapolation (',F6.1,')')

    END SUBROUTINE COFRM
     
    SUBROUTINE DYPOL (COLAT,ELON,VP)
!          Computes parameters for dipole component of geomagnetic field.
!          COFRM must be called before calling DYPOL!
!          940504 A. D. Richmond

!          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed

!          RETURNS:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole
!                    (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude, in T.m, of dipole component of magnetic
!                    potential at geomagnetic pole and geocentric radius
!                    of 6371.2 km
     
    use magcof_mod
    implicit none

    ! Arguments
    real, intent(out) :: COLAT, ELON, VP

    ! Local variables (previously implicit)
    real :: GPL, CTP, STP

    real,PARAMETER :: RTOD = 57.2957795130823, RE = 6371.2
     
!          Compute geographic colatitude and longitude of the north pole of
!          earth centered dipole
    GPL   = SQRT (GB(2)**2 + GB(3)**2 + GB(4)**2)
    CTP   = GB(2) / GPL
    STP   = SQRT (1. - CTP*CTP)
    COLAT = ACOS (CTP) * RTOD
    ELON  = ATAN2 (GB(4),GB(3)) * RTOD
     
!          Compute magnitude of magnetic potential at pole, radius Re.
    VP = .2*GPL*RE
!          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).
     
    RETURN
    END SUBROUTINE DYPOL
     
    SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
!          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
!          COFRM must be called to establish coefficients for correct date
!          prior to calling FELDG.

!          IENTY is an input flag controlling the meaning and direction of the
!                remaining formal arguments:
!          IENTY = 1
!            INPUTS:
!              GLAT = Latitude of point (deg)
!              GLON = Longitude (east=+) of point (deg)
!              ALT  = Ht of point (km)
!            RETURNS:
!              BNRTH  north component of field vector (Gauss)
!              BEAST  east component of field vector  (Gauss)
!              BDOWN  downward component of field vector (Gauss)
!              BABS   magnitude of field vector (Gauss)

!          IENTY = 2
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.2 km )
!              GLON = Y coordinate (in units of earth radii 6371.2 km )
!              ALT  = Z coordinate (in units of earth radii 6371.2 km )
!            RETURNS:
!              BNRTH = X component of field vector (Gauss)
!              BEAST = Y component of field vector (Gauss)
!              BDOWN = Z component of field vector (Gauss)
!              BABS  = Magnitude of field vector (Gauss)
!          IENTY = 3
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.2 km )
!              GLON = Y coordinate (in units of earth radii 6371.2 km )
!              ALT  = Z coordinate (in units of earth radii 6371.2 km )
!            RETURNS:
!              BNRTH = Dummy variable
!              BEAST = Dummy variable
!              BDOWN = Dummy variable
!              BABS  = Magnetic potential (T.m)

!          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed

!          HISTORY:
!          Apr 1983: written by Vincent B. Wickwar (Utah State Univ.).

!          May 1994 (A.D. Richmond): Added magnetic potential calculation

!          Oct 1995 (Barnes): Added ICHG
     
    use magcof_mod
    implicit none

    ! Arguments
    integer, intent(in)  :: IENTY
    real,    intent(in)  :: GLAT, GLON, ALT
    real,    intent(out) :: BNRTH, BEAST, BDOWN, BABS

    ! Local variables (previously implicit)
    integer :: IS, IHMAX, LAST, IMAX, MK, K, I, IH, IL, M, IHM, ILM
    real    :: RLAT, CT, ST, RLON, CP, SP
    real    :: XXX, YYY, ZZZ, RQ, S, T, F, X, Y, Z
    real    :: BXXX, BYYY, BZZZ, BRHO
    integer, save :: IENTYP = -10000

    real,PARAMETER :: DTOR = 0.01745329251994330, RE = 6371.2
    real :: G(255), H(255), XI(3)
    SAVE G
     
    IF (IENTY == 1) THEN
        IS   = 1
        RLAT = GLAT * DTOR
        CT   = SIN (RLAT)
        ST   = COS (RLAT)
        RLON = GLON * DTOR
        CP   = COS (RLON)
        SP   = SIN (RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
    ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
    ENDIF
    RQ    = 1./(XXX**2+YYY**2+ZZZ**2)
    XI(1) = XXX*RQ
    XI(2) = YYY*RQ
    XI(3) = ZZZ*RQ
    IHMAX = NMAX*NMAX+1
    LAST  = IHMAX+NMAX+NMAX
    IMAX  = NMAX+NMAX-1
     
    IF (IENTY /= IENTYP .OR. ICHG == 1) THEN
        IENTYP = IENTY
        ICHG = 0
        IF (IENTY /= 3) THEN
            DO 10 I=1,LAST
                G(I) = GB(I)
            10 END DO
        ELSE
            DO 20 I=1,LAST
                G(I) = GV(I)
            20 END DO
        ENDIF
    ENDIF
     
    DO 30 I=IHMAX,LAST
        H(I) = G(I)
    30 END DO

    MK = 3
    IF (IMAX == 1) MK=1

    DO 100 K=1,MK,2
        I  = IMAX
        IH = IHMAX

        60 IL = IH-I
        F = 2./FLOAT(I-K+2)
        X = XI(1)*F
        Y = XI(2)*F
        Z = XI(3)*(F+F)

        I = I-2
        IF (I < 1) GO TO 90
        IF (I == 1) GO TO 80

        DO 70 M=3,I,2
            IHM = IH+M
            ILM = IL+M
            H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1)) &
            -Y*(H(IHM+2)+H(IHM-2))
            H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2)) &
            +Y*(H(IHM+3)+H(IHM-1))
        70 END DO

        80 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
        H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))

        90 H(IL)   = G(IL)   + Z*H(IH)   + 2.*(X*H(IH+1)+Y*H(IH+2))
        IH = IL
        IF (I >= K) GO TO 60
    100 END DO
     
    S = .5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
    T = (RQ+RQ)*SQRT(RQ)
    BXXX = T*(H(3)-S*XXX)
    BYYY = T*(H(4)-S*YYY)
    BZZZ = T*(H(2)-S*ZZZ)
    BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
    IF (IS == 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
    ELSEIF (IS == 2) THEN        ! (leave in earth centered cartesian)
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
    ENDIF
     
!          Magnetic potential computation makes use of the fact that the
!          calculation of V is identical to that for r*Br, if coefficients
!          in the latter calculation have been divided by (n+1) (coefficients
!          GV).  Factor .1 converts km to m and gauss to tesla.
    IF (IENTY == 3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1
     
    RETURN
    END SUBROUTINE FELDG
     
    SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
!          Convert geodetic to cartesian coordinates by calling CONVRT
!          940503 A. D. Richmond
    implicit none

    ! Arguments
    real, intent(in)  :: GDLAT, GLON, ALT
    real, intent(out) :: X, Y, Z

    ! Local variables (previously implicit)
    real :: RHO, ANG

    real,PARAMETER :: DTOR = 0.01745329251994330
    CALL CONVRT (1,GDLAT,ALT,RHO,Z)
    ANG = GLON*DTOR
    X = RHO*COS(ANG)
    Y = RHO*SIN(ANG)
    RETURN
    END SUBROUTINE GD2CART
     
    SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
!          Convert space point from geodetic to geocentric or vice versa.
!          [See header comments above for full description of I, inputs, outputs]

    implicit none

    ! Arguments
    integer, intent(in)    :: I
    real,    intent(inout) :: GDLAT, ALT
    real,    intent(inout) :: X1, X2

    ! Local variables (previously implicit)
    real :: SINLAT, COSLAT, D, Z, RHO, RKM, GCLAT, SCL
    real :: RI, A2, A4, A6, A8
    real :: CCL, S2CL, C2CL, S4CL, C4CL, S8CL, S6CL, DLTCL, SGL

!          I = 1  (convert from geodetic to cylindrical geocentric)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance above (north of) Earth's equatorial plane (km)

!          I = 2  (convert from geodetic to spherical geocentric)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)

!          I = 3  (convert from cylindrical geocentric to geodetic)
!            INPUTS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance from Earth's equatorial plane (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)

!          I = 4  (convert from spherical geocentric to geodetic)
!            INPUTS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)


!          HISTORY:
!          940503 (A. D. Richmond):  Based on a routine originally written
!          by V. B. Wickwar.

!          Mar 2004: (Barnes) Revise spheroid definition to WGS-1984 to conform
!          with IGRF-9 release (EOS Volume 84 Number 46 November 18 2003).

!          REFERENCE: ASTRON. J. VOL. 66, p. 15-16, 1961
     
!          E2  = square of eccentricity of ellipse
!          REP = earth's polar      radius (km)
!          REQ = earth's equatorial radius (km)
    real,PARAMETER :: RTOD = 57.2957795130823, DTOR = 0.01745329251994330, &
    REP  = 6356.752, REQ = 6378.137, E2 = 1.-(REP/REQ)**2, &
    E4 = E2*E2, E6 = E4*E2, E8 = E4*E4, OME2REQ = (1.-E2)*REQ, &
    A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. , &
    A22 =     (                        E6 +     E8)/  32. , &
    A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. , &
    A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. , &
    A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. , &
    A43 =                                   15.*E8 / 256. , &
    A44 =                                      -E8 /  16. , &
    A61 =  3.*(                     4.*E6 +  5.*E8)/1024. , &
    A62 = -3.*(                        E6 +     E8)/  32. , &
    A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. , &
    A81 =                                   -5.*E8 /2048. , &
    A82 =                                   64.*E8 /2048. , &
    A83 =                                 -252.*E8 /2048. , &
    A84 =                                  320.*E8 /2048. 
     
    IF (I >= 3) GO TO 300
     
!          Geodetic to geocentric
     
!          Compute RHO,Z
    SINLAT = SIN(GDLAT*DTOR)
    COSLAT = SQRT(1.-SINLAT*SINLAT)
    D      = SQRT(1.-E2*SINLAT*SINLAT)
    Z      = (ALT+OME2REQ/D)*SINLAT
    RHO    = (ALT+REQ/D)*COSLAT
    X1 = RHO
    X2 = Z
    IF (I == 1) RETURN
     
!          Compute GCLAT,RKM
    RKM   = SQRT(Z*Z + RHO*RHO)
    GCLAT = RTOD*ATAN2(Z,RHO)
    X1 = GCLAT
    X2 = RKM
    RETURN
     
!          Geocentric to geodetic
    300 IF (I == 3) THEN
        RHO = X1
        Z = X2
        RKM = SQRT(Z*Z+RHO*RHO)
        SCL = Z/RKM
        GCLAT = ASIN(SCL)*RTOD
    ELSEIF (I == 4) THEN
        GCLAT = X1
        RKM = X2
        SCL = SIN(GCLAT*DTOR)
    ELSE
        RETURN
    ENDIF
     
    RI = REQ/RKM
    A2 = RI*(A21+RI*(A22+RI* A23))
    A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
    A6 = RI*(A61+RI*(A62+RI* A63))
    A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
    CCL = SQRT(1.-SCL*SCL)
    S2CL = 2.*SCL*CCL
    C2CL = 2.*CCL*CCL-1.
    S4CL = 2.*S2CL*C2CL
    C4CL = 2.*C2CL*C2CL-1.
    S8CL = 2.*S4CL*C4CL
    S6CL = S2CL*C4CL+C2CL*S4CL
    DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
    GDLAT = DLTCL*RTOD+GCLAT
    SGL = SIN(GDLAT*DTOR)
    ALT = RKM*COS(DLTCL)-REQ*SQRT(1.-E2*SGL*SGL)
    RETURN
    END SUBROUTINE CONVRT
