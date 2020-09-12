//
// Created by cc on 7/10/20.
//

#ifndef PPPLIB_CONSTANT_H
#define PPPLIB_CONSTANT_H

namespace PPPLib {

    #define SQR(x)   ((x)*(x))
    #define SQRT(x)  ((x)<0.0||(x)!=(x)?0.0:sqrt(x))
    #define MIN(x,y) ((x)<=(y)?(x):(y))
    #define MAX(x,y) ((x)>=(y)?(x):(y))

    // common constant definitions
    const static int MAX_BUFF=4096;
    const static double PI=3.14159265358979;
    const static double CLIGHT=299792458.0;
    const static double D2R=(PI/180.0);
    const static double R2D=(180.0/PI);
    const static double AS2R=(D2R/3600);
    const static double MG2MS2=9.80665E-6;
    const static double AU=149597870691.0;
    const static double DTTOL=0.005;

    enum TIME_TYPE {
        TIME_TYPE_NONE = 0,
        TIME_utc,
        TIME_GPS,
        TIME_BDS,
        TIME_GAL,
        TIME_GLO,
        TIME_QZS
    };

    enum TIME_FORMAT {
        TIME_FMT_STR,
        TIME_FMT_WS,
    };

    enum SAT_SYS_TYPE {
        SYS_NONE=0x00,
        SYS_GPS=0x01,
        SYS_BDS=0x02,
        SYS_GAL=0x04,
        SYS_GLO=0x08,
        SYS_QZS=0x10,
        SYS_IRN=0x20,
        SYS_SBS=0x40,
        SYS_ALL=0xFF
    };

    const static double kSysTimeStart[] = {1970, 1, 1, 0, 0, 0};
    const static double kGpsTimeStart[] =  {1980, 1, 6, 0, 0, 0};
    const static double kGalTimeStart[] =  {1999, 8, 22, 0, 0, 0};
    const static double kBdsTimeStart[] =  {2006, 1, 1, 0, 0, 0};
    const static double kUtcLeapSec[65][7]={
            {2017,1,1,0,0,0,-18}, {2015,7,1,0,0,0,-17}, {2012,7,1,0,0,0,-16},
            {2009,1,1,0,0,0,-15}, {2006,1,1,0,0,0,-14}, {1999,1,1,0,0,0,-13},
            {1997,7,1,0,0,0,-12}, {1996,1,1,0,0,0,-11}, {1994,7,1,0,0,0,-10},
            {1993,7,1,0,0,0, -9}, {1992,7,1,0,0,0, -8}, {1991,1,1,0,0,0, -7},
            {1990,1,1,0,0,0, -6}, {1988,1,1,0,0,0, -5}, {1985,7,1,0,0,0, -4},
            {1983,7,1,0,0,0, -3}, {1982,7,1,0,0,0, -2}, {1981,7,1,0,0,0, -1},
            {0}
    };

    const double WGS84_EARTH_LONG_RADIUS=6378137.0;
    const double WGS84_EARTH_SHORT_RADIUS=6356752.3142;
    const double WGS84_EARTH_OBLATEO=1.0/298.257223563;
    const double WGS84_FIRST_E2=0.00669437999013;
    const double WGS84_SECOND_E2=0.006739496742227;
    const double WGS84_GM=3.986005e14;

    const int COORDINATE_ACCURACY=4;  // 空间直角坐标系以及大地坐标系高程方向精确度
    const int MSEC_ACCURACY=3;        // 秒的精确度(ms)
    const int MATRIX_ACCURACY=6;      // 矩阵中double数据与角分秒格式中的秒的精确度
    const int DEGREE_ACCURACY=11;     // 大地坐标中经纬度与小数度的精确度(小数后11位)
    const double LAT_ACCURACY=1.0e-11;// 计算大地纬度B的精度

    enum COORD_FRAME_TYPE {
        WGS84=1,
        CGCS2000,
        ITRF96,
        PZ90
    };

    enum COORDINATE_TYPE {
        COORD_BLH=1,
        COORD_XYZ,
        COORD_ENU,
        COORD_NED
    };

    // gnss related constant definitions
//    const int GNSS_NUM_FREQ=5;
    const int GNSS_NUM_EXOBS=3;
    const int MAX_GNSS_OBS_TYPE=36;
    const int MAX_GNSS_FRQ_NUM=6;
    const int MAX_GNSS_USED_FRQ_NUM=3;
    static const string kSatSysCode="GCERJI";
    static const string kGnssObsCode="CLDS";

    enum GPS{
        GPS_MIN_PRN=1,
        GPS_MAX_PRN=32,
        NUM_GPS_SAT=(GPS_MAX_PRN-GPS_MIN_PRN+1),
        NSYS_GPS=1,
        SYS_INDEX_GPS=0,
    };

    enum BDS{
        BDS_MIN_PRN=1,
        BDS_MAX_PRN=45,
        NUM_BDS_SAT=(BDS_MAX_PRN-BDS_MIN_PRN+1),
        NSYS_BDS=1,
        SYS_INDEX_BDS=1
    };

    enum GAL{
        GAL_MIN_PRN=1,
        GAL_MAX_PRN=36,
        NUM_GAL_SAT=(GAL_MAX_PRN-GAL_MIN_PRN+1),
        NSYS_GAL=1,
        SYS_INDEX_GAL=2
    };

    enum GLO{
        GLO_MIN_PRN=1,
        GLO_MAX_PRN=27,
        NUM_GLO_SAT=(GLO_MAX_PRN-GLO_MIN_PRN+1),
        NSYS_GLO=1,
        SYS_INDEX_GLO=3
    };

    enum QZS{
        QZS_MIN_PRN=193,
        QZS_MAX_PRN=201,
        NUM_QZS_SAT=(QZS_MAX_PRN-QZS_MIN_PRN+1),
        NSYS_QZS=1,
        SYS_INDEX_QZS=4
    };

    enum IRN{
        IRN_MIN_PRN=0,
        IRN_MAX_PRN=10,
        NUM_IRN_SAT=(IRN_MAX_PRN-IRN_MIN_PRN+1),
        NSYS_IRN=1,
        SYS_INDEX_IRN=5
    };
    const int NSYS=NSYS_GPS+NSYS_BDS+NSYS_GAL+NSYS_GLO+NSYS_QZS;
    const int MAX_SAT_NUM=NUM_GPS_SAT+NUM_BDS_SAT+NUM_GAL_SAT+NUM_GLO_SAT+NUM_QZS_SAT;

    enum RECEIVER_INDEX {
        REC_ROVER=0,
        REC_BASE=1
    };

    constexpr double MHZ_TO_HZ=1000000.0;
    constexpr double FREQ_NONE=0.0;
    constexpr double FREQ_GPS_L1=1575.42*MHZ_TO_HZ;
    constexpr double FREQ_GPS_L2=1227.60*MHZ_TO_HZ;
    constexpr double FREQ_GPS_L5=1176.45*MHZ_TO_HZ;
    constexpr double FREQ_BDS_B1=1561.098*MHZ_TO_HZ;
    constexpr double FREQ_BDS_B2=1207.140*MHZ_TO_HZ;
    constexpr double FREQ_BDS_B3=1268.52*MHZ_TO_HZ;
    constexpr double FREQ_BDS_B1C=1575.42*MHZ_TO_HZ;
    constexpr double FREQ_BDS_B2A=1176.45*MHZ_TO_HZ;
    constexpr double FREQ_BDS_B2B=1207.140*MHZ_TO_HZ;
    constexpr double FREQ_GAL_E1=1575.42*MHZ_TO_HZ;
    constexpr double FREQ_GAL_E5A=1176.45*MHZ_TO_HZ;
    constexpr double FREQ_GAL_E5B=1207.140*MHZ_TO_HZ;
    constexpr double FREQ_GAL_E5AB=1191.795*MHZ_TO_HZ;
    constexpr double FREQ_GAL_E6=1278.75*MHZ_TO_HZ;
    constexpr double FREQ_GLO_G1=1602.00*MHZ_TO_HZ;
    constexpr double FREQ_GLO_G2=1246.00*MHZ_TO_HZ;
    constexpr double FREQ_GLO_D1=0.5625*MHZ_TO_HZ;
    constexpr double FREQ_GLO_D2=0.4375*MHZ_TO_HZ;
    constexpr double FREQ_QZS_L1=1575.42*MHZ_TO_HZ;
    constexpr double FREQ_QZS_L2=1227.60*MHZ_TO_HZ;
    constexpr double FREQ_QZS_L5=1176.45*MHZ_TO_HZ;
    constexpr double FREQ_QZS_LEX=1278.75*MHZ_TO_HZ;
    constexpr int GLO_MIN_FREQ=-7;
    constexpr int GLO_MAX_FREQ=13;

    const int GNSS_CODE_NONE=0;
    const int GNSS_CODE_L1C=1;
    const int GNSS_CODE_L1S=2;
    const int GNSS_CODE_L1L=3;
    const int GNSS_CODE_L1X=4;
    const int GNSS_CODE_L1P=5;
    const int GNSS_CODE_L1W=6;
    const int GNSS_CODE_L1Y=7;
    const int GNSS_CODE_L1M=8;
    const int GNSS_CODE_L1A=9;
    const int GNSS_CODE_L1B=10;
    const int GNSS_CODE_L1Z=11;
    const int GNSS_CODE_L1D=12;
    const int GNSS_CODE_L2C=13;
    const int GNSS_CODE_L2D=14;
    const int GNSS_CODE_L2S=15;
    const int GNSS_CODE_L2L=16;
    const int GNSS_CODE_L2X=17;
    const int GNSS_CODE_L2P=18;
    const int GNSS_CODE_L2W=19;
    const int GNSS_CODE_L2Y=20;
    const int GNSS_CODE_L2M=21;
    const int GNSS_CODE_L2I=22;
    const int GNSS_CODE_L2Q=23;
    const int GNSS_CODE_L3I=24;
    const int GNSS_CODE_L3Q=25;
    const int GNSS_CODE_L3X=26;
    const int GNSS_CODE_L4A=27;
    const int GNSS_CODE_L4B=28;
    const int GNSS_CODE_L4X=29;
    const int GNSS_CODE_L5I=30;
    const int GNSS_CODE_L5Q=31;
    const int GNSS_CODE_L5X=32;
    const int GNSS_CODE_L5D=33;
    const int GNSS_CODE_L5P=34;
    const int GNSS_CODE_L5Z=35;
    const int GNSS_CODE_L6A=36;
    const int GNSS_CODE_L6B=37;
    const int GNSS_CODE_L6X=38;
    const int GNSS_CODE_L6C=39;
    const int GNSS_CODE_L6Z=40;
    const int GNSS_CODE_L6S=41;
    const int GNSS_CODE_L6L=42;
    const int GNSS_CODE_L6E=43;
    const int GNSS_CODE_L6I=44;
    const int GNSS_CODE_L6Q=45;
    const int GNSS_CODE_L7I=46;
    const int GNSS_CODE_L7Q=47;
    const int GNSS_CODE_L7X=48;
    const int GNSS_CODE_L7D=49;
    const int GNSS_CODE_L7P=50;
    const int GNSS_CODE_L7Z=51;
    const int GNSS_CODE_L8I=52;
    const int GNSS_CODE_L8Q=53;
    const int GNSS_CODE_L8X=54;
    const int GNSS_CODE_L8D=55;
    const int GNSS_CODE_L8P=56;
    const int GNSS_CODE_L8Z=57;
    const int MAX_GNSS_CODE_TYPE=57;

    const int GPS_C1CC2W=0;
    const int GPS_C1CC5Q=1;
    const int GPS_C1CC5X=2;
    const int GPS_C1WC2W=3;
    const int GPS_C1CC1W=4;
    const int GPS_C2CC2W=5;
    const int GPS_C2WC2S=6;
    const int GPS_C2WC2L=7;
    const int GPS_C2WC2X=8;

    const int BD2_C2IC7I=0;
    const int BD2_C2IC6I=1;
    const int BD3_C1XC5X=2;
    const int BD3_C1PC5P=3;
    const int BD3_C1DC5D=4;
    const int BD3_C1XC6I=5;
    const int BD3_C1PC6I=6;
    const int BD3_C1DC6I=7;
    const int BD3_C2IC6I=8;
    const int BD3_C1XC7Z=9;
    const int BD3_C1XC8X=10;

    const int GAL_C1CC5Q=0;
    const int GAL_C1CC6C=1;
    const int GAL_C1CC7Q=2;
    const int GAL_C1CC8Q=3;
    const int GAL_C1XC5X=4;
    const int GAL_C1XC7X=5;
    const int GAL_C1XC8X=6;

    const int GLO_C1CC2C=0;
    const int GLO_C1CC2P=1;
    const int GLO_C1PC2P=2;
    const int GLO_C1CC1P=3;
    const int GLO_C2CC2P=4;

    const int QZS_C1CC2L=0;
    const int QZS_C1CC5X=1;
    const int QZS_C1CC5Q=2;
    const int QZS_C1XC2X=3;
    const int QZS_C1XC5X=4;
    const int QZS_C1CC1X=5;

    const double kGnssFreqs[NSYS][MAX_GNSS_FRQ_NUM]{
            {FREQ_GPS_L1,FREQ_GPS_L2, FREQ_GPS_L5,  FREQ_NONE,    FREQ_NONE,    FREQ_NONE},
            {FREQ_BDS_B1,FREQ_BDS_B2, FREQ_BDS_B3,  FREQ_BDS_B1C, FREQ_BDS_B2A, FREQ_BDS_B2B},
            {FREQ_GAL_E1,FREQ_GAL_E5A,FREQ_GAL_E5B, FREQ_NONE,    FREQ_NONE,    FREQ_NONE},
            {FREQ_GLO_G1,FREQ_GLO_G2, FREQ_GLO_D1,  FREQ_GLO_D2,  FREQ_NONE,    FREQ_NONE},
            {FREQ_QZS_L1,FREQ_QZS_L2, FREQ_QZS_L5,  FREQ_NONE,    FREQ_NONE,    FREQ_NONE},
    };

    const double MAX_DTOE_GPS=7200.0;
    const double MAX_DTOE_BDS=21600.0;
    const double MAX_DTOE_GAL=10800.0;
    const double MAX_DTOE_GLO=1800.0;
    const double MAX_DTOE_QZS=7200.0;
    const double OMGE_GPS=7.2921151467E-5;
    const double OMGE_BDS=7.292115E-5;
    const double OMGE_GAL=7.2921151467E-5;
    const double OMGE_GLO=7.292115E-5;
    const double MU_GPS=3.9860050E14;
    const double MU_BDS=3.986004418E14;
    const double MU_GAL=3.986004418E14;
    const double MU_GLO=3.9860044E14;

    enum GPS_FRQ {
        GPS_L1,
        GPS_L2,
        GPS_L5
    };

    enum BDS_FRQ {
        BDS_B1I,
        BDS_B2I,
        BDS_B3I,
        BDS_B1C,
        BDS_B2a,
        BDS_B2b
    };

    enum GAL_FRQ {
        GAL_E1,
        GAL_E5a,
        GAL_E5b,
    };

    enum GLO_FRQ {
        GLO_G1,
        GLO_G2
    };

    enum QZS_FRQ {
        QZS_L1,
        QZS_L2,
        QZS_L5
    };

    const string kGnssSignalCodes[]={
            "",  "1C","1S","1L","1X","1P","1W","1Y","1M","1A",
            "1B","1Z","1D","2C","2D","2S","2L","2X","2P","2W",
            "2Y","2M","2I","2Q","3I","3Q","3X","4A","4B","4X",
            "5I","5Q","5X","5D","5P","5Z","6A","6B","6X","6C",
            "6Z","6S","6L","6E","6I","6Q","7I","7Q","7X","7D",
            "7P","7Z","8I","8Q","8X","8D","8P", "8Z" , "" , ""
    };

    // GPS L1(1),L2(2),L5(5)
    const unsigned char kGpsFreqBand[]={
            0, 1, 1, 1, 1, 1, 1, 1, 1, 0,
            0, 0, 0, 2, 2, 2, 2, 2, 2, 2,
            2, 2, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 3, 3, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    //BDS B1(2),B3(6),B2(7),B1C(1),B2a(5),B2b(7)
    const unsigned char kBdsFreqBand[]={
            0, 0, 0, 0, 4, 4, 0, 0, 0, 4,
            0, 0, 4, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 5, 5, 5, 0, 3, 0, 3, 0,
            0, 0, 0, 0, 3, 3, 2, 2, 2, 6,
            6, 6, 0, 0, 0, 0, 0, 0, 0, 0,
    };
    //GAL E1(1) E5a(5) E5b(7) E5ab(8) E6(6)
    const unsigned char kGalFreqBand[]={
            0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
            1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            2, 2, 2, 0, 0, 0, 5, 5, 5, 5,
            5, 0, 0, 0, 0, 0, 3, 3, 3, 0,
            0, 0, 4, 4, 4, 0, 0, 0, 0, 0,
    };

    //GLO
    const unsigned char kGloFreqBand[]={
            0, 1, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 2, 0, 0, 0, 0, 2, 0,
            0, 0, 0, 0, 3, 3, 3, 4, 4, 4,
            0, 0, 0, 0, 0, 0, 5, 5, 5, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    //QZS
    const unsigned char kQzsFreqBand[]={
            0, 1, 1, 1, 1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 2, 2, 2, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            3, 3, 3, 3, 3, 0, 0, 0, 4, 0,
            4, 4, 4, 4, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    };

    const string kGnssSignalPriors[NSYS][MAX_GNSS_FRQ_NUM]{
            {"CWPYMNSL", "CWPYMNDSLX",      "IQX",         "",        "",    ""},
            {     "IQX",      "IQXA",      "IQX",     "DPXA",     "DPX", "DPZ"},
            {   "CABXZ",       "IQX",      "IQX",      "IQX",   "CABXZ",    ""},
            {      "PC",        "PC",      "IQX",      "ABX",    "ABXP",    ""},
            {   "CSLXZ",       "SLX",   "IQXDPZ",    "SLXEZ",        "",    ""},
    };

    const int MAX_GNSS_CODE_BIAS_PAIRS=12;
    const string kGnssCodeBiasPairs[NSYS][MAX_GNSS_CODE_BIAS_PAIRS]{
            {"C1C-C2W", "C1C-C5Q", "C1C-C5X", "C1W-C2W", "C1C-C1W", "C2C-C2W", "C2W-C2S", "C2W-C2L", "C2W-C2X", "", "", ""},
            {"C2I-C7I", "C2I-C6I", "C1X-C5X", "C1P-C5P", "C1D-C5D", "C1X-C6I", "C1P-C6I", "C1D-C6I", "C2I-C6I", "C1X-C7Z", "C1X-C8X", ""},
            {"C1C-C5Q", "C1C-C6C", "C1C-C7Q", "C1C-C8Q", "C1X-C5X", "C1X-C7X", "C1X-C8X", "", "", "", "", ""},
            {"C1C-C2C", "C1C-C2P", "C1P-C2P", "C1C-C1P", "C2C-C2P", "", "", "", "", "", "", ""},
            {"C1C-C2L", "C1C-C5X", "C1C-C5Q", "C1X-C2X", "C1X-C5X", "C1C-C1X", "", "", "", "", "", ""},
    };

    const int NUM_BD2_GEO=5;
    const int NUM_BD2_IGSO=7;
    const int NUM_BD2_MEO=3;
    const int kBD2_GEO[]={1,2,3,4,5};
    const int kBD2_IGSO[]={6,7,8,9,10,13,16};
    const int kBD2_MEO[]={11,12,14};
    const int GLO_FRQ_NUM[]={-7,-6,-4,-3,-2,-1,0,1,2,3,4,5,6};

    enum GNSS_OBS {
        GNSS_OBS_CODE,
        GNSS_OBS_PHASE,
        GNSS_OBS_DOPPLER,
        GNSS_OBS_SNR
    };

    enum GNSS_TRP_OPT {
        TRP_OFF,
        TRP_SAAS,
        TRP_EST_WET,
        TRP_EST_GRAD
    };

    enum GNSS_TID_OPT {
        TID_OFF,
        TID_SOLID,
        TID_OCEAN,
        TID_POLE
    };

    enum GNSS_EPH_OPT {
        EPH_BRD,
        EPH_PRE
    };

    enum GNSS_AC_OPT {
        AC_BRD,
        AC_WUM,
        AC_GBM,
        AC_COM,
    };

    const string kGnssAcStr[]={
            "brd","wum","gbm","com"
    };

    enum GNSS_ION_OPT {
        ION_OFF,
        ION_KLB,
        ION_TEC,
        ION_IF,
        ION_IF_DUAL,
        ION_EST,
        ION_CONST
    };

    const string kIonOptStr[]={
            "OFF","KLB","TEC","IF","IF2","UC","CONST"
    };

    enum GLO_IFCB_OPT {
        GLO_IFCB_OFF,
        GLO_IFCB_LNF,
        GLO_IFCB_QUAD,
        GLO_IFCB_1SAT,
        GLO_IFCB_1FRQ
    };

    enum GNSS_FRQ_OPT {
        FRQ_SINGLE,
        FRQ_DUAL,
        FRQ_TRIPLE
    };

    const string kGnssFrqStr[]={
            "SF","DF","TF"
    };


    enum GNSS_SAT_STAT {
        SAT_NO_USE=-1,
        SAT_USED=0,
        SAT_NO_CODE,
        SAT_NO_PHASE,
        SAT_NO_PROD,
        SAT_LOW_EL,
        SAT_PRI_RES_C,
        SAT_PRI_RES_P,
        SAT_LAR_CP_DIFF,
        SAT_POS_RES_C,
        SAT_POS_RES_P,
    };

    const string kGnssSatStatStr[]{
        "NO_USE",
        "USED",
        "NO_CODE",
        "NO_PHASE",
        "NO_PROD",
        "LOE_EL",
        "PRIOR_RES_C",
        "PRIOR_RES_P",
        "LAR_CP_DIFF",
        "POST_RES_C",
        "POST_RES_P"
    };

    const string kGnssSysStr[]{
        "NONE",
        "GPS",
        "BDS",
        "GAL",
        "GLO",
        "QZS"
    };

    enum GNSS_OBS_COMB {
        COMB_IF,
        COMB_MW,
        COMB_GF,
        COMB_MP,
        COMB_CSC,
        COMB_TDCP,
        COMB_SSD
    };

    enum GNSS_AR_MODE {
        AR_OFF,
        AR_CONT,
        AR_INST,
        AR_FIX_HOLD,
        AR_PPP_AR,
        AR_PPP_AR_ILS,
    };

    enum GLO_AR_MODE {
        GLO_AR_OFF,
        GLO_AR_ON,
        GLO_AR_AUTO,
        GLO_AR_FIXHOLD,
    };

    enum GNSS_RES_QC {
        RES_QC_SIMPLE,
        RES_QC_IGG,
        RES_QC_STEP,
    };

    // ins related constant definitions
    enum IMU_TYPE {
        IMU_UNKNOW=-1,
        IMU_NOVTEL_CPT,
        IMU_NOVTEL_A1,
        IMU_M39,
        IMU_MTI_CSV
    };

    enum IMU_COORD_TYPE {
        IMU_COORD_RFU,
        IMU_COORD_FRD,
    };

    enum IMU_DATA_FORMAT {
        IMU_FORMAT_INCR,
        IMU_FORMAT_RATE,
    };

    enum GYRO_DATA_FORMAT {
        GYRO_FORMAT_DEG,
        GYRO_FORMAT_RAD
    };

    enum SOLVE_ESTIMATOR {
        SOLVE_LSQ,
        SOLVE_KF,
    };

    enum PPPLIB_MODE {
        MODE_SPP=0,
        MODE_PPP,
        MODE_DGNSS,
        MODE_PPK,
        MODE_INS,
        MODE_IGLC,
        MODE_IGTC,
    };

    const string kPpplibModeStr[]={
            "SPP","PPP","DGPS","PPK","INS","IGLC","IGTC"
    };


    enum PPPLIB_MODE_OPT {
        MODE_OPT_STATIC=0,
        MODE_OPT_KINE_SIM,
        MODE_OPT_KINEMATIC,
        MODE_OPT_GSOF,
        MODE_OPT_SPP,
        MODE_OPT_PPP,
        MODE_OPT_PPK,
    };

    const string kPpplibModeOptStr[]={
            "STATIC","KINE_SIM","KINE","GSOF","SPP","PPP","PPK"
    };

    enum SOL_STAT {
        SOL_NONE=0,
        SOL_FIX,
        SOL_FLOAT,
        SOL_SBAS,
        SOL_DGPS,
        SOL_SPP,
        SOL_PPP,
        SOL_REF,
    };

    enum SOL_INS_STAT {
        SOL_INS_NONE,
        SOL_IG_TC,
        SOL_IG_LC,
        SOL_INS_MECH,
        SOL_INS_INIT,
        SOL_INS_REF
    };

    const double kChiSqr[100]={  //chi-sqr(n) (alpha=0.001)
            10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
            31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
            46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
            61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
            74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
            88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
            101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
            113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
            126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
            138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149
    };

    const double UNC_POS=30.0;
    const double UNC_VEL=30.0;
    const double UNC_ATT=(10.0*D2R);
    const double UNC_BA=(1000.0*MG2MS2);
    const double UNC_BG=(10.0*D2R/3600.0);

    #define DIS_FLAG 123456

    enum STATISTICAL_MODE{
        STATISTICAL_CT,
        STATISTICAL_PWC,
        STATISTICAL_RW,
        STATISTICAL_WN,
    };
}

#endif //PPPLIB_CONSTANT_H
