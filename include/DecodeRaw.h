//
// Created by cc on 8/22/20.
//

#ifndef PPPLIB_DECODERAW_H
#define PPPLIB_DECODERAW_H

#include "CmnFunc.h"
#include "GnssFunc.h"
#include "InsFunc.h"

#define MAX_RAW_LEN 4096

namespace PPPLib{
    typedef struct{
        cTime t_tag;

        tGsofUnit gsof_unit;
        tImuDataUnit imu_unit;

        void *strp;
        unsigned char dire; //direction of input raw data(0:forward 1:backward 2:combined)
        unsigned char buff[MAX_RAW_LEN];
        int nbyte;          // number of bytes in message buffer
        int len;            // message length
        int curb;           //index of raw data for backward solution

    }tRaw;

    class cDecodeGsof{
    public:
        cDecodeGsof();
        ~cDecodeGsof();

    public:
        bool DecodeGsof(const string file,vector<tSolInfoUnit>& gsof);

    private:
        void DecodeMessage(tRaw *raw);
        int DecodeMessageForward(tRaw *raw, unsigned char data);
        int DecodeMessageBackward(tRaw *raw, unsigned char data);
        bool InputGsofMessage(tRaw *raw, unsigned char data);

        bool CheckHead(tRaw *raw);
        bool CheckSum(tRaw *raw);
        unsigned short GetU16(unsigned char **data);
        unsigned long GetU32(unsigned char **data);
        float GetFloat(unsigned char **data);
        double GetDouble(unsigned char **data);

        void DecodeType01(tRaw *raw, unsigned char* data);
        void DecodeType02(tRaw *raw, unsigned char* data);
        void DecodeType03(tRaw *raw, unsigned char* data);
        void DecodeType06(tRaw *raw, unsigned char* data);
        void DecodeType08(tRaw *raw, unsigned char* data);
        void DecodeType09(tRaw *raw, unsigned char* data);
        void DecodeType10(tRaw *raw, unsigned char* data);
        void DecodeType11(tRaw *raw, unsigned char* data);
        void DecodeType12(tRaw *raw, unsigned char* data);
        void DecodeType38(tRaw *raw, unsigned char* data);
        void DecodeType41(tRaw *raw, unsigned char* data);

        void GsofUnit2ppplibsol(tGsofUnit& gsof,tSolInfoUnit& ppplib_sol);
        static bool CmpSolInfo(const tSolInfoUnit& p1, const tSolInfoUnit& p2);
    };

    class cDecodeImuM39 {
    public:
        cDecodeImuM39();
        ~cDecodeImuM39();

    public:
        int DecodeM39(string file,tInsConf insC,vector<tImuDataUnit>& imus);

    private:
        int InputM39Message(tRaw *raw, unsigned char data);
        bool DecodeM39Forward(tRaw *raw, unsigned char data);
        bool DecodeM39Backward(tRaw *raw, unsigned char data);
        int DecodeMessage(tRaw *raw);
        void DecodeSowTime(tRaw *raw,double *sow,int *start);
        void AdjustTime(tRaw *raw,const double sowi,double *sowo,double *timu,int *week, unsigned int *dcc);
        bool DecodeM39Odo(tRaw *raw,double dt);

        bool CheckHead(const unsigned char *buff);
        unsigned char CheckSum(const unsigned char *buff,int len);
        static int CmpImuData(const void *p1, const void *p2);

    };
}



#endif //PPPLIB_DECODERAW_H
