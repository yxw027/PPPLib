//
// Created by cc on 8/3/20.
//

#ifndef PPPLIB_OUTSOL_H
#define PPPLIB_OUTSOL_H

#include "CmnFunc.h"
#include "GnssFunc.h"
#include "InsFunc.h"

#define COMMENTH   "%"

namespace PPPLib {
    class cOutSol{
    public:
        cOutSol();
        cOutSol(tPPPLibConf C);
        cOutSol(tPPPLibConf C,vector<tSolInfoUnit>& ref_sols);
        ~cOutSol();

    private:
        int MatchRefSol(cTime sol_time);
        tSolInfoUnit CompareSol(tSolInfoUnit& sol,tSolInfoUnit& ref_sol);
        int OutEcef(unsigned char* buff,const char *s,tSolInfoUnit& sol);
        int OutSolStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos,char *buff);
    public:
        bool InitOutSol(tPPPLibConf C,string file);
        void WriteHead();
        void WriteSol(tSolInfoUnit sol,int epoch);
        void WriteSatStat(tSolInfoUnit *sol,tSatInfoUnit *sat_infos);

        void WriteImuHead();
        void WriteImuObs();

    private:
        int ref_index_=0;
        tPPPLibConf C_;
        FILE *fout_;
        string buff_;
        FILE *fout_stat_;

    public:
        vector<tSolInfoUnit> ref_sols_;
        cImuData *imus= nullptr;
    };
}


#endif //PPPLIB_OUTSOL_H
