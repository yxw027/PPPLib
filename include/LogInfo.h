//
// Created by chenc on 2020/7/15.
//

#ifndef PPPLIB_LOGINFO_H
#define PPPLIB_LOGINFO_H

#include "easylogging++.h"

using namespace std;

namespace PPPLib{

#ifdef WIN32
    #define LOGINI_PATH "..\\conf\\log.ini"
#else
    #define LOGINI_PATH "../conf/log.ini"
#endif

    string SetLogConfPath(string path);
    int SetLogLevel(int level);
    void InitLog(int argc, char** argv, string path, int level);


}

#endif //PPPLIB_LOGINFO_H
