//
// Created by chenc on 2020/7/15.
//

#include "LogInfo.h"

namespace PPPLib {

    string SetLogConfPath(string path) {
        if(!path.empty()) return path;
        else return LOGINI_PATH;
    }

    int SetLogLevel(int level) {
        return level;
    }

    void InitLog(int argc, char** argv,string path, int level){
        START_EASYLOGGINGPP(argc, argv);
        el::Configurations conf(path);
        el::Loggers::reconfigureAllLoggers(conf);
        el::Loggers::addFlag(el::LoggingFlag::StrictLogFileSizeCheck);
        el::Loggers::addFlag(el::LoggingFlag::ColoredTerminalOutput);
        el::Loggers::addFlag(el::LoggingFlag::HierarchicalLogging);
        el::Loggers::setLoggingLevel(static_cast<el::Level>(level));
    }

}