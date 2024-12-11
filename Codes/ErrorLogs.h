#ifndef ERROR_LOGS_H
#define ERROR_LOGS_H

#include <fstream>
#include <iostream>
#include <map>

#include <TString.h>

namespace ErrorHandling
{
  enum class ErrorCodes
  {
    DATA_TYPE = 100,
    RANGE = 101,
    MENU_RANGE = 102,
    FILE_NOT_EXIST = 103,
    TREE_NOT_EXIST = 104,

    DELTA_LT_ZERO = 200,
    DENOM_EQ_ZERO = 201,

    NOT_RECOGNIZED = 666,

    NUM_CODES = 8
  };

  class ErrorLogs
  {
  private:
    ErrorCodes errCodesList[int(ErrorCodes::NUM_CODES)];
    std::map<ErrorCodes, TString> Logs;
    std::map<ErrorCodes, int> errCount;

  public:
    ErrorLogs()
    {
      // Creation of iterable list of error codes
      errCodesList[0] = ErrorCodes::DATA_TYPE;
      errCodesList[1] = ErrorCodes::RANGE;
      errCodesList[2] = ErrorCodes::MENU_RANGE;
      errCodesList[3] = ErrorCodes::FILE_NOT_EXIST;
      errCodesList[4] = ErrorCodes::TREE_NOT_EXIST;
      errCodesList[5] = ErrorCodes::DELTA_LT_ZERO;
      errCodesList[6] = ErrorCodes::DENOM_EQ_ZERO;
      errCodesList[7] = ErrorCodes::NOT_RECOGNIZED;

      // General Logs
      Logs[ErrorCodes::DATA_TYPE] = Form("Invalid input data type. Error code %d", int(ErrorCodes::DATA_TYPE));
      Logs[ErrorCodes::RANGE] = Form("File number outside available range. Check const.h file. Error code %d", int(ErrorCodes::RANGE));
      Logs[ErrorCodes::MENU_RANGE] = Form("Chosen menu option outside range. Error code %d", int(ErrorCodes::MENU_RANGE));
      Logs[ErrorCodes::FILE_NOT_EXIST] = Form("File does not exist. Check the name. Error code %d", int(ErrorCodes::FILE_NOT_EXIST));
      Logs[ErrorCodes::TREE_NOT_EXIST] = Form("Tree does not exist. Check the name. Error code %d", int(ErrorCodes::TREE_NOT_EXIST));

      // Mathematics logs
      Logs[ErrorCodes::DELTA_LT_ZERO] = Form("Quadratic delta less than 0. Error code %d", int(ErrorCodes::DELTA_LT_ZERO));
      Logs[ErrorCodes::DENOM_EQ_ZERO] = Form("Denominator in fraction equal to 0. Error code %d", int(ErrorCodes::DENOM_EQ_ZERO));

      // Not recognized logs
      Logs[ErrorCodes::NOT_RECOGNIZED] = Form("Error not recognized. Error code %d", int(ErrorCodes::NOT_RECOGNIZED));

      // Init of errCount
      for(int i = 0; i < int(ErrorCodes::NUM_CODES); i++)
        errCount[errCodesList[i]] = 0;
    }

    void setErrCount(ErrorCodes errCode)
    {
      errCount[errCode]++;
    }

    void getErrLog(ErrorCodes errCode)
    {
      std::cerr << "--------------------------------------------------" << std::endl;
      std::cerr << Logs[errCode] << std::endl;
      std::cerr << "--------------------------------------------------" << std::endl;
    };

    void getErrLog(ErrorCodes errCode, std::ofstream& LogFile)
    {
      LogFile << "--------------------------------------------------" << std::endl;
      LogFile << Logs[errCode] << std::endl;
      LogFile << "--------------------------------------------------" << std::endl;
    };

    void errLogStats()
    {
      std::cout << "Statistics of errors during the execution (by codes)" << std::endl;
      for(int i = 0; i < int(ErrorCodes::NUM_CODES); i++)
      {
        std::cout << "--------------------------------------------------" << std::endl;
        std::cout << "Error code " << int(errCodesList[i]) << ": " << errCount[errCodesList[i]] << " exceptions" << std::endl;
        std::cout << "--------------------------------------------------" << std::endl;
      }
    };

    void errLogStats(std::ofstream& LogFile)
    {
      LogFile << "Statistics of errors during the execution (by codes)" << std::endl;
      for(int i = 0; i < int(ErrorCodes::NUM_CODES); i++)
      {
        LogFile << "--------------------------------------------------" << std::endl;
        LogFile << "Error code " << int(errCodesList[i]) << ": " << errCount[errCodesList[i]] << " exceptions" << std::endl;
        LogFile << "--------------------------------------------------" << std::endl;
      }
    };
  };
}

#endif