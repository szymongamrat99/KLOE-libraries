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

    DELTA_LT_ZERO = 200,
    DENOM_EQ_ZERO = 201,

    NOT_RECOGNIZED = 666
  };

  class ErrorLogs
  {
  private:
    std::map<ErrorCodes, TString> Logs;

  public:
    ErrorLogs()
    {
      // General Logs
      Logs[ErrorCodes::DATA_TYPE] = Form("Invalid input data type. Error code %d", int(ErrorCodes::DATA_TYPE));
      Logs[ErrorCodes::RANGE] = Form("File number outside available range. Check const.h file. Error code %d", int(ErrorCodes::RANGE));
      Logs[ErrorCodes::MENU_RANGE] = Form("Chosen menu option outside range. Error code %d", int(ErrorCodes::MENU_RANGE));
      Logs[ErrorCodes::FILE_NOT_EXIST] = Form("File does not exist. Error code %d", int(ErrorCodes::FILE_NOT_EXIST));

      // Mathematics logs
      Logs[ErrorCodes::DELTA_LT_ZERO] = Form("Quadratic delta less than 0. Error code %d", int(ErrorCodes::DELTA_LT_ZERO));
      Logs[ErrorCodes::DENOM_EQ_ZERO] = Form("Denominator in fraction equal to 0. Error code %d", int(ErrorCodes::DENOM_EQ_ZERO));

      // Not recognized logs
      Logs[ErrorCodes::NOT_RECOGNIZED] = Form("Error not recognized. Error code %d", int(ErrorCodes::NOT_RECOGNIZED));
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
  };
}

#endif