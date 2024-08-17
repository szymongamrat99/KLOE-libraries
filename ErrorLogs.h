#ifndef ERROR_LOGS_H
#define ERROR_LOGS_H

#include <map>

namespace ErrorHandling
{
  enum class ErrorCodes
  {
    DATA_TYPE = 100,
    RANGE = 101,

    NOT_RECOGNIZED = 666
  };

  class ErrorLogs
  {
  private:
    std::map<ErrorCodes, TString> Logs;

  public:
    ErrorLogs()
    {
      Logs[ErrorCodes::DATA_TYPE] = Form("Invalid input data type. Error code %d", int(ErrorCodes::DATA_TYPE));
      Logs[ErrorCodes::RANGE] = Form("File number outside available range. Check const.h file. Error code %d", int(ErrorCodes::RANGE));
      Logs[ErrorCodes::NOT_RECOGNIZED] = Form("Error not recognized. Error code %d", int(ErrorCodes::NOT_RECOGNIZED));

      std::cout << Logs[ErrorCodes::DATA_TYPE] << std::endl;
    }

    void getErrLog(ErrorCodes errCode)
    {
      std::cerr << "--------------------------------------------------" << std::endl;
      std::cerr << Logs[errCode] << std::endl;
      std::cerr << "--------------------------------------------------" << std::endl;
    };
  };
}

#endif