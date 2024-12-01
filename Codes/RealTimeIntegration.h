#ifndef REALTIMEINTEGRATION_H
#define REALTIMEINTEGRATION_H

#include <curl/curl.h>
#include <json.hpp>
#include <fstream>

#include <TString.h>
#include <iostream>

using json = nlohmann::json;

namespace API
{
  class RTI
  {
  private:
    CURL *
        _curlHandle;
    CURLcode
        _result;
    TString
        _URLschema,
        _interface;

    std::vector<TString>
        _MultiURL;

    std::vector<std::string>
                      _resptable;

    std::string
            _response;

    Int_t
        _pdgidMax;

    json
      _responseJSON;

  public:
    RTI(const char *URLschema);
    RTI(TString URLschema);
    RTI();

    void setURL(const char *URLschema);
    void setURL(TString URLschema);
    void setOpt();

    void setMultiURL(TString propName);

    void getAPICall();
    void getMultiAPICall();
    Double_t getValue(Int_t i);

    static size_t write_to_string(void *ptr, size_t size, size_t count, void *stream);
  };

}

#endif // !REALTIMEINTEGRATION_H