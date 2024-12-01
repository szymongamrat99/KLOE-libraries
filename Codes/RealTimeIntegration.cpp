#include <RealTimeIntegration.h>

API::RTI::RTI(const char *URLschema) : _URLschema(URLschema)
{
  _curlHandle = curl_easy_init();
}

API::RTI::RTI(TString URLschema) : _URLschema(URLschema)
{
  _curlHandle = curl_easy_init();
}

API::RTI::RTI()
{
  _curlHandle = curl_easy_init();
}

void API::RTI::setURL(const char *URLschema)
{
  _URLschema = URLschema;
}

void API::RTI::setURL(TString URLschema)
{
  _URLschema = URLschema;
}

void API::RTI::setOpt()
{
  // curl_easy_setopt(_curlHandle, CURLOPT_URL, _URLschema.Data());
  curl_easy_setopt(_curlHandle, CURLOPT_WRITEFUNCTION, &API::RTI::write_to_string);
  curl_easy_setopt(_curlHandle, CURLOPT_WRITEDATA, &_response);
}

void API::RTI::setMultiURL(TString propName)
{
  std::ifstream fprop(propName);
  json MultiURL = json::parse(fprop);

  _URLschema = (std::string)MultiURL["PDGURL"];
  _interface = (std::string)MultiURL["interfaces"]["summary"];
  _pdgidMax = (Int_t)MultiURL["pdgids"]["length"];

  TString
      pdgID = "",
      combined = "";

  for (Int_t i = 0; i < _pdgidMax; i++)
  {
    pdgID = (std::string)MultiURL["pdgids"]["values"][i];
    combined = _URLschema + _interface + pdgID;
    _MultiURL.push_back(combined);
  }

  // fprop.close();
}

void API::RTI::getAPICall()
{
  _result = curl_easy_perform(_curlHandle);
}

void API::RTI::getMultiAPICall()
{
  for (Int_t i = 0; i < _pdgidMax; i++)
  {
    curl_easy_setopt(_curlHandle, CURLOPT_URL, _MultiURL[i].Data());
    _result = curl_easy_perform(_curlHandle);
    _resptable.push_back(_response);
  }
}

Double_t API::RTI::getValue(Int_t i)
{
  _responseJSON = json::parse(_resptable[i]);

  return _responseJSON["pdg_values"][0]["value"];
}



size_t API::RTI::write_to_string(void *ptr, size_t size, size_t count, void *stream) {
  ((std::string*)stream)->append((char*)ptr, 0, size*count);
  return size*count;
}
