#ifndef READPARAMETERS_H
#define READPARAMETERS_H

#include <string>
#include <fstream>
#include <vector>


class ReadParameters
{
 public:

  ReadParameters (const char* InputFileName, std::string howOpen = "in");
  ~ReadParameters () {fileHandler.close();};
  void ReadFromFile (unsigned int blockNumber, std::vector<std::string>* ParVector);
  void WriteToFile  (std::vector<std::string>* ParVector);


 private:

  std::fstream fileHandler;
};

#endif
