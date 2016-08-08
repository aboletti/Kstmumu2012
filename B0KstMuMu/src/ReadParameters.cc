#include "../interface/ReadParameters.h"

#include <cstdlib>
#include <iostream>

ReadParameters::ReadParameters(const char* InputFileName, std::string howOpen)
{
  if (howOpen == "in")
    {
      fileHandler.open(InputFileName, std::ifstream::in);
      
      if (fileHandler.good() == false)
	{
	  std::cout << "[ReadParameters::ReadFromFile]\tError opening file : " << InputFileName << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else if (howOpen == "out")
    {
      fileHandler.open(InputFileName, std::ofstream::out);
      
      if (fileHandler.good() == false)
	{
	  std::cout << "[ReadParameters::ReadFromFile]\tError opening file : " << InputFileName << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else if (howOpen == "app")
    {
      fileHandler.open(InputFileName, std::ofstream::app | std::ofstream::out);
      
      if (fileHandler.good() == false)
	{
	  std::cout << "[ReadParameters::ReadFromFile]\tError opening file : " << InputFileName << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
}

void ReadParameters::ReadFromFile (unsigned int blockNumber, std::vector<std::string>* ParVector)
// ####################################################################################################################################
// # blockNumber  = 0 --> the file contains just one block of data without the number of lines as header                              #
// # blockNumber != 0 --> the file contains at least one block of data, before each block there must be the number of lines as header #
// ####################################################################################################################################
{
  std::string ReadRows;

  if (blockNumber == 0)
    {
      fileHandler.clear();
      fileHandler.seekg(std::ifstream::beg);      
      while(fileHandler.good() == true)
	{
	  getline(fileHandler,ReadRows);
	    
	  if ((ReadRows.length() != 0) && (ReadRows[0] != '#')) ParVector->push_back(ReadRows);
	}
    }
  else
    {
      std::vector<std::string>* tmpVec = new std::vector<std::string>;

      fileHandler.clear();
      fileHandler.seekg(std::ifstream::beg);
      while (fileHandler.good() == true)
	{
	  getline(fileHandler,ReadRows);

	  if ((ReadRows.length() != 0) && (ReadRows[0] != '#')) tmpVec->push_back(ReadRows);
	}
      fileHandler.clear();
      fileHandler.seekg(std::ifstream::beg);

      int Start = 0;
      int NParam = atoi(tmpVec->operator[](Start).c_str());
      blockNumber--;
      while (blockNumber != 0)
	{
	  Start = Start + 1 + NParam;
	  NParam = atoi(tmpVec->operator[](Start).c_str());
	  blockNumber--;
	}

      for (int i = 0; i < NParam; i++) ParVector->push_back(tmpVec->operator[](Start+1+i));

      tmpVec->clear();
      delete tmpVec;
    }
}

void ReadParameters::WriteToFile (std::vector<std::string>* ParVector)
{
  for (unsigned int i = 0; i < ParVector->size(); i++) fileHandler << ParVector->operator[](i) << std::endl;
}
