#ifndef MonoAlgorithms_MonoEcalCalibReader_h
#define MonoAlgorithms_MonoEcalCalibReader_h

#include <fstream>

#include "FWCore/Utilities/interface/Exception.h"

#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"


namespace Mono {


class MonoEcalCalibReader {


public:
  inline MonoEcalCalibReader() { }

  virtual inline ~MonoEcalCalibReader() { }

  virtual inline void readCalib(const std::string &fileName,MIJType *calibMap, MIJType *avgMap) const {
    assert(calibMap);
    calibMap->clear();

    assert(avgMap);
    avgMap->clear();

    // open file for reading
    std::ifstream inf;
    inf.open(fileName.c_str(),std::ios::in);

    // read the file 
    while ( inf.good() && !inf.eof() ) {
      char line[500];
      inf.getline(line,500);
      char *strCheck = strstr(line,"BEGIN MAP");
      if ( strCheck ){
	std::string str(strCheck); 
 	std::cout << str << std::endl;
        std::size_t space = str.find(' ',10);
        assert( space != std::string::npos );
  	const unsigned length = atoi( str.substr(10,space).c_str() );
        const unsigned width = atoi( str.substr(space+1).c_str() );
  	const unsigned side = length*width;	
	ClustCategorizer cat(length,width);
        std::vector<double> & vec = (*calibMap)[cat];
        vec.resize(side*side);
	for ( unsigned i=0; i != side; i++ ) {
	  for ( unsigned j=0; j != side-1U; j++ ) {
	    inf.getline(line,500,' ');
	    vec[i*side+j] = atof(line);
	  }
	  inf.getline(line,500,'\n');
	  vec[i*side+side-1] = atof(line);
	}
      } 	
      //inf.getline(line,500);
      strCheck = strstr(line,"BEGIN MEAN");
      if ( strCheck ){
	std::string str(strCheck); 
 	std::cout << str << std::endl;
        std::size_t space = str.find(' ',11);
        assert( space != std::string::npos );
  	const unsigned length = atoi( str.substr(11,space).c_str() );
        const unsigned width = atoi( str.substr(space+1).c_str() );
  	const unsigned size = length*width;	
	ClustCategorizer cat(length,width);
        std::vector<double> & vec = (*avgMap)[cat];
        vec.resize(size);
	for ( unsigned i=0; i != width; i++ ) {
	  for ( unsigned j=0; j != length-1U; j++ ) {
	    inf.getline(line,500,' ');
	    vec[i*length+j] = atof(line);
	  }
	  inf.getline(line,500,'\n');
	  vec[i*length+length-1U] = atof(line);
	} 	
	
      }
 
    }

    assert( calibMap->size() == avgMap->size() ); 

  }

  virtual inline void dumpCalib(const std::string &fileName,const MIJType &calibMap, const MIJType &avgMap) const {

    if ( !calibMap.size() )
      throw cms::Exception("MonoEcalCalib calibration map is emtpy");

    if ( !avgMap.size() )
      throw cms::Exception("MonoEcalCalib mean values missing");

    // open file for writing
    std::ofstream ouf;
    ouf.open(fileName.c_str(),std::ios::out);

    // dump map contents  
    MIJType::const_iterator iter = calibMap.begin();
    MIJType::const_iterator iEnd = calibMap.end();
    for ( ; iter != iEnd; iter++ ) {
      const ClustCategorizer & catz = iter->first;
      std::cout << "Dumping Map: " << catz.length << " " << catz.width << std::endl;
      ouf << "BEGIN MAP " << catz.length << " " << catz.width << std::endl;
      const std::vector<double> & data = iter->second;
      const unsigned size = data.size();
      const unsigned side = catz.length*catz.width;
      if ( size == 0 ) continue;
      assert( side*side == size );
      for ( unsigned i=0; i != side; i++ ) {
 	for ( unsigned j=0; j != side; j++ ) {
	  ouf << data[i*side+j] << " ";
	}
	ouf << std::endl;
      }
      ouf << "END MAP" << std::endl; 

    }

    // dump mean contents
    iter = avgMap.begin();
    iEnd = avgMap.end();
    for ( ; iter != iEnd; iter++ ) {
      const ClustCategorizer & catz = iter->first;
      std::cout << "Dumping Mean: " << catz.length << " " << catz.width << std::endl;
      ouf << "BEGIN MEAN " << catz.length << " " << catz.width << std::endl;
      const std::vector<double> & data = iter->second;
      const unsigned size = data.size();
      const unsigned length = catz.length;
      const unsigned width = catz.width;
      assert( length*width == size );
      for ( unsigned j=0; j != width; j++ ) {
   	for ( unsigned i=0; i != length; i++ ) {
	  ouf << data[j*length+i] << " ";
  	}
	ouf << std::endl;
      }
      ouf << "END MEAN" << std::endl;
    }
    

    ouf.close();

  }

private:


};


} // end mono namespace

#endif
