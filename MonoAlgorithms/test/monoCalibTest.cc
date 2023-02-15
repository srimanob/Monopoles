///////////////////////////////////////////////
// Test the MonoEcalCalibReader read and write utilities.
///////////////////////////////////////////////

#include <vector>
#include <iostream>
#include <cassert>
#include <cstring>
#include <string>

#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalCalibReader.h"



int main(int argc, char **argv) {

  // process options
  std::string calibName("testCalib.dat");

  const bool checkEMap = true;
  const bool checkHMap = true;

  std::vector<Mono::ClustCategorizer> categories;
  categories.push_back(Mono::ClustCategorizer(3,5));
  categories.push_back(Mono::ClustCategorizer(5,5));
  categories.push_back(Mono::ClustCategorizer(6,5));

  const unsigned catSize = categories.size();

  // generate mean maps
  Mono::MIJType origEMap;
  for ( unsigned c = 0; c != catSize; c++ ) {
    const Mono::ClustCategorizer & cat = categories[c];
    std::vector<double> & data = origEMap[cat];
    const unsigned length = cat.length;
    const unsigned width = cat.width;
    const unsigned size = length*width;

    data.resize(size);
    for ( unsigned i=0; i != size; i++ ) 
      data[i] = i;
  }
  assert( origEMap.size() == catSize ); 

  // generate H maps
  Mono::MIJType origHMap;
  for ( unsigned c = 0; c != catSize; c++ ) {
    const Mono::ClustCategorizer & cat = categories[c];
    std::vector<double> & data = origHMap[cat];
    const unsigned length = cat.length;
    const unsigned width = cat.width;
    const unsigned side = length*width;
    const unsigned size = side*side;
  
    data.resize(size);
    for( unsigned i=0; i != size; i++ ) 
      data[i] = i;
  }
  assert( origHMap.size() == catSize );

  // dump calibration
  Mono::MonoEcalCalibReader calibReader;
  calibReader.dumpCalib(calibName,origHMap,origEMap); 

  // read calibration
  Mono::MIJType recoEMap;
  Mono::MIJType recoHMap;
  calibReader.readCalib(calibName,&recoHMap,&recoEMap);

  // check maps
  if ( checkEMap ) assert( recoEMap.size() == origEMap.size() );
  if ( checkHMap ) assert( recoHMap.size() == origHMap.size() );

  if ( checkEMap ) {
  for ( unsigned c=0; c != catSize; c++ ) {
    const Mono::ClustCategorizer & cat = categories[c];
    // check EMaps
    Mono::MIJType::const_iterator origIter = origEMap.find(cat);
    Mono::MIJType::const_iterator recoIter = recoEMap.find(cat);

    if ( recoIter == recoEMap.end() ) {
      std::cerr << "Could not find category in recoEMap: " << cat.length << " " << cat.width << std::endl;
      Mono::MIJType::const_iterator it = recoEMap.begin();
      for ( ; it != recoEMap.end(); it++ ) {
        std::cout << "Found categoryin recoEMap: " << it->first.length << " " << it->first.width << std::endl;
      }
      return 1;
    } else if ( origIter == origEMap.end() ) {
      std::cerr << "Could not find category in origEMap: " << cat.length << " " << cat.width << std::endl;
      return 1;
    }
   
    const std::vector<double> & origData = origIter->second;
    const std::vector<double> & recoData = recoIter->second;
  
    const unsigned origSize = origData.size();
    const unsigned recoSize = recoData.size();
    if ( origSize != recoSize ) {
      std::cerr << "Data Sizes do not match to origEMap " << std::endl;
      return 1;
    }

    for ( unsigned i=0; i != origSize; i++ ) {
      if ( origData[i] != recoData[i] ) {
	std::cerr << "Data mismatch in EMap" << std::endl;
   	return 1;
      }
    }

  }
  } // end EMap check



  if ( checkHMap ) {
  for ( unsigned c=0; c != catSize; c++ ) { 
    const Mono::ClustCategorizer & cat = categories[c]; 
    // check HMaps
    Mono::MIJType::const_iterator origIter = origHMap.find(cat);
    Mono::MIJType::const_iterator recoIter = recoHMap.find(cat);

    if ( recoIter == recoHMap.end() ) {
      std::cerr << "Could not find category in recoHMap: " << cat.length << " " << cat.width << std::endl;
      return 1;
    } else if ( origIter == origHMap.end() ) {
      std::cerr << "Could not find category in origHMap: " << cat.length << " " << cat.width << std::endl;
      return 1;
    }
   
    const std::vector<double> & origHData = origIter->second;
    const std::vector<double> & recoHData = recoIter->second;
  
    const unsigned origHSize = origHData.size();
    const unsigned recoHSize = recoHData.size();
    if ( origHSize != recoHSize ) {
      std::cerr << "Data Sizes do not match to origHMap " << std::endl;
      return 1;
    }

    for ( unsigned i=0; i != origHSize; i++ ) {
      if ( origHData[i] != recoHData[i] ) {
	std::cerr << "Data mismatch in HMap" << std::endl;
   	return 1;
      }
    }
  }
  } // end HMap check


  std::cout << "All tests passed successfully!! Have a nice day." << std::endl;

  return 0;
}
