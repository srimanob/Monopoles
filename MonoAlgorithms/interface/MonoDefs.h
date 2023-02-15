#ifndef MONOANALYSIS_MONODEFS_H 
#define MONOANALYSIS_MONODEFS_H

#define MONOID 4110000

#include <map>
#include <vector>
#include <string>

#include "Monopoles/MonoAlgorithms/interface/ClustCategorizer.h"


namespace Mono {

// enumeration to determine monopole from anti_monopole
enum MonoEnum {
  monopole=0
  ,anti_monopole
};

// forward monopole class declarations
class ClustCategorizer;

typedef std::map<ClustCategorizer,std::vector<double> >	          MIJType;


} // end Mono namespace


#endif
