#ifndef MONOALGORITHMS_CLUSTCATEGORIZER_H
#define MONOALGORITHMS_CLUSTCATEGORIZER_H

namespace Mono {

struct ClustCategorizer {

  unsigned length;
  unsigned width;

  inline ClustCategorizer(const unsigned l,const unsigned w) {
    length = l;
    width = w;
  }
  inline bool operator ==(const ClustCategorizer &cat){
    bool res = length == cat.length;
    res = res && width == cat.width;
    return res;
  }

  inline bool operator <(const ClustCategorizer &cat)const {
    return length*width < cat.length*cat.width;
  }
  inline bool operator >(const ClustCategorizer &cat)const {
    return length*width > cat.length*cat.width;
  } 

};

} // end Mono namespace

#endif
