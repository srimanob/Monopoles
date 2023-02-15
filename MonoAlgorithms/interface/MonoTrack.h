#ifndef Monopoles_MonoAlgorithms_MonoTrack_h
#define Monopoles_MonoAlgorithms_MonoTrack_h

//////////////////
// class of monopole tracks
/////////////////

namespace Mono {

class MonoTrack {

public:
  inline MonoTrack()
   :m_xyp0(0.),m_xyp1(0.),m_xyp2(0.)
   ,m_rzp0(0.),m_rzp1(0.),m_rzp2(0.)
   { }

  inline MonoTrack(double xyp0, double xyp1, double xyp2
    ,double rzp0, double rzp1, double rzp2)
    :m_xyp0(xyp0),m_xyp1(xyp1),m_xyp2(xyp2)
    ,m_rzp0(rzp0),m_rzp1(rzp1),m_rzp2(rzp2)
    { }

  inline virtual ~MonoTrack() { }

  // -------- accessor methods ---------------
  inline double xyp0() const { return m_xyp0; }
  inline double xyp1() const { return m_xyp1; }
  inline double xyp2() const { return m_xyp2; }

  inline double rzp0() const { return m_rzp0; }
  inline double rzp1() const { return m_rzp1; }
  inline double rzp2() const { return m_rzp2; }

private:
  double m_xyp0;
  double m_xyp1;
  double m_xyp2;
  
  double m_rzp0;
  double m_rzp1;
  double m_rzp2;

};

} // end mono namespace

#endif
