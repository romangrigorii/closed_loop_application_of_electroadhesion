#include "dsp.h"

long double lat_filter(){
  return (latfb[0]*LATa[0] + latfb[1]*LATa[1] + latfb[2]*LATa[2] - latfa[1]*LATFa[0] - latfa[2]*LATFa[1]);
}
long double e_filter(){
  return (efb[0]*Ea[0] + efb[1]*Ea[1] + efb[2]*Ea[2] + efb[3]*Ea[3] - efa[1]*EFa[0] - efa[2]*EFa[1] - efa[3]*EFa[2]);
}
double nor_filter(){
  return (norfb[0]*NORa[0] + norfb[1]*NORa[1] + norfb[2]*NORa[2] - norfa[1]*NORFa[0] - norfa[2]*NORFa[1]);
}
double pos_filter(){
  return (posfb[0]*POSa[0] + posfb[1]*POSa[1] + posfb[2]*POSa[2] + posfb[3]*POSa[3] - posfa[1]*POSFa[0] - posfa[2]*POSFa[1] - posfa[3]*POSFa[2]);
}

long double d_filter(){
  return (dfb[0]*Da[0] + dfb[1]*Da[1] + dfb[2]*Da[2] + dfb[3]*Da[3] - dfa[1]*DFa[0] - dfa[2]*DFa[1] - dfa[3]*DFa[2]);
}
