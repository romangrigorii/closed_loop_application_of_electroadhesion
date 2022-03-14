#include <xc.h>          // Load the proper header for the processor
#include "constants.h"
#include "interrupts.h"
#include "control.h"
#include "NU32.h"
#include <math.h>
#include "dsp.h"
#include "helpers.h"

void wait(t){ // = 8 cycles + t cycles
  int ii;
  for (ii=0;ii<t;ii++){
  }
}

void array_insert(int type){
  switch(type){
    case 0:
    LATa[2] = LATa[1];
    LATa[1] = LATa[0];
    LATa[0] = lat;
    break;
    case 1:
    LATFa[2] = LATFa[1];
    LATFa[1] = LATFa[0];
    LATFa[0] = latf;
    break;
    case 2:
    NORa[2] = NORa[1];
    NORa[1] = NORa[0];
    NORa[0] = nor;
    break;
    case 3:
    NORFa[2] = NORFa[1];
    NORFa[1] = NORFa[0];
    NORFa[0] = norf;
    break;
    case 4:
    POSa[3] = POSa[2];
    POSa[2] = POSa[1];
    POSa[1] = POSa[0];
    POSa[0] = pos;
    break;
    case 5:
    POSFa[3] = POSFa[2];
    POSFa[2] = POSFa[1];
    POSFa[1] = POSFa[0];
    POSFa[0] = posf;
    break;
    case 6:
    Ea[3] = Ea[2];
    Ea[2] = Ea[1];
    Ea[1] = Ea[0];
    Ea[0] = e;
    break;
    case 7:
    EFa[3] = EFa[2];
    EFa[2] = EFa[1];
    EFa[1] = EFa[0];
    EFa[0] = ef;
    break;
    case 8:
    Da[2] = Da[1];
    Da[1] = Da[0];
    Da[0] = des;
    break;
    case 9:
    DFa[2] = DFa[1];
    DFa[1] = DFa[0];
    DFa[0] = desf;
    break;
  }
}

double sq(double t, double f){
  double sig = sin(2*pi*t*f);
  if (sig>0){
    return 1.0;
  } else {
    return -1.0;
  }
}

int wrap(int p, int p_lim){
  if(p >= p_lim){
    p -= p_lim;
  } else if (p < 0){
    p += p_lim;
  }
  return p;
}

long double absd(long double el){
  if (el<0.0){
    el = -el;
  }
  return el;
}
