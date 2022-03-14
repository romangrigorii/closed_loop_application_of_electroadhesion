#include <xc.h>          // Load the proper header for the processor
#include "constants.h"
#include "interrupts.h"
#include "control.h"
#include "NU32.h"
#include <math.h>
#include "helpers.h"
#include "dsp.h"

double controller(){

  lato = latn;
  latn = alatf;
  dlato = dlatn;
  dlatn = (latn - lato);
  if (latn <= 0.01){ // detects change in swipe direction
    c_flag = 1;
  }
  if (c_flag==1){
    if (dlatn > dlatmax){
      dlatmax = dlatn;
    }
    if (dlatn<(dlatmax*0.5)){ // detects slip
      dlatmax = 0.0;
      c_flag = 2;
      startt = 0.0;
      muk = 0.0;
      mukp = 0;
    }
  }
  if (c_flag==2){
    if (startt<=0.01){
      startt = startt + .0001;
    }
    else if (startt<=0.02){
      startt = startt + .0001;
      muk = muk + latn;
    }
    else {
      c_flag = 3;
      muk = muk/100.0;
      //sprintf(message, "%lf\n\r", muk);
      //NU32_WriteUART3(message);
      e = (muk + des - alat);
      mukp = 1;
      mukf = 0;
      muka[0] = muk;
      mukav = muk;

      mukha[0] = 0.03;
      mukhp = 1;
      mukh = 0.0;
      mukhf = 0;
      mukhav = .03;
      mukpr = 0;

      array_insert(6);
      array_insert(6);
      array_insert(6);
      ef = 0.0;
      array_insert(7);
      array_insert(7);
      array_insert(7);
      return 0.0;
    }
  }
  if (c_flag == 3){
    //if (latn>(mukav*0.6)){
      return feedback_compute();
    //}
  }
  return 0.0;
}

double feedback_compute(){
  /*
  if (absd(ef)<= 0.025){
    if (mukf==0){
      muka[mukp] = muka[0];
    }
    mukav = mukav + ((alat - muka[mukp])/100.0);
    muka[mukp] = alat;
    mukp++;
    if (mukp==100)
    {
      mukp = 0;
      if (mukf == 0){
        mukf = 1;
      }
       sprintf(message,"*%lf\n\r", mukav);
      //NU32_WriteUART3(message);
    }
    muk = mukav;
  }
  */

  e = (muk + des - alat);
  array_insert(6);
  ef = e_filter();
  if (ef>5.0){
    ef = 5.0;
  }
  if (ef<-5.0){
    ef = -5.0;
  }
  //if (absd(ef)<= 0.0153333){
    //ef = 0.0;
  //}
//  if (absd(ef)<= 0.025){
//    ef = 0.0;
//  }
  array_insert(7);
  return ef;
}

double open_loop(){
  return des*24.661;
}
