#include <xc.h>          // Load the proper header for the processor
#include "constants.h"
#include "interrupts.h"
#include "control.h"
#include "NU32.h"
#include "helpers.h"
#include <math.h>
#include <stdio.h>
#include "dsp.h"

void digital_init(){

  TRISBbits.TRISB12 = 0; // indicator pin to synchronize the data
  TRISBbits.TRISB3 = 0; //  CH pin
  TRISBbits.TRISB10 = 0; // ADC spi
  TRISBbits.TRISB11 = 0; // DAC spi
  TRISBbits.TRISB4 = 0; //  CNVST pin

}

void SPI_com_init(){
  // setting up communication with friction control chip_write
  CS1 = 1;
  SPI3CON = 0;
  SPI3BRG = 10; // communication at 640kHz
  SPI3BUF;
  SPI3BUF;
  SPI3STATbits.SPIROV = 0;
  SPI3CONbits.MODE32 = 1;
  SPI3CONbits.MODE16 = 0;
  SPI3CONbits.MSTEN = 1;
  SPI3CONbits.CKE = 1;
  SPI3CONbits.CKP = 0;
  SPI3CONbits.SMP = 0;
  SPI3CONbits.ON = 1;

  // setting up communiction with DAC
  CS2 = 1;
  SPI4CON = 0;
  SPI4BRG = 10;
  SPI4BUF;
  SPI4BUF;
  SPI4STATbits.SPIROV = 0;
  SPI4CONbits.MODE32 = 0;
  SPI4CONbits.MODE16 = 1;
  SPI4CONbits.MSTEN = 1;
  SPI4CONbits.CKE = 1;
  SPI4CONbits.CKP = 1;
  SPI4CONbits.ON = 1;
}

void chip_write_data(int d_sig){
  CS2 = 0;
  SPI4BUF = d_sig;
  while(!SPI4STATbits.SPIRBF){
  }
  SPI4BUF;
  CS2 = 1;
}

int twocompconv(int s){
  if (s>=8192){
    s = s - 8192;
  } else {
    s = s + 8191;
  }
  return s;
}

void chip_read_data(){
  CNVST = 0;
  wait(0);
  CNVST = 1;
  wait(5);
  CS1 = 0;
  SPI3BUF = 0;
  while(!SPI3STATbits.SPIRBF){
  }
  sig = SPI3BUF;
  CS1 = 1;
  sig1 = (long double) twocompconv((sig>>4) & 0x3FFF);
  sig2 = (long double) twocompconv((sig>>18) & 0x3FFF);
}

void interrupt_init(){
  PR4 =  7962; // freq = 80,000,000/(1+3199) = 25kHz
  TMR4 = 0;
  T4CONbits.TCKPS = 0; // = 1
  T4CONbits.ON = 1;
  IPC4bits.T4IP = 5;
  IPC4bits.T4IS = 0;
  IFS0bits.T4IF = 0;
  IEC0bits.T4IE = 1;
}

void sample_lat_nor(){
  chip_read_data();
  lat = (sig1 - sig1av)/adcbits*VtoFlat*20.0;
  nor = (sig2 - sig2av)/adcbits*VtoFnor*20.0;
  array_insert(2);
  norf = nor_filter();
  array_insert(3);
  array_insert(0);
  latf = lat_filter();
  array_insert(1);
  alatf = absd(latf);
  alat = absd(lat);
}

void sample_lat_des(){
  chip_read_data();
  lat = (sig1 - sig1av)/adcbits*20.0*VtoFlat; // converts to lateral force in N
  des = (sig2 - 8192.0)/adcbits/5.0; // des input can span -10..10 V, which corresponds to +-.1N of desired force
  array_insert(0);
  latf = lat_filter();
  array_insert(1);
  array_insert(8);
  desf = d_filter();
  array_insert(9);
  alatf = absd(latf);
  alat = absd(lat);
}

void output_compute(){
  if (feedback_enable){
    chip_out = ((int) ((controller() + 5)*dacbits/10.0));
  } else {
    chip_out = ((int) ((open_loop() + 5)*dacbits/10.0));
  }
  if (chip_out<0){
    chip_out = 0;
  }
  if (chip_out > dacbits){
    chip_out = dacbits;
  }
}

void find_avg(){
  if (avc<avn){
    sig1s = sig1s + sig1;
    sig2s = sig2s + sig2;
    avc++;
    if (avc==avn){
      sig1av = sig1s/((long double) avn);
      sig2av = sig2s/((long double) avn);
      sig1s = 0.0;
      sig2s = 0.0;
    }
  }
}

void des_update(){
  if (counter == 250){
    des = 0.05*sn;
    sn = -sn;
    counter = 1;
  } else {
    counter++;
  }
}

void __ISR(_TIMER_4_VECTOR,IPL5SOFT) Timer4ISR(void){
  // friction sample, filter, and feedback routine
  t = t + .0001;
  sample_lat_des();
  find_avg();
  output_compute();
  chip_write_data(chip_out);
  IFS0bits.T4IF = 0;
}
