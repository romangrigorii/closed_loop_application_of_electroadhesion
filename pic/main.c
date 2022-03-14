#include <xc.h>          // Load the proper header for the processor
#include "constants.h"
#include "interrupts.h"
#include "control.h"
#include "NU32.h"
#include "helpers.h"
#include <math.h>
#include <stdio.h>
#include "dsp.h"

int main(void) {
  int i;  // generic variable
  __builtin_disable_interrupts();
  NU32_Startup();
  digital_init();
  interrupt_init();
  SPI_com_init();
  NU32_ReadUART3(message,100);
  sprintf(message,"%s\n\r","CONNECTED");
  NU32_WriteUART3(message);
  __builtin_enable_interrupts();
  avc = 0;
  init_count = 0;
  SYN = 0;
  CCD = 0;
  CNVST = 1;
  CH = 0;
  while (1){
    NU32_ReadUART3(message,100);
    sscanf(message,"%c",&option);
    sprintf(message,"%c\r\n", option);
    NU32_WriteUART3(message);
    switch (option){
      case 'r':
      pos_vec_i = 0;
      t = 0;
      break;
      case 's':
      for(i=0;i<2000;i++){
        sprintf(message,"%lf\t%lf\r\n",pos_vec[i],vel_vec[i]);
        NU32_WriteUART3(message);
      }
      break;
      case 'p':
      sprintf(message,"%Lf\n\r",des);
      NU32_WriteUART3(message);
      break;
      case 'z':
      // zeros the force readings
      avc = 0;
      break;
      case 'f':
      filter_on = !filter_on;
      break;
      case 'e':
      feedback_enable = !feedback_enable;
      fe[0] = 0;
      ie = 0;
      de = 0;
      break;
      case 'v':
      NU32_ReadUART3(message,100);
      sscanf(message,"%lf",&freq);
      sprintf(message,"%lf\n\r",freq);
      NU32_WriteUART3(message);
      break;
    }

  }
}
