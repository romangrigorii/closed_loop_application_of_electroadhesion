#ifndef control
#define control

void error_compute();
double feedback_compute();
double controller();
double open_loop();

char message[100], option;

// friction control specific paramters

double fKp = 1.0, fKd = 0.0000000, fKi = 2.0, fe[3], de = 0.0, ie = 0.0, ielim = 0.1, startt = 0.0, muk = 0.0, freq = 40.0;
int c_flag = 0;

#endif
