#ifndef interrupts
#define interrupts

#define CS1 LATBbits.LATB10
#define CS2 LATBbits.LATB11
#define SYN LATBbits.LATB12
#define CCD LATBbits.LATB9
#define CH LATBbits.LATB3
#define CNVST LATBbits.LATB4

void chip_read_data();
void chip_write_data(int);
int twocompconv(int);
void digital_init();
void SPI_com_init();
void interrupt_init();
void find_pos();
void find_vel();
int twocompconv(int);
void sample_lat_des();
void sample_lat_nor();
void output_compute();

int p1 = 0, p2 = 0, p3 = 0, t1 = 0, t2 = 0, t3 = 0, t4 = 0, init_count = 1001, th = 0, print_r = 0, pos_vec_i = 2000, wr = 0, avn = 10000, avc = 10000, filter_on = 0, chip_out = 0, sig = 0, feedback_enable = 0;
char option;
double time_vec[2000], pos_vec[2000], vel_vec[2000], st_t = 0, end_t = 0, t_counts = 40000000,t = 0.0, dt = 0;
double TtoP = 0.00785567;
double VtoFlat = 0.25, VtoFnor = -.1526;
long double adcbits = 16383, dacbits = 65535;
double value;
int counter = 1;
long double start_lat;
double sn = 1.0;

#endif
