#ifndef dsp
#define dsp

long double des = 0.0, desf, Da[4] = {0.0,0.0,0.0,0.0}, DFa[4] = {0.0,0.0,0.0,0.0}, sig1 = 0.0, sig1s = 0.0, sig1av = 0.0, sig2 = 0.0, sig2s = 0.0, sig2av = 0.0;
long double LATa[3] = {0.0,0.0,0.0}, LATFa[3] = {0.0,0.0,0.0}, lat = 0.0, latf = 0.0, lato = 0.0, latn = 0.0, latff = 0.0, lats = 0.0, alatf = 0.0, alat = 0.0, dlatmax = 0.0, dlatn = 0.0, dlato = 0.0;
long double Ea[4] = {0.0,0.0,0.0,0.0}, EFa[4] = {0.0,0.0,0.0,0.0}, e = 0.0, ef = 0.0, ec = 1.0;
double NORa[3] = {0.0,0.0,0.0}, NORFa[3] = {0.0,0.0,0.0}, nor = 0.0, norf = 0.0;
int mukp = 0, mukf = 0, muko, mukhp = 0, mukhf = 0, muklp = 0, muklf = 0;

double muka[1000], mukav = 0, mukha[1000], mukhav, mukh, mukl, mukpr = 0;

double POSa[3] = {0.0,0.0,0.0}, POSFa[3] = {0.0,0.0,0.0}, pos = 0.0, posf = 0.0;
double vel;

long double latfa[3] = {1.0, -1.9111970674260727598436915286584, 0.91497583480143340750601055333391};
long double latfb[3] = {0.00094469184384016191557975616888143, 0.0018893836876803238311595123377629, 0.00094469184384016191557975616888143};

long double efa[4] = {1.0, -0.30816050209602025011790260577982, -0.45803063578866187732785419939319, -0.23380886211531795582097004171374};
long double efb[4] = {11.401839036186085962754077627324, 12.602996150186724477748612116557, 2.4933422600854204809195380221354, 1.2921851460847815218357936828397};

long double dfa[3] = {0.3137920502098252018008395225479, 0.62758410041965040360167904509581, 0.3137920502098252018008395225479};
long double dfb[3] = {1.0, -0.078944387242913499624918927111139, 0.33411258808221450111730632670515};

double norfa[3] = {1.0, -1.9111970674260727598, 0.91497583480143340751};
double norfb[3] = {0.00094469184384016191558, 0.0018893836876803238312, 0.00094469184384016191558};

double posfa[4] = {1.0, -2.58186, 2.24667, -0.657275};
double posfb[4] = {0.000941314, 0.00282394, 0.00282394, 0.000941314};


long double lat_filter();
long double e_filter();
double nor_filter();
double pos_filter();
void sample_forces();
void find_avg();

#endif