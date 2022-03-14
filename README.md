# closed_loop_application_of_electroadhesion

Tactile displays based on friction modulation offer wide-bandwidth forces rendered directly on the fingertip. However, due to a number of touch conditions (e.g., normal force, skin hydration) that result in variations in the friction force and the strength of modulation effect, the precision of the force rendering remains limited. In this paper we demonstrate a closed-loop electroadhesion method for precise playback of friction force profiles on a human finger and we apply this method to the tactile rendering of several textiles encountered in everyday life.

This repository contains all of the code used in writing an IEEE publication
titled "Closed loop application of electroadhesion for increased precision in
texture rendering" which can be found here:

DOI: 10.1109/TOH.2020.2972350

IEEE link: https://ieeexplore.ieee.org/document/8986537

/pic contains the C code used to program PIC32 microcontroller. Microcontroller
was responsible for:

1) communication with DAC/ADC chips
2) applying optimal closed loop friction control
3) applying closed loop fingertip state controller

/matlab contains code used to collect, process, and generate graphics of the
performance of the closed loop controller
