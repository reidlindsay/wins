/*  :file: dot11n_channel_model.h
 *         Enumeration of channel models used by Dot11N_Channel.
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-06-17 14:29:04 -0500 (Fri, 17 Jun 2011) $
 *  $LastChangedRevision: 5004 $
 *
 *  :author: Ketan Mandke <kmandke@mail.utexas.edu>
 *  
 *  :copyright:
 *    Copyright 2009-2010 The University of Texas at Austin
 *    
 *    Licensed under the Apache License, Version 2.0 (the "License");
 *    you may not use this file except in compliance with the License.
 *    You may obtain a copy of the License at
 *    
 *      http://www.apache.org/licenses/LICENSE-2.0
 *    
 *    Unless required by applicable law or agreed to in writing, software
 *    distributed under the License is distributed on an "AS IS" BASIS,
 *    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *    See the License for the specific language governing permissions and
 *    limitations under the License.
 */

#ifndef INCLUDED_DOT11N_CHANNEL_MODEL_H
#define INCLUDED_DOT11N_CHANNEL_MODEL_H

#define MODELA_KDB  "0"
#define MODELA_AOA  "45"
#define MODELA_ASA  "40"
#define MODELA_AOD  "45"
#define MODELA_ASD  "40"
#define MODELA_PWR  "0"
#define MODELA_DLY  "0"

#define MODELB_KDB  "0"
#define MODELB_AOA  "4.3    118.4"
#define MODELB_ASA  "14.4   25.2"
#define MODELB_AOD  "225.1  106.5"
#define MODELB_ASD  "14.4   25.4"
#define MODELB_PWR  "0    -5.4  -10.8 -16.2 -21.7 -inf  -inf  -inf  -inf; \
                     -inf -inf  -3.2  -6.3  -9.4  -12.5 -15.6 -18.7 -21.8"
#define MODELB_DLY  "0    10    20    30    40    50    60    70    80"

#define MODELC_KDB  "0"
#define MODELC_AOA  "290.3  332.3"
#define MODELC_ASA  "24.6   22.4"
#define MODELC_AOD  "13.5   56.4"
#define MODELC_ASD  "24.7   22.5"
#define MODELC_PWR  "0    -2.1  -4.3  -6.5  -8.6  -10.8 -13.0 -15.2 -17.3 -19.5 -inf  -inf  -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -5.0  -7.2  -9.3  -11.5 -13.7 -15.8 -18.0 -20.2"
#define MODELC_DLY  "0    10    20    30    40    50    60    70    80    90    110   140   170   200"

#define MODELD_KDB  "3"
#define MODELD_AOA  "158.9  320.2 276.1"
#define MODELD_ASA  "27.7   31.4  37.4"
#define MODELD_AOD  "332.1  49.3  275.9"
#define MODELD_ASD  "27.4   32.1  36.8"
#define MODELD_PWR  "0    -0.9  -1.7  -2.6  -3.5  -4.3  -5.2  -6.1  -6.9  -7.8  -9.0  -11.1 -13.7 -16.3 -19.3 -23.2 -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -6.6  -9.5  -12.1 -14.7 -17.4 -21.9 -25.5 -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -18.8 -23.2 -25.2 -26.7"
#define MODELD_DLY  "0    10    20    30    40    50    60    70    80    90    110   140   170   200   240   290   340   390"

#define MODELE_KDB  "6"
#define MODELE_AOA  "163.7  251.8 80.0  182.0"
#define MODELE_ASA  "35.8   41.6  37.4  40.3"
#define MODELE_AOD  "105.6  293.1 61.9  275.7"
#define MODELE_ASD  "36.1   42.5  38.0  38.7"
#define MODELE_PWR  "-2.6 -3.0  -3.5  -3.9  -4.5  -5.6  -6.9  -8.2  -9.8  -11.7 -13.9 -16.1 -18.3 -20.5 -22.9 -inf  -inf  -inf; \
                     -inf -inf  -inf  -inf  -1.8  -3.2  -4.5  -5.8  -7.1  -9.9  -10.3 -14.3 -14.7 -18.7 -19.9 -22.4 -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -7.9  -9.6  -14.2 -13.8 -18.6 -18.1 -22.8 -inf  -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -20.6 -20.5 -20.7 -24.6"
#define MODELE_DLY  "0    10    20    30    50    80    110   140   180   230   280   330   380   430   490   560   640   730"


#define MODELF_KDB  "6"
#define MODELF_AOA  "315.1  180.4 74.7  251.5  68.5  246.2"
#define MODELF_ASA  "48.0   55.0  42.0  28.6   30.7  38.2"
#define MODELF_AOD  "56.2   183.7 153.0 112.5  291.0 62.3"
#define MODELF_ASD  "41.6   55.2  47.4  27.2   33.0  38.0"
#define MODELF_PWR  "-3.3 -3.6  -3.9  -4.2  -4.6  -5.3  -6.2  -7.1  -8.2  -9.5  -11.0 -12.5 -14.3 -16.7 -19.9 -inf  -inf  -inf; \
                     -inf -inf  -inf  -inf  -1.8  -2.8  -3.5  -4.4  -5.3  -7.4  -7.0  -10.3 -10.4 -13.8 -15.7 -19.9 -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -5.7  -6.7  -10.4 -9.6  -14.1 -12.7 -18.5 -inf  -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -8.8  -13.3 -18.7 -inf  -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -12.9 -14.2 -inf  -inf; \
                     -inf -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -inf  -16.3 -21.2"
#define MODELF_DLY  "0    10    20    30    50    80    110   140   180   230   280   330   400   490   600   730   880   1050"



enum {
  DOT11N_LOS_MODEL,
  DOT11N_TGN_MODEL_A,
  DOT11N_TGN_MODEL_B,
  DOT11N_TGN_MODEL_C,
  DOT11N_TGN_MODEL_D,
  DOT11N_TGN_MODEL_E,
  DOT11N_TGN_MODEL_F
};

/* Flags for configuring Dot11N_Channel */
enum {
  DOT11N_USE_LOS=(1<<0),                    /* LOS component on */
  DOT11N_USE_NLOS=(1<<1),                   /* NLOS component on */
  DOT11N_LOS_RICE=(1<<2),                   /* LOS  - Rice steering matrix [default=identity] */
  DOT11N_NLOS_FADING=(1<<3),                /* NLOS - Rayleigh fading [default=identity] */
  DOT11N_NLOS_EQUALGAIN=(1<<4),             /* NLOS - Equal gain (random phase) */
  DOT11N_NLOS_DOPPLER=(1<<5),               /* NLOS - Doppler filter Rayleigh fading */
  DOT11N_TGN_DEFAULT=(1<<0|1<<1|1<<2|1<<3)
};

#endif  /* INCLUDED_DOT11N_CHANNEL_MODEL_H */
