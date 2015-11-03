/*  :file: dot11n_channel.h
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-10-23 19:18:52 -0500 (Sun, 23 Oct 2011) $
 *  $LastChangedRevision: 5287 $
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

#ifndef INCLUDED_DOT11N_CHANNEL_H
#define INCLUDED_DOT11N_CHANNEL_H

#include <wins-80211n/transmission_interface_11n.h>
#include <wins-80211n/waveform.h>
#include <wins-80211n/txvector.h>
#include <wins-80211n/channel.h>
#include <dot11n_channel_model.h>
#include <string>

/*  :class: Dot11N_Channel
 *  :description: Interface to IEEE 802.11n Channel.
 */
class Dot11N_Channel : public Channel {
  private:
    int d_model;
    int d_ntx;
    int d_nrx;
    int d_flag;
    // parameters for doppler filter
    int d_nfft;           // FFT size needed for computing Rayleigh fading samples
    double d_tcorr;       // Coherence time
    double d_tsamp;       // Sample time
    double d_toffset;     // Current time offset into fading sample buffer
    double d_phi_last;    // Phase of last update on LOS component
    double d_env_speed;   // Environmental speed
    // other precomputed values
    Channel d_corr_sqrt;  // R = (Rtx x Rrx)^(1/2) (Kronecker product)
    Channel d_rice;       // Rice steering matrix (depends on AoA/AoD)
    Channel d_hs;         // [IID] Rayleigh fading samples used to calculate Hw
    itpp::mat d_doppler;  // Doppler filter

    /* internal helper methods */
    void calc_rice();                // initialize Rice steering matrix for LOS
    void calc_correlation();         // compute correlation matrix for NLOS
    void calc_doppler();             // initialize doppler filter
    void update_rayleigh(double, unsigned int); // update rayleigh component Hw
    itpp::cmat get_rayleigh(unsigned int);      // get rayleigh component Hw

    static void get_Kdb(double &, int);      // K-factor (in dB)
    static void get_Pdb(itpp::mat &, int);   // Power of taps (in dB)
    static void get_dly(itpp::ivec &, int);  // excess delay (in ns)
    static void get_aoa(itpp::vec &, int);   // AoA (in degrees)
    static void get_aod(itpp::vec &, int);   // AoD
    static void get_asa(itpp::vec &, int);   // AS (receiver)
    static void get_asd(itpp::vec &, int);   // AS (transmitter)

  public:
    Dot11N_Channel(int model=DOT11N_TGN_MODEL_A, int ntx=-1, int nrx=-1, int flag=DOT11N_TGN_DEFAULT);
    virtual ~Dot11N_Channel() {}
    bool setup(int model, int ntx, int nrx, int flag=DOT11N_TGN_DEFAULT);
    void set_nrx(int ntx,   bool update=true);
    void set_ntx(int nrx,   bool update=true);
    void set_flags(int flag, bool update=true);
    void set_model(int model,   bool update=true);
    void set_environmentspeed(double v, bool=true);
    inline std::string model() const { return model_string(); }
    inline int ntx() const { return d_ntx; }
    inline int nrx() const { return d_nrx; }
    inline int maxdelay() { return size()-1; }
    inline double coherencetime() const { return d_tcorr; }
    inline double samplefrequency() const { return (d_tsamp>0)?(1.0/d_tsamp):-1; }
    inline double environmentspeed() const { return d_env_speed; }
    inline int flags() const { return d_flag; }
    inline int nfft() const { return d_nfft; }
    void update(double tupdate=-1.0);

    Dot11N_Channel & operator=(const Channel &);

    /* compute RMS delay spread of current channel */
    double rmsdelay();
    /* calculate expected RMS delay for model */
    static double calc_rmsdelay(int model, int ntx, int nrx);
    /* compute median power of subcarriers */
    double medianpower();
    inline double medianpowerdb() { return itpp::dB(medianpower()); }

    /* internal method to update of channel */
    bool configure();
    std::string model_string() const;
};

#endif  /* INCLUDED_DOT11N_CHANNEL_H */
