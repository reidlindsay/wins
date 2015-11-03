/*  :file: dot11n_transmitter.h
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-03-24 11:29:29 -0500 (Thu, 24 Mar 2011) $
 *  $LastChangedRevision: 4942 $
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

#ifndef INCLUDED_DOT11N_TRANSMITTER_H
#define INCLUDED_DOT11N_TRANSMITTER_H

#include <wins-80211n/transmission_interface_11n.h>
#include <wins-80211n/waveform.h>
#include <wins-80211n/txvector.h>
#include <string>

/*  :class: Dot11N_Transmitter
 *  :description: Interface to IEEE 802.11n Transmitter.
 */
class Dot11N_Transmitter {
  private:
    Transmission_Interface_11n d_txi;
    txvector d_txvector;
    /* parameters for txvector */
    int  d_ntx;
    int  d_rate;
    bool d_is_channel_smoothing_on;
    bool d_is_full_estimation_on;

  public:
    Dot11N_Transmitter();
    virtual ~Dot11N_Transmitter() {}
    Waveform encode(const std::string & p="");
    std::string to_string();

    void reset();
    static bool is_ntx_valid(int);
    static bool is_rate_valid(int rate, int ntx);
    static inline bool is_mcs_valid(int mcs, int ntx) { return Dot11N_Transmitter::is_rate_valid(mcs, ntx); }

    void set_ntx(int);
    void set_rate(int);
    inline void set_mcs(int m) { set_rate(m);}
    void set_channel_smoothing_on();
    void set_channel_smoothing_off();
    void set_full_estimation_on();
    void set_full_estimation_off();

    inline int  get_ntx() { return d_ntx;}
    inline int  get_rate() { return d_rate;}
    inline int  get_mcs() { return get_rate(); }
    inline bool is_channel_smoothing_on()  { return d_is_channel_smoothing_on;}
    inline bool is_channel_smoothing_off() { return !is_channel_smoothing_on();}
    inline bool is_full_estimation_on()    { return d_is_full_estimation_on;}
    inline bool is_full_estimation_off()   { return !is_full_estimation_on();}
};

#endif  /* INCLUDED_DOT11N_TRANSMITTER_H */
