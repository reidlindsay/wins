/*  :file: dot11n_receiver.h
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-07-15 16:49:08 -0500 (Fri, 15 Jul 2011) $
 *  $LastChangedRevision: 5057 $
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

#ifndef INCLUDED_DOT11N_RECEIVER_H
#define INCLUDED_DOT11N_RECEIVER_H

#include <wins-80211n/header_decoding_11n.h>
#include <wins-80211n/data_decoding_11n.h>
#include <wins-80211n/rxvector.h>
#include <wins-80211n/waveform.h>
#include <string>

/*  :class: Dot11N_Receiver
 *  :description: Interface to IEEE 802.11n Receiver.
 */
class Dot11N_Receiver {
  private:
    Header_Decoding_11n d_hed;
    Data_Decoding_11n   d_rdp;
    rxvector            d_rx_param;
    unsigned int d_start_index;
    unsigned int d_symbols_left;
    double       d_sampling_frequency;
    bool         d_use_cfo_correction;

  public:
    Dot11N_Receiver(double=20e6);
    virtual ~Dot11N_Receiver() {}
    bool decode_header(const Waveform & w, bool checkcrc=false);
    std::string decode_data(const Waveform & w);
    std::string to_string();

    inline int mcs() { return d_rx_param.get_MCS(); }
    inline int length() { return d_rx_param.get_LENGTH(); }
    inline double sampling_frequency() const { return d_sampling_frequency; }
    inline double get_avg_snr() { return 10*log10(d_hed.get_avg_snr()); }

    inline void set_sampling_frequency(const double f) { d_sampling_frequency = f; }

    inline unsigned int start_index() const { return d_start_index; }
    inline void   set_start_index(const unsigned int s) { d_start_index = s; }
    inline double get_cfo() { return d_hed.get_coarse_cfo() + d_hed.get_fine_cfo(); }

    inline void enable_cfo_correction()  { d_use_cfo_correction = true; }
    inline void disable_cfo_correction() { d_use_cfo_correction = false; }
    inline bool using_cfo_correction() const { return d_use_cfo_correction; }
};

#endif  /* INCLUDED_DOT11N_RECEIVER_H */
