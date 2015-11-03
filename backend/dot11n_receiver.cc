/*  :file: dot11n_receiver.cc
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-07-14 16:12:04 -0500 (Thu, 14 Jul 2011) $
 *  $LastChangedRevision: 5055 $
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

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include <dot11n_receiver.h>
#include <wins-80211n/itpp.h>
#include <iostream>
#include <sstream>

/* MACROS and DEFINITIONS */
#define DOT11N_CORRELATION_THRESHOLD  0.75
//#define DOT11N_COARSE_WINDOW_SIZE     48
//#define DOT11N_EXACT_WINDOW_SIZE      24

Dot11N_Receiver::Dot11N_Receiver(double f)
    : d_hed(), d_rdp(), d_start_index(0), d_symbols_left(0), 
      d_sampling_frequency(f), d_use_cfo_correction(true)
{
}

bool
Dot11N_Receiver::decode_header(const Waveform & w, bool checkcrc)
{
  double correlation_threshold = DOT11N_CORRELATION_THRESHOLD;
  d_hed.set_sampling_frequency(sampling_frequency()); // set to 20 MHz
  d_hed.set_correlation_threshold(correlation_threshold);
  d_hed.set_timing_tolerance(5);

#ifdef DOT11N_COARSE_WINDOW_SIZE
  unsigned int coarse_window   = DOT11N_COARSE_WINDOW_SIZE;
  d_hed.set_coarse_packet_window(coarse_window);
#endif
  
#ifdef DOT11N_EXACT_WINDOW_SIZE
  unsigned int exact_window   = DOT11N_EXACT_WINDOW_SIZE;
  d_hed.set_exact_packet_window(exact_window);
#endif
  /* Set up CFO correction parameter */
  if (using_cfo_correction()) {
    d_hed.enable_cfo_correction();
  } else {
    d_hed.disable_cfo_correction();
  }

  //unsigned int start_index=0;
  //unsigned int symbols_left=0;
  d_hed.egc_decode(w,d_start_index,d_symbols_left,d_rx_param);
#if 0
  printf("Dot11N_Receiver::decode_header() called \n");
  printf("--> packet detected? = %s\n", (d_hed.is_packet_detected()) ? "true":"false");
  printf("--> waveform size    = %d\n", w.cols() );
  printf("--> start index      = %d\n", d_start_index);
  printf("--> symbols left     = %d\n", d_symbols_left);
  if (d_hed.is_packet_detected()) {
    double cfo = d_hed.get_coarse_cfo() + d_hed.get_fine_cfo();
    printf("--> Coarse CFO       = %.2f kHz \n", d_hed.get_coarse_cfo()*1e-3 );
    printf("--> Fine CFO         = %.2f kHz \n", d_hed.get_fine_cfo()*1e-3 );
    printf("--> Total CFO        = %.2f kHz \n", cfo*1e-3);
    printf("--> Average SNR      = %.2f dB\n", 10*log10(d_hed.get_avg_snr()) );
    printf("--> Length           = %d Bytes \n", d_rx_param.get_LENGTH() );
    printf("--> MCS              = %d \n", d_rx_param.get_MCS() );
    printf("--> CRC OK?          = %s \n", (d_hed.is_crc_ok()) ? "true": "false" );
  }
#endif
  bool rval = d_hed.is_packet_detected();
  if (checkcrc) rval = rval && d_hed.is_crc_ok();
  return rval;
}

std::string
Dot11N_Receiver::decode_data(const Waveform & w)
{
  itpp::bvec output_data_bits;

  // shut off rx timing
  d_rdp.trigger_rx_timing(false);
  /* Set up CFO correction parameter */
  if (using_cfo_correction()) {
    d_rdp.enable_cfo_correction();
  } else {
    d_rdp.disable_cfo_correction();
  }

  // use the d_rx_param computed from previous call to decode_header()
  double cfo = d_hed.get_coarse_cfo() + d_hed.get_fine_cfo();
  d_rdp.set_frequency_offset(cfo, sampling_frequency());
  unsigned int nrx = w.rows();
  if (w.cols()>=d_start_index+d_symbols_left) {
    // check size of input before passing data to recover_bits()
    d_rdp.recover_bits(d_rx_param,
                       w.get(0, nrx-1, d_start_index, d_symbols_left+d_start_index-1),
                       d_hed.get_avg_snr(),
                       output_data_bits);
  } else {
    output_data_bits.clear();
    output_data_bits.set_size(0);
  }
  std::string data = to_bytes(output_data_bits);
  return data;
}


std::string
Dot11N_Receiver::to_string()
{
  std::ostringstream s;
  s << "Dot11N_Receiver";

  std::string rstr = s.str();
  return rstr;
}
