/*  :file: dot11n_transmitter.cc
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

#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>

#include <dot11n_transmitter.h>
#include <wins-80211n/itpp.h>
#include <iostream>
#include <sstream>

Dot11N_Transmitter::Dot11N_Transmitter()
    : d_txi(), d_txvector(), d_ntx(1), d_rate(0), 
      d_is_channel_smoothing_on(false), d_is_full_estimation_on (true)
{
  reset();
}

Waveform
Dot11N_Transmitter::encode(const std::string & p)
{
  int plen = p.length();
  itpp::bvec inbits;
  itpp::cmat txwaveform;

  // set parameters needed for encoding
  to_bvec(p, inbits);
  d_txvector.set_LENGTH(plen);
  d_txvector.calculate_parameters(true); // calculate validity of TXVECTOR
  d_txi.format_frame(d_txvector, inbits, txwaveform);
  Waveform w(txwaveform);
#if 0
  std::cout << "len(p) = " << p.length() << std::endl;
  std::cout << "len(inbits) = " << inbits.length() << std::endl;
  std::cout << "txwaveform ~ " << txwaveform.rows() << " x " << txwaveform.cols() << std::endl;
#endif
  return w;
}

std::string
Dot11N_Transmitter::to_string()
{
  std::ostringstream s;
  s << "Dot11N_Transmitter";

  std::string rstr = s.str();
  return rstr;
}


void
Dot11N_Transmitter::reset()
{
  // set default parameters and txvector
  d_ntx = 1;
  d_rate = 0;
  d_is_channel_smoothing_on = false;
  d_is_full_estimation_on = true;
  // set up txvector
  d_txvector.set_TXPWR_LEVEL(1);
  d_txvector.set_ANTENNA_SET("1 0 0 0 0 0 0 0");
  d_txvector.set_N_TX(d_ntx);
  set_ntx(d_ntx);
  set_rate(d_rate);
  set_channel_smoothing_off();
  if (d_is_channel_smoothing_on) set_channel_smoothing_on();
}

bool
Dot11N_Transmitter::is_ntx_valid(int ntx)
{
  bool ntx_is_valid = true;
  if ((ntx<1) || (ntx>4) ) ntx_is_valid = false;
  return ntx_is_valid;
}

bool
Dot11N_Transmitter::is_rate_valid(int rate, int ntx)
{
  // check for valid rate
  bool rate_is_valid = true;
  if (!is_ntx_valid(ntx) ) rate_is_valid = false;
  if ((rate<0) || (rate>8*ntx) ) rate_is_valid = false;
  return rate_is_valid;
}

void
Dot11N_Transmitter::set_ntx(int ntx)
{
  // check for valid ntx
  assert(is_ntx_valid(ntx) );
  
  // set parameters for txvector
  int num_exten_ss = ntx - 1;
  if (is_full_estimation_off() ) num_exten_ss = 0;
  d_ntx = ntx;
  d_txvector.set_N_TX(ntx);
  d_txvector.set_NUM_EXTEN_SS(num_exten_ss);

  switch (ntx)
  {
    case 2:
      d_txvector.set_ANTENNA_SET("1 1 0 0 0 0 0 0");
      break;
    case 3:
      d_txvector.set_ANTENNA_SET("1 1 1 0 0 0 0 0");
      break;
    case 4:
      d_txvector.set_ANTENNA_SET("1 1 1 1 0 0 0 0");
      break;
    default:
      d_txvector.set_ANTENNA_SET("1 0 0 0 0 0 0 0");
  }
}

void
Dot11N_Transmitter::set_rate(int m)
{
  // check for valid rate
  assert(is_rate_valid(m, d_ntx) );

  // set parameters for txvector
  d_txvector.set_MCS(m);
}

void
Dot11N_Transmitter::set_channel_smoothing_on()
{
  d_is_channel_smoothing_on = true;
  d_txvector.set_SMOOTHING(SMOOTHING_REC);
}

void
Dot11N_Transmitter::set_channel_smoothing_off()
{
  d_is_channel_smoothing_on = false;
  d_txvector.set_SMOOTHING(SMOOTHING_NOT_REC);
}

void
Dot11N_Transmitter::set_full_estimation_on()
{
  d_is_full_estimation_on = true;
  set_ntx(d_ntx);
}

void
Dot11N_Transmitter::set_full_estimation_off()
{
  d_is_full_estimation_on = false;
  set_ntx(d_ntx);
}
