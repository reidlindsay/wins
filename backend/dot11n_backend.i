/*  :file: dot11n_backend.i
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-06-01 15:50:22 -0500 (Wed, 01 Jun 2011) $
 *  $LastChangedRevision: 4991 $
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

%module(docstring="Interface to IEEE 802.11n Physical Layer") dot11n_backend

%{
#include <itpp/base/random.h>
%}

%include <dot11n_transmitter.i>
%include <dot11n_receiver.i>
%include <dot11n_channel.i>

%inline %{
  void RNG_init(int seed=-1) {
    if (seed<0) {
      itpp::RNG_randomize();
    } else {
      itpp::RNG_reset((unsigned int) seed);
    }
  }
%}
