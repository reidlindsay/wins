/*  :file: dot11n_transmitter.i
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

/* %module(docstring="Interface to IEEE 802.11n Transmitter") dot11n_transmitter */

%include <std_string.i>

%feature("autodoc","1");
%{
#include <dot11n_transmitter.h>
%}

%inline %{
  std::string dot11n_transmitter_py_dot11n_transmitter__repr(Dot11N_Transmitter *w) {
    std::string s;
    s = w->to_string();
    return s;
  }
%}

class Dot11N_Transmitter {
  public:
    Dot11N_Transmitter();
    virtual ~Dot11N_Transmitter() {}
    Waveform encode(const std::string & p="");

    void reset();
    static bool is_ntx_valid(int);
    static bool is_rate_valid(int rate, int ntx);
    static inline bool is_mcs_valid(int mcs, int ntx);

    void set_ntx(int);
    void set_rate(int);
    inline void set_mcs(int);
    void set_channel_smoothing_on();
    void set_channel_smoothing_off();
    void set_full_estimation_on();
    void set_full_estimation_off();

    inline int  get_ntx();
    inline int  get_rate();
    inline int  get_mcs();
    inline bool is_channel_smoothing_on();
    inline bool is_channel_smoothing_off();
    inline bool is_full_estimation_on();
    inline bool is_full_estimation_off();
};

%pythoncode %{
Dot11N_Transmitter.__repr__ = dot11n_transmitter_py_dot11n_transmitter__repr
%}
