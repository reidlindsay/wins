/*  :file: dot11n_receiver.i
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

/* %module(docstring="Interface to IEEE 802.11n Receiver") dot11n_receiver */

%include <std_string.i>

%feature("autodoc","1");
%{
#include <dot11n_receiver.h>
%}

%inline %{
  std::string dot11n_receiver_py_dot11n_receiver__repr(Dot11N_Receiver *w) {
    std::string s;
    s = w->to_string();
    return s;
  }
%}

class Dot11N_Receiver {
  public:
    Dot11N_Receiver(double=20e6);
    virtual ~Dot11N_Receiver() {}
    bool decode_header(const Waveform & w, bool checkcrc=false);
    std::string decode_data(const Waveform & w);
    inline int mcs();
    inline int length();

    inline double sampling_frequency() const;
    inline void set_sampling_frequency(const double f);
    inline double get_avg_snr();

    inline unsigned int start_index() const;
    inline void   set_start_index(const unsigned int s);
    inline double get_cfo();

    inline void enable_cfo_correction();
    inline void disable_cfo_correction();
    inline bool using_cfo_correction() const;
};

%pythoncode %{
Dot11N_Receiver.__repr__ = dot11n_receiver_py_dot11n_receiver__repr
%}
