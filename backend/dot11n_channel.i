/*  :file: dot11n_channel.i
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-10-01 14:15:29 -0500 (Sat, 01 Oct 2011) $
 *  $LastChangedRevision: 5172 $
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

/* %module(docstring="Interface to IEEE 802.11n Channel") dot11n_channel */

%include <std_string.i>
%include <channel.i>

%feature("autodoc","1");
%{
#include <dot11n_channel_model.h>
#include <dot11n_channel.h>
%}

%include <dot11n_channel_model.h>

%inline %{
  std::string dot11n_channel_py_dot11n_channel__repr(Dot11N_Channel *H) {
    std::string s;
    s = H->to_string();
    return s;
  }

  Dot11N_Channel dot11n_channel_py_channel__add(Dot11N_Channel *A, Dot11N_Channel *B) {
    Dot11N_Channel H = (*A);
    for (int k=0; k< A->size(); k++) H(k) = (*A)(k) + (*B)(k);
    return H;
  }

  Dot11N_Channel dot11n_channel_py_channel__sub(Dot11N_Channel *A, Dot11N_Channel *B) {
    Dot11N_Channel H = (*A);
    for (int k=0; k< A->size(); k++) H(k) = (*A)(k) - (*B)(k);
    return H;
  }

  Dot11N_Channel dot11n_channel_py_channel__add(Dot11N_Channel *H, std::complex<double> x) {
    Dot11N_Channel R = (*H);
    for (int k=0; k<R.size(); k++)
      R(k) = (*H)(k) + x;
    return R;
  }

  Dot11N_Channel dot11n_channel_py_channel__sub(Dot11N_Channel *H, std::complex<double> x) {
    Dot11N_Channel R = (*H);
    for (int k=0; k<R.size(); k++)
      R(k) = (*H)(k) - x;
    return R;
  }

  Dot11N_Channel dot11n_channel_py_channel__mul(Dot11N_Channel *H, std::complex<double> t) {
    Dot11N_Channel R = (*H);
    R = mult(*H, t);
    return R;
  }
%}

class Dot11N_Channel : public Channel {
  public:
    Dot11N_Channel(int model=DOT11N_TGN_MODEL_A, int ntx=-1, int nrx=-1, int flag=DOT11N_TGN_DEFAULT);
    virtual ~Dot11N_Channel() {}
    bool setup(int model, int ntx, int nrx, int flag=DOT11N_TGN_DEFAULT);
    void set_nrx(int ntx,   bool update=true);
    void set_ntx(int nrx,   bool update=true);
    void set_flags(int flag, bool update=true);
    void set_model(int model,   bool update=true);
    void set_environmentspeed(double v, bool=true);
    inline std::string model() const;
    inline int ntx() const;
    inline int nrx() const;
    inline int maxdelay();
    inline double coherencetime() const;
    inline double samplefrequency() const;
    inline double environmentspeed() const;
    inline int flags() const;
    inline int nfft() const;
    void update(double tupdate=-1.0);      /* create next channel using set parameters */

    /* compute RMS delay spread of current channel */
    double rmsdelay();
    /* calculate expected RMS delay for model */
    static double calc_rmsdelay(int model, int ntx, int nrx);
    /* compute median power of subcarriers */
    double medianpower();
    inline double medianpowerdb();
};

%pythoncode %{
Dot11N_Channel.__repr__ = dot11n_channel_py_dot11n_channel__repr
Dot11N_Channel.__add__  = dot11n_channel_py_channel__add
Dot11N_Channel.__radd__ = dot11n_channel_py_channel__add
Dot11N_Channel.__sub__  = dot11n_channel_py_channel__sub
Dot11N_Channel.__rsub__ = dot11n_channel_py_channel__sub
Dot11N_Channel.__mul__  = dot11n_channel_py_channel__mul
Dot11N_Channel.__rmul__ = dot11n_channel_py_channel__mul
%}
