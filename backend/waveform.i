/*  :file: waveform.i
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-06-08 16:41:49 -0500 (Wed, 08 Jun 2011) $
 *  $LastChangedRevision: 4998 $
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

%module(docstring="Wrapper for IT++ Complex Matrices") waveform

%include <std_string.i>
%include <std_vector.i>
%include <std_complex.i>

%feature("autodoc","1");
%{
#include <wins-80211n/waveform.h>
#include <itpp/base/random.h>
%}

%template(vectori) std::vector<int>;
%template(vectorc) std::vector<std::complex<double> >;

%inline %{
  std::string waveform_py_waveform__repr(Waveform *w) {
    std::string s;
    s = w->to_string();
    return s;
  }

  Waveform waveform_py_waveform__add(Waveform *a, Waveform *b) {
    Waveform w;
    int arows = a->rows(), acols = a->cols();
    int brows = b->rows(), bcols = b->cols();
    if (acols<bcols) {
      w = itpp::concat_horizontal(*a, itpp::zeros_c(arows, bcols-acols)) + *b;
    } else if (acols>bcols) {
      w = *a + itpp::concat_horizontal(*b, itpp::zeros_c(brows, acols-bcols) );
    } else {
      w = *a + *b;
    }
    return w;
  }

  Waveform waveform_py_waveform__add(Waveform *x, std::complex<double> a) {
    Waveform w;
    w = (*x) + a;
    return w;
  }

  Waveform waveform_py_waveform__mul(Waveform *a, Waveform *b) {
    Waveform w;
    w = itpp::elem_mult(*a,*b);
    return w;
  }

  Waveform waveform_py_waveform__mul(Waveform *x, std::complex<double> a) {
    Waveform w;
    w = (*x) * a;
    return w;
  }

  std::complex<double> waveform_py_waveform__get(Waveform *x, int r, int c) {
    return (*x)(r,c);
  }

  std::complex<double> waveform_py_waveform__get(Waveform *x, const std::vector<int> & idx) {
    return x->get(idx[0], idx[1]);
  }

  void waveform_py_waveform__set(Waveform *x, const std::vector<int> & idx, const std::complex<double> & v) {
    x->set(idx[0], idx[1],v);
  }

  Waveform waveform_py_waveform__get(Waveform *x, int r1, int r2,  int c1, int c2) {
    return (*x)(r1,r2,c1,c2);
  }

  Waveform waveform_py_waveform__get_cols(Waveform *x, int c1, int c2) {
    return x->get_cols(c1,c2);
  }

  Waveform waveform_py_waveform__concat_horizontal(Waveform *x, Waveform *y) {
    Waveform w = itpp::concat_horizontal(*x, *y);
    return w;
  }

  Waveform waveform_py_waveform__concat_vertical(Waveform *x, Waveform *y) {
    Waveform w = itpp::concat_vertical(*x, *y);
    return w;
  }
  
  void waveform_py_waveform__apply_offset(Waveform *w, double fco, double fs=1.0) {
    int ncols = w->cols();
    double epsilon = fco/fs;
    std::complex<double> iorj(0,1);
    itpp::cmat phi = itpp::exp(iorj*(2*epsilon*itpp::pi*itpp::linspace(0,ncols-1,ncols) ) );
    phi = itpp::repmat(phi.transpose(), w->rows(), 1);
    w->set_rows(0, itpp::elem_mult(phi, *w) );
  }
%}

class Waveform {
  public:
    Waveform();
    virtual ~Waveform();
    int rows();
    int cols();
    double energy();            // average energy per symbol (column)
    inline double power();      // total power of waveform
    inline double powerdb();
};

%extend Waveform {
  static Waveform ones(int nrows, int ncols) {
    Waveform w = itpp::ones_c(nrows, ncols);
    return w;
  }

  static Waveform zeros(int nrows, int ncols) {
    Waveform w = itpp::zeros_c(nrows, ncols);
    return w;
  }

  static Waveform randn(int nrows, int ncols, double mean=0.0, double variance=1.0) {
    double v = 0.5*variance/(nrows);
    itpp::Normal_RNG rng(mean, v);
    Waveform w = itpp::to_cmat(rng(nrows,ncols), rng(nrows,ncols) );
    return w;
  }
};

%pythoncode %{
Waveform.__repr__ = waveform_py_waveform__repr
Waveform.__add__  = waveform_py_waveform__add
Waveform.__radd__ = waveform_py_waveform__add
Waveform.__mul__  = waveform_py_waveform__mul
Waveform.__rmul__ = waveform_py_waveform__mul
Waveform.__call__ = waveform_py_waveform__get
Waveform.__getitem__ = waveform_py_waveform__get
Waveform.__setitem__ = waveform_py_waveform__set
Waveform.get_cols = waveform_py_waveform__get_cols
Waveform.concat   = waveform_py_waveform__concat_horizontal
Waveform.concat_horizontal = waveform_py_waveform__concat_horizontal
Waveform.concat_vertical   = waveform_py_waveform__concat_vertical
Waveform.apply_offset = waveform_py_waveform__apply_offset
%}
