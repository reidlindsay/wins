/*  :file: dot11n_channel.cc
 *
 *  Revision Info
 *  =============
 *  $LastChangedBy: mandke $
 *  $LastChangedDate: 2011-11-13 17:45:52 -0600 (Sun, 13 Nov 2011) $
 *  $LastChangedRevision: 5337 $
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

#define DOT11N_DEBUG  if (0) std::cerr

#include <dot11n_channel.h>
#include <wins-80211n/itpp.h>
#include <iostream>
#include <sstream>

#include <itpp/signal/filter.h>
#include <itpp/signal/transforms.h>
#include <itpp/stat/misc_stat.h>

/* Macros and constants */
#define CHANNEL_BANDWIDTH   100e6
#define DOT11N_BANDWIDTH    20e6
#define FIR_LENGTH          120
#define IS_TGN_MODEL(x)     (((x)>=DOT11N_TGN_MODEL_A) && ((x)<=DOT11N_TGN_MODEL_F))

/* Doppler parameters */
#define ENVIRONMENTAL_SPEED   1.2     // (km/hr)
#define SPEED_OF_LIGHT        3e8     // (m/s)
#define CARRIER_FREQUENCY     5.180e9 // (Hz)   # Channel 36
#define DOPPLER_MAX           5       // multiplier for fd -> maximum doppler spread
#define DOPPLER_A_NORM        9       // normalization constant
#define CORRELATION_FLOOR     0.5     // below this assume uncorrelated samples
#define DOPPLER_SAMPLE_FREQ   0.5e3   // doppler sampling frequency
#define DOPPLER_NFFT_MIN      128     // minimum FFT size for doppler filter
#define DOPPLER_NFFT_MAX      (DOPPLER_NFFT_MIN * (1<<4)) // max FFT size

using itpp::pi;
using itpp::sqr;
using itpp::exp;

std::complex<double> iorj(0,1);

Dot11N_Channel::Dot11N_Channel(int model, int ntx, int nrx, int flag)
  : d_model(model), d_ntx(ntx), d_nrx(nrx), d_flag(flag),
    d_nfft(DOPPLER_NFFT_MIN), d_tcorr(0), d_tsamp(1), 
    d_toffset(-1), d_phi_last(-1), d_env_speed(ENVIRONMENTAL_SPEED)
{
#if 0
  fprintf(stdout, "Dot11N_Channel(%s) created, %d x %d\n", model.data(), nrx, ntx);
#endif
  set_size(0);
  d_corr_sqrt.set_size(0);
  d_rice.set_size(0);
  d_hs.set_size(0);
  d_doppler.set_size(0,0);
  configure();
}

bool
Dot11N_Channel::setup(int model, int ntx, int nrx, int flag)
{
  d_model = model;
  set_ntx(ntx, false);
  set_nrx(nrx, false);
  set_flags(flag, false);
  set_environmentspeed(d_env_speed, false);
  return configure();
}

void
Dot11N_Channel::set_nrx(int nrx, bool update)
{
  if (d_nrx == nrx) return;   // do nothing
  d_nrx = nrx;
  if (update) configure();
}

void
Dot11N_Channel::set_ntx(int ntx, bool update)
{
  if (d_ntx == ntx) return;   // do nothing
  d_ntx = ntx;
  if (update) configure();
}

void
Dot11N_Channel::set_flags(int flag, bool update)
{
  if (d_flag == flag) return;   // do nothing
  d_flag= flag;
  if (update) configure();
}

void
Dot11N_Channel::set_model(int model, bool update)
{
  d_model= model;
  if (update) configure();
}

void
Dot11N_Channel::set_environmentspeed(double v, bool up)
{
  double ev = v;
  if (v<=0) ev = ENVIRONMENTAL_SPEED;   // use default speed?
  d_env_speed = ev;
  if (up) {
    int nfft = d_nfft;
    calc_doppler();
    if (d_nfft!=nfft) update();         // need to reinitialize?
  }
}

double
Dot11N_Channel::rmsdelay()
{
  double tau=0;
  int ntaps = taps();
  if (ntaps<2) return tau;   // single tap -> tau = 0

  /* set up parameters */
  assert(size()==ntaps);
  int ntx = d_ntx;
  int nrx = d_nrx;
  double Tsamp = 1e9/DOT11N_BANDWIDTH;    // sample period in ns
  itpp::vec dly = Tsamp*itpp::linspace(0,ntaps-1,ntaps);
  /* Calculate RMS delay spread from current channel */
  assert(ntx>0);
  assert(nrx>0);
  assert(ntaps>1);
  assert(dly.size()==ntaps);
  double den=0;
  double num_1=0;   // first moment of  delay spread = num_1/den
  double num_2=0;   // second moment of delay spread = num_2/den
  double p=0;
  for (int l=0; l<ntaps; l++) {
    p = sqr(itpp::norm((*this)(l), "fro"));
    num_1 += dly(l)*p;
    num_2 += sqr(dly(l))*p;
    den   += p;
  }
  assert(den>0);
  tau = sqrt(num_2/den - sqr(num_1/den));   // rms delay in ns
  return tau*1e-9;
}

double
Dot11N_Channel::calc_rmsdelay(int model, int ntx, int nrx)
{
  double tau = 0.0;
  if (!IS_TGN_MODEL(model)) return tau;

  /* Get power delay profile parameters */
  itpp::mat Pdb;                      // power of tap in dB
  itpp::ivec dly;                     // excess delay in ns
  /* Get model parameters */
  get_Pdb(Pdb, model);
  get_dly(dly, model);
  /* set up parameters */
  int ncluster = Pdb.rows();
  int ntaps    = Pdb.cols();
  itpp::mat pwr = itpp::inv_dB(Pdb);
  if (ntaps<2) return 0;    // single tap -> tau = 0
  /* Calculate RMS delay spread */
  assert(ntx>0);
  assert(nrx>0);
  assert(ntaps>1);
  assert(dly.size()==ntaps);
  double den=0;
  double num_1=0;   // first moment of  delay spread = num_1/den
  double num_2=0;   // second moment of delay spread = num_2/den
  double p=0;
  for (int l=0; l<ntaps; l++) {
    for (int c=0; c<ncluster; c++) {
      p = pwr(c,l)*ntx*nrx;
      num_1 += dly(l)*p;
      num_2 += sqr(dly(l))*p;
      den   += p;
    }
  }
  assert(den>0);
  tau = sqrt(num_2/den - sqr(num_1/den));   // rms delay in ns
  return tau*1e-9;
}

double
Dot11N_Channel::medianpower()
{
  /* set up channel matrix - h */
  itpp::cmat h(d_nrx*d_ntx, d_nfft);
  h.zeros();
  for (int l=0; l<size(); l++) {
    assert((*this)(l).size() == h.rows());
    h.set_col(l, itpp::cvectorize((*this)(l)) );
  }
  /* compute FFT - H = FFT(h) */
  itpp::cmat H(d_nrx*d_ntx, d_nfft);
  for (int n=0; n<h.rows(); n++) {
    H.set_row(n, itpp::fft(h.get_row(n)) );
  }
  /* compute power of each subcarrier */
  itpp::vec P; P.set_size(H.cols());
  for (int k=0; k<H.cols(); k++) {
    P(k) = sqr(itpp::norm(H.get_col(k), "fro"));
  }
  /* sort subcarrier powers and find median */
  return itpp::median(P);
}

void
Dot11N_Channel::update(double tupdate)
{
  bool isvalid = true ;   // are parameters valid?
  bool isnonstd = false;  // is model non-standard?
  // check if parameters are valid
  isvalid = isvalid && (d_nrx>0);
  isvalid = isvalid && (d_ntx>0);
  assert(isvalid);    // channel parameters are valid?
  // declare channel parameters
  double Kdb = itpp::vec("+inf")[0];  // K-factor (in dB)
  itpp::mat Pdb;                      // power of tap in dB
  itpp::ivec dly;                     // excess delay in ns
  /* Get model parameters */
  if (IS_TGN_MODEL(d_model) ) {
    get_Kdb(Kdb, d_model);
    get_Pdb(Pdb, d_model);
    get_dly(dly, d_model);
  } else {
    // non-std channel
    assert(d_nrx==d_ntx);
    isnonstd = true; set_size(1);
    (*this)(0) = itpp::eye_c(d_nrx);
    return;
  }
  /* Verify that correlation and Rice steering matrices are set */
  int ncluster = Pdb.rows();
  assert(d_rice.size()==ncluster);        // Rice steering matrices not set
  assert(d_corr_sqrt.size()==ncluster);   // Correlation matrices not set
  itpp::Uniform_RNG randu(0, 1);
  /* Implement one of the TGn standard models:
   *    1. Calculate Hf and Hv for each channel tap
   *    2. Combine Hf, Hv to create oversampled channel
   *    3. Apply decimation filter and downsample (if needed)
   */
  assert(!isnonstd);   // check for standard model flag
  if (d_flag & DOT11N_USE_NLOS) {
    Kdb = (d_flag & DOT11N_USE_LOS) ? Kdb : itpp::dB(0); 
  } else if (d_flag & DOT11N_USE_LOS) {
    Kdb = -itpp::dB(0);
  }
  // initialize channel parameters
  Channel H; itpp::cmat Hf, Hv, Hw;
  double K = itpp::inv_dB(Kdb);
  itpp::mat pwr = itpp::inv_dB(Pdb);
  double Kv, Kf;
  double A = 0; int ntaps = 1;
  int maxdelay = 1; int d = 0;
  // assemble oversampled channel
  int Tsamp = int(1e9/CHANNEL_BANDWIDTH); // channel sample period in ns (10 ns)
  assert(d_nrx<=4);
  assert(d_ntx<=4);
  if (d_flag & DOT11N_USE_NLOS) {
    // set NLOS parameters
    maxdelay = itpp::max(dly)/Tsamp + 1;
    ntaps = dly.size();
    assert(maxdelay>=ntaps);
  }

  H.set_size(maxdelay);
  assert(pwr.rows()>=ncluster);
  assert(pwr.cols()>=ntaps);
  assert(H.size()==maxdelay);
  for (d=0; d<H.size(); d++) H(d) = itpp::zeros_c(d_nrx, d_ntx);

  /* Compute phase of LOS component */
  double phi = 0;
  double vo = environmentspeed() * (1e3/3600);      // (km/hr) -> (m/s)
  double speed_of_light = SPEED_OF_LIGHT;           // (m/s)
  double carrier_frequency = CARRIER_FREQUENCY;     // (Hz)
  if ((d_toffset<0)||(tupdate<0)) d_phi_last = -1;  // reset phi
  if (d_phi_last<0) {
    phi = 2.0*pi*randu();   // uninitialized phi, so just use random phase
  } else {
    assert(!(tupdate<0));
    phi = d_phi_last;
    phi += 2.0*pi*(vo/(speed_of_light/carrier_frequency))*tupdate*sqrt(0.5);
  }
  while (phi>2*pi) phi -= 2*pi*int(phi/(2*pi)); // modulo 2*pi
  d_phi_last = phi;
  /* Calculate channel taps */
  assert(ntaps*ncluster>0);
  DOT11N_DEBUG << "Calling update_rayleigh() directly, nelem = "<< ntaps*ncluster;
  DOT11N_DEBUG << ", ntaps = "<< ntaps << ", ncluster = " << ncluster <<std::endl;
  update_rayleigh(tupdate, ntaps*ncluster);   // update Rayleigh components
  for (int l=0; l<ntaps; l++) {
    A = sqrt(pwr(0,l));
    d = dly(l)/Tsamp;   // convert tap index to delay index
    assert(d<=H.size());
    // create LOS channel
    Hf = itpp::zeros_c(d_nrx, d_ntx);
    if (d_flag & DOT11N_USE_LOS) {
      assert(d_rice.size()==ncluster);
      if (d==0) { // only use LOS channel for first tap
        for (int c=0; c<ncluster; c++) {
          A = sqrt(pwr(c,d));
          // set Hf according to flags
          if (d_flag & DOT11N_LOS_RICE) {
            // use Rice steering matrix
            if (A>0) Hf += A*exp(iorj*phi)*d_rice(c);
          } else {
            // use identity matrix - [DEFAULT]
            if (A>0) Hf += A*itpp::eye_c(d_nrx);
          }
        }
      }
    }
    // create NLOS channel
    Hv = itpp::zeros_c(d_nrx, d_ntx);
    if (d_flag & DOT11N_USE_NLOS) {
      assert(d_corr_sqrt.size()==ncluster);
      for (int c=0; c<ncluster; c++) {
        A = sqrt(pwr(c,l));
        if (A>0) {
          // set Hw according to flags
          if (d_flag & DOT11N_NLOS_FADING) {
            // use Rayleigh fading for cluster c, tap l
            DOT11N_DEBUG << "Calling get_rayleigh on tap "<< l << ", cluster " << c;
            DOT11N_DEBUG << ", toffset = "<< d_toffset<<std::endl;
            Hw = get_rayleigh(c+l*ncluster);
          } else if (d_flag & DOT11N_NLOS_EQUALGAIN) {
            // use matrix with random phases
            Hw = itpp::elem_mult(exp(iorj*2.0*pi*randu(d_nrx,d_ntx)), itpp::ones_c(d_nrx, d_ntx) );
          } else {
            // use identity matrix - [DEFAULT]
            Hw = itpp::eye_c(d_nrx);
          }
          Hv += itpp::reshape(A*d_corr_sqrt(c)*itpp::reshape(Hw, d_nrx*d_ntx, 1), d_nrx, d_ntx);
        }
      }
    }
    // combine LOS and NLOS
    if (d) Kv = 1;        // NLOS; Keff = -inf (dB)
    else   Kv = 1/(K+1);  // NLOS; Keff = K (dB)
    Kf = 1 - Kv;          // LOS = 1 - NLOS
    H(d) = sqrt(Kf)*Hf + sqrt(Kv)*Hv;
  }   // end for (... ntaps ...)

  /* Apply decimation filter and downsample channel */
  if (maxdelay>1) {
    // filter parameters
    int upsample = int(CHANNEL_BANDWIDTH/DOT11N_BANDWIDTH); // upsample factor
    itpp::vec fir = itpp::fir1(FIR_LENGTH, 1.0/upsample);   // low pass filter
    assert(fir.size()==FIR_LENGTH+1);
    // initialize downsampled channel
    int numtaps = (maxdelay-1)/upsample + 1;
    set_size(numtaps);
    for (int l=0; l<size(); l++) (*this)(l) = itpp::zeros_c(d_nrx, d_ntx);
    // filter setup
    itpp::MA_Filter<std::complex<double>,double,std::complex<double> > filt(fir);
    itpp::cvec h, y;
    // downsampling setup
    itpp::bvec idx;
    idx.set_size(maxdelay); idx.zeros();
    for (int l=0; l<idx.size(); l++) idx(l) = (l%upsample==0) ? 1:0;
    // Apply filter
    for (int m=0; m<d_nrx; m++) {
      for (int n=0; n<d_ntx; n++) {
        // apply filter to each antenna pair
        h.set_size(maxdelay);
        assert(h.size()==H.size());
        for (int l=0; l<h.size(); l++) h(l) = H(l)(m,n);
        // filter - only if needed
        if (upsample>1) {
          filt.clear();
          y = filt(itpp::concat(h, itpp::zeros_c(FIR_LENGTH)));
          assert(y.size()==(h.size()+FIR_LENGTH) );
          y = upsample*y.get(FIR_LENGTH/2, FIR_LENGTH/2 +  h.size() - 1);
        } else {
          y = h;
        }
        assert(y.size()==idx.size());
        // downsample
        y = y.get(idx);
        int s = 0; for (int z=0; z<idx.size(); z++) if (idx(z)) s++;
        assert(y.size()==size());
        for (int l=0; l<size(); l++) (*this)(l)(m,n) = y(l);
      }
    }
  } else {
    // no decimation filtering needed - only use first channel tap
    assert(ntaps==1);
    set_size(ntaps);
    (*this)(0) = H(0);
  }
}

/* XXX FIXME XXX:
 *    Since the doppler filter is currently implemented using an IFFT, there is
 *    a "block-like" behavior. After Nfft*Ts samples, the filtered output will
 *    run out and need to be reinitialized. Thus causing, a step-wise change in
 *    the filtered channel.
 *
 *    This might be fixed by implementing the doppler filter as an IIR or FIR
 *    (i.e. time-domain) filter.
 */
void
Dot11N_Channel::update_rayleigh(double tupdate, unsigned int nelem)
{
  /* Only apply filter if DOT11N_NLOS_DOPPLER flag is set */
  bool usefading = (d_flag & DOT11N_NLOS_FADING);
  bool usenlos   = (d_flag & DOT11N_USE_NLOS);
  bool usefilter = (d_flag & DOT11N_NLOS_DOPPLER);
  if (!(usefilter && usefading && usenlos)) return;

  itpp::cmat hw;  // IID Gaussian samples
  itpp::cmat Hs;  // filtered frequency response of channel
  itpp::Complex_Normal_RNG randn(0, 1);

  // d_hs -> buffer for filtered time domain channel
  if (d_hs.size()==nelem) {
    if ((d_hs.rows()!=(d_nrx*d_ntx)) || (d_hs.cols()!=d_nfft)) tupdate = -1;
  } else {
    tupdate = -1;
  }

  // verify dimensions of doppler filter
  assert(d_doppler.rows() == d_nrx*d_ntx);
  assert(d_doppler.cols() == d_nfft);

  if (tupdate<0) {
    // Reset all components in d_hs
    assert(nelem>0);
    d_hs.set_size(nelem);
    d_toffset = 0;
    DOT11N_DEBUG << "Reinitializing rayleigh, nelem = " << nelem << std::endl;
    for (int d=0; d<nelem; d++) {
      // Reinitialize Gaussian samples
      hw = randn(d_nrx*d_ntx, d_nfft);
      // Fitler Gaussian samples
      Hs = itpp::to_cmat(itpp::elem_mult(itpp::real(hw), d_doppler), \
                         itpp::elem_mult(itpp::imag(hw), d_doppler));
      d_hs(d) = itpp::zeros_c(Hs.rows(), Hs.cols());
      // Apply IFFT to filtered channel response
      for (int k=0; k<d_hs(d).rows(); k++) {
        d_hs(d).set_row(k, itpp::ifft(Hs.get_row(k)) );
      }
      assert(d_hs(d).rows()==d_nrx*d_ntx);
      assert(d_hs(d).cols()==d_nfft);
    }
  } else {
    assert(d_hs.size()==nelem);
    // Reinitialize Gaussian samples if 4*Tcorr has passed since last update
    if ((tupdate>4*d_tcorr)||(d_toffset<0)||(tupdate+d_toffset>=d_tsamp*d_nfft)) {
      DOT11N_DEBUG << " -> out of range: calling update_rayleigh()" <<std::endl;
      return update_rayleigh(-1, nelem);
    }
    // tupdate is in buffer
    assert(tupdate+d_toffset<d_tsamp*d_nfft);
    // Update offset into buffer of fading samples
    d_toffset += tupdate;
  }
  // Verify dimensions of d_hs
  assert(d_hs.size()==nelem);
  assert(d_hs.cols()==d_nfft);
  assert(d_hs.rows()==d_nrx*d_ntx);
  DOT11N_DEBUG << "-> update finished: d_hs.size = " << d_hs.size() << std::endl;
}

itpp::cmat
Dot11N_Channel::get_rayleigh(unsigned int c)
{
  itpp::cmat Hw;  // return filtered channel
  itpp::cmat hw;  // IID Gaussian samples
  itpp::cmat Hs;  // filtered frequency response of channel
  itpp::Complex_Normal_RNG randn(0, 1);

  /* Only apply filter if DOT11N_NLOS_DOPPLER flag is set */
  bool usefilter = (d_flag & DOT11N_NLOS_DOPPLER);
  if (!usefilter) {
    /* create Rayleigh channel sample */
    Hw = randn(d_nrx, d_nrx);
    return Hw;
  }

  // d_hs -> buffer for filtered time domain channel
  bool needupdate = false;
  if ((d_hs.size()>c) && (d_toffset>=0)) {
    if ((d_hs(c).rows()!=(d_nrx*d_ntx)) || (d_hs(c).cols()!=d_nfft)) needupdate = true;
  } else {
    needupdate = true;
  }
  if (needupdate){
    DOT11N_DEBUG << " -> need update, d_hs.size = "<<d_hs.size();
    DOT11N_DEBUG << ", nrows = " <<d_hs.rows();
    DOT11N_DEBUG << ", ncols = " <<d_hs.cols() << std::endl;
    DOT11N_DEBUG << "                 c = " << c << ", nfft  = " <<d_nfft;
    DOT11N_DEBUG << ", toffset = " << d_toffset <<std::endl;
    update_rayleigh(-1, c+1);
  }

  // verify dimensions of doppler filter
  assert(d_doppler.rows() == d_nrx*d_ntx);
  assert(d_doppler.cols() == d_nfft);

  // Get column corresponding to offset (interpolate between samples)
  assert(d_toffset>=0);
  assert(d_hs.size()>c);
  itpp::cmat Ha, Hb;
  int offset; double gamma=0.0;
  offset = int(d_toffset/d_tsamp);
  assert(offset<d_hs.cols());
  if (offset+1<d_hs(c).cols()) {
    Ha = itpp::reshape(d_hs(c).get_col(offset), d_nrx, d_ntx);
    Hb = itpp::reshape(d_hs(c).get_col(offset+1), d_nrx, d_ntx);
    gamma = (d_toffset - d_tsamp*offset)/d_tsamp;
    assert(gamma<1.0); assert(gamma>=0.0);
    Hw = (1.0 - gamma)*Ha + gamma*Hb;
  } else {
    Hw = itpp::reshape(d_hs(c).get_col(offset), d_nrx, d_ntx);
  }
  return Hw;
}

bool
Dot11N_Channel::configure()
{
  bool isvalid = true ;   // are parameters valid?
  // check if parameters are valid
  isvalid = isvalid && (d_nrx>0);
  isvalid = isvalid && (d_ntx>0);
  // Set up channel if valid
  if (isvalid) {
    // calculate correlation and Rice steering matrices
    calc_rice();
    calc_correlation();
    calc_doppler();
    update();
  } else {
    /* create default channel (empty) */
    set_size(0);
  }
  return isvalid;
}

void
Dot11N_Channel::calc_rice() {
  bool isvalid = true ;   // are parameters valid?
  // check if parameters are valid
  isvalid = isvalid && (d_nrx>0);
  isvalid = isvalid && (d_ntx>0);
  assert(isvalid);    // channel parameters are valid?
  // declare channel parameters
  itpp::mat Pdb;        // power of tap in dB
  itpp::vec aoa, aod;   // AoA, AoD
  if (IS_TGN_MODEL(d_model) ) {
    get_Pdb(Pdb, d_model);
    get_aoa(aoa, d_model);
    get_aod(aod, d_model);
  } else {
    // non-std channel
    d_rice.set_size(0);  // reset and return
    return;
  }
  /* Compute Rice steering matrices */
  Channel S;
  itpp::vec idx;
  itpp::cmat Srx, Stx;
  int ncluster = Pdb.rows();
  double phi,D;
  assert(ncluster>0);             // too few clusters in calc_rice()
  assert(aoa.size()==ncluster);   // wrong AoA size   in calc_rice()
  assert(aod.size()==ncluster);   // wrong AoD size   in calc_rice()
  // set S for each cluster - d = lambda/2 (antenna separation)
  S.set_size(ncluster);
  for (int k=0; k<S.size(); k++) {
    // create Srx
    D = pi; phi = aoa(k)*(pi/180);    // receive parameters (AoA)
    idx = itpp::linspace(0,d_nrx-1,d_nrx);
    Srx = exp(iorj*D*sin(phi*idx));
    // create Stx
    D = pi; phi = aod(k)*(pi/180);    // transmit parameters (AoD)
    idx = itpp::linspace(0,d_ntx-1,d_ntx);
    Stx = exp(iorj*D*sin(phi*idx));
    // S =  Srx * (Stx)'
    S(k) = Srx*Stx.T();
  }
  /* Copy to local member */
  d_rice = S;
}

void
Dot11N_Channel::calc_correlation() {
  bool isvalid = true ;   // are parameters valid?
  // check if parameters are valid
  isvalid = isvalid && (d_nrx>0);
  isvalid = isvalid && (d_ntx>0);
  assert(isvalid);    // channel parameters are valid?
  // declare channel parameters
  itpp::mat Pdb;                    // power of tap in dB
  itpp::vec aoa, aod, asa, asd;     // AoA, AoD, AS (receiver), AS (transmitter)
  if (IS_TGN_MODEL(d_model) ) {
    get_Pdb(Pdb, d_model);
    get_aoa(aoa, d_model); get_aod(aod, d_model);
    get_asa(asa, d_model); get_asd(asd, d_model);
  } else {
    // non-std channel
    d_corr_sqrt.set_size(0);  // reset and return
    return;
  }
  /* Compute RX correlation */
  Channel Rrx, Rtx; itpp::cmat a_ula;
  itpp::vec idx; itpp::mat diff;
  itpp::cmat rnum; itpp::mat rden;
  int ncluster = Pdb.rows();
  double phi, dev, beta, D;
  assert(ncluster>0);             // too few clusters in calc_correlation()
  assert(aoa.size()==ncluster);   // wrong AoA size   in calc_correlation()
  assert(asa.size()==ncluster);   // wrong AS size    in calc_correlation()
  // calculate Rrx for each cluster - d = lambda/2 (antenna separation)
  Rrx.set_size(ncluster);
  for (int k=0; k<Rrx.size(); k++) {
    D = pi;       // normalized antenna separation
    phi = aoa(k)*(pi/180);
    dev = asa(k)*(pi/180);
    idx = itpp::linspace(0,d_nrx-1,d_nrx);
    beta = 1.0/(1 - exp(-sqrt(2)*pi/dev));
    diff = itpp::toeplitz(idx, -idx);
    rnum = beta*exp(iorj*D*sin(phi)*diff);
    rden = 1 + 0.5*sqr(dev*D*cos(phi)*diff);
    Rrx(k) = itpp::elem_div(rnum, itpp::to_cmat(rden));
  }
  /* Compute TX correlation */
  assert(aod.size()==ncluster);   // wrong AoD size   in nlos_channel()
  assert(asd.size()==ncluster);   // wrong AS size    in nlos_channel()
  // calculate Rtx for each cluster
  Rtx.set_size(ncluster);
  for (int k=0; k<Rtx.size(); k++) {
    D = pi;       // normalized antenna separation
    phi = aod(k)*(pi/180);
    dev = asd(k)*(pi/180);
    idx = itpp::linspace(0,d_ntx-1,d_ntx);
    beta = 1.0/(1 - exp(-sqrt(2)*pi/dev));
    diff = itpp::toeplitz(idx, -idx);
    rnum = beta*exp(iorj*D*sin(phi)*diff);
    rden = 1 + 0.5*sqr(dev*D*cos(phi)*diff);
    Rtx(k) = itpp::elem_div(rnum, itpp::to_cmat(rden));
  }
  /* Combine Rtx and Rrx */
  d_corr_sqrt.set_size(ncluster);
  assert(Rtx.size()==ncluster);
  assert(Rrx.size()==ncluster);
  for (int k=0; k<d_corr_sqrt.size(); k++) {
    d_corr_sqrt(k) = itpp::sqrtm(itpp::kron(Rtx(k), Rrx(k)));
  }
}

void
Dot11N_Channel::calc_doppler() {
  /* Doppler parameters */
  double fc = CARRIER_FREQUENCY;                  // (Hz)
  double c = SPEED_OF_LIGHT;                      // (m/s)
  double vo = environmentspeed() * (1e3/3600);    // (km/hr) -> (m/s)
  double fd = vo / (c/fc);                        // doppler spread
  double fmax = DOPPLER_MAX*fd;                   // maximum doppler spread
  double A = DOPPLER_A_NORM;                      // normalization constant
  /* Coherence time parameter */
  double xcorr = CORRELATION_FLOOR;                       // corrleation floor
  double Tcorr = std::log(1.0/xcorr)*sqrt(A)/(2*pi*fd);   // coherence time
  /* Filter & IFFT parameters */
  int    Nfft = 1;                  // initial FFT size
  int    Nmin = DOPPLER_NFFT_MIN;   // Minimum FFT size
  int    Nmax = DOPPLER_NFFT_MAX;   // Maximum FFT size
  double fs = DOPPLER_SAMPLE_FREQ;  // sample frequency
  /* Adjust Nfft and fs as needed */
  bool done = false;
  double df = fs/Nfft, falpha = 0.10;
  while (!done) {
    Nfft = 1;
    while (Nfft<2*fs*Tcorr) Nfft*=2;  // ensure Nfft*Ts > 2*Tcorr
    while (Nfft<5*fs/fd)    Nfft*=2;  // ensure 5*df < fd
    if (Nfft<Nmin) {
      fs = (1+falpha)*fs;       // increase fs
    } else if (Nfft>Nmax) {
      fs = (1-falpha)*fs;       // decrease fs
    } else {
      done = true;
    }
    df = fs/Nfft;   // update df
  }
  assert(Nfft<=Nmax);
  assert(Nfft>=Nmin);
  /* Compute doppler filter */
  double K = Nfft*df*sqrt(A)/(pi*fd);
  itpp::vec permk = itpp::concat(itpp::linspace(0, Nfft/2-1, Nfft/2), \
                                 itpp::linspace(-Nfft/2, -1, Nfft/2));
  itpp::vec box;
  // create truncated doppler filter
  box.set_size(permk.size());
  for (int k=0; k<box.size(); k++) {
    box(k) = (abs(permk(k)*df)<fmax)?1:0;
  }
  // filter is normalized by sqrt(Nfft*K)
  d_doppler = sqrt(Nfft*K*itpp::elem_div(box, 1 + A*sqr(permk*df/fd)));
  assert(d_doppler.size() == Nfft);
  d_doppler = itpp::repmat(itpp::reshape(d_doppler, 1, Nfft), d_nrx*d_ntx, 1);
  assert(d_doppler.rows() == d_nrx*d_ntx);
  assert(d_doppler.cols() == Nfft);
  // set other parameters
  d_nfft = Nfft;
  d_tcorr = Tcorr;
  d_tsamp = 1.0/fs;
  assert(d_tcorr<d_tsamp*d_nfft);   // verify that Tcorr<Ts*Nfft
}

void
Dot11N_Channel::get_Kdb(double & Kdb, int model) {
  switch (model)
  {
    case DOT11N_TGN_MODEL_A:
      Kdb = itpp::vec(MODELA_KDB)[0]; break;
    case DOT11N_TGN_MODEL_B:
      Kdb = itpp::vec(MODELB_KDB)[0]; break;
    case DOT11N_TGN_MODEL_C:
      Kdb = itpp::vec(MODELC_KDB)[0]; break;
    case DOT11N_TGN_MODEL_D:
      Kdb = itpp::vec(MODELD_KDB)[0]; break;
    case DOT11N_TGN_MODEL_E:
      Kdb = itpp::vec(MODELE_KDB)[0]; break;
    case DOT11N_TGN_MODEL_F:
      Kdb = itpp::vec(MODELF_KDB)[0]; break;
    default:
      // non-std channel
      Kdb = itpp::vec("+inf")[0];  // K-factor (in dB)
  }
}

void
Dot11N_Channel::get_Pdb(itpp::mat & Pdb, int model) {
  switch (model)
  {
    case DOT11N_TGN_MODEL_A:
      Pdb = MODELA_PWR; break;
    case DOT11N_TGN_MODEL_B:
      Pdb = MODELB_PWR; break;
    case DOT11N_TGN_MODEL_C:
      Pdb = MODELC_PWR; break;
    case DOT11N_TGN_MODEL_D:
      Pdb = MODELD_PWR; break;
    case DOT11N_TGN_MODEL_E:
      Pdb = MODELE_PWR; break;
    case DOT11N_TGN_MODEL_F:
      Pdb = MODELF_PWR; break;
    default:
      // non-std channel
      Pdb = "0";      // power of taps (in dB)
  }
}

void
Dot11N_Channel::get_dly(itpp::ivec & dly, int model) {
  switch (model)
  {
    case DOT11N_TGN_MODEL_A:
      dly = MODELA_DLY; break;
    case DOT11N_TGN_MODEL_B:
      dly = MODELB_DLY; break;
    case DOT11N_TGN_MODEL_C:
      dly = MODELC_DLY; break;
    case DOT11N_TGN_MODEL_D:
      dly = MODELD_DLY; break;
    case DOT11N_TGN_MODEL_E:
      dly = MODELE_DLY; break;
    case DOT11N_TGN_MODEL_F:
      dly = MODELF_DLY; break;
    default:
      // non-std channel
      dly = "0";      // excess delay in ns
  }
}

void
Dot11N_Channel::get_aoa(itpp::vec & aoa, int model) {
  switch (model)
  {
    case  DOT11N_TGN_MODEL_A:
      aoa = itpp::vec(MODELA_AOA); break;
    case  DOT11N_TGN_MODEL_B:
      aoa = itpp::vec(MODELB_AOA); break;
    case  DOT11N_TGN_MODEL_C:
      aoa = itpp::vec(MODELC_AOA); break;
    case  DOT11N_TGN_MODEL_D:
      aoa = itpp::vec(MODELD_AOA); break;
    case  DOT11N_TGN_MODEL_E:
      aoa = itpp::vec(MODELE_AOA); break;
    case  DOT11N_TGN_MODEL_F:
      aoa = itpp::vec(MODELF_AOA); break;
    default:
      // non-std channel
      aoa = "0"; // AoA
  }
}

void
Dot11N_Channel::get_aod(itpp::vec & aod, int model) {
  switch (model)
  {
    case  DOT11N_TGN_MODEL_A:
      aod = itpp::vec(MODELA_AOD); break;
    case  DOT11N_TGN_MODEL_B:
      aod = itpp::vec(MODELB_AOD); break;
    case  DOT11N_TGN_MODEL_C:
      aod = itpp::vec(MODELC_AOD); break;
    case  DOT11N_TGN_MODEL_D:
      aod = itpp::vec(MODELD_AOD); break;
    case  DOT11N_TGN_MODEL_E:
      aod = itpp::vec(MODELE_AOD); break;
    case  DOT11N_TGN_MODEL_F:
      aod = itpp::vec(MODELF_AOD); break;
    default:
      // non-std channel
      aod = "0"; // AoD
  }
}

void
Dot11N_Channel::get_asa(itpp::vec & asa, int model) {
  switch (model)
  {
    case  DOT11N_TGN_MODEL_A:
      asa = itpp::vec(MODELA_ASA); break;
    case  DOT11N_TGN_MODEL_B:
      asa = itpp::vec(MODELB_ASA); break;
    case  DOT11N_TGN_MODEL_C:
      asa = itpp::vec(MODELC_ASA); break;
    case  DOT11N_TGN_MODEL_D:
      asa = itpp::vec(MODELD_ASA); break;
    case  DOT11N_TGN_MODEL_E:
      asa = itpp::vec(MODELE_ASA); break;
    case  DOT11N_TGN_MODEL_F:
      asa = itpp::vec(MODELF_ASA); break;
    default:
      // non-std channel
      asa = "0"; // AS (receiver)
  }
}

void
Dot11N_Channel::get_asd(itpp::vec & asd, int model) {
  switch (model)
  {
    case  DOT11N_TGN_MODEL_A:
      asd = itpp::vec(MODELA_ASD); break;
    case  DOT11N_TGN_MODEL_B:
      asd = itpp::vec(MODELB_ASD); break;
    case  DOT11N_TGN_MODEL_C:
      asd = itpp::vec(MODELC_ASD); break;
    case  DOT11N_TGN_MODEL_D:
      asd = itpp::vec(MODELD_ASD); break;
    case  DOT11N_TGN_MODEL_E:
      asd = itpp::vec(MODELE_ASD); break;
    case  DOT11N_TGN_MODEL_F:
      asd = itpp::vec(MODELF_ASD); break;
    default:
      // non-std channel
      asd = "0"; // AS (transmitter)
  }
}

std::string
Dot11N_Channel::model_string() const
{
  std::string mstr = "";
  std::string fstr = "";
  bool ismodel = true;
  bool hasflag = false;

  /* Convert model to string */
  switch(d_model)
  {
    case DOT11N_TGN_MODEL_A:
      mstr += "TGN Model A"; break;
    case DOT11N_TGN_MODEL_B:
      mstr += "TGN Model B"; break;
    case DOT11N_TGN_MODEL_C:
      mstr += "TGN Model C"; break;
    case DOT11N_TGN_MODEL_D:
      mstr += "TGN Model D"; break;
    case DOT11N_TGN_MODEL_E:
      mstr += "TGN Model E"; break;
    case DOT11N_TGN_MODEL_F:
      mstr += "TGN Model F"; break;
    default:
      // do nothing;
      mstr += "NULL";
      ismodel = false;
  }

  if (!ismodel) return mstr;

  /* Identify TGN model flags */
#define NFLAGS  6
  assert(NFLAGS==6);  // only supports flags below
  int tgnflags[NFLAGS] = {DOT11N_USE_LOS, DOT11N_USE_NLOS, \
                          DOT11N_LOS_RICE, DOT11N_NLOS_FADING, \
                          DOT11N_NLOS_EQUALGAIN, DOT11N_NLOS_DOPPLER};
  std::string vals[NFLAGS] = {"LOS", "NLOS", "RICE", "FADING", \
                              "EQUALGAIN", "DOPPLER"};
  for (int k=0; k<NFLAGS; k++) {
    // check if flag matches and set string
    if (d_flag & tgnflags[k]) {
      if (hasflag) fstr += "|";
      hasflag = true; fstr += vals[k];
    }
  }

  if (hasflag) {
    mstr += " ("; mstr += fstr; mstr += ")";
  }
  return mstr;
}

Dot11N_Channel &
Dot11N_Channel::operator=(const Channel & H) {
  set_size(H.size());
  for (int k=0; k<size(); k++) {
    (*this)(k) = H(k);
  }
  return (*this);
}
