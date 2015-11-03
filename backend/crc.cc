/*  :file: crc.cc
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

#include <crc.h>
#include <itpp/itbase.h>


CRC::CRC() : coder("CCITT-32")
{
}

CRC::CRC(const std::string &code) : coder(code)
{
}

void
CRC::set_code(const std::string &code)
{
  coder.set_code(code);
}

int
CRC::parity(const std::string &input, int nbits)
{
  int blen = nbits;
  int slen = input.length();

  /* check length parameters */
  if (blen>8*slen)
      throw "[CRC] Error: Input string not long enough!";
  if (blen<0)
      blen = 8*slen;

  /* set remaining parameters */
  unsigned int stop = blen/8;
  unsigned int remain = blen%8;

  /* convert string to bit-string */
  assert (blen>=0);
  itpp::bvec inbits(blen); inbits.zeros();
  for (unsigned int k=0; k<stop; k++)
      inbits.set_subvector(k*8, itpp::dec2bin(8, input[k]) );
  if (remain>0) {
      unsigned char mask=0x00;
      for (unsigned int k=0; k<8-remain; k++) mask |= 0x01<<k;
      assert(stop<slen);
      inbits.set_subvector(stop*8, itpp::dec2bin(remain, input[stop]&(~mask)) );
  }

  /* calculate parity bits */
  itpp::bvec par;
  coder.parity(inbits, par);

  /* return parity bits */
  int crc = itpp::bin2dec(par);
  return crc;
}
