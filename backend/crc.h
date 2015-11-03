/*  :file: crc.h
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

#ifndef INCLUDED_CRC_H
#define INCLUDED_CRC_H

#include <itpp/comm/crc.h>
#include <string>

/*  :class: CRC
 *  :description: Wrapper for IT++ CRC_code.
 */
class CRC {
  private:
    itpp::CRC_Code coder;

  public:
    CRC();
    CRC(const std::string &code);
    virtual ~CRC() {}
    void set_code(const std::string &code);
    int parity(const std::string &input, int nbits=-1);
};

#endif  /* INCLUDED_CRC_H */
