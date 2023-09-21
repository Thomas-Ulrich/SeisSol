/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Common kernel-level functions
 **/

#ifndef KERNElS_DENSEMATRIXOPS_HPP_
#define KERNElS_DENSEMATRIXOPS_HPP_

#if defined(__SSE3__) || defined(__MIC__)
#include <immintrin.h>
#endif

#if defined(__AVX512F__)
#include "denseMatrixOps_AVX512.hpp"
#elif defined(__MIC__)
#include "denseMatrixOps_MIC.hpp"
#elif defined(__AVX2__)
#include "denseMatrixOps_AVX2.hpp"
#elif defined(__AVX__)
#include "denseMatrixOps_AVX.hpp"
#elif defined(__SSE3__)
#include "denseMatrixOps_SSE3.hpp"
#elif defined(__ARM_FEATURE_SVE)
#include "denseMatrixOps_SVE.hpp"
#elif defined(__aarch64__)
#include "denseMatrixOps_AARCH64.hpp"
#else
#include "denseMatrixOps_noarch.hpp"
#endif

#include <cassert>

namespace seissol::kernels {
    /** Stores X in Y with non-temporal hint.
     * 
     * @param numberOfReals The size of X and Y.
     * @param X
     * @param Y
     */
    inline void streamstore( unsigned numberOfReals,
                             real const* X,
                             real* Y )
    {
      assert(numberOfReals % DMO_INCREMENT == 0);
      
      for (unsigned i = 0; i < numberOfReals; i += DMO_INCREMENT) {
        DMO_STREAM(&X[i], &Y[i])
      }
    }
} // namespace seissol::kernels

#endif
