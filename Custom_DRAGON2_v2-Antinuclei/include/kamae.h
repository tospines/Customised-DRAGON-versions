/**
 * @file cparammodel.h
 * @brief Header file for cparamlib
 *
 * It implements Kamae et al model for secondary electron/positron production in pp interactions.
 */

#ifndef _KAMAE_H_
#define _KAMAE_H_

#include "cparamlib/cparamlib.h"

namespace KamaeYields {
    double GetSigma(double, double, PARTICLE_IDS);
};

#endif
