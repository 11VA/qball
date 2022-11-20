////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2013, Lawrence Livermore National Security, LLC.
// qb@ll:  Qbox at Lawrence Livermore
//
// This file is part of qb@ll.
//
// Produced at the Lawrence Livermore National Laboratory.
// Written by Erik Draeger (draeger1@llnl.gov) and Francois Gygi (fgygi@ucdavis.edu).
// Based on the Qbox code by Francois Gygi Copyright (c) 2008
// LLNL-CODE-635376. All rights reserved.
//
// qb@ll is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details, in the file COPYING in the
// root directory of this distribution or <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////
//
// MRRKWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: MRRKWavefunctionStepper.h,v 1.5 2011-06-02 15:56:19 schleife Exp $

#include <config.h>

#ifndef MRRKWAVEFUNCTIONSTEPPER_H
#define MRRKWAVEFUNCTIONSTEPPER_H

#include "EnergyFunctional.h"
#include "Wavefunction.h"
#include "WavefunctionStepper.h"

using namespace std;

class MRRKWavefunctionStepper : public WavefunctionStepper {
private:

    double tddt_;
    void updateV();
    void copy(const Wavefunction& oldwf);
    void copy(Wavefunction& newwf, const Wavefunction& oldwf);
    void update_mrrk(Wavefunction& dwf);
    void update_fast(Wavefunction& dwf);
    void update_slow(Wavefunction& dwf);

protected:

    EnergyFunctional & ef_;
    Sample & s_;

public:
    void update(Wavefunction& dwf);

    MRRKWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s);
    ~MRRKWavefunctionStepper() {};
};
#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:

// Local Variables:
// mode: c++
// End:
