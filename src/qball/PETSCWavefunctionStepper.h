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
// PETSCWavefunctionStepper.h
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ExponentialWavefunctionStepper.h,v 1.5 2011-06-02 15:56:19 schleife Exp $

#include <config.h>

#ifndef PETSCWAVEFUNCTIONSTEPPER_H
#define PETSCWAVEFUNCTIONSTEPPER_H

#include "EnergyFunctional.h"
#include "SelfConsistentPotential.h"
#include "Wavefunction.h"
#include "WavefunctionStepper.h"
#include <deque>

#include <petscts.h>

using namespace std;

// AK: context for rhs function
struct PETSC_CTX {
    EnergyFunctional & ef_;
    Sample & s_;
    Wavefunction & dwf_;
    Vec & dwf_vec_;
    vector<vector<double>> & fion_;
    valarray<double> & sigma_eks_;
    TimerMap & tmap_;
    Vec & hamil_wf_vec_;
    // AK: needed to initialize
    PETSC_CTX(EnergyFunctional & ef, Sample & s, Wavefunction & dwf, Vec & dwf_vec, vector<vector<double>> & fion, valarray<double> & sigma_eks, TimerMap & tmap,Vec& hamil_wf_vec)
        : ef_(ef), s_(s), dwf_(dwf), dwf_vec_(dwf_vec), fion_(fion), sigma_eks_(sigma_eks), tmap_(tmap), hamil_wf_vec_(hamil_wf_vec) {}
};

class PETSCWavefunctionStepper : public WavefunctionStepper {
private:

    double tddt_;
    std::vector<SelfConsistentPotential> potential_;
    Wavefunction newwf_;
    Wavefunction dwf_;

    PetscErrorCode ierr;
    TS petsc_ts;
    Vec petsc_wf_vec;        // holds flattened wf array
    Vec petsc_dwf_vec;

    Vec hamil_wf_vec;

    PETSC_CTX * petsc_ctx;

protected:

    EnergyFunctional & ef_;
    Sample & s_;

public:
    void update(Wavefunction& dwf);
    static PetscErrorCode dummy_RHS(TS ts, PetscReal t, Vec wf_vec, Vec rhs, void *ctx_);
    static PetscErrorCode RHS(TS ts, PetscReal t, Vec wf_vec, Vec rhs, void *ctx_ );
    static PetscErrorCode RegisterMyRKC2(void);

    //static PetscErrorCode unit_RHS(TS ts, PetscReal t, Vec wf_vec, Vec rhs, void *ctx_); 
    //static PetscErrorCode unit_evolve_RHS(TS ts, PetscReal t, Vec wf_vec, Vec rhs, void *ctx_);

    PETSCWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s, vector<vector<double>> & fion, valarray<double> & sigma_eks);
    ~PETSCWavefunctionStepper();
};
#endif

// Local Variables:
// mode: c++
// coding: utf-8
// End:

// Local Variables:
// mode: c++
// End:
