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
// SPOWavefunctionStepper.cc
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>
#include "SPOWavefunctionStepper.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include "FourierTransform.h"
using namespace std;

////////////////////////////////////////////////////////////////////////////////
SPOWavefunctionStepper::SPOWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s)
    : tddt_(tddt), WavefunctionStepper(wf,tmap), ef_(ef), s_(s) {
}

////////////////////////////////////////////////////////////////////////////////
void SPOWavefunctionStepper::expVr(Wavefunction& newwf,const double dt) {
    // transform states to real space, multiply states by exp(-iv[r]) in real space
    // transform back to reciprocal space
    for ( int ispin = 0; ispin < newwf.nspin(); ispin++ ) {
        if (newwf.spinactive(ispin)) {
            for ( int ikp=0; ikp<newwf.nkp(); ikp++) {
                if (newwf.kptactive(ikp)) {
                    assert(wf_.sd(ispin,ikp) != 0);
                    SlaterDet& sdp = *(newwf.sd(ispin,ikp));
                    const Basis& wfbasis = sdp.basis();
                    FourierTransform ft(wfbasis,wfbasis.np(0),wfbasis.np(1),wfbasis.np(2));
                    if(s_.ctrl.petsc_Vr) {
                        ComplexMatrix& cp = sdp.c();
                        vector<complex<double> > tmp(ft.np012loc());
                        vector<complex<double> > ctmp(2*cp.mloc());
                        const int np012loc = ft.np012loc();
                        const int mloc = cp.mloc();
                        for ( int n = 0; n < sdp.nstloc(); n++ ) {
                            ft.backward(cp.cvalptr(n*mloc),&tmp[0]);

                            if(n==0 && s_.ctxt_.oncoutpe()){
                                cout<<"real space before expVr*psi ";
                                for (int ic=0;ic<10;ic++){
                                    cout<<tmp[ic]<<" ";
                                }
                                cout<<endl;
                            }

                            #pragma omp parallel for
                            for ( int i = 0; i < np012loc; i++ ) {
                                tmp[i]=exp(-complex<double>(0,1)*ef_.v_r[ispin][i]*dt)*tmp[i];
                            }

                            if(n==0&&s_.ctxt_.oncoutpe()){
                                cout<<"real space after expVr*psi ";
                                for (int ic=0;ic<10;ic++){
                                    cout<<tmp[ic]<<" ";
                                }
                                cout<<endl;
                            }
                            ft.forward(&tmp[0],cp.valptr(n*mloc));
                        }
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void SPOWavefunctionStepper::expKE(Wavefunction& newwf,const double dt) {
    for ( int ispin = 0; ispin < newwf.nspin(); ispin++ ) {
        if (newwf.spinactive(ispin)) {
            for ( int ikp=0; ikp<newwf.nkp(); ikp++) {
                if (newwf.kptactive(ikp)) {
                    assert(wf_.sd(ispin,ikp) != 0);
                    SlaterDet& sdp = *(newwf.sd(ispin,ikp));
                    const Basis& wfbasis = sdp.basis();
                    ComplexMatrix& cp = sdp.c();//newwf.sd(ispin,ikp)->c();
                    const int mloc = cp.mloc();
                    const double* kpg2 = wfbasis.kpg2_ptr();
                    const int ngwloc = wfbasis.localsize();
                    if(ef_.vp) kpg2 = ef_.vp->get_kpgpa2(wfbasis);
                    FourierTransform ft(wfbasis,wfbasis.np(0),wfbasis.np(1),wfbasis.np(2));
                    for ( int n = 0; n < sdp.nstloc(); n++ ) {
                        for ( int ig = 0; ig < ngwloc; ig++ ) {
                            cp[ig+mloc*n] = exp(complex<double>(0,-0.5*kpg2[ig]*dt))*cp[ig+mloc*n];
                        }
                    }
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
void SPOWavefunctionStepper::copy(const Wavefunction& newwf) {
    for ( int ispin = 0; ispin < newwf.nspin(); ispin++ ) {
        if (newwf.spinactive(ispin)) {
            for ( int ikp=0; ikp<newwf.nkp(); ikp++) {
                if (newwf.kptactive(ikp)) {
                    assert(newwf.sd(ispin,ikp) != 0);
                    wf_.sd(ispin,ikp)->c()=newwf.sd(ispin, ikp)->c();
                    if(s_.ctrl.petsc_tdH)
                        s_.hamil_wf->sd(ispin, ikp)->c() =newwf.sd(ispin, ikp)->c();
                }
            }
        }
    }

}

////////////////////////////////////////////////////////////////////////////////
void SPOWavefunctionStepper::update(Wavefunction& newwf) {
    for ( int ispin = 0; ispin < newwf.nspin(); ispin++ ) {
        if (newwf.spinactive(ispin)) {
            for ( int ikp=0; ikp<newwf.nkp(); ikp++) {
                if (newwf.kptactive(ikp)) {
                    assert(newwf.sd(ispin,ikp) != 0);
                    newwf.sd(ispin,ikp)->c()=wf_.sd(ispin,ikp)->c();
                }
            }
        }
    }
    if(s_.ctrl.petsc_KE) {
        expKE(newwf,tddt_/2);
    }

    if(s_.ctrl.petsc_tdH) {
        copy(newwf);
        ef_.hamil_cd()->update_density();
        ef_.update_hamiltonian();
        ef_.update_vhxc();
    }

    if(s_.ctrl.petsc_Vr) {
        expVr(newwf,tddt_);
    }

    if(s_.ctrl.petsc_KE) {
        expKE(newwf,tddt_/2);
    }
    copy(newwf);
}
