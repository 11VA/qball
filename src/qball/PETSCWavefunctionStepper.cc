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
// PETSCWavefunctionStepper.C
//
////////////////////////////////////////////////////////////////////////////////
// $Id: ExponentialWavefunctionStepper.C,v 1.8 2011-06-02 15:56:19 schleife Exp $

#include <config.h>

#include "PETSCWavefunctionStepper.h"
#include "SelfConsistentPotential.h"
#include "SlaterDet.h"
#include "Sample.h"
#include <iostream>
#include <deque>
#include <sys/stat.h>
#include <petscsys.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
PETSCWavefunctionStepper::PETSCWavefunctionStepper(Wavefunction& wf, double tddt, TimerMap& tmap, EnergyFunctional & ef, Sample & s, vector<vector<double> > & fion, valarray<double> & sigma_eks)
    : tddt_(tddt), WavefunctionStepper(wf,tmap), ef_(ef), s_(s), newwf_(s.wf), dwf_(s.wf) {
    TSAdapt adapt;
    // AK: context for RHS function
    petsc_ctx=new PETSC_CTX(ef_,s_,dwf_,petsc_dwf_vec,fion,sigma_eks,tmap,hamil_wf_vec);

//  PetscInitialize(0,0,NULL,NULL);             // AK: no command line options for now
    PetscInitializeNoArguments();
    RegisterMyRKC2();


    // EC: Force the nonlinear solver to use finite differences to approximate the Jacobian
    // expensive but almost as accurate as providing the anlytical Jacobian
    //PetscOptionsSetValue(NULL,"-snes_fd","true");
    // cheaper but less accurate
    PetscOptionsSetValue(NULL,"-snes_mf","true");
    PetscOptionsSetValue(NULL,"-snes_max_it","1000000");
    PetscOptionsSetValue(NULL,"-snes_max_funcs","10000000");

    // AK: flatten wf and create petsc vector that points to flat array
    wf_.allocate_flat();
    wf_.flatten();
    // AK: need to check communicator
    VecCreateMPIWithArray(s_.ctxt_.comm(), 1, wf_.flatarrsize, PETSC_DECIDE, wf_.flatarr, &petsc_wf_vec);

    const complex<double> * wf_vec_arr;
    VecGetArrayRead(petsc_wf_vec, &wf_vec_arr);
    if (s_.ctxt_.oncoutpe()) {
        cout << "wf_vec at beginning of RHS: " << endl;
        for (int i=0; i<10; i++)
            cout << wf_vec_arr[i] << " ";
        cout << endl;
    }
    VecRestoreArrayRead(petsc_wf_vec, &wf_vec_arr);

    //if (unit_evolve) {
    VecDuplicate(petsc_wf_vec,&hamil_wf_vec);
    VecCopy(petsc_wf_vec,hamil_wf_vec);
    const complex<double> * hamil_wf_vec_arr;
    VecGetArrayRead(hamil_wf_vec, &hamil_wf_vec_arr);
    if (s_.ctxt_.oncoutpe()) {
        cout << "hamil_wf_vec at beginning of RHS: " << endl;
        for (int i=0; i<10; i++)
            cout << hamil_wf_vec_arr[i] << " ";
        cout << endl;
    }
    VecRestoreArrayRead(hamil_wf_vec, &hamil_wf_vec_arr);
    //}

    // AK: also initialize dwf aka the working rhs vector
    dwf_.allocate_flat();
    dwf_.flatten();
    VecCreateMPIWithArray(s_.ctxt_.comm(), 1, dwf_.flatarrsize, PETSC_DECIDE, dwf_.flatarr, &petsc_dwf_vec);

    s_.hamil_wf->allocate_flat();

    // AK: set up time stepper
    TSCreate(PETSC_COMM_WORLD, &petsc_ts);      // AK: create petsc time stepper
    TSSetFromOptions(petsc_ts);
    TSSetProblemType(petsc_ts, TS_NONLINEAR);   // AK: problem has form U_t - A(U,t) U = 0
    TSSetSolution(petsc_ts, petsc_wf_vec);      // AK: where to put solution
    TSSetType(petsc_ts, s_.ctrl.petsc_ts_type.c_str()); // AK: input file provides time stepper type

    if (s_.ctrl.petsc_ts_type == "rk") {        // AK: set RK type; RK4 by default
        if (s_.ctrl.petsc_ts_subtype == "")
            TSRKSetType(petsc_ts,"4");
        else
            TSRKSetType(petsc_ts,s_.ctrl.petsc_ts_subtype.c_str());
    }

    if (s_.ctrl.petsc_ts_type == "ssp") {      // AK: set SSP type; RK104 by default
        if (s_.ctrl.petsc_ts_subtype == "")
            TSSSPSetType(petsc_ts,"rk104");
        else {
            string ssptype=s_.ctrl.petsc_ts_subtype;
            if (ssptype.back() == '2') {           // AK: 2nd order, variable stages
                TSSSPSetType(petsc_ts,"rks2");
                TSSSPSetNumStages(petsc_ts, stoi(ssptype.substr(2,ssptype.size()-3)));
            }
            else if (ssptype.back() == '3') {           // AK: 3rd order, variable stages
                TSSSPSetType(petsc_ts,"rks3");
                TSSSPSetNumStages(petsc_ts, stoi(ssptype.substr(2,ssptype.size()-3)));
            }
            else
                TSSSPSetType(petsc_ts,"rk104");      // AK: 4th order, 10 stages
        }
    }

    TSSetTimeStep(petsc_ts, tddt_);            // AK: input file provides time step
    TSSetExactFinalTime(petsc_ts, TS_EXACTFINALTIME_STEPOVER);  // Don't do anything if final time is exceeded
    TSGetAdapt(petsc_ts, &adapt);
    TSAdaptSetType(adapt, s_.ctrl.petsc_adapt.c_str());
    if (s_.ctrl.petsc_adapt != "none") {
        TSAdaptSetMonitor(adapt, PETSC_TRUE);
        if (s_.ctrl.petsc_adapt_atol != -1)
            TSSetTolerances(petsc_ts, s_.ctrl.petsc_adapt_atol, NULL, s_.ctrl.petsc_adapt_rtol, NULL);
    }
//  TSSetMaxSteps(petsc_ts,4);
//  TSSetMaxTime(petsc_ts, tddt_);	     // AK: do single time step

    // AK: output information
    if (s_.ctxt_.oncoutpe()) {
        cout << "initialized TS" << endl;
        cout << "TS type: " << s_.ctrl.petsc_ts_type.c_str() << endl;
        if (s_.ctrl.petsc_ts_type == "ssp") {
            int nstages;
            TSSSPGetNumStages(petsc_ts, &nstages);
            cout << "stages: " << nstages << endl;
        }
    }
    TSView(petsc_ts, PETSC_VIEWER_STDOUT_WORLD);

    // AK: not sure whether to give it NULL or petsc_dwf_vec
    //TSSetRHSFunction(petsc_ts, NULL, unit_RHS, petsc_ctx);
    //TSSetRHSFunction(petsc_ts, NULL, RHS, petsc_ctx);
    {
        if (s_.ctxt_.oncoutpe() && s_.ctrl.petsc_unit) {
            struct stat info;
            int rc=stat("evolve",&info);
            if (rc==-1) {
                mkdir("evolve", 0775);
                //rc = stat("evolve", &info);
            }
        }
        //    TSSetRHSFunction(petsc_ts, NULL, unit_evolve_RHS, petsc_ctx);
        TSSetRHSFunction(petsc_ts, NULL, RHS, petsc_ctx);
    }
    // AK: may read/write wfs between constructor and update, but probably don't need this
    VecGetArray(petsc_wf_vec, &wf_.flatarr);
    VecGetArray(petsc_dwf_vec, &dwf_.flatarr);


}


////////////////////////////////////////////////////////////////////////////////
PetscErrorCode PETSCWavefunctionStepper::dummy_RHS(TS ts, PetscReal t, Vec wf_vec, Vec rhs, void *ctx_) {
    cout << "inside RHS" << endl;

    PETSC_CTX * ctx=(PETSC_CTX*) ctx_;

    if (ctx->s_.ctxt_.oncoutpe())
        cout << "converted petsc_ctx, doing VecSet" << endl;

    // AK: U_t = 0 for testing
    VecSet(rhs,0);

    if (ctx->s_.ctxt_.oncoutpe())
        cout << "finished RHS" << endl;

    return 0;
}

////////////////////////////////////////////////////////////////////////////////

PetscErrorCode PETSCWavefunctionStepper::RHS(TS ts, PetscReal t, Vec wf_vec, Vec rhs, void *ctx_ ) {
    // AK: compute rhs=-i*H(t)*|psi> where |psi> is given by wf_vec

    PETSC_CTX * ctx=(PETSC_CTX*) ctx_;
    PetscInt step;
    TSGetStepNumber(ts,&step);
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout<<"current step number "<<step<<endl;
    }

#if PETSC_DEBUG
    if (ctx->s_.ctxt_.oncoutpe())
        cout << "inside petsc RHS" << endl;

    // AK: print wf_vec
    const complex<double> * wf_vec_arr;
    VecGetArrayRead(wf_vec, &wf_vec_arr);
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout << "wf_vec at beginning of RHS: " << endl;
        for (int i=0; i<10; i++)
            cout << wf_vec_arr[i] << " ";
        cout << endl;
    }
    VecRestoreArrayRead(wf_vec, &wf_vec_arr);
#endif

    // AK: set hamil_wf and wf to contents of wf_vec. turns out both are used in H*psi
    // AK: would be a problem if the RHS function is ever called on different vectors at the same time
    // AK: set_from_vec already contains VecGetArrayRead and VecRestoreArrayRead
    ctx->s_.hamil_wf->set_from_vec(wf_vec);
    ctx->s_.wf.set_from_vec(wf_vec);

#if PETSC_DEBUG
    // AK: print hamil_wf
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout << "hamil_wf in RHS: " << endl;
        for (int i=0; i<10; i++)
            cout << ctx->s_.hamil_wf->sd(0,0)->c()[i] << " ";
        cout << endl;
    }
#endif

    // calculate hamil_cd from hamil_wf and update hamil
    ctx->tmap_["charge"].start();
    ( ctx->ef_.hamil_cd() )->update_density();
    // do I need cd_.update_density(); here???
    ctx->tmap_["charge"].stop();

#if PETSC_DEBUG
    // AK: print hamil_cd
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout << "hamil_cd in RHS: " << endl;
        for (int i=0; i<10; i++)
            cout << ctx->ef_.hamil_cd()->rhog[0][i] << "  ";
        cout << endl;
    }
#endif

    ctx->tmap_["efn"].start();
    ctx->ef_.update_hamiltonian();
    ctx->ef_.update_vhxc();
    ctx->tmap_["efn"].stop();

#if PETSC_DEBUG
    // AK output fion and sigm_eks
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout << "fion in RHS:" << endl;
        for (int i=0; i<3; i++)
            cout << ctx->fion_[0][i] << " ";
        cout << endl;

        cout << "sigma_eks in RHS:" << endl;
        for (int i=0; i<6; i++)
            cout << ctx->sigma_eks_[i] << " ";
        cout << endl;
    }
#endif

    // AK: output v_r
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout << "v_r in RHS:" << endl;
        for (int i=0; i<10; i++)
            cout << ctx->ef_.v_r[0][i] << " ";
        cout << endl;
    }

    // AK: calculate dwf=H*psi
    // EnergyFunctional::energy(bool compute_hpsi, Wavefunction& dwf, bool compute_forces, vector<vector<double> >& fion, bool compute_stress, valarray<double>& sigma)
    ctx->ef_.energy(true,ctx->dwf_,false,ctx->fion_,false,ctx->sigma_eks_);	// fion, sigma_eks are passed in from EhrenSampleStepper

#if PETSC_DEBUG
    if (ctx->s_.ctxt_.oncoutpe())
        cout << "computed Hpsi" << endl;
#endif

    // AK: set contents of dwf to rhs

    VecGetArray(ctx->dwf_vec_, &(ctx->dwf_.flatarr));
    ctx->dwf_.flatten();
    VecRestoreArray(ctx->dwf_vec_, &(ctx->dwf_.flatarr));

    VecCopy(ctx->dwf_vec_, rhs);
    VecScale(rhs, -1.0*PETSC_i);

#if PETSC_DEBUG
    const complex<double> * rhs_arr;
    VecGetArrayRead(rhs, &rhs_arr);
    if (ctx->s_.ctxt_.oncoutpe()) {
        cout << "rhs at end of RHS:" << endl;
        for (int i=0; i<10; i++)
            cout << rhs_arr[i] << " ";
        cout << endl;
    }
    VecRestoreArrayRead(rhs, &rhs_arr);

    if (ctx->s_.ctxt_.oncoutpe())
        cout << "finished RHS" << endl;
#endif

    return 0;
}


/////////////////////////////////////////////////////////////////////////////////
PETSCWavefunctionStepper::~PETSCWavefunctionStepper() {
//  cout << "Destructor is called" << endl;
    VecDestroy(&petsc_wf_vec);
    VecDestroy(&petsc_dwf_vec);
    PetscFinalize();
}

/////////////////////////////////////////////////////////////////////////////////

void PETSCWavefunctionStepper::update(Wavefunction& dwf) {
//   const bool oncoutpe = s_.ctxt_.oncoutpe();

    TSView(petsc_ts, PETSC_VIEWER_STDOUT_WORLD);

    PetscInt step;
    TSGetStepNumber(petsc_ts,&step);
    if (s_.ctxt_.oncoutpe()) {
        cout<<"update current step number "<<step<<endl;
    }
    PetscInt N;
    VecGetSize(petsc_wf_vec,&N);
    if (s_.ctxt_.oncoutpe()) {
        cout << "doing petsc update" << endl;
//     TSView(petsc_ts, PETSC_VIEWER_STDOUT_SELF);

//   tmap_["petsc"].start();   // timing

#if PETSC_DEBUG

        cout << "wf_ before update: ";
        for (int i=0; i<10; i++)
            cout << wf_.flatarr[i] << "  ";
        cout<<"...";
        for (int i=N-10; i<N; i++)
            cout << wf_.flatarr[i] << "  ";
        cout << endl;
        cout << "dwf_ before update: ";
        for (int i=0; i<10; i++)
            cout << dwf_.flatarr[i] << "  ";
        cout<<"...";
        for (int i=N-10; i<N; i++)
            cout << dwf_.flatarr[i] << "  ";
        cout << endl;
#endif
    }
    if (s_.ctrl.petsc_unit) {
        const complex<double> * vec_arr;
        VecGetArrayRead(petsc_wf_vec, &vec_arr);
        if (s_.ctxt_.oncoutpe()) {
            cout<<"petsc_wf_vec before update:"<<endl;
            for (int i=0; i<10; i++)
                cout << vec_arr[i] << " ";
            cout<<"...";
            for (int i=N-10; i<N; i++)
                cout << vec_arr[i] << " ";
            cout << endl;
        }
        VecRestoreArrayRead(petsc_wf_vec, &vec_arr);
        Vec tmp;
        PetscErrorCode ierr=VecDuplicate(petsc_wf_vec,&tmp);
        if (s_.ctxt_.oncoutpe())
            cout<<"error code: "<<ierr<<endl;
        VecCopy(petsc_wf_vec,tmp);
        if(step==0) {
            PetscInt j=s_.ctrl.petsc_unit_k;
            if (j>=N) {
                if (s_.ctxt_.oncoutpe())
                    cout<<"index out of range. max index="<<N<<" current index="<<j<<endl;
                s_.ctxt_.abort(2);
            }
            PetscScalar eps=s_.ctrl.petsc_unit_eps;
            VecSet(tmp,0);
            VecSetValues(tmp,1,&j,&eps,INSERT_VALUES);
            VecAssemblyBegin(tmp);
            VecAssemblyEnd(tmp);
            const complex<double> * tmp_arr;
            VecGetArrayRead(tmp, &tmp_arr);
            if (s_.ctxt_.oncoutpe()) {
                cout << "tmp:" << endl;
                for (int i=0; i<10; i++)
                    cout << tmp_arr[i] << " ";
                cout<<"...";
                //for (int i=N-10; i<N; i++)
                //    cout << tmp_arr[i] << " ";
                cout << endl;
            }
            VecRestoreArrayRead(tmp, &tmp_arr);
            PetscScalar sum=0;
            VecSum(tmp,&sum);
            if (s_.ctxt_.oncoutpe())
                cout<<"index at "<<j<<" sum tmp: "<<sum<<endl;
            VecAXPY(petsc_wf_vec,1,tmp);

            PetscViewer   output; /* file to output data to */
            char name[50];
            sprintf(name,"evolve/gs_%d.dat",s_.ctrl.petsc_unit_k);
            PetscViewerBinaryOpen(PETSC_COMM_WORLD,name,FILE_MODE_WRITE,&output);
            VecView(petsc_wf_vec, output);
            if (s_.ctxt_.oncoutpe())
                cout<<"save tmp"<<endl;
            PetscViewerDestroy(&output);
        }
    }

    // AK: qball may read/write wfs between updates. probably don't need this
    VecRestoreArray(petsc_wf_vec, &wf_.flatarr);
    VecRestoreArray(petsc_dwf_vec, &dwf_.flatarr);

#if PETSC_DEBUG
    PetscReal time;
    TSGetTime(petsc_ts, &time);

    if (s_.ctxt_.oncoutpe()) {
        cout << "current time " << time << endl;
        cout << "next time " << time+tddt_ << endl;
        cout << "starting TSStep" << endl;
    }
#endif
    TSSetMaxTime(petsc_ts, time+tddt_);
//   TSSolve(petsc_ts, petsc_wf_vec);

    if (s_.ctxt_.oncoutpe())
        cout << "woc TSStep" << endl;
    // AK: advance once time step
    TSStep(petsc_ts);

#if PETSC_DEBUG
    if (s_.ctxt_.oncoutpe())
        cout << "finished TSStep" << endl;
#endif

    // AK: qball may read/write wfs between update
    VecGetArray(petsc_wf_vec, &wf_.flatarr);
    VecGetArray(petsc_dwf_vec, &dwf_.flatarr);

#if PETSC_DEBUG
    if (s_.ctxt_.oncoutpe()) {
        cout << "wf_ after update: ";
        for (int i=0; i<10; i++)
            cout << wf_.flatarr[i] << "  ";
        cout << endl;
    }
#endif

    Vec petsc_sol;
    TSGetSolution(petsc_ts, &petsc_sol);

    wf_.set_from_vec(petsc_sol);

    Vec rhs;
    VecDuplicate(petsc_dwf_vec,&rhs);
    VecCopy(petsc_dwf_vec, rhs);
    VecScale(rhs, -1.0*PETSC_i);
    // output the rhs
    if(s_.ctrl.petsc_unit && step%s_.ctrl.petsc_savefreq==0) {
        if (s_.ctxt_.oncoutpe()) 
            cout<<"save unit rhs at "<<s_.ctrl.mditer<<endl;
        PetscViewer   outputfile; /* file to output data to */
        char filename[50];
        PetscViewerCreate(PETSC_COMM_WORLD, &outputfile);
        PetscViewerSetType(outputfile, PETSCVIEWERBINARY);
        sprintf(filename, "evolve/unit_%d_%d.dat",s_.ctrl.petsc_unit_k,step);
        PetscViewerFileSetMode(outputfile, FILE_MODE_WRITE);
        PetscViewerFileSetName(outputfile, filename);
        VecView(rhs, outputfile);
        PetscViewerDestroy(&outputfile);
    }

#if PETSC_DEBUG
    if (s_.ctxt_.oncoutpe()) {
        cout << "wf_ after update, set from petsc_sol: ";
        for (int i=0; i<10; i++)
            cout << wf_.flatarr[i] << "  ";
        cout << endl;

        cout << "dwf_ after update: ";
        for (int i=0; i<10; i++)
            cout << dwf_.flatarr[i] << "  ";
        cout << endl;
    }
#endif

//   wf_.unflatten();	// copy wf_.flatarr (aka petsc_wf_vec) into wf_.sd_

    /*
       // AK: test methods to flatten and unflatten wf
       wf_.set_from_vec(petsc_wf_vec);		// copy vec into wf_.sd_
       wf_.flatten();				// copy wf_.sd_ into wf_.flatarr, also array held by petsc_wf_vec
       wf_.unflatten();				// copy wf_.flatarr back into wf_.sd_
    */

//   tmap_["petsc"].stop();   // timing
}
//

PetscErrorCode PETSCWavefunctionStepper::RegisterMyRKC2(void) {
    {
        const PetscReal A[5][5] = {{0, 0, 0, 0, 0}, {0.0315862022077577455912627, 0, 0, 0,0}, {-0.128686897689645502255899, 0.255809213036559752204734, 0, 0,0}, {-0.562734368899226551546441, 0.601462745045115244568001,0.298892036910798860068196, 0, 0}, {-1.19932558300501514367783,0.942407757213538464197643, 0.621901167048249099059611,0.264501915438166696455688, 0}};
        const PetscReal B[5]= {-2.00026220580022632018238, 1.27250824704230238438162,0.940925234639918312902377, 0.531421930235961355176752, 0.255406793882044267721627};
        PetscCall(TSRKRegister("myrkc2",2,5,&A[0][0],B,NULL,NULL,0,NULL));
    }
    PetscFunctionReturn(0);
}
