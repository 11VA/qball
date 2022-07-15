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
// PETScAdapt.h
//
////////////////////////////////////////////////////////////////////////////////

#include <config.h>

#ifndef PETSCADAPT_H
#define PETSCADAPT_H

#include<iostream>
#include<iomanip>
#include<sstream>
#include<stdlib.h>

#include <qball/SlaterDet.h>

class PETScAdapt : public Var {
    Sample *s;

public:

    char const*name ( void ) const {
        return "petsc_adapt";
    };

    int set ( int argc, char **argv ) {
        if ( argc < 2 || argc > 4 ) {
            if ( ui->oncoutpe() )
                cout << " <ERROR> petsc_adapt takes 1-3 values </ERROR>" << endl;
            return 1;
        }

        string v = argv[1];
        if ( !(v == "basic" || v == "none") ) {
            if ( ui->oncoutpe() )
                cout << " <ERROR> petsc_adapt must be in [basic,none] </ERROR>" << endl;
            return 1;
        }

        s->ctrl.petsc_adapt = v;

        // AK: set adapter tolerances
        // AK: maybe don't allow this for "none" adapter?
        if (argc == 3) {
            s->ctrl.petsc_adapt_atol = atof(argv[2]);
            s->ctrl.petsc_adapt_rtol = atof(argv[2]);
        }

        if (argc == 4) {
            s->ctrl.petsc_adapt_atol = atof(argv[2]);
            s->ctrl.petsc_adapt_rtol = atof(argv[3]);
        }

        return 0;
    }

    string print (void) const {
        ostringstream st;
        st.setf(ios::left,ios::adjustfield);
        st << setw(10) << name() << " = ";
        st.setf(ios::right,ios::adjustfield);
        st << setw(10) << s->ctrl.petsc_adapt;
        return st.str();
    }

    PETScAdapt(Sample *sample) : s(sample) {
        s->ctrl.petsc_adapt = "none";
    };
};
#endif

// Local Variables:
// mode: c++
// End:
