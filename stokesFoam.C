/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    simpleFoam

Description
    Steady-state solver for incompressible, turbulent flow, using the SIMPLE
    algorithm.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    double eqnResidual = 1.0;
    
    while (eqnResidual > 0.001)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        p.storePrevIter();

        tmp<fvVectorMatrix> tUEqn
        (
           -fvm::laplacian(nu, U)
        );

        fvVectorMatrix& UEqn = tUEqn.ref();
        
        UEqn.relax();

        solve(UEqn== - fvc::grad(p));

        p.correctBoundaryConditions();

        volScalarField AU = UEqn.A();
        U = UEqn.H()/AU;
        U.correctBoundaryConditions();
        // UEqn.clear();

        phi = fvc::interpolate(U) & mesh.Sf();
        //adjustPhi(phi, U, p);

        fvScalarMatrix pEqn
        (
           fvm::laplacian(1.0/AU, p) == fvc::div(phi)
        );
        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve();

        phi -= pEqn.flux();

        # include "continuityErrs.H"

        p.relax();
        U -= fvc::grad(p)/AU;
        U.correctBoundaryConditions();
        
        runTime.write();

        eqnResidual = pEqn.solve().initialResidual();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
