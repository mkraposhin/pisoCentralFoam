/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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
    pisoCentralFoam

Description
    Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "cellQuality.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    
    #include "createTime.H"
    #include "createMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "createFields.H"
    #include "createFvOptions.H"
    #include "createMRF.H"
    #include "initContinuityErrs.H"
    #include "readCourantType.H"
    
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    
    #include "createSurfaceFields.H"
    surfaceScalarField rhoHat_pos
    (
        "rhoHat_pos",
        fvc::interpolate(rhoHat, pos, "reconstruct(rhoHat)")
    );
    surfaceScalarField rhoHat_neg
    (
         "rhoHat_neg",
         fvc::interpolate(rhoHat, neg, "reconstruct(rhoHat)")
    );
    surfaceScalarField rhoHatPhi_pos
    (
        "rhoHatPhi_pos",
        phiPos*0.0
    );
    surfaceScalarField rhoHatPhi_neg
    (
         "rhoHatPhi_neg",
         phiNeg*0.0
    );
    surfaceScalarField rhof
    (
        "rhof",
        p_pos*psi_pos*a_pos + rhoHat_pos*a_pos + p_neg*psi_neg*a_neg + rhoHat_neg*a_neg
    );

    #include "markBadQualityCells.H"
    
    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    Info<< "\nStarting time loop\n" << endl;
    
    #include "initKappaField.H"
    
    while (runTime.run())
    {
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "acousticCourantNo.H"
            #include "compressibleCourantNo.H"
            #include "readTimeControls.H"
            #include "setDeltaT.H"
        }
        
        runTime++;
        
        rho.oldTime();
        rhoHat.oldTime();
        p.oldTime();
        U.oldTime();
        h.oldTime();
        Ek.oldTime();
        psi.oldTime();
        
        Info<< "Time = " << runTime.timeName() << nl << endl;
        
        // --- Solve density
        #include "systemDensity.H"
        
        // --- Solve momentum
        #include "UEqn.H"
        
        // --- Solve energy
        #include "hEqn.H"
        
        // --- Solve pressure (PISO)
        {
            while (pimple.correct())
            {
                #include "pEqn.H"
            }
            
            #define PISOCENTRALFOAM_LTS
            #include "updateKappa.H"
        }
        
        // --- Solve turbulence
        turbulence->correct();
        
        Ek = 0.5*magSqr(U);
        EkChange = fvc::ddt(rho,Ek) + fvc::div(phiPos,Ek) + fvc::div(phiNeg,Ek);
        dpdt = fvc::ddt(p);
        
        runTime.write();
        
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
