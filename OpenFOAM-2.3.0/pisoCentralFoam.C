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
#include "psiThermo.H"
#include "pimpleControl.H"
#include "turbulenceModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "cellQuality.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    
    #include "createTime.H"
    #include "createMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createFields.H"
    #include "initContinuityErrs.H"
    #include "readTimeControls.H"
    
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    #include "readCourantType.H"
    
    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);
    
    Info<< "\nStarting time loop\n" << endl;
    
    #include "createSurfaceFields.H"
    #include "markBadQualityCells.H"
    
    while (runTime.run())
    {
	#include "acousticCourantNo.H"
	#include "compressibleCourantNo.H"
	#include "readTimeControls.H"
	#include "setDeltaT.H"
	
	runTime++;
        
	psi.oldTime();
	rho.oldTime();
	p.oldTime();
	U.oldTime();
	h.oldTime();
	
	Info<< "Time = " << runTime.timeName() << nl << endl;
	
	// --- Solve density
	solve
	(
	    fvm::ddt(rho) + fvc::div(phi)
	);
	
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
