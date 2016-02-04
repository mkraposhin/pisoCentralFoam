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
    pisoCentralDyMFoam

Description
    Pressure-based semi implicit compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor with support of moving meshes

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "pimpleControl.H"
#include "turbulentFluidThermoModel.H"
#include "zeroGradientFvPatchFields.H"
#include "coupledFvsPatchFields.H"
#include "cellQuality.H"
#include "fvIOoptionList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createFields.H"
    #include "createMRF.H"
    #include "createTimeControls.H"
    #include "createFvOptions.H"
    bool checkMeshCourantNo =
            readBool(pimple.dict().lookup("checkMeshCourantNo"));
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
    #include "initContinuityErrs.H"
    #include "readCourantType.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;
    
    #include "createSurfaceFields.H"
    #include "markBadQualityCells.H"
    
    surfaceScalarField meshPhi
    (
	"volMeshPhi",
	phiv_pos * 0.0
    );
    
    surfaceScalarField rel_phiv_pos
    (
        "rel_phiv_pos",
        phiv_pos - meshPhi
    );
    
    surfaceScalarField rel_phiv_neg
    (
        "rel_phiv_neg",
        phiv_neg - meshPhi
    );
    
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
        
        // --- Move mesh and update fluxes
        {
           // Do any mesh changes
           mesh.update();
           
            if (mesh.changing())
            {
                meshPhi = fvc::meshPhi(rho,U)();
                
                if (runTime.timeIndex() > 1)
                {
                    surfaceScalarField amNew = min(min(phiv_pos - meshPhi - cSf_pos, phiv_neg - meshPhi - cSf_neg), v_zero);
                    phiNeg += kappa*(amNew - am)*p_neg*psi_neg;
                    phiPos += (1.0 - kappa)*(amNew - am)*p_neg*psi_neg;
                }
                else
                {
                    phiNeg -= meshPhi * fvc::interpolate(rho);
                }
                
                phi = phiPos + phiNeg;
                
                if (checkMeshCourantNo)
                {
                    #include "customMeshCourantNo.H"
                }
                
                #include "markBadQualityCells.H"
            }
        }
        
        // --- Solve density
        {
            fvScalarMatrix rhoEqn
            (
                fvm::ddt(rho) + fvc::div(phi)
                ==
                fvOptions(rho)
            );
            
            fvOptions.constrain(rhoEqn);
            
            rhoEqn.solve();
            
            fvOptions.correct(rho);

            Info<< "rho max/min : " << max(rho).value()
            << " / " << min(rho).value() << endl;
        }
        
        
        // --- Solve momentum
        #include "UEqn.H"
        
        // --- Solve energy
        #include "hEqn.H"
        
        // --- Solve pressure (PISO)
        {
            while (pimple.correct())
            {
                #include "pEqnDyM.H"
            }
            #include "updateKappa.H"
        }
        
        // --- Solve turbulence
        turbulence->correct();
        
        Ek = 0.5*magSqr(U);
        EkChange = fvc::ddt(rho,Ek) + fvc::div(phiPos,Ek) + fvc::div(phiNeg,Ek);
        dpdt = fvc::ddt(p) - fvc::div(meshPhi, p);
        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
