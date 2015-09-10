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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    
    pimpleControl pimple(mesh);
    
    #include "createTimeControls.H"
    #include "createFields.H"
    #include "readTimeControls.H"
    bool checkMeshCourantNo =
            readBool(pimple.dict().lookup("checkMeshCourantNo"));
    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    #include "readFluxScheme.H"

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    Info<< "\nStarting time loop\n" << endl;
    
    #include "createSurfaceFields.H"
    
    surfaceScalarField kappa
    (
	IOobject
	(
	    "kappa",
	    runTime.timeName(),
	    mesh,
	    IOobject::NO_READ,
	    IOobject::NO_WRITE
	),
	mesh,
	dimless
    );
    
    forAll(kappa, iFace)
    {
	kappa[iFace] = 0.0;
    }
    forAll(kappa.boundaryField(), iPatch)
    {
	kappa.boundaryField()[iPatch] = 0.0;
    }
    
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
		surfaceScalarField amNew = min(min(phiv_pos - fvc::meshPhi(rho,U) - cSf_pos, phiv_neg - fvc::meshPhi(rho,U) - cSf_neg), v_zero);
		phiNeg += kappa*(amNew - am)*p_neg*psi_neg;
		phiPos += (1.0 - kappa)*(amNew - am)*p_neg*psi_neg;
		
		phi = phiPos + phiNeg;
		
		if (checkMeshCourantNo)
		{
		    #include "meshCourantNo.H"
		}
	    }
	}
	
	// --- Solve density
	solve
	(
	    fvm::ddt(rho) + fvc::div(phi)
	);
	Info<< "rhoEqn max/min : " << max(rho).value()
	<< " " << min(rho).value() << endl;


	// --- Solve momentum
	#include "UEqn.H"
	
	// --- Solve energy
	#include "hEqn.H"
	
	// --- Solve pressure (PISO)
	{
	    while (pimple.correct())
	    {
		volVectorField HbyA ("HbyA", U);
		volScalarField rAU ("rAU", 1.0 / UEqn.A());
		HbyA = UEqn.H() * rAU;

		psi_pos = fvc::interpolate(psi, pos, "reconstruct(psi)");
		psi_neg = fvc::interpolate(psi, neg, "reconstruct(psi)");

		psiU_pos= fvc::interpolate(psi*HbyA, pos, "reconstruct(U)");
		psiU_neg= fvc::interpolate(psi*HbyA, neg, "reconstruct(U)");
		
		phiv_pos= (psiU_pos / psi_pos) & mesh.Sf();
		phiv_neg= (psiU_neg / psi_neg) & mesh.Sf();
        
		c = sqrt(thermo.Cp()/thermo.Cv() / psi);

		cSf_pos = fvc::interpolate(c, pos, "reconstruct(psi)")*mesh.magSf();
		cSf_neg = fvc::interpolate(c, neg, "reconstruct(psi)")*mesh.magSf();

		ap = max(max(phiv_pos - fvc::meshPhi(rho,U) + cSf_pos, phiv_neg - fvc::meshPhi(rho,U) + cSf_neg), v_zero);
		am = min(min(phiv_pos - fvc::meshPhi(rho,U) - cSf_pos, phiv_neg - fvc::meshPhi(rho,U) - cSf_neg), v_zero);

		a_pos = ap/(ap - am);
		aSf = am*a_pos;
		a_neg = 1.0 - a_pos;

		aphiv_pos = a_pos * (phiv_pos - fvc::meshPhi(rho,U)) - aSf;
		aphiv_neg = a_neg * (phiv_neg - fvc::meshPhi(rho,U)) + aSf;
	
		phid_pos = aphiv_pos * psi_pos;
		phid_neg = aphiv_neg * psi_neg;
	
		surfaceScalarField Dp_pos
		(
		    "Dp_pos",
		    fvc::interpolate(rho*rAU, pos, "reconstruct(Dp)")
		);
		surfaceScalarField Dp_neg
		(
		    "Dp_neg",
		    fvc::interpolate(rho*rAU, neg, "reconstruct(Dp)")
		);

		while (pimple.correctNonOrthogonal())
		{
		    fvScalarMatrix pEqn_pos
		    (
			fvm::div(phid_pos,p) - fvm::laplacian(Dp_pos*a_pos,p)
		    );
		
		    fvScalarMatrix pEqn_neg
		    (
			fvm::div(phid_neg,p) - fvm::laplacian(Dp_neg*a_neg,p)
		    );
		    
		    solve
		    (
			fvm::ddt(psi,p)
			+
			pEqn_pos
			+
			pEqn_neg,
			mesh.solver(p.select(pimple.finalInnerIter()))
		    );
		    
		    if (pimple.finalNonOrthogonalIter())
		    {
			phiPos = pEqn_pos.flux();
			phiNeg = pEqn_neg.flux();
		    }
		}
		
		p_pos = fvc::interpolate(p, pos, "reconstruct(p)");
		p_neg = fvc::interpolate(p, neg, "reconstruct(p)");
		
		phiv_pos= phiPos / (p_pos*psi_pos);
		phiv_neg= phiNeg / (p_neg*psi_neg);
		
		ap = max(max(phiv_pos + cSf_pos, phiv_neg + cSf_neg), v_zero);
		am = min(min(phiv_pos - cSf_pos, phiv_neg - cSf_neg), v_zero);
		
		a_pos = ap/(ap - am);
		a_neg = 1.0 - a_pos;
		
		//gradp = fvc::grad(p);
		gradp = fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf());
		U = HbyA - rAU * gradp;
		U.correctBoundaryConditions();
		
		Info << "max(U): " << max(U).value() << endl;
		
		rho = thermo.rho();
	    }
	    
	    aSf = am*a_pos;
	    phiv_pos *= a_pos;
	    phiv_neg *= a_neg;
	    aphiv_pos = phiv_pos - aSf;
	    aphiv_neg = phiv_neg + aSf;
	    amaxSf = max(mag(aphiv_pos), mag(aphiv_neg));
	
	    surfaceScalarField amaxSfbyDelta
	    (
		mesh.surfaceInterpolation::deltaCoeffs()*amaxSf
	    );
	    
	    surfaceScalarField Maf
	    (
		
		mag(phi) / (psi_pos*p_pos*a_pos + psi_neg*p_neg*a_neg)
		/ (cSf_pos*a_pos + cSf_neg*a_neg)
	    );
	    
	    Info << "max/min Maf: " << max(Maf).value() << "/" << min(Maf).value() << endl;
	    
	    kappa = 
		min
		(
		    Maf / (amaxSfbyDelta/mesh.magSf() * runTime.deltaT()),
		    scalar(1.0)
		);
	    
	    forAll(kappa.boundaryField(), iPatch)
	    {
		fvsPatchField<scalar>& kappapf = kappa.boundaryField()[iPatch];
		if (isA<coupledFvsPatchField<scalar> > (kappapf))
		{
		    forAll (kappapf, iFace)
		    {
			kappapf[iFace] = 0.0;
		    }
		}
	    }
	    
	    Info << "max / min kappa: " << max(kappa).value() << "/" << min(kappa).value() << endl;
	    
	    phiPos = phiPos + (1.0 - kappa) * phiNeg;
	    phiNeg = kappa * phiNeg;
	    phi = phiPos + phiNeg;
	}
	
	// --- Solve turbulence
	turbulence->correct();
	
	Ek = 0.5*magSqr(U);
	EkChange = fvc::ddt(rho,Ek) + fvc::div(phiPos,Ek) + fvc::div(phiNeg,Ek);
	dpdt = fvc::ddt(p) - fvc::div(fvc::meshPhi(rho,U), p);
	
	runTime.write();

	Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
	    << "  ClockTime = " << runTime.elapsedClockTime() << " s"
	    << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
