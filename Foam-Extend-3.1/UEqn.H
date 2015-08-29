volVectorField gradP
(
    "gradP",
    fvc::div((a_pos*p_pos + a_neg*p_neg)*mesh.Sf())
);

fvVectorMatrix UEqn
(
    fvm::ddt(rho,U)
    +
    fvm::div(phiPos,U)
    +
    fvm::div(phiNeg,U)
    +
    turbulence->divDevRhoReff(U)
);

if (momentumPredictor)
{
    solve
    (
	UEqn == -gradP
    );
}
