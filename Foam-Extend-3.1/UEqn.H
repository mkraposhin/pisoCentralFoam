
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
	UEqn == -gradp
    );
}
