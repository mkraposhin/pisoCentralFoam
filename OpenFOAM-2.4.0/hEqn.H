fvScalarMatrix hEqn
(
    fvm::ddt(rho,h)
    +
    fvm::div(phiPos,h)
    +
    fvm::div(phiNeg,h)
    -
    fvm::laplacian(turbulence->alphaEff(), h)
    +
    EkChange
    ==
    dpdt
    +
    fvc::div( ((-turbulence->devRhoReff()) & U) )
);

hEqn.solve();

thermo.correct();

Info << "max/min T:" << max(T).value() << "/" << min(T).value() << endl;
