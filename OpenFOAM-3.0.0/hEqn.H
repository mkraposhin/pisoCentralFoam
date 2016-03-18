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
    dissip
    +
    fvOptions(rho,h)
);

fvOptions.constrain(hEqn);

hEqn.solve();

fvOptions.correct(h);

thermo.correct();

Info << "max/min T:" << max(T).value() << "/" << min(T).value() << endl;
