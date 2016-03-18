Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> pThermo
(
    rhoThermo::New(mesh)
);
rhoThermo& thermo = pThermo();
thermo.correct();

volScalarField& p = thermo.p();
volScalarField& h = thermo.he();

const volScalarField& psi = thermo.psi();
const volScalarField& T = thermo.T();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField Ek
(
    "Ek",
    0.5*magSqr(U)
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField rhoHat
(
    IOobject
    (
        "rhoHat",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    rho - psi*p
);

volScalarField EkChange
(
    "EkChange",
    rho*Ek / runTime.deltaT() * 0.0
);

volScalarField dpdt
(
    "dpdt",
    p / runTime.deltaT() * 0.0
);

volScalarField dissip
(
    "dissip",
    p / runTime.deltaT() * 0.0
);

volVectorField gradp
(
    "gradp",
    fvc::grad(p)
);
volVectorField limGradP 
(
    "limGradP",
    fvc::grad(p)
);

DynamicList<label> badQualityCells(0);

surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi
(
    IOobject
    (
	"phi",
	runTime.timeName(),
	mesh,
	IOobject::READ_IF_PRESENT,
	IOobject::AUTO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);


