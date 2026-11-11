Info<< "Reading materialProperties" << nl << endl;

IOdictionary materialProperties
(
    IOobject
    (
        "materialProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    )
);

const dimensionedScalar rho("rho", dimDensity, materialProperties);
const dimensionedScalar E("E", dimPressure, materialProperties);
const dimensionedScalar nu("nu", dimless, materialProperties);

// Bulk and shear moduli, and the elastic P-wave speed
const dimensionedScalar K(E/(3.0*(1.0 - 2.0*nu)));
const dimensionedScalar G(E/(2.0*(1.0 + nu)));
const dimensionedScalar waveSpeed(sqrt((K + 4.0/3.0*G)/rho));

const scalar maxCo =
    runTime.controlDict().getOrDefault<scalar>("maxCo", 0.5);

Info<< "Reading field U" << nl << endl;

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
