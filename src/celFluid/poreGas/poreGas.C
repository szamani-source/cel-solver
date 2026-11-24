/*---------------------------------------------------------------------------*\
    poreGas: gas fractions and permeability (WP F.1, decisions 0125, 0126)
\*---------------------------------------------------------------------------*/

#include "poreGas.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Local helpers * * * * * * * * * * * * * * //

namespace
{

Foam::IOobject gasPropertiesIO(const Foam::fvMesh& mesh)
{
    return Foam::IOobject
    (
        "gasProperties",
        mesh.time().constant(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
    );
}


Foam::IOobject fractionIO(const Foam::word& name, const Foam::fvMesh& mesh)
{
    // Computed from rho and alpha, never read (a restart recomputes them)
    return Foam::IOobject
    (
        name,
        mesh.time().timeName(),
        mesh,
        Foam::IOobject::NO_READ,
        Foam::IOobject::AUTO_WRITE
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * Static functions  * * * * * * * * * * * * * //

Foam::scalar Foam::poreGas::sutherland(const scalar T)
{
    // tools/gas_reference.py: mu0 1.716e-5 Pa s, T0 273.15 K, S 110.4 K
    const scalar mu0 = 1.716e-5, T0 = 273.15, S = 110.4;
    return mu0*Foam::pow(T/T0, 1.5)*(T0 + S)/(T + S);
}


Foam::scalar Foam::poreGas::permeability(const scalar e, const scalar d)
{
    if (e <= 0)
    {
        return 0;
    }
    if (e < eSwitch)
    {
        // Ergun's viscous term
        return e*e*e*d*d/(150*(1 - e)*(1 - e));
    }
    if (e >= 1)
    {
        return GREAT;
    }
    // Wen and Yu with Schiller-Naumann Cd at Re -> 0:
    // K = 18 mu (1-e) e^-2.65 / d^2, k = mu e^2 / K
    return d*d*Foam::pow(e, 4.65)/(18*(1 - e));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::poreGas::poreGas(const eulerianSolid& solid, const IOobject& io)
:
    solid_(solid),
    dict_(io),
    rhoP_(dict_.get<scalar>("rhoP")),
    d_(dict_.get<scalar>("d")),
    T_(dict_.getOrDefault<scalar>("T", 293.15)),
    p0_(dict_.getOrDefault<scalar>("p0", 101325)),
    drag_(dict_.getOrDefault<word>("drag", "ergunWenYu")),
    rhoG_(p0_/(Rair*T_)),
    mu_(sutherland(T_)),
    phiS_
    (
        fractionIO("phiS", solid.mesh()),
        solid.mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    phiG_
    (
        fractionIO("phiG", solid.mesh()),
        solid.mesh(),
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    k_
    (
        fractionIO("k", solid.mesh()),
        solid.mesh(),
        dimensionedScalar(dimArea, Zero),
        fvPatchFieldBase::zeroGradientType()
    )
{
    if (rhoP_ <= 0 || d_ <= 0 || T_ <= 0 || p0_ <= 0)
    {
        FatalIOErrorInFunction(dict_)
            << "rhoP, d, T and p0 must be positive" << exit(FatalIOError);
    }
    if (drag_ != "ergunWenYu")
    {
        // Decision 0125 item 7: F.2 adds the run-time selector, F.7 the
        // other laws (ergun, wenYu, diFelice, syamlalOBrien, beetstra)
        FatalIOErrorInFunction(dict_)
            << "drag " << drag_ << ": only ergunWenYu is available yet"
            << exit(FatalIOError);
    }
    if (solid.rho0() > rhoP_)
    {
        WarningInFunction
            << "materialProperties rho " << solid.rho0()
            << " is above the grain density rhoP " << rhoP_
            << ": the packing would exceed 1" << endl;
    }

    const label oldPrecision = Info().precision(8);
    Info<< "Pore gas (WP F.1): air at T = " << T_ << " K, p0 = " << p0_
        << " Pa: rho_g = " << rhoG_ << " kg/m^3, mu = " << mu_ << " Pa s;"
        << nl << "    grains rhoP = " << rhoP_ << " kg/m^3, d = " << d_
        << " m, drag " << drag_ << "; set packing rho/rhoP = "
        << solid.rho0()/rhoP_ << nl << endl;
    Info().precision(oldPrecision);

    correctFractions();
}


Foam::autoPtr<Foam::poreGas> Foam::poreGas::New(const eulerianSolid& solid)
{
    IOobject io(gasPropertiesIO(solid.mesh()));
    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return nullptr;
    }

    const word model = IOdictionary(io).getOrDefault<word>("model", "poreGas");
    if (model == "none")
    {
        Info<< "Pore gas: model none (as without gasProperties)" << nl << endl;
        return nullptr;
    }
    if (model != "poreGas")
    {
        FatalErrorInFunction
            << "gasProperties model " << model
            << " unknown; use none or poreGas" << exit(FatalError);
    }
    return autoPtr<poreGas>::New(solid, io);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::poreGas::correctFractions()
{
    const scalarField& rho = solid_.rho().primitiveField();
    const scalarField& alpha = solid_.alpha().primitiveField();
    const tmp<scalarField> tOpen(solid_.cellOpen());
    const scalarField& open = tOpen();

    scalarField& phiS = phiS_.primitiveFieldRef();
    scalarField& phiG = phiG_.primitiveFieldRef();
    scalarField& k = k_.primitiveFieldRef();

    forAll(rho, celli)
    {
        if (open[celli] < eulerianSolid::coveredOpen)
        {
            // Inside a body: no grains and no gas
            phiS[celli] = 0;
            phiG[celli] = 0;
            k[celli] = 0;
            continue;
        }

        // Fractions of the open volume: alpha rho_m / rhoP = rho / rhoP
        phiS[celli] = rho[celli]/rhoP_;
        phiG[celli] = 1 - phiS[celli];

        // Permeability of the filled part (decision 0125 item 4)
        k[celli] =
            alpha[celli] < alphaMin
          ? 0
          : permeability(1 - rho[celli]/(alpha[celli]*rhoP_), d_);
    }

    phiS_.correctBoundaryConditions();
    phiG_.correctBoundaryConditions();
    k_.correctBoundaryConditions();
}


// ************************************************************************* //
