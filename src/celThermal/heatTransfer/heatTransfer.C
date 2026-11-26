/*---------------------------------------------------------------------------*\
    Temperature in the granular material (WP T.1, decision 0135).
    See heatTransfer.H.
\*---------------------------------------------------------------------------*/

#include "heatTransfer.H"
#include "fvMesh.H"
#include "Time.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

namespace
{

Foam::IOobject thermalPropertiesIO(const Foam::fvMesh& mesh)
{
    return Foam::IOobject
    (
        "thermalProperties",
        mesh.time().constant(),
        mesh,
        Foam::IOobject::MUST_READ,
        Foam::IOobject::NO_WRITE
    );
}

} // End anonymous namespace


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::heatTransfer::heatTransfer
(
    eulerianSolid& solid,
    const IOobject& io
)
:
    solid_(solid),
    mesh_(solid.mesh()),
    dict_(io),
    model_(dict_.get<word>("model")),
    c_(dict_.get<scalar>("c")),
    conductivity_(dict_.getOrDefault<word>("conductivity", "constant")),
    k_(dict_.get<scalar>("k")),
    T_
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    heatIn_(0),
    energy0_(0),
    carriedIn_(0),
    TOld_(),
    limiter_(dict_.getOrDefault<word>("limiter", "superbee")),
    voidT_(dict_.getOrDefault<scalar>("voidT", 293.15)),
    nSub_(1)
{
    if (solid_.immersed())
    {
        FatalIOErrorInFunction(dict_)
            << "temperature with immersed (STL) walls or bodies is not yet "
            << "supported (WP T.2b)" << exit(FatalIOError);
    }
    solid_.addPassenger(T_, limiter_);
    solid_.storeRemapMass();

    if (model_ != "equilibrium")
    {
        FatalIOErrorInFunction(dict_)
            << "thermalProperties model " << model_ << " unknown; use none or "
            << "equilibrium (twoTemperature comes with WP T.4)"
            << exit(FatalIOError);
    }
    if (conductivity_ != "constant")
    {
        FatalIOErrorInFunction(dict_)
            << "thermalProperties conductivity " << conductivity_
            << " unknown; use constant" << exit(FatalIOError);
    }
    if (c_ <= 0 || k_ < 0)
    {
        FatalIOErrorInFunction(dict_)
            << "c must be positive and k not negative" << exit(FatalIOError);
    }

    energy0_ = energy();

    const label oldPrecision = Info().precision(16);
    Info<< "Temperature (WP T.1): model " << model_ << ", c = " << c_
        << " J/(kg K), conductivity " << conductivity_ << " k = " << k_
        << " W/(m K)" << nl
        << "Thermal energy at start = " << energy0_ << " J" << nl
        << "Conduction step limit   = " << stableDeltaT() << " s" << nl << endl;
    Info().precision(oldPrecision);
}


Foam::autoPtr<Foam::heatTransfer> Foam::heatTransfer::New
(
    eulerianSolid& solid
)
{
    IOobject io(thermalPropertiesIO(solid.mesh()));
    if (!io.typeHeaderOk<IOdictionary>(true))
    {
        return nullptr;
    }

    const word model = IOdictionary(io).getOrDefault<word>("model", "none");
    if (model == "none")
    {
        Info<< "Temperature: model none (as without thermalProperties)" << nl
            << endl;
        return nullptr;
    }
    return autoPtr<heatTransfer>::New(solid, io);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::heatTransfer::fill() const
{
    if (solid_.freeSurface())
    {
        return tmp<scalarField>::New
        (
            max(min(solid_.alpha().primitiveField(), scalar(1)), scalar(0))
        );
    }
    return tmp<scalarField>::New(mesh_.nCells(), scalar(1));
}


Foam::tmp<Foam::volScalarField> Foam::heatTransfer::fillField() const
{
    tmp<volScalarField> tf
    (
        new volScalarField
        (
            IOobject
            (
                "thermalFill",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh_,
            dimensionedScalar(dimless, Zero),
            fvPatchFieldBase::zeroGradientType()
        )
    );
    tf.ref().primitiveFieldRef() = fill();
    tf.ref().correctBoundaryConditions();
    return tf;
}


Foam::tmp<Foam::scalarField> Foam::heatTransfer::capacity() const
{
    return tmp<scalarField>::New
    (
        solid_.rho().primitiveField()*c_*mesh_.V().field()
    );
}


Foam::scalar Foam::heatTransfer::energy() const
{
    return gSum(capacity()*T_.primitiveField());
}


Foam::scalar Foam::heatTransfer::stableDeltaT() const
{
    // Explicit limit per cell: C_P / sum_f k_f |S_f| delta_f
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceScalarField& magSf = mesh_.magSf();
    const surfaceScalarField& ndc = mesh_.nonOrthDeltaCoeffs();
    const tmp<volScalarField> tAVol(fillField());
    const volScalarField& aVol = tAVol();
    const scalarField& a = aVol.primitiveField();

    scalarField sumG(mesh_.nCells(), Zero);
    forAll(nei, facei)
    {
        const scalar g =
            k_*min(a[own[facei]], a[nei[facei]])*magSf[facei]*ndc[facei];
        sumG[own[facei]] += g;
        sumG[nei[facei]] += g;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (isA<emptyFvPatch>(p))
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const scalarField& pMag = magSf.boundaryField()[patchi];
        const scalarField& pNdc = ndc.boundaryField()[patchi];
        const scalarField aN
        (
            p.coupled()
          ? aVol.boundaryField()[patchi].patchNeighbourField()
          : tmp<scalarField>::New(fc.size(), scalar(1))
        );
        forAll(fc, i)
        {
            sumG[fc[i]] += k_*min(a[fc[i]], aN[i])*pMag[i]*pNdc[i];
        }
    }

    const tmp<scalarField> tC(capacity());
    const scalarField& C = tC();
    scalar dtMin = GREAT;
    forAll(C, celli)
    {
        if (sumG[celli] > VSMALL && C[celli] > VSMALL)
        {
            dtMin = min(dtMin, C[celli]/sumG[celli]);
        }
    }
    return returnReduce(dtMin, minOp<scalar>());
}


void Foam::heatTransfer::conduct(const scalar dt)
{
    T_.correctBoundaryConditions();

    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();
    const scalarField& w = mesh_.weights();
    const surfaceScalarField& ndc = mesh_.nonOrthDeltaCoeffs();
    const surfaceVectorField& corr = mesh_.nonOrthCorrectionVectors();
    const scalarField& V = mesh_.V();
    const tmp<volScalarField> tAVol(fillField());
    const volScalarField& aVol = tAVol();
    const scalarField& a = aVol.primitiveField();
    scalarField& T = T_.primitiveFieldRef();

    // Gauss gradient with linear face values
    volVectorField gradT
    (
        IOobject
        (
            "gradT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            IOobject::NO_REGISTER
        ),
        mesh_,
        dimensionedVector(dimTemperature/dimLength, Zero),
        fvPatchFieldBase::calculatedType()
    );
    vectorField& g = gradT.primitiveFieldRef();
    forAll(nei, facei)
    {
        // a face to an empty cell takes the material side's value
        const label P = own[facei];
        const label N = nei[facei];
        const scalar wf =
            a[P] > 0 && a[N] > 0 ? w[facei] : (a[P] > 0 ? 1 : 0);
        const scalar Tf = wf*T[P] + (1 - wf)*T[N];
        g[P] += Sf[facei]*Tf;
        g[N] -= Sf[facei]*Tf;
    }
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField& pSf = Sf.boundaryField()[patchi];
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        if (pT.coupled())
        {
            // a processor patch holds the neighbour cell values
            const scalarField TN(pT.patchNeighbourField());
            const scalarField aN
            (
                aVol.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pw = p.weights();
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar wf =
                    a[P] > 0 && aN[i] > 0 ? pw[i] : (a[P] > 0 ? 1 : 0);
                g[P] += pSf[i]*(wf*T[P] + (1 - wf)*TN[i]);
            }
        }
        else
        {
            forAll(fc, i)
            {
                g[fc[i]] += pSf[i]*pT[i];
            }
        }
    }
    g /= V;
    gradT.correctBoundaryConditions();

    // Heat flow into each cell over dt
    scalarField dE(mesh_.nCells(), Zero);
    const scalarField& ndcI = ndc.primitiveField();
    const vectorField& corrI = corr.primitiveField();
    forAll(nei, facei)
    {
        const label P = own[facei];
        const label N = nei[facei];
        const scalar kf = k_*min(a[P], a[N]);
        if (kf <= 0)
        {
            continue;
        }
        const vector gf = w[facei]*g[P] + (1 - w[facei])*g[N];
        const scalar snGrad =
            (T[N] - T[P])*ndcI[facei] + (corrI[facei] & gf);
        const scalar q = dt*kf*magSf[facei]*snGrad;    // into P
        dE[P] += q;
        dE[N] -= q;
    }

    scalar inflow = 0;
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (isA<emptyFvPatch>(p))
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        const scalarField& pMag = magSf.boundaryField()[patchi];
        const scalarField& pNdc = ndc.boundaryField()[patchi];

        if (pT.coupled())
        {
            const scalarField TN(pT.patchNeighbourField());
            const scalarField aN(aVol.boundaryField()[patchi].patchNeighbourField());
            const vectorField gN
            (
                gradT.boundaryField()[patchi].patchNeighbourField()
            );
            const scalarField& pw = p.weights();
            const vectorField& pCorr = corr.boundaryField()[patchi];
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar kf = k_*min(a[P], aN[i]);
                if (kf <= 0)
                {
                    continue;
                }
                const vector gf = pw[i]*g[P] + (1 - pw[i])*gN[i];
                const scalar snGrad =
                    (TN[i] - T[P])*pNdc[i] + (pCorr[i] & gf);
                dE[P] += dt*kf*pMag[i]*snGrad;
            }
        }
        else if (pT.fixesValue())
        {
            // Dirichlet: (T_b - T_P) delta + (n - delta d) . grad T_P
            // Cf - C (fvPatch::delta() is the normal part only for
            // non-coupled patches)
            const vectorField d(p.Cf() - p.Cn());
            const vectorField nf(p.nf());
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar kf = k_*a[P];
                const scalar dl = 1/max(nf[i] & d[i], 0.05*mag(d[i]));
                const scalar snGrad =
                    (pT[i] - T[P])*dl + ((nf[i] - dl*d[i]) & g[P]);
                const scalar q = dt*kf*pMag[i]*snGrad;
                dE[P] += q;
                inflow += q;
            }
        }
        else
        {
            // zeroGradient, fixedGradient (heat flux), symmetry, wedge
            const scalarField sn(pT.snGrad());
            forAll(fc, i)
            {
                const label P = fc[i];
                const scalar q = dt*k_*a[P]*pMag[i]*sn[i];
                dE[P] += q;
                inflow += q;
            }
        }
    }

    const tmp<scalarField> tC(capacity());
    const scalarField& C = tC();
    forAll(T, celli)
    {
        if (C[celli] > VSMALL)
        {
            T[celli] += dE[celli]/C[celli];
        }
    }
    T_.correctBoundaryConditions();

    heatIn_ += returnReduce(inflow, sumOp<scalar>());
}


void Foam::heatTransfer::beforeRemap()
{
    TOld_ = T_.primitiveField();
}


void Foam::heatTransfer::afterRemap()
{
    const surfaceScalarField& dm = solid_.remapMass();
    scalar in = 0;
    forAll(mesh_.boundary(), patchi)
    {
        const fvPatch& p = mesh_.boundary()[patchi];
        if (p.coupled() || isA<emptyFvPatch>(p))
        {
            continue;
        }
        const labelUList& fc = p.faceCells();
        const scalarField& pdm = dm.boundaryField()[patchi];
        const fvPatchScalarField& pT = T_.boundaryField()[patchi];
        forAll(fc, i)
        {
            if (pdm[i] > 0)
            {
                in -= c_*pdm[i]*TOld_[fc[i]];       // out with the cell's T
            }
            else if (pdm[i] < 0)
            {
                in -= c_*pdm[i]*pT[i];              // in at the inlet's T
            }
        }
    }
    carriedIn_ += returnReduce(in, sumOp<scalar>());

    // Empty cells hold no heat; write voidT there
    if (solid_.freeSurface())
    {
        const scalarField& rho = solid_.rho().primitiveField();
        scalarField& T = T_.primitiveFieldRef();
        forAll(T, celli)
        {
            if (rho[celli] <= 0)
            {
                T[celli] = voidT_;
            }
        }
    }
    T_.correctBoundaryConditions();
}


void Foam::heatTransfer::solve(const scalar dt)
{
    const scalar dtc = stableDeltaT();
    nSub_ = 1;
    if (dt > 0.5*dtc)
    {
        nSub_ = label(std::ceil(dt/(0.5*dtc)));
    }
    const scalar h = dt/nSub_;
    for (label n = 0; n < nSub_; ++n)
    {
        conduct(h);
    }
}


void Foam::heatTransfer::report() const
{
    const scalar E = energy();
    const label oldPrecision = Info().precision(16);
    Info<< "Thermal energy  = " << E << nl
        << "Heat in         = " << heatIn_ << nl
        << "Heat carried in = " << carriedIn_ << nl
        << "Energy balance  = " << E - energy0_ - heatIn_ - carriedIn_
        << " (energy - start - heat in - carried in)" << endl;
    Info().precision(oldPrecision);
}


// ************************************************************************* //
