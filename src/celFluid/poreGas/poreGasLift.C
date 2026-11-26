/*---------------------------------------------------------------------------*\
    celFoam: explicit Coupled Eulerian-Lagrangian solver on OpenFOAM
    Licence: GPL-3.0-or-later

    Range-of-validity check (decision 0135; gasProperties liftCheck, on by
    default).

    The gas coupling holds the grains at their packed (stressed) density:
    a bed that stretches opens void inside its cells, it does not loosen
    (decision 0021 rule 4), and the gas sees the permeability of the
    packed part (decision 0125 item 4). That is right while the grains'
    own stress carries their weight, less the gas's share. Once the gas
    carries the whole weight the bed lifts or expands, which a two-fluid
    model represents and this one does not (decision 0135). This check
    measures how much of the grain mass is in that state.

    Each cell with grains (not void) gets the force the solid momentum
    uses, F = -alpha sum(p_f S_f) over its open faces (eulerianSolid,
    WP F.3; the face of an immersed wall is left out), and its weight
    m |g| with m = rho V_open. The cell is carried when the upward part of
    F, F.(-g)/|g|, is at least liftRatio m |g|. The carried share is the
    carried mass over the grain mass. The step's share is kept
    (liftShare()), printed at write times, and a warning is given the
    first time it has stayed above liftWarn for liftTime (a lifting bed
    stays carried; the opening of an outlet gives spikes of a few steps).
\*---------------------------------------------------------------------------*/

#include "poreGas.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::poreGas::checkLift()
{
    liftShare_ = 0;
    const vector g(solid_.g());
    const scalar gMag = mag(g);
    if (!liftCheck_ || !coupled_ || gMag <= 0)
    {
        return;
    }
    const vector up(-g/gMag);

    const fvMesh& mesh = solid_.mesh();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();
    const scalarField& rho = solid_.rho().primitiveField();
    const scalarField& alpha = solid_.alpha().primitiveField();
    const scalarField& eM = eM_.primitiveField();
    const scalarField& vd = void_.primitiveField();
    const tmp<scalarField> tVol(solid_.openVolume());
    const scalarField& V = tVol();

    tmp<surfaceVectorField> tSf(mesh.Sf());
    if (const surfaceScalarField* faceOpen = solid_.faceOpenPtr())
    {
        tSf = mesh.Sf()*(*faceOpen);
    }
    const surfaceVectorField& Sf = tSf();

    // Upward gas force per cell
    scalarField Fup(mesh.nCells(), Zero);
    const scalarField& pfi = pf_.primitiveField();
    forAll(own, facei)
    {
        const scalar f = pfi[facei]*(Sf[facei] & up);
        Fup[own[facei]] -= f;
        Fup[nei[facei]] += f;
    }
    forAll(mesh.boundary(), patchi)
    {
        const labelUList& fc = mesh.boundary()[patchi].faceCells();
        const vectorField& pS = Sf.boundaryField()[patchi];
        const scalarField& pp = pf_.boundaryField()[patchi];
        forAll(fc, i)
        {
            Fup[fc[i]] -= pp[i]*(pS[i] & up);
        }
    }

    const bool fs = solid_.freeSurface();
    scalar mAll = 0, mCarried = 0;
    forAll(Fup, celli)
    {
        // Cells with grains: a filled part with a gas fraction (eM > 0),
        // not a void cell
        if (eM[celli] <= 0 || (fs && vd[celli] > 0.5))
        {
            continue;
        }
        const scalar m = rho[celli]*V[celli];
        if (m <= 0)
        {
            continue;
        }
        mAll += m;
        if ((fs ? alpha[celli] : 1)*Fup[celli] >= liftRatio_*m*gMag)
        {
            mCarried += m;
        }
    }
    reduce(mAll, sumOp<scalar>());
    reduce(mCarried, sumOp<scalar>());
    liftShare_ = mAll > 0 ? mCarried/mAll : 0;

    if (mesh.time().writeTime())
    {
        Info<< "Gas lift: the gas carries the weight of " << 100*liftShare_
            << " % of the grain mass (decision 0135)" << nl;
    }
    const scalar t = mesh.time().value();
    if (liftShare_ <= liftWarn_)
    {
        liftSince_ = -1;
    }
    else if (liftSince_ < 0)
    {
        liftSince_ = t;
    }
    if (!liftWarned_ && liftSince_ >= 0 && t - liftSince_ >= liftTime_)
    {
        liftWarned_ = true;
        WarningInFunction
            << "since t = " << liftSince_ << " s the gas has carried the "
            << "weight of more than " << 100*liftWarn_ << " % of the grain "
            << "mass" << nl << "    (" << 100*liftShare_ << " % now)." << nl
            << "    The bed is lifting or expanding: from here on the "
            << "results are outside the" << nl
            << "    range of the gas coupling (grains held at their packed "
            << "density, decision 0135)." << nl
            << "    A fluidised or expanded bed needs a two-fluid model. "
            << "(liftCheck no in" << nl
            << "    constant/gasProperties turns this check off.)" << endl;
    }
}


// ************************************************************************* //
