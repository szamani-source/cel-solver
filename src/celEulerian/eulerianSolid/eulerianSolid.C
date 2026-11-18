/*---------------------------------------------------------------------------*\
    eulerianSolid: see eulerianSolid.H
\*---------------------------------------------------------------------------*/

#include "eulerianSolid.H"
#include "unitConversion.H"
#include "surfaceFields.H"
#include "zeroGradientFvPatchFields.H"
#include "uniformDimensionedFields.H"
#include "mathematicalConstants.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{

//- IOobject for a field read if present and written at output times
static IOobject fieldIO(const word& name, const fvMesh& mesh)
{
    return IOobject
    (
        name,
        mesh.time().timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    );
}

//- Gravity from constant/g if present, else zero
static vector readGravity(const fvMesh& mesh)
{
    IOobject io
    (
        "g",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        IOobject::NO_REGISTER
    );
    if (!io.typeHeaderOk<uniformDimensionedVectorField>(true))
    {
        return Zero;
    }
    return uniformDimensionedVectorField(io).value();
}


//- Wall patches and tan(friction) from materialProperties walls
static HashTable<scalar> readWalls(const fvMesh& mesh, const dictionary& dict)
{
    HashTable<scalar> walls;
    const dictionary* wd = dict.findDict("walls");
    if (!wd)
    {
        return walls;
    }
    for (const entry& e : *wd)
    {
        if (mesh.boundaryMesh().findPatchID(e.keyword()) < 0)
        {
            FatalIOErrorInFunction(*wd)
                << "walls: no patch " << e.keyword() << exit(FatalIOError);
        }
        const scalar deg = e.dict().get<scalar>("friction");
        if (deg < 0 || deg >= 90)
        {
            FatalIOErrorInFunction(*wd)
                << "walls: friction of " << e.keyword()
                << " must be in [0, 90) degrees, got " << deg
                << exit(FatalIOError);
        }
        walls.set(e.keyword(), Foam::tan(deg*constant::mathematical::pi/180));
    }
    return walls;
}


//- Outlet patches from materialProperties outlets (WP 3.3): opening time,
//  and their closed-wall friction added to walls
static HashTable<scalar> readOutlets
(
    const fvMesh& mesh,
    const dictionary& dict,
    HashTable<scalar>& walls
)
{
    HashTable<scalar> outlets;
    const dictionary* od = dict.findDict("outlets");
    if (!od)
    {
        return outlets;
    }
    for (const entry& e : *od)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(e.keyword());
        if (patchi < 0)
        {
            FatalIOErrorInFunction(*od)
                << "outlets: no patch " << e.keyword() << exit(FatalIOError);
        }
        if (walls.found(e.keyword()))
        {
            FatalIOErrorInFunction(*od)
                << "outlets: " << e.keyword() << " is also in walls"
                << exit(FatalIOError);
        }
        const scalar deg = e.dict().get<scalar>("friction");
        if (deg < 0 || deg >= 90)
        {
            FatalIOErrorInFunction(*od)
                << "outlets: friction of " << e.keyword()
                << " must be in [0, 90) degrees, got " << deg
                << exit(FatalIOError);
        }
        walls.set(e.keyword(), Foam::tan(deg*constant::mathematical::pi/180));
        outlets.set
        (
            e.keyword(),
            e.dict().getOrDefault<scalar>("openTime", VGREAT)
        );
    }
    return outlets;
}

} // End namespace Foam


Foam::eulerianSolid::eulerianSolid(const fvMesh& mesh)
:
    mesh_(mesh),
    materialProperties_
    (
        IOobject
        (
            "materialProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    ),
    material_(constitutiveModel::New(materialProperties_)),
    maxCo_(mesh.time().controlDict().getOrDefault<scalar>("maxCo", 0.4)),
    solveMomentum_
    (
        mesh.time().controlDict().getOrDefault<Switch>("solveMomentum", true)
    ),
    secondOrder_
    (
        mesh.time().controlDict().getOrDefault<Switch>("secondOrder", true)
    ),
    rho_
    (
        fieldIO("rho", mesh),
        mesh,
        dimensionedScalar("rho", dimDensity, materialProperties_),
        fvPatchFieldBase::zeroGradientType()
    ),
    U_
    (
        fieldIO("U", mesh),
        mesh,
        dimensionedVector(dimVelocity, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    sigma_
    (
        fieldIO("sigma", mesh),
        mesh,
        dimensionedSymmTensor(dimPressure, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    epsP_
    (
        fieldIO("epsP", mesh),
        mesh,
        dimensionedSymmTensor(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    epsPEq_
    (
        fieldIO("epsPEq", mesh),
        mesh,
        dimensionedScalar(dimless, Zero),
        fvPatchFieldBase::zeroGradientType()
    ),
    phiL_
    (
        IOobject
        (
            "phiL",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimVolume/dimTime, Zero)
    ),
    remapLimiter_(readRemapLimiters(materialProperties_)),
    g_(readGravity(mesh)),
    gravityRamp_
    (
        mesh.time().controlDict().getOrDefault<scalar>("gravityRamp", 0)
    ),
    damping_(mesh.time().controlDict().getOrDefault<scalar>("damping", 0)),
    wallFriction_(readWalls(mesh, materialProperties_)),
    vStar_
    (
        IOobject
        (
            "vStar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimVelocity, Zero)
    ),
    Fstar_
    (
        IOobject
        (
            "Fstar",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedVector(dimForce, Zero)
    ),
    freeSurface_(materialProperties_.getOrDefault<Switch>("freeSurface", false)),
    alpha_
    (
        fieldIO("alpha", mesh),
        mesh,
        dimensionedScalar(dimless, 1),
        fvPatchFieldBase::zeroGradientType()
    ),
    rho0_(materialProperties_.get<scalar>("rho")),
    noTension_
    (
        (materialProperties_.found("c") && materialProperties_.get<scalar>("c") == 0)
     || (
            materialProperties_.found("sigmaT")
         && materialProperties_.get<scalar>("sigmaT") == 0
        )
    ),
    tensileStrength_(GREAT),
    shearCap_(materialProperties_.getOrDefault<Switch>("shearCap", false)),
    normalCap_(materialProperties_.getOrDefault<Switch>("normalCap", false)),
    capSinPhi_
    (
        (shearCap_ || normalCap_)
      ? Foam::sin(degToRad(materialProperties_.get<scalar>("phi")))
      : 0
    ),
    capCCos_
    (
        (shearCap_ || normalCap_)
      ? materialProperties_.getOrDefault<scalar>("c", 0)
       *Foam::cos(degToRad(materialProperties_.get<scalar>("phi")))
      : 0
    ),
    flowImpedance_
    (
        materialProperties_.getOrDefault<Switch>("flowImpedance", false)
    ),
    fullBlend_
    (
        materialProperties_.getOrDefault<Switch>("fullBlend", false)
    ),
    blendVoid_
    (
        materialProperties_.getOrDefault<FixedList<scalar, 2>>
        (
            "fullBlendVoid", FixedList<scalar, 2>({1e-4, 1e-2})
        )
    ),
    flowLaw_(Zero),
    spreadLayers_
    (
        materialProperties_.getOrDefault<Switch>("spreadLayers", false)
    ),
    nCapped_(0),
    nNormalCapped_(0),
    outletOpenTime_(readOutlets(mesh, materialProperties_, wallFriction_)),
    discharged_(0),
    inletOfPatch_(mesh.boundary().size(), -1),
    charged_(0),
    immersed_(false),
    nUnmixed_(0),
    movingWalls_(false),
    lostMass_(0),
    nHandOver_(0)
{
    if (outletOpenTime_.size() && !freeSurface_)
    {
        FatalIOErrorInFunction(materialProperties_)
            << "outlets need freeSurface yes" << exit(FatalIOError);
    }
    forAllConstIters(outletOpenTime_, iter)
    {
        const label patchi = mesh.boundaryMesh().findPatchID(iter.key());
        const fvPatch& p = mesh.boundary()[patchi];
        if
        (
            p.coupled() || isA<emptyFvPatch>(p)
         || U_.boundaryField()[patchi].fixesValue()
        )
        {
            FatalErrorInFunction
                << "Outlet " << iter.key() << " must be a plain patch "
                << "without a fixed velocity" << exit(FatalError);
        }
        Info<< "Outlet " << iter.key() << ": closed (friction coefficient "
            << wallFriction_[iter.key()] << "), opens at t = "
            << iter.val() << " s" << nl;
    }
    readInlets();
    if (spreadLayers_ && !freeSurface_)
    {
        FatalIOErrorInFunction(materialProperties_)
            << "spreadLayers needs freeSurface yes" << exit(FatalIOError);
    }
    if (shearCap_)
    {
        Info<< "Shear cap at Mohr-Coulomb yield (rule 8): sin(phi) "
            << capSinPhi_ << ", c cos(phi) " << capCCos_ << " Pa" << nl;
    }
    if (normalCap_)
    {
        Info<< "Normal jump cap at Mohr-Coulomb yield (rule 12)" << nl;
    }
    if (fullBlend_)
    {
        if (!(blendVoid_[0] > 0 && blendVoid_[1] > blendVoid_[0]))
        {
            FatalIOErrorInFunction(materialProperties_)
                << "fullBlendVoid needs 0 < v1 < v2" << exit(FatalIOError);
        }
        Info<< "Full / partly filled rules blended over void "
            << blendVoid_[0] << " to " << blendVoid_[1]
            << " (no nearly-full snap)" << nl;
    }
    if (flowImpedance_)
    {
        const dictionary& d = materialProperties_;
        if (d.get<word>("type") != "mohrCoulomb" || !d.found("mu2"))
        {
            FatalIOErrorInFunction(d)
                << "flowImpedance needs type mohrCoulomb with mu(I) (mu2)"
                << exit(FatalIOError);
        }
        flowLaw_[0] = Foam::tan(degToRad(d.get<scalar>("phi")));
        flowLaw_[1] = d.get<scalar>("mu2");
        flowLaw_[2] = d.get<scalar>("I0");
        flowLaw_[3] = d.get<scalar>("d");
        flowLaw_[4] = d.get<scalar>("rhoS");
        if (flowLaw_[1] <= flowLaw_[0])
        {
            FatalIOErrorInFunction(d)
                << "flowImpedance needs mu2 > tan(phi)" << exit(FatalIOError);
        }
        if (capSinPhi_ == 0)
        {
            capSinPhi_ = Foam::sin(degToRad(d.get<scalar>("phi")));
            capCCos_ =
                d.getOrDefault<scalar>("c", 0)
               *Foam::cos(degToRad(d.get<scalar>("phi")));
        }
        Info<< "Flow impedance (WP 3.6c): mu(I) tangent viscosity at faces "
            << "between flowing cells, plastic Hancock predictor" << nl;
    }
    if (spreadLayers_)
    {
        Info<< "Top layers spread sideways (rule 9)" << nl;
    }

    if (maxCo_ <= 0)
    {
        FatalIOErrorInFunction(mesh.time().controlDict())
            << "maxCo must be positive, got " << maxCo_
            << exit(FatalIOError);
    }

    // Tensile strength (WP 3.5, decision 0054)
    if (noTension_)
    {
        tensileStrength_ = 0;
    }
    else if (materialProperties_.found("sigmaT"))
    {
        tensileStrength_ = materialProperties_.get<scalar>("sigmaT");
    }
    else if
    (
        materialProperties_.found("c") && materialProperties_.found("phi")
     && materialProperties_.get<scalar>("phi") > 0
    )
    {
        tensileStrength_ =
            materialProperties_.get<scalar>("c")
           /Foam::tan(degToRad(materialProperties_.get<scalar>("phi")));
    }

    if (freeSurface_)
    {
        setCartesianGeometry();
        lean_.setSize(mesh.nCells(), Zero);
        alphaL_ = alpha_.primitiveField();
        smooth_.setSize(mesh.nCells(), Zero);
        Info<< "Free surface: material plus void, alpha in ["
            << gMin(alpha_.primitiveField()) << ", "
            << gMax(alpha_.primitiveField()) << "], "
            << (
                   noTension_ ? "no tensile strength"
                 : tensileStrength_ < GREAT ? "tensile strength "
                 : "no tensile limit"
               );
        if (!noTension_ && tensileStrength_ < GREAT)
        {
            Info<< tensileStrength_ << " Pa";
        }
        Info<< nl;
    }
    else
    {
        alpha_.writeOpt(IOobject::NO_WRITE);
    }

    const scalar minRho = gMin(rho_.primitiveField());
    if (minRho < 0 || (minRho == 0 && !freeSurface_))
    {
        FatalErrorInFunction
            << "rho must be positive in every cell, smallest is " << minRho
            << exit(FatalError);
    }

    Info<< "Eulerian solid: rho in [" << minRho << ", "
        << gMax(rho_.primitiveField()) << "] kg/m^3, maxCo = " << maxCo_
        << ", solveMomentum " << solveMomentum_
        << ", secondOrder " << secondOrder_ << nl;
    if (mag(g_) > 0)
    {
        Info<< "Gravity " << g_ << " m/s^2, ramp " << gravityRamp_ << " s"
            << nl;
    }
    if (damping_ > 0)
    {
        Info<< "Velocity damping " << damping_ << " 1/s" << nl;
    }
    forAllConstIters(wallFriction_, iter)
    {
        Info<< "Wall " << iter.key() << ": friction coefficient "
            << iter.val() << nl;
    }
    Info<< endl;
}


void Foam::eulerianSolid::openOutlets()
{
    forAllConstIters(outletOpenTime_, iter)
    {
        if (wallFriction_.erase(iter.key()))
        {
            Info<< "Outlet " << iter.key() << " opened at t = "
                << mesh_.time().value() << " s; damping off" << endl;
            damping_ = 0;
        }
    }
}


Foam::scalar Foam::eulerianSolid::discharged() const
{
    return returnReduce(discharged_, sumOp<scalar>());
}


Foam::scalar Foam::eulerianSolid::charged() const
{
    return returnReduce(charged_, sumOp<scalar>());
}


void Foam::eulerianSolid::readInlets()
{
    const dictionary* id = materialProperties_.findDict("inlets");
    if (!id)
    {
        return;
    }
    if (!freeSurface_)
    {
        FatalIOErrorInFunction(materialProperties_)
            << "inlets need freeSurface yes" << exit(FatalIOError);
    }
    for (const entry& e : *id)
    {
        const word& name = e.keyword();
        if (!e.isDict())
        {
            FatalIOErrorInFunction(*id)
                << "inlets: " << name << " must be a dictionary"
                << exit(FatalIOError);
        }
        const dictionary& d = e.dict();
        const label patchi = mesh_.boundaryMesh().findPatchID(name);
        if (patchi < 0)
        {
            FatalIOErrorInFunction(*id)
                << "inlets: no patch " << name << exit(FatalIOError);
        }
        const fvPatch& p = mesh_.boundary()[patchi];
        if
        (
            p.coupled() || isA<emptyFvPatch>(p)
         || U_.boundaryField()[patchi].fixesValue()
        )
        {
            FatalErrorInFunction
                << "Inlet " << name << " must be a plain patch without a "
                << "fixed velocity" << exit(FatalError);
        }
        if (wallFriction_.found(name) || outletOpenTime_.found(name))
        {
            FatalIOErrorInFunction(*id)
                << "inlets: " << name << " is also a wall or an outlet"
                << exit(FatalIOError);
        }
        const scalar deg = d.getOrDefault<scalar>("friction", 0);
        if (deg < 0 || deg >= 90)
        {
            FatalIOErrorInFunction(d)
                << "inlets: friction of " << name
                << " must be in [0, 90) degrees, got " << deg
                << exit(FatalIOError);
        }
        const scalar rhoIn = d.getOrDefault<scalar>("density", rho0_);
        if (rhoIn <= 0)
        {
            FatalIOErrorInFunction(d)
                << "inlets: density of " << name << " must be positive"
                << exit(FatalIOError);
        }
        const label k = inletName_.size();
        inletName_.append(name);
        inletPatch_.append(patchi);
        inletSpeed_.append(Function1<scalar>::New("speed", d, &mesh_));
        inletStart_.append(d.getOrDefault<scalar>("startTime", 0));
        inletEnd_.append(d.getOrDefault<scalar>("endTime", VGREAT));
        inletFriction_.append(Foam::tan(degToRad(deg)));
        inletRho_.append(rhoIn);
        inletP_.append(d.getOrDefault<scalar>("pressure", 0));
        inletOpen_.append(false);
        inletNow_.append(0);
        inletOfPatch_[patchi] = k;

        // Closed until its start time
        wallFriction_.set(name, inletFriction_[k]);
        Info<< "Inlet " << name << ": feeds material at " << rhoIn
            << " kg/m^3, pressure " << inletP_[k] << " Pa, from t = "
            << inletStart_[k] << " to " << inletEnd_[k]
            << " s; closed wall otherwise (friction coefficient "
            << inletFriction_[k] << ")" << nl;
    }
}


void Foam::eulerianSolid::updateInlets(const scalar t0, const scalar dt)
{
    // An inlet is open for the steps that start in [startTime, endTime)
    // (WP 3.7a, decision 0057); its speed is taken at mid-step
    forAll(inletName_, k)
    {
        const word& name = inletName_[k];
        const bool open =
            t0 >= inletStart_[k] - 1e-9*dt && t0 < inletEnd_[k] - 1e-9*dt;
        if (open)
        {
            const scalar s = inletSpeed_[k].value(t0 + 0.5*dt);
            if (s < 0)
            {
                FatalErrorInFunction
                    << "Inlet " << name << ": negative speed " << s
                    << " at t = " << t0 + 0.5*dt << " s (material only flows "
                    << "in; use an outlet to let it out)" << exit(FatalError);
            }
            inletNow_[k] = s;
            if (wallFriction_.erase(name))
            {
                Info<< "Inlet " << name << " opened at t = " << t0 << " s"
                    << endl;
            }
        }
        else
        {
            inletNow_[k] = 0;
            if (inletOpen_[k])
            {
                Info<< "Inlet " << name << " closed at t = " << t0 << " s"
                    << endl;
            }
            wallFriction_.set(name, inletFriction_[k]);
        }
        inletOpen_[k] = open;
    }
}


Foam::vector Foam::eulerianSolid::gravityAt(const scalar t) const
{
    if (gravityRamp_ <= 0 || t >= gravityRamp_)
    {
        return g_;
    }
    return
        0.5*(1 - Foam::cos(constant::mathematical::pi*max(t, 0)/gravityRamp_))
       *g_;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::eulerianSolid::setCartesianGeometry()
{
    // h_a = V / A_a with A_a half the summed |S_a| of the cell's faces;
    // every face must be normal to an axis (decision 0031)
    const vectorField& Sf = mesh_.Sf().primitiveField();
    const scalarField& magSf = mesh_.magSf().primitiveField();
    const labelUList& own = mesh_.owner();
    const labelUList& nei = mesh_.neighbour();
    vectorField A(mesh_.nCells(), Zero);

    auto check = [&](const vector& S, const scalar mS)
    {
        if (cmptMax(cmptMag(S)) < (1 - 1e-6)*mS)
        {
            FatalErrorInFunction
                << "freeSurface needs a Cartesian hex mesh: face normal "
                << S/mS << " is not along an axis" << exit(FatalError);
        }
    };
    forAll(own, facei)
    {
        check(Sf[facei], magSf[facei]);
        A[own[facei]] += 0.5*cmptMag(Sf[facei]);
        A[nei[facei]] += 0.5*cmptMag(Sf[facei]);
    }
    // Poly patches: empty patches (2D front and back) count here too
    forAll(mesh_.boundaryMesh(), patchi)
    {
        const polyPatch& p = mesh_.boundaryMesh()[patchi];
        const labelUList& fc = p.faceCells();
        const vectorField pS(p.faceAreas());
        forAll(fc, i)
        {
            check(pS[i], mag(pS[i]));
            A[fc[i]] += 0.5*cmptMag(pS[i]);
        }
    }
    const scalarField& V = mesh_.V().field();
    hCell_.setSize(mesh_.nCells());
    forAll(hCell_, celli)
    {
        for (direction a = 0; a < vector::nComponents; ++a)
        {
            hCell_[celli][a] = V[celli]/A[celli][a];
        }
    }
}


Foam::scalar Foam::eulerianSolid::stableDeltaT() const
{
    if (freeSurface_)
    {
        // Wave speed of the material (density rho/alpha, bounded as in the
        // Lagrange step); empty cells take the reference density
        const scalarField& r = rho_.primitiveField();
        const scalarField& a = alpha_.primitiveField();
        const scalarField delta(cbrt(mesh_.V().field()));
        scalar dt = GREAT;
        forAll(r, celli)
        {
            const scalar rm =
                a[celli] > 1e-12
              ? min(max(r[celli]/a[celli], 0.5*rho0_), 2*rho0_)
              : rho0_;
            dt = min(dt, delta[celli]/sqrt(material_->waveModulus()/rm));
        }
        return maxCo_*returnReduce(dt, minOp<scalar>());
    }

    // Cell size cbrt(V) over P-wave speed sqrt(M/rho), smallest over the mesh
    const scalarField delta(cbrt(mesh_.V().field()));
    const scalarField c(sqrt(material_->waveModulus()/rho_.primitiveField()));

    return maxCo_*gMin(delta/c);
}


// remap(dt): eulerianSolidRemap.C (WP 2.4)


Foam::scalar Foam::eulerianSolid::totalMass() const
{
    return gSum(rho_.primitiveField()*openVolume());
}


Foam::vector Foam::eulerianSolid::totalMomentum() const
{
    return gSum
    (
        rho_.primitiveField()*openVolume()*U_.primitiveField()
    );
}


// ************************************************************************* //
