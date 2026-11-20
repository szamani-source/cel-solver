/*---------------------------------------------------------------------------*\
Namespace
    cel::constitutive

Description
    Plain C++ constitutive kernel (WP 2.2). No OpenFOAM headers, so it builds
    and unit-tests with g++ alone, and the Lagrangian FE library can call it
    directly. The OpenFOAM wrapper (constitutiveModel) converts to and from
    symmTensor and adds run-time selection.

    Conventions (decision 0004, kept by decision 0006):
    - Voigt order xx, yy, zz, xy, yz, zx.
    - Strain increments use engineering shear (gamma = 2 eps); stress uses
      tensor components. Stress positive in tension. SI units.
    - One contract for every model: strain increment, old stress, old history
      and parameters in; new stress and new history out. Inputs are never
      changed, so a model may be called with sigNew aliasing nothing else.

    Reference implementation: prototype/constitutive/models.py.

SourceFiles
    celConstitutiveKernel.C
\*---------------------------------------------------------------------------*/

#ifndef celConstitutiveKernel_H
#define celConstitutiveKernel_H

#include <map>
#include <memory>
#include <string>

namespace cel
{
namespace constitutive
{

class FrictionLaw;
class Nonlocal;

//- Six Voigt components, order xx, yy, zz, xy, yz, zx
typedef double Voigt[6];

//- Per-point history, stored per cell or integration point and remapped
struct History
{
    //- Plastic strain (engineering shear), Voigt order
    double epsP[6] = {0, 0, 0, 0, 0, 0};

    //- Equivalent plastic strain, sqrt(2/3 e:e) accumulated
    double epsPEq = 0;

    //- True if the last increment returned to a cone apex
    bool apex = false;
};

//- Relative tolerance on the trial yield function below which a step is
//  elastic (YIELD_TOL in models.py)
constexpr double yieldTol = 1e-12;

//- Material parameters by name (SI units, angles in degrees)
typedef std::map<std::string, double> Parameters;

//- Non-numeric model options by name (e.g. Drucker-Prager "fit"), kept apart
//  from Parameters so parameters stay numbers only (WP 2.5a)
typedef std::map<std::string, std::string> Options;


//- Elastic constants from E and nu
struct Elastic
{
    double K;   // bulk modulus
    double G;   // shear modulus

    Elastic(double E, double nu);

    //- sigma + D deps (D maps engineering strain to stress)
    void addStress(const Voigt sigma, const Voigt deps, Voigt out) const;

    //- C dsigma (engineering strain from a stress increment)
    void strain(const Voigt dsigma, Voigt out) const;

    //- P-wave modulus K + 4G/3, for the stable time step
    double waveModulus() const
    {
        return K + 4.0*G/3.0;
    }
};


// Invariants (stress positive in tension)

double meanStress(const Voigt s);

void deviator(const Voigt s, Voigt out);

//- J2 = 1/2 s:s of the deviator
double J2(const Voigt s);



//- Abstract material model. Parameters are fixed at construction.
class Model
{
public:

    virtual ~Model() = default;

    //- Model name as used in dictionaries ("linearElastic", "vonMises",
    //  "druckerPrager", "mohrCoulomb")
    virtual const char* name() const = 0;

    //- One increment: (deps, sigOld, hOld) -> (sigNew, hNew)
    virtual void update
    (
        const Voigt deps,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const = 0;

    //- One increment over a time step dt [s]. Rate-dependent models (mu(I)
    //  friction, WP 3.6a) need dt; the others ignore it and call update.
    virtual void updateRate
    (
        const Voigt deps,
        double dt,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const
    {
        (void)dt;
        update(deps, sigOld, hOld, sigNew, hNew);
    }

    //- True if updateRate depends on dt
    virtual bool rateDependent() const
    {
        return false;
    }

    //- The mu(I) friction law of a rate-dependent model, or null
    virtual const FrictionLaw* frictionLaw() const
    {
        return nullptr;
    }

    //- The nonlocal fluidity law (WP N.1), or null
    virtual const Nonlocal* nonlocal() const
    {
        return nullptr;
    }

    //- One increment over dt at fluidity g [1/s] (nonlocal models); the
    //  others ignore g and call updateRate
    virtual void updateFluidity
    (
        const Voigt deps,
        double dt,
        double g,
        const Voigt sigOld,
        const History& hOld,
        Voigt sigNew,
        History& hNew
    ) const
    {
        (void)g;
        updateRate(deps, dt, sigOld, hOld, sigNew, hNew);
    }

    //- Elastic P-wave modulus K + 4G/3 (stable time step)
    virtual double waveModulus() const = 0;

    //- Elastic shear modulus G (shear impedance of the Lagrange step)
    virtual double shearModulus() const = 0;
};


//- Build a model by name. Throws std::invalid_argument for an unknown name,
//  a missing parameter or an unknown option value.
std::unique_ptr<Model> makeModel
(
    const std::string& name,
    const Parameters& p,
    const Options& opts = Options()
);

//- Read parameter key from p; throws std::invalid_argument if absent
double parameter(const Parameters& p, const std::string& key);

//- Read parameter key from p, or def if absent
double parameter(const Parameters& p, const std::string& key, double def);

//- Read option key from opts, or def if absent
std::string option
(
    const Options& opts,
    const std::string& key,
    const std::string& def
);

//- Finish a plastic step from its corrector (_finish in models.py):
//  hNew.epsP = hOld.epsP + C (trial - sigNew), hNew.epsPEq = hOld.epsPEq +
//  sqrt(2/3 e:e) of the deviatoric plastic strain increment, hNew.apex = apex.
void plasticCorrector
(
    const Elastic& el,
    const Voigt trial,
    const Voigt sigNew,
    const History& hOld,
    bool apex,
    History& hNew
);


} // End namespace constitutive
} // End namespace cel

#endif
