#include "constitutiveModel.H"
#include "error.H"

#include <exception>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(constitutiveModel, 0);
    defineRunTimeSelectionTable(constitutiveModel, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constitutiveModel::constitutiveModel
(
    const dictionary& dict,
    const word& modelType,
    const wordList& required,
    const wordList& optional,
    const wordList& wordOptions
)
{
    cel::constitutive::Parameters p;

    for (const word& key : required)
    {
        p[key] = dict.get<scalar>(key);
    }
    for (const word& key : optional)
    {
        scalar value = 0;
        if (dict.readIfPresent(key, value))
        {
            p[key] = value;
        }
    }

    cel::constitutive::Options opts;
    for (const word& key : wordOptions)
    {
        word value;
        if (dict.readIfPresent(key, value))
        {
            opts[key] = value;
        }
    }

    try
    {
        kernel_ = cel::constitutive::makeModel(modelType, p, opts);
    }
    catch (const std::exception& e)
    {
        FatalIOErrorInFunction(dict)
            << e.what() << exit(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::constitutiveModel::toVoigt
(
    const symmTensor& t,
    double v[6],
    const scalar shearFactor
)
{
    v[0] = t.xx();
    v[1] = t.yy();
    v[2] = t.zz();
    v[3] = shearFactor*t.xy();
    v[4] = shearFactor*t.yz();
    v[5] = shearFactor*t.xz();
}


Foam::symmTensor Foam::constitutiveModel::fromVoigt(const double v[6])
{
    // symmTensor(xx, xy, xz, yy, yz, zz)
    return symmTensor(v[0], v[3], v[5], v[1], v[4], v[2]);
}


void Foam::constitutiveModel::update
(
    const symmTensor& deps,
    const symmTensor& sigOld,
    const history& hOld,
    symmTensor& sigNew,
    history& hNew
) const
{
    double d[6], s[6], sn[6];
    toVoigt(deps, d, 2);    // engineering shear for the kernel
    toVoigt(sigOld, s);
    kernel_->update(d, s, hOld, sn, hNew);
    sigNew = fromVoigt(sn);
}


// ************************************************************************* //
