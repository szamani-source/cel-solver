/*---------------------------------------------------------------------------*\
Application
    celConstitutiveWrapperCheck

Description
    WP 2.2: checks the OpenFOAM wrapper of the constitutive kernel (needs
    OpenFOAM; run by tests/test_constitutive_cpp.sh).

    - Run-time selection builds both models from a dictionary.
    - symmTensor <-> Voigt conversion: each tensor shear strain component eps_ij
      gives sigma_ij = 2 G eps_ij in the same component and nothing elsewhere;
      a normal strain gives lambda and lambda + 2G.
    - A large shear in xz takes vonMises to the yield surface: q = sigma_y.
\*---------------------------------------------------------------------------*/

#include "constitutiveModel.H"
#include "IOstreams.H"

using namespace Foam;

int main()
{
    const scalar E = 20e6, nu = 0.3, sigmaY = 200e3;
    const scalar G = E/(2*(1 + nu));
    const scalar lambda = E*nu/((1 + nu)*(1 - 2*nu));

    label nFail = 0;
    auto check = [&](const char* what, bool ok)
    {
        Info<< (ok ? "PASS  " : "FAIL  ") << what << nl;
        if (!ok) ++nFail;
    };

    dictionary elDict;
    elDict.add("type", word("linearElastic"));
    elDict.add("E", E);
    elDict.add("nu", nu);

    dictionary vmDict(elDict);
    vmDict.set("type", word("vonMises"));
    vmDict.add("sigma_y", sigmaY);

    autoPtr<constitutiveModel> el(constitutiveModel::New(elDict));
    autoPtr<constitutiveModel> vm(constitutiveModel::New(vmDict));

    check("selected linearElastic", word(el->kernel().name()) == "linearElastic");
    check("selected vonMises", word(vm->kernel().name()) == "vonMises");
    check
    (
        "wave modulus K + 4G/3",
        mag(el->waveModulus() - (lambda + 2*G)) < 1e-9*(lambda + 2*G)
    );

    const scalar e = 1e-4;
    const direction comps[6] =
    {
        symmTensor::XX, symmTensor::XY, symmTensor::XZ,
        symmTensor::YY, symmTensor::YZ, symmTensor::ZZ
    };
    for (const direction c : comps)
    {
        symmTensor deps(Zero);
        deps.component(c) = e;

        symmTensor sigNew;
        constitutiveModel::history h, hNew;
        el->update(deps, symmTensor(Zero), h, sigNew, hNew);

        symmTensor expected(Zero);
        const bool normal =
            (c == symmTensor::XX || c == symmTensor::YY || c == symmTensor::ZZ);
        if (normal)
        {
            expected.xx() = expected.yy() = expected.zz() = lambda*e;
        }
        expected.component(c) += 2*G*e;

        const scalar err = mag(sigNew - expected)/mag(expected);
        Info<< "    component " << symmTensor::componentNames[c]
            << "  sigma = " << sigNew << "  rel error " << err << nl;
        check("component order and shear factor", err < 1e-12);
    }

    // Plastic: large xz shear from an isotropic compressive state
    {
        symmTensor deps(Zero);
        deps.xz() = 0.01;
        const symmTensor sig0(-1e5*symmTensor(I));
        symmTensor sigNew;
        constitutiveModel::history h, hNew;
        vm->update(deps, sig0, h, sigNew, hNew);

        const scalar q = sqrt(1.5)*mag(dev(sigNew));
        Info<< "    vonMises xz shear: sigma = " << sigNew << "  q = " << q
            << "  epsPEq = " << hNew.epsPEq << nl;
        check("vonMises returns to q = sigma_y", mag(q - sigmaY) < 1e-9*sigmaY);
        check("vonMises keeps the mean stress", mag(tr(sigNew)/3 + 1e5) < 1e-6);
        check("vonMises plastic strain grows", hNew.epsPEq > 0);
        check
        (
            "plastic shear lands in xz (Voigt slot zx)",
            hNew.epsP[5] > 0 && hNew.epsP[3] == 0 && hNew.epsP[4] == 0
        );
    }

    Info<< nl << (nFail ? "FAILED" : "All wrapper checks passed") << endl;
    return nFail ? 1 : 0;
}
