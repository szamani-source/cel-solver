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
    - WP 2.5a: druckerPrager is selected with its word option fit; a large
      xz shear returns to the cone, a large volumetric extension to the apex.
    - WP 2.5b: run-time selection of mohrCoulomb; an isochoric axial compression
      from an isotropic state returns to the triaxial compression edge
      (sxx = syy) on the yield surface f = 0.
    - WP 2.5c: mohrCoulomb reads the optional sigmaT; a hydrostatic stretch
      ends at the tension apex sigma = sigmaT I.
\*---------------------------------------------------------------------------*/

#include "constitutiveModel.H"
#include "celDruckerPrager.H"
#include "celMohrCoulomb.H"
#include "IOstreams.H"
#include "unitConversion.H"

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

        const scalar q = Foam::sqrt(scalar(1.5))*mag(dev(sigNew));
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

    // Drucker-Prager (WP 2.5a): phi 30, psi 5, c 10 kPa, plane-strain fit
    {
        dictionary dpDict(elDict);
        dpDict.set("type", word("druckerPrager"));
        dpDict.add("phi", 30.0);
        dpDict.add("c", 10e3);
        dpDict.add("psi", 5.0);
        dpDict.add("fit", word("plane_strain"));

        autoPtr<constitutiveModel> dp(constitutiveModel::New(dpDict));
        check("selected druckerPrager", word(dp->kernel().name()) == "druckerPrager");

        const auto* k =
            dynamic_cast<const cel::constitutive::DruckerPrager*>(&dp->kernel());
        check
        (
            "fit read as plane_strain",
            k && k->fit() == cel::constitutive::DruckerPrager::Fit::planeStrain
        );

        dictionary defDict(dpDict);
        defDict.remove("fit");
        autoPtr<constitutiveModel> dpDef(constitutiveModel::New(defDict));
        const auto* kDef =
            dynamic_cast<const cel::constitutive::DruckerPrager*>(&dpDef->kernel());
        check
        (
            "fit defaults to compression",
            kDef && kDef->fit() == cel::constitutive::DruckerPrager::Fit::compression
        );

        if (k)
        {
            const symmTensor sig0(-1e5*symmTensor(I));
            constitutiveModel::history h, hNew;
            symmTensor sigNew;

            // Large xz shear: cone return
            symmTensor deps(Zero);
            deps.xz() = 0.01;
            dp->update(deps, sig0, h, sigNew, hNew);
            double v[6];
            constitutiveModel::toVoigt(sigNew, v);
            Info<< "    druckerPrager xz shear: sigma = " << sigNew
                << "  f = " << k->yield(v) << nl;
            check("druckerPrager shear lands on the cone", mag(k->yield(v)) < 1e-9*1e5);
            check("druckerPrager shear is not an apex return", !hNew.apex);
            // psi > 0: plastic dilation; at fixed total volume that raises
            // the confinement (mean stress more compressive)
            check
            (
                "druckerPrager dilates (plastic volume grows)",
                hNew.epsP[0] + hNew.epsP[1] + hNew.epsP[2] > 0
            );
            check("druckerPrager confinement grows", tr(sigNew)/3 < -1e5);

            // Large volumetric extension: apex return
            deps = symmTensor(0.01*symmTensor(I));
            dp->update(deps, sig0, h, sigNew, hNew);
            const scalar pApex = k->xi()*10e3/k->eta();
            Info<< "    druckerPrager extension: sigma = " << sigNew
                << "  apex p = " << pApex << nl;
            check("druckerPrager extension returns to the apex", hNew.apex);
            check
            (
                "druckerPrager apex stress is p = xi c/eta, zero deviator",
                mag(sigNew - pApex*symmTensor(I)) < 1e-9*pApex
            );
        }
    }

    // Mohr-Coulomb (WP 2.5b): selection and a compression edge return
    {
        const scalar phi = 30, c = 10e3, psi = 5;
        dictionary mcDict(elDict);
        mcDict.set("type", word("mohrCoulomb"));
        mcDict.add("phi", phi);
        mcDict.add("c", c);
        mcDict.add("psi", psi);
        autoPtr<constitutiveModel> mc(constitutiveModel::New(mcDict));
        check("selected mohrCoulomb", word(mc->kernel().name()) == "mohrCoulomb");

        // Isochoric axial compression (uniaxial strain alone stays elastic)
        symmTensor deps(Zero);
        deps.xx() = deps.yy() = 0.005;
        deps.zz() = -0.01;
        const symmTensor sig0(-1e5*symmTensor(I));
        symmTensor sigNew;
        constitutiveModel::history h, hNew;
        mc->update(deps, sig0, h, sigNew, hNew);

        const vector w = eigenValues(sigNew);   // ascending
        const scalar sp = Foam::sin(degToRad(phi));
        const scalar f =
            (w.z() - w.x()) + (w.z() + w.x())*sp
          - 2*c*Foam::cos(degToRad(phi));
        Info<< "    mohrCoulomb axial compression: sigma = " << sigNew
            << "  f = " << f << "  epsPEq = " << hNew.epsPEq << nl;
        check("mohrCoulomb returns to f = 0", mag(f) < 1e-9*mag(w.x()));
        check
        (
            "mohrCoulomb on the compression edge sxx = syy",
            mag(sigNew.xx() - sigNew.yy()) < 1e-9*mag(w.x())
        );
        check("mohrCoulomb plastic strain grows", hNew.epsPEq > 0);
    }

    // Mohr-Coulomb tension cut-off (WP 2.5c): optional sigmaT keyword;
    // hydrostatic stretch from 100 kPa compression ends at the tension apex
    {
        const scalar phi = 30, c = 10e3, psi = 5, sigmaT = 5e3;
        dictionary mcDict(elDict);
        mcDict.set("type", word("mohrCoulomb"));
        mcDict.add("phi", phi);
        mcDict.add("c", c);
        mcDict.add("psi", psi);
        mcDict.add("sigmaT", sigmaT);
        autoPtr<constitutiveModel> mc(constitutiveModel::New(mcDict));
        const auto* k =
            dynamic_cast<const cel::constitutive::MohrCoulomb*>(&mc->kernel());
        check("mohrCoulomb reads sigmaT", k && k->hasTensionCutOff());

        const symmTensor deps(0.01*symmTensor(I));
        const symmTensor sig0(-1e5*symmTensor(I));
        symmTensor sigNew;
        constitutiveModel::history h, hNew;
        mc->update(deps, sig0, h, sigNew, hNew);
        Info<< "    mohrCoulomb hydrostatic tension: sigma = " << sigNew
            << "  apex = " << hNew.apex << nl;
        check
        (
            "mohrCoulomb hydrostatic tension capped at sigmaT",
            mag(sigNew - sigmaT*symmTensor(I)) < 1e-9*sigmaT && hNew.apex
        );
    }

    Info<< nl << (nFail ? "FAILED" : "All wrapper checks passed") << endl;
    return nFail ? 1 : 0;
}
