/*---------------------------------------------------------------------------*\
Application
    compareReference

Description
    WP 2.2 kernel test. Reads the reference paths written by
    tests/gen_constitutive_reference.py (Python models of WP 1.4), replays each
    strain increment through the C++ kernel and compares stress and history
    step by step.

    Error measure for each case: max over steps of |C++ - Python|, divided by
    the largest reference value on the path (stress: max |sigma|; plastic
    strain: max of |epsP| and epsPEq). Passes if every case is below 1e-10,
    the apex flags agree on every step, every plastic case has plastic steps,
    and every Drucker-Prager random path has apex returns (WP 2.5a).

    Mohr-Coulomb cases (WP 2.5b) also count which return each step took and
    require the returns the path is there to test: triaxial compression the
    s1 = s2 edge, triaxial extension the s2 = s3 edge, uniaxial strain the
    apex, and the random paths at least one edge and one apex return.
    Paths with a tension cut-off (WP 2.5c) must reach the cut-off and never
    the Mohr-Coulomb apex; triaxial tension and uniaxial strain need the
    tension plane, hydrostatic tension the tension apex, tension-then-shear a mixed Mohr-Coulomb/tension corner, and
    the random paths every tension return (plane, edge, apex, corner).

    WP 3.6a: a case with a "dt <value>" line (tests/gen_mu_i_reference.py) is
    replayed with updateRate(dt) and must be rate dependent; return counting
    is skipped for it. Its apex flag may differ from Python's on a step where
    the cone return lands at the apex to round-off (the root of the mu(I)
    return decides which side); stress and history must still match.

    compareReference --mu2-equal-mus <file> (WP 3.6a): replays only the
    Drucker-Prager and Mohr-Coulomb cases, with mu2 = tan(phi) (and I0, d,
    rhoS) added and updateRate(dt = 1e-4): the model must not be rate
    dependent and must give the rate-independent C++ model's path along the
    file's increments bit for bit (tolerance 0).

    Build and run: tests/test_constitutive_cpp.sh, tests/test_mu_i.sh
\*---------------------------------------------------------------------------*/

#include "celConstitutiveKernel.H"
#include "celMohrCoulomb.H"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace cel::constitutive;

namespace
{

struct Step
{
    double deps[6];
    double sigma[6];
    double epsP[6];
    double epsPEq;
    int apex;
};

struct Case
{
    std::string name;
    std::string model;
    Parameters params;
    Options options;
    double dt = 0;
    double sigma0[6];
    std::vector<Step> steps;
};


std::vector<Case> readCases(const char* file)
{
    std::ifstream in(file);
    if (!in)
    {
        throw std::runtime_error(std::string("cannot open ") + file);
    }

    std::vector<Case> cases;
    std::string line;
    while (std::getline(in, line))
    {
        std::istringstream is(line);
        std::string key;
        is >> key;
        if (key == "case")
        {
            Case c;
            int n;
            is >> c.name >> c.model >> n;
            for (int i = 0; i < n; ++i)
            {
                std::string k;
                double v;
                is >> k >> v;
                c.params[k] = v;
            }
            cases.push_back(c);
        }
        else if (key == "dt")
        {
            is >> cases.back().dt;
        }
        else if (key == "option")
        {
            std::string k, v;
            is >> k >> v;
            cases.back().options[k] = v;
        }
        else if (key == "init")
        {
            for (double& v : cases.back().sigma0) is >> v;
        }
        else if (key == "step")
        {
            Step s;
            for (double& v : s.deps) is >> v;
            for (double& v : s.sigma) is >> v;
            for (double& v : s.epsP) is >> v;
            is >> s.epsPEq >> s.apex;
            cases.back().steps.push_back(s);
        }
        if (!is && key != "end")
        {
            throw std::runtime_error("bad line: " + line);
        }
    }
    return cases;
}

} // End anonymous namespace


int main(int argc, char* argv[])
{
    const bool muSame =
        argc == 3 && std::string(argv[1]) == "--mu2-equal-mus";
    if (argc != 2 && !muSame)
    {
        std::cerr
            << "usage: compareReference [--mu2-equal-mus] <reference file>\n";
        return 2;
    }

    const double tol = muSame ? 0 : 1e-10;
    int failures = 0;

    std::vector<Case> cases = readCases(argv[argc - 1]);
    if (muSame)
    {
        std::vector<Case> kept;
        for (Case c : cases)
        {
            if (c.model == "druckerPrager" || c.model == "mohrCoulomb")
            {
                // The expected path is the rate-independent C++ model's own
                const auto plain = makeModel(c.model, c.params, c.options);
                double sig[6];
                std::copy(c.sigma0, c.sigma0 + 6, sig);
                History h;
                for (Step& st : c.steps)
                {
                    History hNew;
                    plain->update(st.deps, sig, h, st.sigma, hNew);
                    std::copy(hNew.epsP, hNew.epsP + 6, st.epsP);
                    st.epsPEq = hNew.epsPEq;
                    st.apex = hNew.apex;
                    std::copy(st.sigma, st.sigma + 6, sig);
                    h = hNew;
                }
                c.params["mu2"] =
                    std::tan(c.params["phi"]*3.14159265358979323846/180.0);
                c.params["I0"] = 0.279;
                c.params["d"] = 1e-3;
                c.params["rhoS"] = 2650;
                c.dt = 1e-4;
                kept.push_back(c);
            }
        }
        cases = kept;
    }
    std::cout << "Read " << cases.size() << " reference paths\n";

    for (const Case& c : cases)
    {
        const auto model = makeModel(c.model, c.params, c.options);
        if (muSame && c.dt > 0 && model->rateDependent())
        {
            std::cout << "FAIL  " << c.name
                << ": mu2 = tan(phi) gave a rate-dependent model\n";
            ++failures;
            continue;
        }

        double sig[6];
        std::copy(c.sigma0, c.sigma0 + 6, sig);
        History h;

        double sigScale = 0, epsScale = 0, sigErr = 0, epsErr = 0;
        int apexMismatch = 0, nPlastic = 0, nApex = 0;

        const bool rate = c.dt > 0;
        // Rate-dependent cases must give a rate-dependent model, except in
        // the --mu2-equal-mus mode, where none may
        const bool wantRate = rate && !muSame;
        const auto* mc =
            rate ? nullptr : dynamic_cast<const MohrCoulomb*>(model.get());
        int nRet[MohrCoulomb::nReturns] = {};

        for (const Step& s : c.steps)
        {
            double sigNew[6];
            History hNew;
            if (rate)
            {
                model->updateRate(s.deps, c.dt, sig, h, sigNew, hNew);
            }
            else if (mc)
            {
                ++nRet[int(mc->updateReturn(s.deps, sig, h, sigNew, hNew))];
            }
            else
            {
                model->update(s.deps, sig, h, sigNew, hNew);
            }

            for (int i = 0; i < 6; ++i)
            {
                sigScale = std::max(sigScale, std::abs(s.sigma[i]));
                epsScale = std::max(epsScale, std::abs(s.epsP[i]));
                sigErr = std::max(sigErr, std::abs(sigNew[i] - s.sigma[i]));
                epsErr = std::max(epsErr, std::abs(hNew.epsP[i] - s.epsP[i]));
            }
            epsScale = std::max(epsScale, std::abs(s.epsPEq));
            epsErr = std::max(epsErr, std::abs(hNew.epsPEq - s.epsPEq));
            apexMismatch += (int(hNew.apex) != s.apex);
            nPlastic += (hNew.epsPEq != h.epsPEq);
            nApex += hNew.apex;

            std::copy(sigNew, sigNew + 6, sig);
            h = hNew;
        }

        const double relSig = sigScale > 0 ? sigErr/sigScale : sigErr;
        // Elastic paths have no plastic strain: then require exact zeros
        const double relEps = epsScale > 0 ? epsErr/epsScale : epsErr;
        const bool needPlastic = (c.model != "linearElastic");
        // Drucker-Prager random paths must exercise the apex return (WP 2.5a)
        const bool needApex =
            c.model == "druckerPrager"
         && c.name.find("random") != std::string::npos;

        // Returns each Mohr-Coulomb path must exercise
        bool returnsOk = true;
        std::string returns;
        if (mc)
        {
            typedef MohrCoulomb::Return R;
            const int plane = nRet[int(R::plane)];
            const int comp = nRet[int(R::compressionEdge)];
            const int ext = nRet[int(R::extensionEdge)];
            const int apx = nRet[int(R::apex)];
            const int tPlane = nRet[int(R::tensionPlane)];
            const int tEdge = nRet[int(R::tensionEdge)];
            const int tApex = nRet[int(R::tensionApex)];
            const int tCorner = nRet[int(R::tensionCorner)];
            const int tAll = tPlane + tEdge + tApex + tCorner;
            auto has = [&](const char* tag)
            {
                return c.name.find(tag) != std::string::npos;
            };
            if (has("triaxialCompression")) returnsOk = comp > 0;
            if (has("triaxialExtension")) returnsOk = ext > 0;
            if (has("uniaxialStrain")) returnsOk = apx > 0;
            if (has("random")) returnsOk = (comp + ext) > 0 && apx > 0;
            if (mc->hasTensionCutOff())
            {
                // WP 2.5c paths: each must reach the cut-off, and the
                // apex flag must only come from the tension apex
                returnsOk = tAll > 0 && apx == 0;
                if (has("hydrostaticTension")) returnsOk = returnsOk && tApex > 0;
                if (has("triaxialTension")) returnsOk = returnsOk && tPlane > 0;
                if (has("uniaxialStrain")) returnsOk = returnsOk && tPlane > 0;
                if (has("tensionShear")) returnsOk = returnsOk && tCorner > 0;
                if (has("random"))
                {
                    returnsOk = returnsOk
                     && tPlane > 0 && tEdge > 0 && tApex > 0 && tCorner > 0;
                }
            }

            char rb[256];
            std::snprintf
            (
                rb, sizeof(rb),
                "      returns: plane %d, compression edge %d, "
                "extension edge %d, apex %d; tension plane %d, edge %d, "
                "apex %d, corner %d%s\n",
                plane, comp, ext, apx, tPlane, tEdge, tApex, tCorner,
                returnsOk ? "" : "  (MISSING)"
            );
            returns = rb;
        }

        const bool ok =
            (muSame ? relSig == 0 : relSig < tol)
         && (epsScale > 0 && !muSame ? relEps < tol : epsErr == 0)
         && ((rate && !muSame) || apexMismatch == 0)
         && (!needPlastic || nPlastic > 0)
         && (!needApex || nApex > 0)
         && returnsOk
         && model->rateDependent() == wantRate;

        char buf[256];
        std::snprintf
        (
            buf, sizeof(buf),
            "%s  %-42s steps %4zu  plastic %4d  apex %3d  stress err %.2e  "
            "history err %.2e  apex mismatches %d\n",
            ok ? "PASS" : "FAIL", c.name.c_str(), c.steps.size(), nPlastic,
            nApex, relSig, relEps, apexMismatch
        );
        std::cout << buf << returns;
        failures += !ok;
    }

    if (cases.empty())
    {
        std::cout << "FAIL  no reference paths read\n";
        return 1;
    }

    std::cout << "\n" << (failures ? "FAILED" : "All paths match")
        << " (tolerance " << tol << ")\n";
    return failures ? 1 : 0;
}
