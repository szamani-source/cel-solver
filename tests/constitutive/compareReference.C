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

    Build and run: tests/test_constitutive_cpp.sh
\*---------------------------------------------------------------------------*/

#include "celConstitutiveKernel.H"

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
    if (argc != 2)
    {
        std::cerr << "usage: compareReference <reference file>\n";
        return 2;
    }

    const double tol = 1e-10;
    int failures = 0;

    const std::vector<Case> cases = readCases(argv[1]);
    std::cout << "Read " << cases.size() << " reference paths\n";

    for (const Case& c : cases)
    {
        const auto model = makeModel(c.model, c.params, c.options);

        double sig[6];
        std::copy(c.sigma0, c.sigma0 + 6, sig);
        History h;

        double sigScale = 0, epsScale = 0, sigErr = 0, epsErr = 0;
        int apexMismatch = 0, nPlastic = 0, nApex = 0;

        for (const Step& s : c.steps)
        {
            double sigNew[6];
            History hNew;
            model->update(s.deps, sig, h, sigNew, hNew);

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

        const bool ok =
            relSig < tol
         && (epsScale > 0 ? relEps < tol : epsErr == 0)
         && apexMismatch == 0
         && (!needPlastic || nPlastic > 0)
         && (!needApex || nApex > 0);

        char buf[256];
        std::snprintf
        (
            buf, sizeof(buf),
            "%s  %-42s steps %4zu  plastic %4d  apex %3d  stress err %.2e  "
            "history err %.2e  apex mismatches %d\n",
            ok ? "PASS" : "FAIL", c.name.c_str(), c.steps.size(), nPlastic,
            nApex, relSig, relEps, apexMismatch
        );
        std::cout << buf;
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
