/*
 * Copyright (c) 2026 Chris Giles
 *
 * Permission to use, copy, modify, distribute and sell this software
 * and its documentation for any purpose is hereby granted without fee,
 * provided that the above copyright notice appear in all copies.
 * Chris Giles makes no representations about the suitability
 * of this software for any purpose.
 * It is provided "as is" without express or implied warranty.
 */

#include "benchmark.h"

#include <algorithm>
#include <chrono>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "scenes.h"
#include "solver.h"

namespace
{
std::string normalizeToken(const std::string &value)
{
    std::string out;
    out.reserve(value.size());

    for (char c : value)
    {
        if (c == ' ' || c == '_' || c == '-')
            continue;
        out.push_back((char)std::tolower((unsigned char)c));
    }

    return out;
}

bool parseFloat(const std::string &value, float &parsed)
{
    char *end = nullptr;
    parsed = std::strtof(value.c_str(), &end);
    return end && *end == '\0' && std::isfinite(parsed);
}

bool parseInt(const std::string &value, int &parsed)
{
    char *end = nullptr;
    long v = std::strtol(value.c_str(), &end, 10);
    if (!(end && *end == '\0'))
        return false;
    if (v < (long)std::numeric_limits<int>::min() || v > (long)std::numeric_limits<int>::max())
        return false;
    parsed = (int)v;
    return true;
}

bool parseFloatList(const std::string &value, std::vector<float> &out)
{
    out.clear();
    std::stringstream ss(value);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        float v = 0.0f;
        if (!parseFloat(item, v))
            return false;
        out.push_back(v);
    }

    return !out.empty();
}

bool parseIntList(const std::string &value, std::vector<int> &out)
{
    out.clear();
    std::stringstream ss(value);
    std::string item;

    while (std::getline(ss, item, ','))
    {
        int v = 0;
        if (!parseInt(item, v))
            return false;
        out.push_back(v);
    }

    return !out.empty();
}

int resolveSceneIndex(const std::string &token)
{
    int sceneIndex = -1;
    if (parseInt(token, sceneIndex))
    {
        if (sceneIndex >= 0 && sceneIndex < sceneCount)
            return sceneIndex;
        return -1;
    }

    std::string query = normalizeToken(token);
    for (int i = 0; i < sceneCount; ++i)
    {
        if (normalizeToken(sceneNames[i]) == query)
            return i;
    }

    return -1;
}

template <typename T>
std::vector<T> listOrBase(const std::vector<T> &list, T base)
{
    if (!list.empty())
        return list;
    return {base};
}

void applyRunSpec(Solver &solver, float dt, int iterations, float alpha, float betaLin, float betaAng, float gamma,
                  bool globalCorrectionEnabled, int globalCorrectionIterations,
                  float globalCorrectionDamping, float globalCorrectionScale)
{
    solver.dt = dt;
    solver.iterations = iterations;
    solver.alpha = alpha;
    solver.betaLin = betaLin;
    solver.betaAng = betaAng;
    solver.gamma = gamma;
    solver.globalCorrectionEnabled = globalCorrectionEnabled;
    solver.globalCorrectionIterations = globalCorrectionIterations;
    solver.globalCorrectionDamping = globalCorrectionDamping;
    solver.globalCorrectionScale = globalCorrectionScale;
}
} // namespace

void printBenchmarkHelp()
{
    std::cout << "Usage: avbd_benchmark [options]\n"
              << "\n"
              << "Options:\n"
              << "  --scene <index|name>            Scene index or scene name\n"
              << "  --steps <int>                   Number of simulation steps\n"
              << "  --dt <float>                    Timestep override\n"
              << "  --iterations <int>              Solver iterations override\n"
              << "  --alpha <float>                 Stabilization parameter override\n"
              << "  --betaLin <float>               Linear penalty ramp override\n"
              << "  --betaAng <float>               Angular penalty ramp override\n"
              << "  --gamma <float>                 Warmstart decay override\n"
              << "  --global-correction <0|1>       Enable or disable global correction\n"
              << "  --global-correction-iters <int> Global correction inner sweeps\n"
              << "  --global-correction-damping <f> Global correction damping\n"
              << "  --global-correction-scale <f>   Global correction step scale\n"
              << "  --output <path>                 Output CSV path\n"
              << "  --sweep-dt <v1,v2,...>          Sweep dt values\n"
              << "  --sweep-iterations <v1,v2,...>  Sweep iteration values\n"
              << "  --sweep-alpha <v1,v2,...>       Sweep alpha values\n"
              << "  --sweep-betaLin <v1,v2,...>     Sweep betaLin values\n"
              << "  --sweep-betaAng <v1,v2,...>     Sweep betaAng values\n"
              << "  --sweep-gamma <v1,v2,...>       Sweep gamma values\n"
              << "  --help                          Show this help\n"
              << "\n"
              << "Available scenes:\n";

    for (int i = 0; i < sceneCount; ++i)
        std::cout << "  " << i << ": " << sceneNames[i] << "\n";
}

bool parseBenchmarkArgs(int argc, char **argv, BenchmarkConfig &config, BenchmarkSweepConfig &sweeps, std::string &errorMessage)
{
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        auto needValue = [&](const char *option) -> const char * {
            if (i + 1 >= argc)
            {
                errorMessage = std::string("Missing value for option: ") + option;
                return nullptr;
            }
            ++i;
            return argv[i];
        };

        if (arg == "--help")
        {
            errorMessage = "help";
            return false;
        }
        else if (arg == "--scene")
        {
            const char *value = needValue("--scene");
            if (!value)
                return false;

            int sceneIndex = resolveSceneIndex(value);
            if (sceneIndex < 0)
            {
                errorMessage = std::string("Invalid scene: ") + value;
                return false;
            }
            config.sceneIndex = sceneIndex;
        }
        else if (arg == "--steps")
        {
            const char *value = needValue("--steps");
            if (!value)
                return false;

            if (!parseInt(value, config.steps) || config.steps < 0)
            {
                errorMessage = std::string("Invalid steps value: ") + value;
                return false;
            }
        }
        else if (arg == "--dt")
        {
            const char *value = needValue("--dt");
            if (!value)
                return false;

            if (!parseFloat(value, config.dt) || config.dt <= 0.0f)
            {
                errorMessage = std::string("Invalid dt value: ") + value;
                return false;
            }
        }
        else if (arg == "--iterations")
        {
            const char *value = needValue("--iterations");
            if (!value)
                return false;

            if (!parseInt(value, config.iterations) || config.iterations <= 0)
            {
                errorMessage = std::string("Invalid iterations value: ") + value;
                return false;
            }
        }
        else if (arg == "--alpha")
        {
            const char *value = needValue("--alpha");
            if (!value)
                return false;

            if (!parseFloat(value, config.alpha) || config.alpha < 0.0f)
            {
                errorMessage = std::string("Invalid alpha value: ") + value;
                return false;
            }
        }
        else if (arg == "--betaLin")
        {
            const char *value = needValue("--betaLin");
            if (!value)
                return false;

            if (!parseFloat(value, config.betaLin) || config.betaLin < 0.0f)
            {
                errorMessage = std::string("Invalid betaLin value: ") + value;
                return false;
            }
        }
        else if (arg == "--betaAng")
        {
            const char *value = needValue("--betaAng");
            if (!value)
                return false;

            if (!parseFloat(value, config.betaAng) || config.betaAng < 0.0f)
            {
                errorMessage = std::string("Invalid betaAng value: ") + value;
                return false;
            }
        }
        else if (arg == "--gamma")
        {
            const char *value = needValue("--gamma");
            if (!value)
                return false;

            if (!parseFloat(value, config.gamma) || config.gamma < 0.0f)
            {
                errorMessage = std::string("Invalid gamma value: ") + value;
                return false;
            }
        }
        else if (arg == "--output")
        {
            const char *value = needValue("--output");
            if (!value)
                return false;
            config.outputCsv = value;
        }
        else if (arg == "--global-correction")
        {
            const char *value = needValue("--global-correction");
            if (!value)
                return false;

            int enabled = -1;
            if (!parseInt(value, enabled) || (enabled != 0 && enabled != 1))
            {
                errorMessage = std::string("Invalid global-correction value: ") + value;
                return false;
            }
            config.globalCorrectionEnabled = enabled;
        }
        else if (arg == "--global-correction-iters")
        {
            const char *value = needValue("--global-correction-iters");
            if (!value)
                return false;

            if (!parseInt(value, config.globalCorrectionIterations) || config.globalCorrectionIterations < 0)
            {
                errorMessage = std::string("Invalid global-correction-iters value: ") + value;
                return false;
            }
        }
        else if (arg == "--global-correction-damping")
        {
            const char *value = needValue("--global-correction-damping");
            if (!value)
                return false;

            if (!parseFloat(value, config.globalCorrectionDamping) || config.globalCorrectionDamping < 0.0f)
            {
                errorMessage = std::string("Invalid global-correction-damping value: ") + value;
                return false;
            }
        }
        else if (arg == "--global-correction-scale")
        {
            const char *value = needValue("--global-correction-scale");
            if (!value)
                return false;

            if (!parseFloat(value, config.globalCorrectionScale) || config.globalCorrectionScale < 0.0f)
            {
                errorMessage = std::string("Invalid global-correction-scale value: ") + value;
                return false;
            }
        }
        else if (arg == "--sweep-dt")
        {
            const char *value = needValue("--sweep-dt");
            if (!value)
                return false;
            if (!parseFloatList(value, sweeps.dt))
            {
                errorMessage = std::string("Invalid sweep-dt list: ") + value;
                return false;
            }
        }
        else if (arg == "--sweep-iterations")
        {
            const char *value = needValue("--sweep-iterations");
            if (!value)
                return false;
            if (!parseIntList(value, sweeps.iterations))
            {
                errorMessage = std::string("Invalid sweep-iterations list: ") + value;
                return false;
            }
        }
        else if (arg == "--sweep-alpha")
        {
            const char *value = needValue("--sweep-alpha");
            if (!value)
                return false;
            if (!parseFloatList(value, sweeps.alpha))
            {
                errorMessage = std::string("Invalid sweep-alpha list: ") + value;
                return false;
            }
        }
        else if (arg == "--sweep-betaLin")
        {
            const char *value = needValue("--sweep-betaLin");
            if (!value)
                return false;
            if (!parseFloatList(value, sweeps.betaLin))
            {
                errorMessage = std::string("Invalid sweep-betaLin list: ") + value;
                return false;
            }
        }
        else if (arg == "--sweep-betaAng")
        {
            const char *value = needValue("--sweep-betaAng");
            if (!value)
                return false;
            if (!parseFloatList(value, sweeps.betaAng))
            {
                errorMessage = std::string("Invalid sweep-betaAng list: ") + value;
                return false;
            }
        }
        else if (arg == "--sweep-gamma")
        {
            const char *value = needValue("--sweep-gamma");
            if (!value)
                return false;
            if (!parseFloatList(value, sweeps.gamma))
            {
                errorMessage = std::string("Invalid sweep-gamma list: ") + value;
                return false;
            }
        }
        else
        {
            errorMessage = std::string("Unknown option: ") + arg;
            return false;
        }
    }

    return true;
}

int main(int argc, char **argv)
{
    BenchmarkConfig config;
    BenchmarkSweepConfig sweeps;
    std::string parseError;

    if (!parseBenchmarkArgs(argc, argv, config, sweeps, parseError))
    {
        if (parseError != "help")
            std::cerr << parseError << "\n\n";
        printBenchmarkHelp();
        return parseError == "help" ? 0 : 1;
    }

    Solver defaults;
    float baseDt = config.dt > 0.0f ? config.dt : defaults.dt;
    int baseIterations = config.iterations > 0 ? config.iterations : defaults.iterations;
    float baseAlpha = config.alpha >= 0.0f ? config.alpha : defaults.alpha;
    float baseBetaLin = config.betaLin >= 0.0f ? config.betaLin : defaults.betaLin;
    float baseBetaAng = config.betaAng >= 0.0f ? config.betaAng : defaults.betaAng;
    float baseGamma = config.gamma >= 0.0f ? config.gamma : defaults.gamma;
    bool baseGlobalCorrectionEnabled = config.globalCorrectionEnabled >= 0 ? (config.globalCorrectionEnabled == 1) : defaults.globalCorrectionEnabled;
    int baseGlobalCorrectionIterations = config.globalCorrectionIterations >= 0 ? config.globalCorrectionIterations : defaults.globalCorrectionIterations;
    float baseGlobalCorrectionDamping = config.globalCorrectionDamping >= 0.0f ? config.globalCorrectionDamping : defaults.globalCorrectionDamping;
    float baseGlobalCorrectionScale = config.globalCorrectionScale >= 0.0f ? config.globalCorrectionScale : defaults.globalCorrectionScale;

    std::vector<float> dtValues = listOrBase(sweeps.dt, baseDt);
    std::vector<int> iterationValues = listOrBase(sweeps.iterations, baseIterations);
    std::vector<float> alphaValues = listOrBase(sweeps.alpha, baseAlpha);
    std::vector<float> betaLinValues = listOrBase(sweeps.betaLin, baseBetaLin);
    std::vector<float> betaAngValues = listOrBase(sweeps.betaAng, baseBetaAng);
    std::vector<float> gammaValues = listOrBase(sweeps.gamma, baseGamma);

    std::filesystem::path outputPath(config.outputCsv);
    if (!outputPath.parent_path().empty())
        std::filesystem::create_directories(outputPath.parent_path());

    std::ofstream csv(outputPath, std::ios::out | std::ios::trunc);
    if (!csv.is_open())
    {
        std::cerr << "Failed to open output CSV: " << config.outputCsv << "\n";
        return 1;
    }

    csv << std::setprecision(9);
    csv << "runIndex,sceneIndex,sceneName,stepIndex,dt,iterations,alpha,betaLin,betaAng,gamma,"
        << "globalCorrectionEnabled,globalCorrectionIterations,globalCorrectionDamping,globalCorrectionScale,"
        << "activeBodies,activeForces,activeContacts,maxJointResidual,avgJointResidual,"
        << "maxPenetration,avgPenetration,maxJointForce,avgJointForce,maxJointTorque,avgJointTorque,"
        << "maxContactNormalForce,avgContactNormalForce,maxContactTangentialForce,avgContactTangentialForce,"
        << "avgContactSlip,stickingContacts,linearMomentumNorm,angularMomentumNorm,"
        << "kineticEnergy,potentialEnergy,totalEnergy,wallClockMs\n";

    int runIndex = 0;
    for (float dt : dtValues)
    {
        for (int iterations : iterationValues)
        {
            for (float alpha : alphaValues)
            {
                for (float betaLin : betaLinValues)
                {
                    for (float betaAng : betaAngValues)
                    {
                        for (float gamma : gammaValues)
                        {
                            Solver solver;
                            applyRunSpec(solver, dt, iterations, alpha, betaLin, betaAng, gamma,
                                         baseGlobalCorrectionEnabled, baseGlobalCorrectionIterations,
                                         baseGlobalCorrectionDamping, baseGlobalCorrectionScale);
                            scenes[config.sceneIndex](&solver);
                            solver.diagnosticsEnabled = true;
                            solver.stepIndex = 0;
                            solver.diagnostics.history.clear();

                            std::cout << "Run " << runIndex << ": scene=" << sceneNames[config.sceneIndex]
                                      << " steps=" << config.steps
                                      << " dt=" << dt
                                      << " iterations=" << iterations
                                      << " alpha=" << alpha
                                      << " betaLin=" << betaLin
                                      << " betaAng=" << betaAng
                                      << " gamma=" << gamma
                                      << " globalCorrection=" << (solver.globalCorrectionEnabled ? 1 : 0)
                                      << " globalCorrectionIters=" << solver.globalCorrectionIterations
                                      << "\n";

                            for (int step = 0; step < config.steps; ++step)
                            {
                                auto t0 = std::chrono::high_resolution_clock::now();
                                solver.step();
                                auto t1 = std::chrono::high_resolution_clock::now();

                                double wallClockMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
                                const StepDiagnostics &d = solver.diagnostics.history.back();

                                csv << runIndex << ','
                                    << config.sceneIndex << ','
                                    << '"' << sceneNames[config.sceneIndex] << '"' << ','
                                    << d.stepIndex << ','
                                    << dt << ','
                                    << iterations << ','
                                    << alpha << ','
                                    << betaLin << ','
                                    << betaAng << ','
                                    << gamma << ','
                                    << (solver.globalCorrectionEnabled ? 1 : 0) << ','
                                    << solver.globalCorrectionIterations << ','
                                    << solver.globalCorrectionDamping << ','
                                    << solver.globalCorrectionScale << ','
                                    << d.activeBodies << ','
                                    << d.activeForces << ','
                                    << d.activeContacts << ','
                                    << d.maxJointResidual << ','
                                    << d.avgJointResidual << ','
                                    << d.maxPenetration << ','
                                    << d.avgPenetration << ','
                                    << d.maxJointForce << ','
                                    << d.avgJointForce << ','
                                    << d.maxJointTorque << ','
                                    << d.avgJointTorque << ','
                                    << d.maxContactNormalForce << ','
                                    << d.avgContactNormalForce << ','
                                    << d.maxContactTangentialForce << ','
                                    << d.avgContactTangentialForce << ','
                                    << d.avgContactSlip << ','
                                    << d.stickingContacts << ','
                                    << d.linearMomentumNorm << ','
                                    << d.angularMomentumNorm << ','
                                    << d.kineticEnergy << ','
                                    << d.potentialEnergy << ','
                                    << d.totalEnergy << ','
                                    << wallClockMs << '\n';
                            }

                            csv.flush();
                            runIndex++;
                        }
                    }
                }
            }
        }
    }

    std::cout << "Benchmark complete. Wrote " << runIndex << " run(s) to " << config.outputCsv << "\n";
    return 0;
}
