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

#pragma once

#include <string>
#include <vector>

struct BenchmarkConfig
{
    int sceneIndex = 4;
    int steps = 600;
    float dt = -1.0f;
    int iterations = -1;
    float alpha = -1.0f;
    float betaLin = -1.0f;
    float betaAng = -1.0f;
    float gamma = -1.0f;
    int globalCorrectionEnabled = -1;
    int globalCorrectionIterations = -1;
    float globalCorrectionDamping = -1.0f;
    float globalCorrectionScale = -1.0f;
    std::string outputCsv = "test/benchmark_results.csv";
};

struct BenchmarkSweepConfig
{
    std::vector<float> dt;
    std::vector<int> iterations;
    std::vector<float> alpha;
    std::vector<float> betaLin;
    std::vector<float> betaAng;
    std::vector<float> gamma;
};

bool parseBenchmarkArgs(int argc, char **argv, BenchmarkConfig &config, BenchmarkSweepConfig &sweeps, std::string &errorMessage);
void printBenchmarkHelp();
