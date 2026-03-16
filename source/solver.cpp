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

#include "solver.h"

#include <unordered_map>

namespace
{
struct CorrectionBodyState
{
    Rigid *body;
    float invMass;
    float3x3 invInertiaWorld;
};

struct JointConstraint
{
    int bodyA;
    int bodyB;
    float3 residual;
    float3x3 jALin;
    float3x3 jAAng;
    float3x3 jBLin;
    float3x3 jBAng;
};

inline bool solve3x3(float3x3 A, float3 b, float3 &x)
{
    float m[3][4] = {
        {A[0][0], A[0][1], A[0][2], b[0]},
        {A[1][0], A[1][1], A[1][2], b[1]},
        {A[2][0], A[2][1], A[2][2], b[2]}};

    for (int col = 0; col < 3; ++col)
    {
        int pivot = col;
        float pivotAbs = fabsf(m[col][col]);
        for (int row = col + 1; row < 3; ++row)
        {
            float v = fabsf(m[row][col]);
            if (v > pivotAbs)
            {
                pivot = row;
                pivotAbs = v;
            }
        }

        if (pivotAbs < 1.0e-8f)
            return false;

        if (pivot != col)
        {
            for (int c = col; c < 4; ++c)
            {
                float tmp = m[col][c];
                m[col][c] = m[pivot][c];
                m[pivot][c] = tmp;
            }
        }

        for (int row = col + 1; row < 3; ++row)
        {
            float factor = m[row][col] / m[col][col];
            for (int c = col; c < 4; ++c)
                m[row][c] -= factor * m[col][c];
        }
    }

    x[2] = m[2][3] / m[2][2];
    x[1] = (m[1][3] - m[1][2] * x[2]) / m[1][1];
    x[0] = (m[0][3] - m[0][2] * x[2] - m[0][1] * x[1]) / m[0][0];
    return true;
}
} // namespace

Solver::Solver()
    : bodies(0), forces(0), diagnosticsEnabled(false), stepIndex(0)
{
    defaultParams();
}

Solver::~Solver()
{
    clear();
}

Rigid *Solver::pick(float3 origin, float3 dir, float3 &local)
{
    const float epsilon = 1.0e-6f;
    float bestT = INFINITY;
    Rigid *bestBody = 0;
    float3 bestLocal = {0, 0, 0};

    // Ray-cast against each OBB by transforming the ray into body local space.
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if (body->mass <= 0.0f)
            continue;

        quat invRot = conjugate(body->positionAng);
        float3 o = rotate(invRot, origin - body->positionLin);
        float3 d = rotate(invRot, dir);
        float3 half = body->size * 0.5f;

        float tEnter = 0.0f;
        float tExit = INFINITY;
        bool hit = true;

        for (int i = 0; i < 3; ++i)
        {
            if (fabsf(d[i]) < epsilon)
            {
                if (o[i] < -half[i] || o[i] > half[i])
                {
                    hit = false;
                    break;
                }
                continue;
            }

            float invD = 1.0f / d[i];
            float t0 = (-half[i] - o[i]) * invD;
            float t1 = (half[i] - o[i]) * invD;
            if (t0 > t1)
            {
                float tmp = t0;
                t0 = t1;
                t1 = tmp;
            }

            tEnter = max(tEnter, t0);
            tExit = min(tExit, t1);
            if (tEnter > tExit)
            {
                hit = false;
                break;
            }
        }

        if (!hit)
            continue;

        float tHit = tEnter >= 0.0f ? tEnter : tExit;
        if (tHit < 0.0f)
            continue;

        if (tHit < bestT)
        {
            bestT = tHit;
            bestBody = body;
            bestLocal = o + d * tHit;
        }
    }

    if (!bestBody)
        return 0;

    local = bestLocal;
    return bestBody;
}

void Solver::clear()
{
    while (forces)
        delete forces;

    while (bodies)
        delete bodies;

    stepIndex = 0;
    diagnostics.history.clear();
}

void Solver::defaultParams()
{
    dt = 1.0f / 60.0f;
    gravity = -10.0f;
    iterations = 10;

    // Note: in the paper, beta is suggested to be [1, 1000]. Technically, the best choice will
    // depend on the length, mass, and constraint function scales (ie units) of your simulation,
    // along with your strategy for incrementing the penalty parameters.
    // If the value is not in the right range, you may see slower convergance for complex scenes.
    // A minor upgrade from the paper is using separate betas for constraints of different units (eg linear vs angular).
    betaLin = 10000.0f;
    betaAng = 100.0f;

    // Alpha controls how much stabilization is applied. Higher values give slower and smoother
    // error correction, and lower values are more responsive and energetic. Tune this depending
    // on your desired constraint error response.
    alpha = 0.99f;

    // Gamma controls how much the penalty and lambda values are decayed each step during warmstarting.
    // This should always be < 1 so that the penalty values can decrease (unless you use a different
    // penalty parameter strategy which does not require decay).
    gamma = 0.999f;

    globalCorrectionEnabled = true;
    globalCorrectionIterations = 2;
    globalCorrectionDamping = 1.0e-4f;
    globalCorrectionScale = 0.5f;

    diagnosticsEnabled = false;
    stepIndex = 0;
    diagnostics.history.clear();
}

void Solver::globalConstraintCorrection()
{
    if (!globalCorrectionEnabled || globalCorrectionIterations <= 0)
        return;

    std::vector<CorrectionBodyState> states;
    states.reserve(64);
    std::unordered_map<Rigid *, int> bodyToState;

    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if (body->mass <= 0.0f)
            continue;

        float invMass = body->mass > 0.0f ? 1.0f / body->mass : 0.0f;
        float3 invMoment = {
            body->moment.x > 0.0f ? 1.0f / body->moment.x : 0.0f,
            body->moment.y > 0.0f ? 1.0f / body->moment.y : 0.0f,
            body->moment.z > 0.0f ? 1.0f / body->moment.z : 0.0f};

        float3x3 invIBody = diagonal(invMoment.x, invMoment.y, invMoment.z);
        float3x3 R = rotation(body->positionAng);

        CorrectionBodyState s;
        s.body = body;
        s.invMass = invMass;
        s.invInertiaWorld = R * invIBody * transpose(R);
        int index = (int)states.size();
        states.push_back(s);
        bodyToState[body] = index;
    }

    if (states.empty())
        return;

    std::vector<JointConstraint> constraints;
    constraints.reserve(128);

    const float3x3 zero = float3x3{0, 0, 0, 0, 0, 0, 0, 0, 0};
    const float3x3 I = float3x3{1, 0, 0, 0, 1, 0, 0, 0, 1};

    for (Force *force = forces; force != 0; force = force->next)
    {
        Joint *joint = dynamic_cast<Joint *>(force);
        if (!joint || joint->broken)
            continue;

        int idxA = -1;
        int idxB = -1;

        if (joint->bodyA)
        {
            auto itA = bodyToState.find(joint->bodyA);
            if (itA != bodyToState.end())
                idxA = itA->second;
        }

        if (joint->bodyB)
        {
            auto itB = bodyToState.find(joint->bodyB);
            if (itB != bodyToState.end())
                idxB = itB->second;
        }

        if (idxA < 0 && idxB < 0)
            continue;

        if (isinf(joint->stiffnessLin))
        {
            JointConstraint c = {};
            c.bodyA = idxA;
            c.bodyB = idxB;
            c.jALin = I;
            c.jAAng = joint->bodyA ? skew(-rotate(joint->bodyA->positionAng, joint->rA)) : zero;
            c.jBLin = -I;
            c.jBAng = joint->bodyB ? skew(rotate(joint->bodyB->positionAng, joint->rB)) : zero;

            float3 worldA = joint->bodyA ? transform(joint->bodyA->positionLin, joint->bodyA->positionAng, joint->rA) : joint->rA;
            float3 worldB = transform(joint->bodyB->positionLin, joint->bodyB->positionAng, joint->rB);
            c.residual = worldA - worldB - joint->C0Lin * alpha;

            constraints.push_back(c);
        }

        if (isinf(joint->stiffnessAng) && joint->torqueArm > 1.0e-6f)
        {
            JointConstraint c = {};
            c.bodyA = idxA;
            c.bodyB = idxB;
            c.jALin = zero;
            c.jAAng = (joint->bodyA ? I : zero) * joint->torqueArm;
            c.jBLin = zero;
            c.jBAng = -I * joint->torqueArm;

            quat qa = joint->bodyA ? joint->bodyA->positionAng : quat{0, 0, 0, 1};
            c.residual = (qa - joint->bodyB->positionAng) * joint->torqueArm - joint->C0Ang * alpha;

            constraints.push_back(c);
        }
    }

    if (constraints.empty())
        return;

    struct BodyDelta
    {
        float3 lin;
        float3 ang;
    };

    std::vector<BodyDelta> delta(states.size(), BodyDelta{{0, 0, 0}, {0, 0, 0}});

    for (int sweep = 0; sweep < globalCorrectionIterations; ++sweep)
    {
        for (const JointConstraint &c : constraints)
        {
            float3 jDelta = {0, 0, 0};

            if (c.bodyA >= 0)
                jDelta += c.jALin * delta[c.bodyA].lin + c.jAAng * delta[c.bodyA].ang;

            if (c.bodyB >= 0)
                jDelta += c.jBLin * delta[c.bodyB].lin + c.jBAng * delta[c.bodyB].ang;

            float3 rhs = -(c.residual + jDelta);
            float3x3 K = diagonal(globalCorrectionDamping, globalCorrectionDamping, globalCorrectionDamping);

            if (c.bodyA >= 0)
            {
                const CorrectionBodyState &sA = states[c.bodyA];
                float3x3 invMA = diagonal(sA.invMass, sA.invMass, sA.invMass);
                K += c.jALin * invMA * transpose(c.jALin);
                K += c.jAAng * sA.invInertiaWorld * transpose(c.jAAng);
            }

            if (c.bodyB >= 0)
            {
                const CorrectionBodyState &sB = states[c.bodyB];
                float3x3 invMB = diagonal(sB.invMass, sB.invMass, sB.invMass);
                K += c.jBLin * invMB * transpose(c.jBLin);
                K += c.jBAng * sB.invInertiaWorld * transpose(c.jBAng);
            }

            float3 dLambda = {0, 0, 0};
            if (!solve3x3(K, rhs, dLambda))
                continue;

            if (c.bodyA >= 0)
            {
                CorrectionBodyState &sA = states[c.bodyA];
                delta[c.bodyA].lin += (transpose(c.jALin) * dLambda) * sA.invMass;
                delta[c.bodyA].ang += sA.invInertiaWorld * (transpose(c.jAAng) * dLambda);
            }

            if (c.bodyB >= 0)
            {
                CorrectionBodyState &sB = states[c.bodyB];
                delta[c.bodyB].lin += (transpose(c.jBLin) * dLambda) * sB.invMass;
                delta[c.bodyB].ang += sB.invInertiaWorld * (transpose(c.jBAng) * dLambda);
            }
        }
    }

    for (size_t i = 0; i < states.size(); ++i)
    {
        Rigid *body = states[i].body;
        float3 dLin = delta[i].lin * globalCorrectionScale;
        float3 dAng = delta[i].ang * globalCorrectionScale;

        body->positionLin += dLin;
        body->positionAng = body->positionAng + dAng;
    }
}

int Solver::countContacts() const
{
    int total = 0;
    for (Force *force = forces; force != 0; force = force->next)
    {
        Manifold *manifold = dynamic_cast<Manifold *>(force);
        if (manifold)
            total += manifold->numContacts;
    }
    return total;
}

void Solver::computeMomentum(float3 &linearMomentum, float3 &angularMomentum) const
{
    linearMomentum = {0, 0, 0};
    angularMomentum = {0, 0, 0};

    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if (body->mass <= 0.0f)
            continue;

        float3 p = body->velocityLin * body->mass;
        linearMomentum += p;

        float3x3 R = rotation(body->positionAng);
        float3x3 IBody = diagonal(body->moment.x, body->moment.y, body->moment.z);
        float3x3 IWorld = R * IBody * transpose(R);
        float3 spin = IWorld * body->velocityAng;
        angularMomentum += spin + cross(body->positionLin, p);
    }
}

void Solver::computeSystemEnergy(float &kineticEnergy, float &potentialEnergy) const
{
    kineticEnergy = 0.0f;
    potentialEnergy = 0.0f;

    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        if (body->mass <= 0.0f)
            continue;

        float v2 = lengthSq(body->velocityLin);
        kineticEnergy += 0.5f * body->mass * v2;

        float3x3 R = rotation(body->positionAng);
        float3x3 IBody = diagonal(body->moment.x, body->moment.y, body->moment.z);
        float3x3 IWorld = R * IBody * transpose(R);
        float3 Iw = IWorld * body->velocityAng;
        kineticEnergy += 0.5f * dot(body->velocityAng, Iw);

        potentialEnergy += -body->mass * gravity * body->positionLin.z;
    }
}

void Solver::computeConstraintResiduals(float &maxJointResidual, float &avgJointResidual, float &maxPenetration, float &avgPenetration,
                                       float &maxJointForce, float &avgJointForce, float &maxJointTorque, float &avgJointTorque,
                                       float &maxContactNormalForce, float &avgContactNormalForce, float &maxContactTangentialForce,
                                       float &avgContactTangentialForce, float &avgContactSlip, int &stickingContacts) const
{
    maxJointResidual = 0.0f;
    avgJointResidual = 0.0f;
    maxPenetration = 0.0f;
    avgPenetration = 0.0f;
    maxJointForce = 0.0f;
    avgJointForce = 0.0f;
    maxJointTorque = 0.0f;
    avgJointTorque = 0.0f;
    maxContactNormalForce = 0.0f;
    avgContactNormalForce = 0.0f;
    maxContactTangentialForce = 0.0f;
    avgContactTangentialForce = 0.0f;
    avgContactSlip = 0.0f;
    stickingContacts = 0;

    int jointCount = 0;
    int contactCount = 0;

    for (Force *force = forces; force != 0; force = force->next)
    {
        Joint *joint = dynamic_cast<Joint *>(force);
        if (joint && !joint->broken)
        {
            float linResidual = joint->constraintResidualLin();
            float angResidual = joint->constraintResidualAng();

            float residual = max(linResidual, angResidual);
            maxJointResidual = max(maxJointResidual, residual);
            avgJointResidual += residual;

            Joint::JointWrench wrench = joint->computeWrench();
            float forceNorm = length(wrench.forceWorld);
            float torqueNorm = length(wrench.torqueWorld);
            maxJointForce = max(maxJointForce, forceNorm);
            avgJointForce += forceNorm;
            maxJointTorque = max(maxJointTorque, torqueNorm);
            avgJointTorque += torqueNorm;
            jointCount++;
        }

        Manifold *manifold = dynamic_cast<Manifold *>(force);
        if (manifold)
        {
            for (int i = 0; i < manifold->numContacts; ++i)
            {
                Manifold::ContactDiagnostics contact = manifold->computeContactDiagnostics(i);
                float penetration = contact.penetration;

                maxPenetration = max(maxPenetration, penetration);
                avgPenetration += penetration;

                maxContactNormalForce = max(maxContactNormalForce, contact.normalForce);
                avgContactNormalForce += contact.normalForce;
                maxContactTangentialForce = max(maxContactTangentialForce, contact.tangentialForce);
                avgContactTangentialForce += contact.tangentialForce;
                avgContactSlip += contact.slip;
                stickingContacts += contact.stick ? 1 : 0;
                contactCount++;
            }
        }
    }

    if (jointCount > 0)
    {
        avgJointResidual /= (float)jointCount;
        avgJointForce /= (float)jointCount;
        avgJointTorque /= (float)jointCount;
    }

    if (contactCount > 0)
    {
        avgPenetration /= (float)contactCount;
        avgContactNormalForce /= (float)contactCount;
        avgContactTangentialForce /= (float)contactCount;
        avgContactSlip /= (float)contactCount;
    }
}

void Solver::collectDiagnostics()
{
    StepDiagnostics stepDiagnostics = {};
    stepDiagnostics.stepIndex = stepIndex;

    for (Rigid *body = bodies; body != 0; body = body->next)
        stepDiagnostics.activeBodies++;

    for (Force *force = forces; force != 0; force = force->next)
        stepDiagnostics.activeForces++;

    stepDiagnostics.activeContacts = countContacts();

    computeConstraintResiduals(
        stepDiagnostics.maxJointResidual,
        stepDiagnostics.avgJointResidual,
        stepDiagnostics.maxPenetration,
        stepDiagnostics.avgPenetration,
        stepDiagnostics.maxJointForce,
        stepDiagnostics.avgJointForce,
        stepDiagnostics.maxJointTorque,
        stepDiagnostics.avgJointTorque,
        stepDiagnostics.maxContactNormalForce,
        stepDiagnostics.avgContactNormalForce,
        stepDiagnostics.maxContactTangentialForce,
        stepDiagnostics.avgContactTangentialForce,
        stepDiagnostics.avgContactSlip,
        stepDiagnostics.stickingContacts);

    float3 linearMomentum = {0, 0, 0};
    float3 angularMomentum = {0, 0, 0};
    computeMomentum(linearMomentum, angularMomentum);
    stepDiagnostics.linearMomentumNorm = length(linearMomentum);
    stepDiagnostics.angularMomentumNorm = length(angularMomentum);

    computeSystemEnergy(stepDiagnostics.kineticEnergy, stepDiagnostics.potentialEnergy);
    stepDiagnostics.totalEnergy = stepDiagnostics.kineticEnergy + stepDiagnostics.potentialEnergy;

    diagnostics.history.push_back(stepDiagnostics);
}

void Solver::step()
{
    // Perform broadphase collision detection
    // This is a naive O(n^2) approach, but it is sufficient for small numbers of bodies in this sample.
    for (Rigid *bodyA = bodies; bodyA != 0; bodyA = bodyA->next)
    {
        for (Rigid *bodyB = bodyA->next; bodyB != 0; bodyB = bodyB->next)
        {
            float3 dp = bodyA->positionLin - bodyB->positionLin;
            float r = bodyA->radius + bodyB->radius;
            if (dot(dp, dp) <= r * r && !bodyA->constrainedTo(bodyB))
                new Manifold(this, bodyA, bodyB);
        }
    }

    // Initialize and warmstart forces
    for (Force *force = forces; force != 0;)
    {
        // Initialization can including caching anything that is constant over the step
        if (!force->initialize())
        {
            // Force has returned false meaning it is inactive, so remove it from the solver
            Force *next = force->next;
            delete force;
            force = next;
        }
        else
            force = force->next;
    }

    // Initialize and warmstart bodies (ie primal variables)
    for (Rigid *body = bodies; body != 0; body = body->next)
    {
        // Compute inertial position (Eq 2)
        body->inertialLin = body->positionLin + body->velocityLin * dt;
        if (body->mass > 0)
            body->inertialLin += float3{0, 0, gravity} * (dt * dt);
        body->inertialAng = body->positionAng + body->velocityAng * dt;

        // Adaptive warmstart (See original VBD paper)
        float3 accel = (body->velocityLin - body->prevVelocityLin) / dt;
        float accelExt = accel.z * sign(gravity);
        float accelWeight = clamp(accelExt / abs(gravity), 0.0f, 1.0f);
        if (!isfinite(accelWeight))
            accelWeight = 0.0f;

        // Save initial position (x-) and compute warmstarted position (See original VBD paper)
        body->initialLin = body->positionLin;
        body->initialAng = body->positionAng;
        if (body->mass > 0)
        {
            body->positionLin = body->positionLin + body->velocityLin * dt + float3{0, 0, gravity} * (accelWeight * dt * dt);
            body->positionAng = body->positionAng + body->velocityAng * dt;
        }
    }

    // Main solver loop
    for (int it = 0; it < iterations; it++)
    {
        // Primal update
        for (Rigid *body = bodies; body != 0; body = body->next)
        {
            // Skip static / kinematic bodies
            if (body->mass <= 0)
                continue;

            // Initialize left and right hand sides of the linear system (Eqs. 5, 6)
            float3x3 MLin = diagonal(body->mass, body->mass, body->mass);
            float3x3 MAng = diagonal(body->moment.x, body->moment.y, body->moment.z);

            float3x3 lhsLin = MLin / (dt * dt);
            float3x3 lhsAng = MAng / (dt * dt);
            float3x3 lhsCross = float3x3{0, 0, 0, 0, 0, 0, 0, 0, 0};

            float3 rhsLin = MLin / (dt * dt) * (body->positionLin - body->inertialLin);
            float3 rhsAng = MAng / (dt * dt) * (body->positionAng - body->inertialAng);

            // Iterate over all forces acting on the body
            for (Force *force = body->forces; force != 0; force = (force->bodyA == body) ? force->nextA : force->nextB)
            {
                // Stamp the force and hessian into the linear system
                force->updatePrimal(body, alpha, lhsLin, lhsAng, lhsCross, rhsLin, rhsAng);
            }

            // Solve the SPD linear system using LDL and apply the update (Eq. 4)
            float3 dxLin, dxAng;
            solve(lhsLin, lhsAng, lhsCross, -rhsLin, -rhsAng, dxLin, dxAng);
            body->positionLin = body->positionLin + dxLin;
            body->positionAng = body->positionAng + dxAng;
        }

        globalConstraintCorrection();

        // Dual update
        for (Force *force = forces; force != 0; force = force->next)
        {
            force->updateDual(alpha);
        }
    }

    // Compute velocities (BDF1) after the final iteration
    for (Rigid* body = bodies; body != 0; body = body->next)
    {
        body->prevVelocityLin = body->velocityLin;
        if (body->mass > 0)
        {
            body->velocityLin = (body->positionLin - body->initialLin) / dt;
            body->velocityAng = (body->positionAng - body->initialAng) / dt;
        }
    }

    if (diagnosticsEnabled)
        collectDiagnostics();

    stepIndex++;
}
