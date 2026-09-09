// Research-only arithmetic, not a production operator or a visual experiment.
// node idea-stage/conductance_star_arithmetic_check.js
"use strict";
const assert = require("node:assert/strict");
const fs = require("node:fs");
const relationEnergyOnly = process.argv.includes("--relation-energy-check");
const relationEnergySelectivityOnly = process.argv.includes("--relation-energy-selectivity-check");
const noWriteLegacy = process.argv.includes("--no-write-legacy");
const selectiveCorrectionStageA = process.argv.includes("--selective-correction-stage-a");
const structureComputationOnly = process.argv.includes("--structure-computation-only") || selectiveCorrectionStageA;
const constructiveCouplingOnly = process.argv.includes("--constructive-coupling-only") || structureComputationOnly;
const associationCreditKillOnly = process.argv.includes("--association-credit-kill-only") || constructiveCouplingOnly;
const visualCreditKillOnly = process.argv.includes("--visual-credit-kill-only") || associationCreditKillOnly;
const noveltyKillOnly = process.argv.includes("--novelty-kill-only") || visualCreditKillOnly;
const repairKillOnly = process.argv.includes("--repair-kill-only") || noveltyKillOnly;
const continuityKillOnly = process.argv.includes("--continuity-kill-only") || repairKillOnly;
const creditKillOnly = process.argv.includes("--credit-kill-only") || continuityKillOnly;
const p = {leak: 0.8, somaLeak: 0.8, axial: 0.4};
const path = [{a: 1, b: 0, t: 1}, {a: 0, b: 1, t: 1}, {a: -1, b: 0, t: 1}, {a: 0, b: -1, t: 1}];
const checks = [];
function check(name, condition) { assert.ok(condition, name); checks.push(name); }
function norm(v) { return Math.max(...v.map(Math.abs)); }
function gap(a, b) { return norm(a.map((x, i) => x - b[i])); }
function rk4(y, h, f) {
  const add = (a, b, c) => a.map((v, i) => v + c * b[i]);
  const k1 = f(y), k2 = f(add(y, k1, h / 2)), k3 = f(add(y, k2, h / 2)), k4 = f(add(y, k3, h));
  return y.map((v, i) => v + h * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) / 6);
}
function integrate(segments, initial, step, rhs) {
  let state = [...initial];
  for (const segment of segments) {
    const count = Math.ceil(segment.t / step), h = segment.t / count;
    for (let j = 0; j < count; j++) state = rk4(state, h, y => rhs(y, segment.a, segment.b));
  }
  return state;
}
const shuntModulation = (a, b) => [b, -b, a, -a];
function teachingCurrent(v, a, b, epsilon) {
  return shuntModulation(a, b).map((g, j) => -epsilon * g * v[j]);
}
function starFlow(v, a, b, epsilon, axial = p.axial, mismatch = 0, tonic = null, axialCubic = null) {
  const lambda = p.leak + axial, somaLambda = p.somaLeak + 4 * axial;
  const drive = [a, -a, -b, b], modulation = shuntModulation(a, b);
  const result = drive.map((input, j) => -(lambda + epsilon * modulation[j] + (j === 0 ? mismatch : 0)) * v[j] + axial * v[4] + input);
  if (tonic !== null) {
    assert.ok(tonic.length === 4 && tonic.every(g => Number.isFinite(g) && g > Math.abs(epsilon) * Math.max(Math.abs(a), Math.abs(b))), "tonic synaptic conductance must dominate modulation");
    result.forEach((_, j) => { result[j] -= (tonic[j] - p.leak / 2) * v[j]; });
  }
  result.push(-somaLambda * v[4] + axial * v.slice(0, 4).reduce((s, x) => s + x, 0));
  if (axialCubic !== null) {
    assert.ok(tonic === null, "choose one plastic carrier, not tonic plus axial");
    assert.ok(axialCubic.length === 4 && axialCubic.every(g => Number.isFinite(g) && g >= 0), "nonlinear axial coefficients must be finite and nonnegative");
    axialCubic.forEach((g, j) => {
      const current = g * (v[j] - v[4]) ** 3;
      result[j] -= current;
      result[4] += current;
    });
  }
  return result;
}

// Shared implicit-cell primitives.  The legacy R/S checks and the relation-
// energy API both use this same five-state free dynamics and arrowhead solve.
// The optional weightContext exists only for the frozen legacy Stage A path;
// the relation-energy path always leaves it null (no history gate).
const sharedConnectionSigns = [1, -1, -1, 1];
const sharedZeroState = [0, 0, 0, 0, 0];
function sharedIncomingWeightData(incoming, Cstar, mode = "full") {
  assert.ok(mode === "full" || mode === "difference" || mode === "common", "unknown selective weight mode");
  assert.ok(Number.isFinite(Cstar) && Cstar > 0, "selective contrast scale must be positive");
  const q = incoming.slice(0, 4).map(value => (value - incoming[4]) ** 2), P = q.reduce((sum, value) => sum + value, 0) / 4;
  const e = q.map(value => (value - P) / Cstar);
  const weights = mode === "full"
    ? e.map((value, j) => 1 + sharedConnectionSigns[j] * value)
    : mode === "difference"
      ? e.map((value, j) => sharedConnectionSigns[j] * value)
      : e.map(() => 1);
  return {q, P, Cstar, e, weights, mode, connectionSigns: [...sharedConnectionSigns]};
}
function sharedResponse(diagonal, edge, rhs) {
  assert.ok(diagonal.length === 5 && edge.length === 4 && rhs.length === 5, "five-state response shape mismatch");
  const denominator = diagonal[4] - edge.reduce((sum, g, j) => sum + g * g / diagonal[j], 0);
  assert.ok(denominator > 0 && diagonal.every(d => d > 0), "response must be positive definite");
  const soma = (rhs[4] + edge.reduce((sum, g, j) => sum + g * rhs[j] / diagonal[j], 0)) / denominator;
  return [...edge.map((g, j) => (rhs[j] + g * soma) / diagonal[j]), soma];
}
function sharedImplicitCell(th, input, carrier = "axial", incoming = sharedZeroState, tolerance = 1e-12, eps = 0, audit = null, weightContext = null, step = 0.4, maxIterations = 40) {
  const h = step;
  assert.ok(carrier === "axial" || carrier === "tonic", "unknown plastic carrier");
  if (audit !== null) audit.forwardCalls += 1;
  const weightData = weightContext === null ? null : sharedIncomingWeightData(incoming, weightContext.Cstar, weightContext.mode);
  const coefficient = th.map((v, j) => (carrier === "axial" ? p.axial : p.leak / 2) * Math.exp(weightData === null ? v : v * weightData.weights[j]));
  if (!coefficient.every(value => Number.isFinite(value) && value >= 0)) throw new Error(`nonlinear coefficient out of domain: theta=${JSON.stringify(th)} weights=${JSON.stringify(weightData === null ? null : weightData.weights)} coefficient=${JSON.stringify(coefficient)}`);
  const tonic = carrier === "tonic" ? coefficient : null, cubic = carrier === "axial" ? coefficient : null;
  const modulation = shuntModulation(input[0], input[1]);
  const leak = [...modulation.map((m, j) => p.leak + eps * m + (tonic ? tonic[j] - p.leak / 2 : 0)), p.somaLeak];
  assert.ok(leak.every(g => g > 0), "strictly dissipative leakage required");
  function system(state) {
    const drop = state.slice(0, 4).map(v => v - state[4]);
    const edge = drop.map((d, j) => h * (p.axial + (cubic ? 3 * cubic[j] * d * d : 0)));
    const diagonal = [...edge.map((g, j) => 1 + h * leak[j] + g), 1 + h * leak[4] + edge.reduce((sum, g) => sum + g, 0)];
    return {drop, edge, diagonal};
  }
  const residual = state => {
    if (audit !== null) audit.residualEvaluations += 1;
    const flow = starFlow(state, input[0], input[1], eps, p.axial, 0, tonic, cubic);
    return state.map((v, j) => v - incoming[j] - h * flow[j]);
  };
  let state = [...incoming], iterations = 0;
  for (; iterations < maxIterations && norm(residual(state)) > tolerance; iterations++) {
    const r = residual(state), sys = system(state), direction = sharedResponse(sys.diagonal, sys.edge, r);
    if (audit !== null) {
      audit.responseCalls += 1;
      audit.newtonResponseCalls += 1;
    }
    let accepted = false;
    for (let backtrack = 0; backtrack < 24; backtrack++) {
      const scale = 2 ** (-backtrack), next = state.map((v, j) => v - scale * direction[j]);
      if (norm(residual(next)) < norm(r)) {
        state = next;
        accepted = true;
        if (audit !== null) audit.backtrackSteps += backtrack;
        break;
      }
    }
    assert.ok(accepted, "implicit cell Newton step failed; no silent fallback");
  }
  if (audit !== null) audit.newtonIterations += iterations;
  assert.ok(norm(residual(state)) <= tolerance, "implicit cell failed to converge");
  const result = {state, ...system(state), coefficient, carrier, iterations, residual: norm(residual(state))};
  return weightData === null ? result : {...result, weights: weightData.weights, weightData};
}

// ---------------------------------------------------------------------------
// Relation-energy API (the current primary arithmetic operator).
// ---------------------------------------------------------------------------
const relationStateDimension = 5, relationCompartmentCount = 4;
const relationZero = () => [0, 0, 0, 0, 0];
const relationDot = (a, b) => a.reduce((sum, value, i) => sum + value * b[i], 0);
const relationMaxAbs = values => Math.max(...values.map(value => Math.abs(value)));
const relationGap = (a, b) => relationMaxAbs(a.map((value, i) => value - b[i]));
const relationAdd = (a, b) => a.map((value, i) => value + b[i]);
const relationSub = (a, b) => a.map((value, i) => value - b[i]);
const relationScale = (a, scale) => a.map(value => scale * value);
const relationFinite = value => Number.isFinite(value);
const relationFiniteVector = (value, length) => Array.isArray(value) && value.length === length && value.every(relationFinite);
const relationEuclidean = values => Math.sqrt(relationDot(values, values));
function relationAssertVector(value, length, name) {
  assert.ok(relationFiniteVector(value, length), `${name} must be a finite length-${length} vector`);
}
function relationAssertTheta(theta) {
  relationAssertVector(theta, relationCompartmentCount, "theta");
}
function relationRecord(record, name, requireY = false) {
  assert.ok(record !== null && typeof record === "object", `${name} record is required`);
  relationAssertVector(record.x, 2, `${name}.x`);
  relationAssertVector(record.incoming, relationStateDimension, `${name}.incoming`);
  if (requireY) relationAssertVector(record.y, 2, `${name}.y`);
  if (record.action !== undefined) assert.ok(typeof record.action === "string" && record.action.length > 0, `${name}.action must be explicit`);
  return record;
}
function relationOptions(options = {}) {
  const h = options.h === undefined ? 0.4 : options.h;
  const epsilon = options.epsilon === undefined ? 0 : options.epsilon;
  const gamma = options.gamma === undefined ? 1 : options.gamma;
  const tau = options.tau === undefined ? 1 : options.tau;
  const tolerance = options.tolerance === undefined ? 1e-12 : options.tolerance;
  const maxIterations = options.maxIterations === undefined ? 80 : options.maxIterations;
  assert.ok(relationFinite(h) && h > 0, "h must be positive and finite");
  assert.ok(relationFinite(epsilon), "epsilon must be finite");
  assert.ok(relationFinite(gamma) && gamma >= 0, "gamma must be nonnegative and finite");
  assert.ok(relationFinite(tau) && tau > 0, "tau must be positive and finite");
  assert.ok(relationFinite(tolerance) && tolerance > 0, "solver tolerance must be positive and finite");
  assert.ok(Number.isInteger(maxIterations) && maxIterations >= 0, "maxIterations must be a nonnegative integer");
  const theta = options.theta === undefined ? [0, 0, 0, 0] : [...options.theta];
  relationAssertTheta(theta);
  assert.ok(p.leak > 0 && p.somaLeak > 0 && p.axial > 0 && [p.leak, p.somaLeak, p.axial].every(relationFinite), "base conductances must be positive and finite");
  return {h, epsilon, gamma, tau, tolerance, maxIterations, theta, energyBias: options.energyBias};
}
function relationCoefficients(theta) {
  relationAssertTheta(theta);
  const alpha = theta.map(value => p.axial * Math.exp(value));
  assert.ok(alpha.every(value => relationFinite(value) && value > 0), `alpha=k*exp(theta) must be positive and finite: ${JSON.stringify({theta, alpha})}`);
  return alpha;
}
function relationInput(record) {
  return [record.x[0], -record.x[0], -record.x[1], record.x[1], 0];
}
function relationLeak(record, epsilon) {
  const modulation = shuntModulation(record.x[0], record.x[1]);
  const leak = [...modulation.map((value, j) => p.leak + epsilon * value), p.somaLeak];
  assert.ok(leak.every(value => relationFinite(value) && value > 0), `strictly positive leakage required: ${JSON.stringify({x: record.x, epsilon, leak})}`);
  return leak;
}
function relationStrongConvexityLowerBound(record, options) {
  return Math.min(...relationLeak(record, options.epsilon).map(value => 1 / options.h + value));
}
function relationBias(record, theta, options) {
  if (options.energyBias === undefined) return 0;
  const value = typeof options.energyBias === "function" ? options.energyBias(record, theta) : options.energyBias;
  assert.ok(relationFinite(value), "state-independent energy bias must be finite");
  return value;
}
function relationPotential(record, theta, state, options) {
  relationRecord(record, "record");
  relationAssertVector(state, relationStateDimension, "state");
  const alpha = relationCoefficients(theta), leak = relationLeak(record, options.epsilon), drive = relationInput(record);
  let value = 0.5 * leak.slice(0, 4).reduce((sum, conductance, j) => sum + conductance * state[j] ** 2, 0) + 0.5 * leak[4] * state[4] ** 2 - relationDot(drive, state);
  for (let j = 0; j < 4; j++) {
    const d = state[j] - state[4];
    value += p.axial * d ** 2 / 2 + alpha[j] * d ** 4 / 4;
  }
  return value;
}
function relationConditionalEnergy(record, theta, state, incoming, options) {
  relationRecord(record, "record");
  relationAssertVector(incoming, relationStateDimension, "incoming");
  relationAssertVector(state, relationStateDimension, "state");
  const value = relationDot(relationSub(state, incoming), relationSub(state, incoming)) / (2 * options.h)
    + relationPotential(record, theta, state, options) + relationBias(record, theta, options);
  assert.ok(relationFinite(value), "conditional free energy is not finite");
  return value;
}
function relationGradient(record, theta, state, incoming, options) {
  relationRecord(record, "record");
  relationAssertVector(incoming, relationStateDimension, "incoming");
  relationAssertVector(state, relationStateDimension, "state");
  const alpha = relationCoefficients(theta), leak = relationLeak(record, options.epsilon), drive = relationInput(record);
  const gradient = state.map((value, i) => (value - incoming[i]) / options.h + leak[i] * value - drive[i]);
  for (let j = 0; j < 4; j++) {
    const d = state[j] - state[4], current = p.axial * d + alpha[j] * d ** 3;
    gradient[j] += current;
    gradient[4] -= current;
  }
  assert.ok(relationFiniteVector(gradient, relationStateDimension), "free-energy gradient is not finite");
  return gradient;
}
function relationHessian(record, theta, state, incoming, options) {
  relationRecord(record, "record");
  relationAssertVector(incoming, relationStateDimension, "incoming");
  relationAssertVector(state, relationStateDimension, "state");
  const alpha = relationCoefficients(theta), leak = relationLeak(record, options.epsilon);
  const hessian = Array.from({length: relationStateDimension}, () => Array(relationStateDimension).fill(0));
  for (let j = 0; j < relationStateDimension; j++) hessian[j][j] = 1 / options.h + leak[j];
  for (let j = 0; j < 4; j++) {
    const d = state[j] - state[4], conductance = p.axial + 3 * alpha[j] * d ** 2;
    hessian[j][j] += conductance;
    hessian[j][4] -= conductance;
    hessian[4][j] -= conductance;
    hessian[4][4] += conductance;
  }
  assert.ok(hessian.flat().every(relationFinite), "free-energy Hessian is not finite");
  return hessian;
}
function relationDenseSolve(matrixInput, rhs) {
  const size = rhs.length;
  assert.ok(matrixInput.length === size && matrixInput.every(row => Array.isArray(row) && row.length === size), "dense solve shape mismatch");
  const rows = matrixInput.map((row, i) => [...row, rhs[i]]);
  for (let column = 0; column < size; column++) {
    let pivot = column;
    for (let row = column + 1; row < size; row++) if (Math.abs(rows[row][column]) > Math.abs(rows[pivot][column])) pivot = row;
    if (!relationFinite(rows[pivot][column]) || Math.abs(rows[pivot][column]) < 1e-14) throw new Error("singular relation-energy dense solve");
    [rows[column], rows[pivot]] = [rows[pivot], rows[column]];
    const scale = rows[column][column];
    rows[column] = rows[column].map(value => value / scale);
    for (let row = 0; row < size; row++) if (row !== column) {
      const factor = rows[row][column];
      rows[row] = rows[row].map((value, index) => value - factor * rows[column][index]);
    }
  }
  const solution = rows.map(row => row[size]);
  assert.ok(solution.every(relationFinite), "dense relation-energy solve returned a nonfinite value");
  return solution;
}
function relationPairGradient(query, candidate, theta, state, queryIncoming, candidateIncoming, options) {
  const queryGradient = relationGradient(query, theta, state.slice(0, 5), queryIncoming, options);
  const candidateGradient = relationGradient(candidate, theta, state.slice(5), candidateIncoming, options);
  for (let j = 0; j < 4; j++) {
    const difference = state[j] - state[5 + j];
    queryGradient[j] += options.gamma * difference;
    candidateGradient[j] -= options.gamma * difference;
  }
  return [...queryGradient, ...candidateGradient];
}
function relationPairHessian(query, candidate, theta, state, queryIncoming, candidateIncoming, options) {
  const queryHessian = relationHessian(query, theta, state.slice(0, 5), queryIncoming, options);
  const candidateHessian = relationHessian(candidate, theta, state.slice(5), candidateIncoming, options);
  const hessian = Array.from({length: 10}, () => Array(10).fill(0));
  for (let i = 0; i < 5; i++) for (let j = 0; j < 5; j++) {
    hessian[i][j] = queryHessian[i][j];
    hessian[5 + i][5 + j] = candidateHessian[i][j];
  }
  for (let j = 0; j < 4; j++) {
    hessian[j][j] += options.gamma;
    hessian[5 + j][5 + j] += options.gamma;
    hessian[j][5 + j] -= options.gamma;
    hessian[5 + j][j] -= options.gamma;
  }
  assert.ok(hessian.flat().every(relationFinite), "paired Hessian is not finite");
  return hessian;
}
function relationPairEnergy(query, candidate, theta, state, queryIncoming, candidateIncoming, options) {
  const queryState = state.slice(0, 5), candidateState = state.slice(5);
  const portDifference = relationSub(queryState.slice(0, 4), candidateState.slice(0, 4));
  const value = relationConditionalEnergy(query, theta, queryState, queryIncoming, options)
    + relationConditionalEnergy(candidate, theta, candidateState, candidateIncoming, options)
    + options.gamma * relationDot(portDifference, portDifference) / 2;
  assert.ok(relationFinite(value), "paired energy is not finite");
  return {value, portDifference, couplingCost: options.gamma * relationDot(portDifference, portDifference) / 2};
}
function relationFreeToken(record, theta, options = {}) {
  const resolved = relationOptions({...options, theta});
  relationRecord(record, "free");
  const audit = options.audit || null;
  const cellTolerance = resolved.tolerance * resolved.h / Math.sqrt(relationStateDimension);
  const cell = sharedImplicitCell(theta, record.x, "axial", record.incoming, cellTolerance, resolved.epsilon, audit, null, resolved.h, resolved.maxIterations);
  const alpha = relationCoefficients(theta);
  assert.ok(relationGap(cell.coefficient, alpha) <= 1e-14, "shared free dynamics did not realize alpha=k*exp(theta)");
  const energy = relationConditionalEnergy(record, theta, cell.state, record.incoming, resolved);
  const energyGradient = relationGradient(record, theta, cell.state, record.incoming, resolved);
  const gradientResidual = relationEuclidean(energyGradient), m = relationStrongConvexityLowerBound(record, resolved), stateErrorUpperBound = gradientResidual / m, energyErrorUpperBound = gradientResidual ** 2 / (2 * m);
  assert.ok(gradientResidual <= resolved.tolerance, `free state energy-gradient residual exceeds tolerance: ${gradientResidual}`);
  return {record: relationRecordCopy(record), theta: [...theta], incoming: [...record.incoming], state: [...cell.state], alpha: [...alpha], energy, energyGradient: [...energyGradient], stationarityResidual: gradientResidual, gradientResidual, strongConvexityLowerBound: m, stateErrorUpperBound, energyErrorUpperBound, residual: cell.residual, iterations: cell.iterations, audit, epsilon: resolved.epsilon, h: resolved.h};
}
function relationPairTokenWithFree(query, candidate, theta, resolved, queryFree, candidateFree) {
  relationRecord(query, "query");
  relationRecord(candidate, "candidate");
  if (resolved.gamma === 0) {
    const state = [...queryFree.state, ...candidateFree.state];
    const gradient = relationPairGradient(query, candidate, theta, state, query.incoming, candidate.incoming, resolved), residual = relationMaxAbs(gradient), gradientResidual = relationEuclidean(gradient), m = Math.min(queryFree.strongConvexityLowerBound, candidateFree.strongConvexityLowerBound);
    return {query: relationRecordCopy(query), candidate: relationRecordCopy(candidate), state, queryState: [...queryFree.state], candidateState: [...candidateFree.state], energy: queryFree.energy + candidateFree.energy, baseEnergy: queryFree.energy + candidateFree.energy, portDifference: relationSub(queryFree.state.slice(0, 4), candidateFree.state.slice(0, 4)), couplingCost: 0, residual, gradientResidual, strongConvexityLowerBound: m, stateErrorUpperBound: gradientResidual / m, energyErrorUpperBound: gradientResidual ** 2 / (2 * m), iterations: 0, dimension: {query: 5, candidate: 5, total: 10}, unique: true, gamma: 0, queryFree, candidateFree};
  }
  let state = [...queryFree.state, ...candidateFree.state], iterations = 0;
  for (; iterations < resolved.maxIterations; iterations++) {
    const gradient = relationPairGradient(query, candidate, theta, state, query.incoming, candidate.incoming, resolved), residual = relationMaxAbs(gradient), gradientResidual = relationEuclidean(gradient);
    if (gradientResidual <= resolved.tolerance) break;
    const direction = relationDenseSolve(relationPairHessian(query, candidate, theta, state, query.incoming, candidate.incoming, resolved), gradient);
    const current = relationPairEnergy(query, candidate, theta, state, query.incoming, candidate.incoming, resolved);
    let accepted = false;
    for (let backtrack = 0; backtrack < 32; backtrack++) {
      const scale = 2 ** (-backtrack), next = state.map((value, i) => value - scale * direction[i]);
      const nextEnergy = relationPairEnergy(query, candidate, theta, next, query.incoming, candidate.incoming, resolved), nextGradient = relationPairGradient(query, candidate, theta, next, query.incoming, candidate.incoming, resolved), nextResidual = relationMaxAbs(nextGradient), nextGradientResidual = relationEuclidean(nextGradient);
      const energyRoundingTolerance = 64 * Number.EPSILON * Math.max(1, Math.abs(current.value), Math.abs(nextEnergy.value));
      if (nextGradientResidual < gradientResidual && nextResidual < residual && (nextEnergy.value < current.value || nextEnergy.value <= current.value + energyRoundingTolerance)) {
        state = next;
        accepted = true;
        break;
      }
    }
    assert.ok(accepted, "paired Newton step failed; no silent fallback");
  }
  const gradient = relationPairGradient(query, candidate, theta, state, query.incoming, candidate.incoming, resolved), residual = relationMaxAbs(gradient), gradientResidual = relationEuclidean(gradient), m = Math.min(relationStrongConvexityLowerBound(query, resolved), relationStrongConvexityLowerBound(candidate, resolved));
  assert.ok(gradientResidual <= resolved.tolerance, `paired state energy-gradient residual exceeds tolerance: ${gradientResidual}`);
  const pairEnergy = relationPairEnergy(query, candidate, theta, state, query.incoming, candidate.incoming, resolved);
  assert.ok(pairEnergy.couplingCost >= 0 && relationFinite(pairEnergy.couplingCost), "paired coupling cost must be nonnegative and finite");
  return {query: relationRecordCopy(query), candidate: relationRecordCopy(candidate), state: [...state], queryState: state.slice(0, 5), candidateState: state.slice(5), energy: pairEnergy.value, baseEnergy: pairEnergy.value - relationBias(query, theta, resolved) - relationBias(candidate, theta, resolved), portDifference: [...pairEnergy.portDifference], couplingCost: pairEnergy.couplingCost, residual, gradientResidual, strongConvexityLowerBound: m, stateErrorUpperBound: gradientResidual / m, energyErrorUpperBound: gradientResidual ** 2 / (2 * m), iterations, dimension: {query: 5, candidate: 5, total: 10}, unique: true, gamma: resolved.gamma, queryFree, candidateFree};
}
function relationPairToken(query, candidate, theta, options = {}) {
  assert.ok(!Object.prototype.hasOwnProperty.call(options, "queryFree") && !Object.prototype.hasOwnProperty.call(options, "candidateFree"), "public pairToken does not accept injected free-state caches");
  const resolved = relationOptions({...options, theta});
  relationRecord(query, "query");
  relationRecord(candidate, "candidate");
  return relationPairTokenWithFree(query, candidate, theta, resolved, relationFreeToken(query, theta, resolved), relationFreeToken(candidate, theta, resolved));
}
function relationLocalTerms(pair, queryFree, candidateFree, alpha) {
  return alpha.map((coefficient, j) => {
    const pairedQuery = coefficient * (pair.queryState[j] - pair.queryState[4]) ** 4 / 4;
    const pairedCandidate = coefficient * (pair.candidateState[j] - pair.candidateState[4]) ** 4 / 4;
    const freeQuery = coefficient * (queryFree.state[j] - queryFree.state[4]) ** 4 / 4;
    const freeCandidate = coefficient * (candidateFree.state[j] - candidateFree.state[4]) ** 4 / 4;
    return {pairedQuery, pairedCandidate, freeQuery, freeCandidate, value: pairedQuery + pairedCandidate - freeQuery - freeCandidate};
  });
}
function relationEvaluateAssociation(query, candidate, theta, options, queryFree = null) {
  const resolved = relationOptions({...options, theta}), freeQuery = queryFree || relationFreeToken(query, theta, resolved), freeCandidate = relationFreeToken(candidate, theta, resolved), pair = relationPairTokenWithFree(query, candidate, theta, resolved, freeQuery, freeCandidate), alpha = relationCoefficients(theta), localTerms = relationLocalTerms(pair, freeQuery, freeCandidate, alpha), value = pair.energy - freeQuery.energy - freeCandidate.energy, e = localTerms.map(term => term.value);
  assert.ok(relationFinite(value) && relationFiniteVector(e, 4), "association energy or local credit is not finite");
  assert.ok(value >= -1e-9, `net association energy is negative: ${value}`);
  return {value, energy: value, e, localTerms, pair, freeQuery, freeCandidate, alpha, gamma: resolved.gamma, tau: resolved.tau};
}
function relationAssociationEnergy(query, candidate, theta, options = {}) {
  relationRecord(query, "query");
  relationRecord(candidate, "candidate");
  return relationEvaluateAssociation(query, candidate, theta, options);
}
function relationAssociationMargin(query, candidateA, candidateB, theta, options = {}) {
  relationRecord(query, "query");
  relationRecord(candidateA, "candidateA");
  relationRecord(candidateB, "candidateB");
  const resolved = relationOptions({...options, theta}), freeQuery = relationFreeToken(query, theta, resolved), associationA = relationEvaluateAssociation(query, candidateA, theta, resolved, freeQuery), associationB = relationEvaluateAssociation(query, candidateB, theta, resolved, freeQuery);
  const value = (associationB.value - associationA.value) / resolved.tau, gradient = associationB.e.map((entry, j) => (entry - associationA.e[j]) / resolved.tau);
  assert.ok(relationFinite(value) && relationFiniteVector(gradient, 4), "association margin is not finite");
  return {value, margin: value, gradient, eA: associationA.e, eB: associationB.e, associationA, associationB, queryFree: freeQuery, tau: resolved.tau, gamma: resolved.gamma};
}
function relationArrowheadParts(hessian) {
  assert.ok(Array.isArray(hessian) && hessian.length === 5 && hessian.every(row => Array.isArray(row) && row.length === 5), "arrowhead Hessian shape mismatch");
  const diagonal = hessian.map((row, j) => row[j]), edge = hessian.slice(0, 4).map(row => -row[4]);
  assert.ok(diagonal.every(value => relationFinite(value)) && edge.every(value => relationFinite(value)), "arrowhead Hessian is nonfinite");
  return {diagonal, edge};
}
function relationFreeThetaJacobian(record, theta, free, options) {
  const hessian = relationHessian(record, theta, free.state, record.incoming, options), {diagonal, edge} = relationArrowheadParts(hessian), jacobian = Array.from({length: relationStateDimension}, () => Array(relationCompartmentCount).fill(0));
  for (let parameter = 0; parameter < relationCompartmentCount; parameter++) {
    const d = free.state[parameter] - free.state[4], source = Array(relationStateDimension).fill(0), partialGradient = free.alpha[parameter] * d ** 3;
    source[parameter] = -partialGradient;
    source[4] = partialGradient;
    const derivative = sharedResponse(diagonal, edge, source);
    for (let stateIndex = 0; stateIndex < relationStateDimension; stateIndex++) jacobian[stateIndex][parameter] = derivative[stateIndex];
  }
  assert.ok(jacobian.flat().every(relationFinite), "free-state theta Jacobian is nonfinite");
  return {jacobian, responseSolves: relationCompartmentCount};
}
function relationDistanceFreeBundle(record, theta, options) {
  const free = relationFreeToken(record, theta, options), derivative = relationFreeThetaJacobian(record, theta, free, options);
  return {free, jacobian: derivative.jacobian, responseSolves: derivative.responseSolves};
}
function relationDistanceAssociationFromBundles(queryBundle, candidateBundle, options) {
  const portDifference = relationSub(queryBundle.free.state.slice(0, 4), candidateBundle.free.state.slice(0, 4)), value = options.gamma * relationDot(portDifference, portDifference) / 2, gradient = Array(relationCompartmentCount).fill(0);
  for (let parameter = 0; parameter < relationCompartmentCount; parameter++) {
    const portDerivative = queryBundle.jacobian.slice(0, 4).map((row, j) => row[parameter] - candidateBundle.jacobian[j][parameter]);
    gradient[parameter] = options.gamma * relationDot(portDifference, portDerivative);
  }
  assert.ok(relationFinite(value) && relationFiniteVector(gradient, relationCompartmentCount), "distance association is nonfinite");
  return {value, gradient, portDifference};
}
function relationDistanceMargin(query, candidateA, candidateB, theta, options = {}) {
  const resolved = relationOptions({...options, theta}), queryBundle = relationDistanceFreeBundle(query, theta, resolved), candidateABundle = relationDistanceFreeBundle(candidateA, theta, resolved), candidateBBundle = relationDistanceFreeBundle(candidateB, theta, resolved), associationA = relationDistanceAssociationFromBundles(queryBundle, candidateABundle, resolved), associationB = relationDistanceAssociationFromBundles(queryBundle, candidateBBundle, resolved), value = (associationB.value - associationA.value) / resolved.tau, gradient = associationB.gradient.map((entry, j) => (entry - associationA.gradient[j]) / resolved.tau), newtonIterations = queryBundle.free.iterations + candidateABundle.free.iterations + candidateBBundle.free.iterations;
  assert.ok(relationFinite(value) && relationFiniteVector(gradient, relationCompartmentCount), "distance margin is nonfinite");
  return {value, margin: value, gradient, associationA, associationB, queryFree: queryBundle.free, solveSummary: {freeSolves: 3, pairSolves: 0, linearResponseSolves: queryBundle.responseSolves + candidateABundle.responseSolves + candidateBBundle.responseSolves, totalSolverCalls: 3, newtonIterations}};
}
function relationEnergyMarginValueOnly(query, candidateA, candidateB, theta, options = {}) {
  const resolved = relationOptions({...options, theta}), queryFree = relationFreeToken(query, theta, resolved), candidateAFree = relationFreeToken(candidateA, theta, resolved), candidateBFree = relationFreeToken(candidateB, theta, resolved), pairA = relationPairTokenWithFree(query, candidateA, theta, resolved, queryFree, candidateAFree), pairB = relationPairTokenWithFree(query, candidateB, theta, resolved, queryFree, candidateBFree), valueA = pairA.energy - queryFree.energy - candidateAFree.energy, valueB = pairB.energy - queryFree.energy - candidateBFree.energy, value = (valueB - valueA) / resolved.tau;
  assert.ok(relationFinite(value), "value-only relation margin is nonfinite");
  return {value, solveSummary: {freeSolves: 3, pairSolves: 2, linearResponseSolves: 0, totalSolverCalls: 5, newtonIterations: queryFree.iterations + candidateAFree.iterations + candidateBFree.iterations + pairA.iterations + pairB.iterations}};
}
function relationDistanceMarginValueOnly(query, candidateA, candidateB, theta, options = {}) {
  const resolved = relationOptions({...options, theta}), queryFree = relationFreeToken(query, theta, resolved), candidateAFree = relationFreeToken(candidateA, theta, resolved), candidateBFree = relationFreeToken(candidateB, theta, resolved), portA = relationSub(queryFree.state.slice(0, 4), candidateAFree.state.slice(0, 4)), portB = relationSub(queryFree.state.slice(0, 4), candidateBFree.state.slice(0, 4)), valueA = resolved.gamma * relationDot(portA, portA) / 2, valueB = resolved.gamma * relationDot(portB, portB) / 2, value = (valueB - valueA) / resolved.tau;
  assert.ok(relationFinite(value), "value-only distance margin is nonfinite");
  return {value, solveSummary: {freeSolves: 3, pairSolves: 0, linearResponseSolves: 0, totalSolverCalls: 3, newtonIterations: queryFree.iterations + candidateAFree.iterations + candidateBFree.iterations}};
}
function relationLogSumExp(values) {
  assert.ok(values.length > 0 && values.every(relationFinite), "logSumExp inputs must be finite");
  const maximum = Math.max(...values), sum = values.reduce((total, value) => total + Math.exp(value - maximum), 0);
  const result = maximum + Math.log(sum);
  assert.ok(relationFinite(result), "logSumExp returned a nonfinite value");
  return result;
}
function relationSigmoid(value) {
  assert.ok(relationFinite(value), "sigmoid input must be finite");
  return value >= 0 ? 1 / (1 + Math.exp(-value)) : Math.exp(value) / (1 + Math.exp(value));
}
function relationLogSigmoid(value) {
  assert.ok(relationFinite(value), "log-sigmoid input must be finite");
  return value >= 0 ? -Math.log1p(Math.exp(-value)) : value - Math.log1p(Math.exp(value));
}
function relationSigmoidDelta(value, increment) {
  assert.ok(relationFinite(value) && relationFinite(increment), "sigmoid delta inputs must be finite");
  if (increment === 0) return 0;
  const shifted = value + increment, logBefore = relationLogSigmoid(value), logAfter = relationLogSigmoid(shifted);
  if (increment > 0) return Math.exp(logAfter) * (-Math.expm1(logBefore - logAfter));
  return Math.exp(logBefore) * Math.expm1(logAfter - logBefore);
}
function relationLogLikelihood(y, mean, variance) {
  const residual = relationSub(y, mean), squared = relationDot(residual, residual);
  assert.ok(relationFinite(squared), "teacher residual square is not finite");
  const value = -0.5 * (y.length * Math.log(2 * Math.PI * variance) + squared / variance);
  assert.ok(relationFinite(value), "teacher log likelihood is not finite");
  return {residual, squared, logLikelihood: value};
}
function relationFreezeTeacherSpec(teacherSpec, requireGaussianOutside = false) {
  assert.ok(teacherSpec !== null && typeof teacherSpec === "object", "coverage teacher specification is required");
  assert.ok(typeof teacherSpec.version === "string" && teacherSpec.version.length > 0, "teacher version must be a nonempty string");
  const pi0 = teacherSpec.pi0, variance = teacherSpec.variance === undefined ? 1 : teacherSpec.variance;
  assert.ok(relationFinite(pi0) && pi0 > 0 && pi0 < 1, "pi0 must be in (0,1)");
  assert.ok(relationFinite(variance) && variance > 0, "teacher variance must be positive and finite");
  let emptyMean = teacherSpec.emptyMean;
  let emptyVariance = teacherSpec.emptyVariance;
  if (teacherSpec.outside !== undefined) {
    assert.ok(teacherSpec.outside !== null && typeof teacherSpec.outside === "object", "outside teacher specification must be an object");
    emptyMean = teacherSpec.outside.mean;
    emptyVariance = teacherSpec.outside.variance;
  }
  if (requireGaussianOutside || emptyMean !== undefined || emptyVariance !== undefined) {
    relationAssertVector(emptyMean, 2, "emptyMean");
    assert.ok(relationFinite(emptyVariance) && emptyVariance > 0, "emptyVariance must be positive and finite");
  }
  let logEmptyLikelihood = teacherSpec.logEmptyLikelihood;
  if (logEmptyLikelihood !== undefined) {
    assert.ok(typeof logEmptyLikelihood === "number" && relationFinite(logEmptyLikelihood), "logEmptyLikelihood must be a finite number");
  }
  assert.ok(teacherSpec.emptyLikelihood === undefined, "emptyLikelihood is not an accepted teacher input; use finite logEmptyLikelihood or Gaussian outside");
  if (requireGaussianOutside) assert.ok(logEmptyLikelihood === undefined, "online Gaussian outside model cannot also provide logEmptyLikelihood");
  if (requireGaussianOutside) assert.ok(emptyMean !== undefined && emptyVariance !== undefined, "online pending requires a frozen Gaussian outside model");
  return Object.freeze({version: teacherSpec.version, pi0, variance, emptyMean: emptyMean === undefined ? undefined : Object.freeze([...emptyMean]), emptyVariance, logEmptyLikelihood});
}
function relationCoverageTeacher(pending, arrival, teacherSpec) {
  assert.ok(pending !== null && typeof pending === "object", "pending event is required");
  assert.ok(arrival !== null && typeof arrival === "object", "arrival record is required");
  relationAssertVector(arrival.y, 2, "arrival.y");
  const frozen = relationFreezeTeacherSpec(teacherSpec, false), pi0 = frozen.pi0, variance = frozen.variance;
  const w = relationSub(arrival.y, pending.queryRecord.x), likelihoodA = relationLogLikelihood(arrival.y, pending.muA, variance), likelihoodB = relationLogLikelihood(arrival.y, pending.muB, variance);
  let logEmpty;
  if (frozen.logEmptyLikelihood !== undefined) logEmpty = frozen.logEmptyLikelihood;
  else if (frozen.emptyMean !== undefined) logEmpty = relationLogLikelihood(arrival.y, frozen.emptyMean, frozen.emptyVariance).logLikelihood;
  else throw new Error("coverage teacher requires an explicit finite log density or frozen Gaussian outside model");
  assert.ok(relationFinite(logEmpty), "logEmptyLikelihood must be finite");
  const f = pending.margin, pBefore = relationSigmoid(f), lambda = likelihoodA.logLikelihood - likelihoodB.logLikelihood, pPlus = relationSigmoid(f + lambda);
  const logPBefore = relationLogSigmoid(f), logOneMinusPBefore = relationLogSigmoid(-f), logMixture = relationLogSumExp([logPBefore + likelihoodA.logLikelihood, logOneMinusPBefore + likelihoodB.logLikelihood]);
  const logEvidence = relationLogSumExp([Math.log1p(-pi0) + logMixture, Math.log(pi0) + logEmpty]), logCoverageMass = Math.log1p(-pi0) + logMixture - logEvidence, omega = Math.exp(logCoverageMass), delta = omega * relationSigmoidDelta(f, lambda);
  const result = {y: [...arrival.y], w, residualA: likelihoodA.residual, residualB: likelihoodB.residual, squaredResidualA: likelihoodA.squared, squaredResidualB: likelihoodB.squared, logEllA: likelihoodA.logLikelihood, logEllB: likelihoodB.logLikelihood, logEmptyLikelihood: logEmpty, lambda, pBefore, pPlus, logMixture, logEvidence, loss: -logEvidence, omega, delta, derivativeOfLossWrtMargin: -delta, pi0, variance, teacherVersion: frozen.version, deltaNumerics: "stable log-sigmoid difference using expm1"};
  assert.ok(Object.values(result).filter(value => typeof value === "number").every(relationFinite), "coverage teacher returned a nonfinite scalar");
  return result;
}
function relationRecordCopy(record, includeY = true) {
  const copy = {...record, x: [...record.x], incoming: [...record.incoming]};
  if (includeY && record.y !== undefined) copy.y = [...record.y];
  if (!includeY) delete copy.y;
  return copy;
}
function relationTransition(record, query) {
  relationRecord(record, "candidate", true);
  relationRecord(query, "query");
  assert.ok(typeof record.id === "string" && record.id.length > 0, "candidate id must be a nonempty string");
  assert.ok(typeof record.recordVersion === "string" && record.recordVersion.length > 0, "candidate recordVersion must be explicit");
  assert.ok(record.completed === true && Number.isInteger(record.arrival) && record.arrival >= 0, "candidate record must be a completed record with an arrival sequence");
  assert.ok(typeof query.action === "string" && query.action.length > 0, "query action must be explicit and nonempty");
  assert.ok(record.action === query.action, "candidate and query actions must be explicit and equal");
  const increment = relationSub(record.y, record.x);
  return {candidateId: record.id, action: record.action, increment, mu: relationAdd(query.x, increment)};
}
function relationCreatePending(spec) {
  assert.ok(spec !== null && typeof spec === "object", "pending specification is required");
  relationRecord(spec.queryRecord, "queryRecord");
  relationRecord(spec.candidateA, "candidateA", true);
  relationRecord(spec.candidateB, "candidateB", true);
  assert.ok(typeof spec.eventId === "string" && spec.eventId.length > 0, "eventId must be nonempty");
  assert.ok(Number.isInteger(spec.thetaVersion) && spec.thetaVersion >= 0, "thetaVersion must be a nonnegative integer");
  assert.ok(spec.queryRecord.completed === false && spec.queryRecord.arrival === null, "query record must be the current incomplete record with no arrival sequence");
  assert.ok(typeof spec.queryRecord.recordVersion === "string" && spec.queryRecord.recordVersion.length > 0, "query recordVersion must be explicit");
  assert.ok(Number.isInteger(spec.queryRecord.observedAt) && spec.queryRecord.observedAt >= 0, "query observedAt must be a nonnegative integer");
  assert.ok(spec.candidateA.arrival <= spec.queryRecord.observedAt && spec.candidateB.arrival <= spec.queryRecord.observedAt, "candidate completion cannot be after query observedAt");
  assert.ok(spec.candidateA.id !== spec.candidateB.id, "candidate ids must be distinct");
  assert.ok(spec.options === undefined || typeof spec.options.energyBias !== "function", "pending relation does not accept mutable energy-bias callbacks");
  const teacherSpec = relationFreezeTeacherSpec(spec.teacherSpec, true), theta = [...spec.theta];
  relationAssertTheta(theta);
  const options = relationOptions({...spec.options, theta}), computed = relationAssociationMargin(spec.queryRecord, spec.candidateA, spec.candidateB, theta, options), predictionA = relationTransition(spec.candidateA, spec.queryRecord), predictionB = relationTransition(spec.candidateB, spec.queryRecord);
  const frozenQuery = relationFrozenRecordCopy(spec.queryRecord, false), frozenA = relationFrozenRecordCopy(spec.candidateA), frozenB = relationFrozenRecordCopy(spec.candidateB);
  const pending = {
    eventId: spec.eventId,
    thetaVersion: spec.thetaVersion,
    thetaSnapshot: Object.freeze([...theta]),
    action: spec.queryRecord.action,
    teacherVersion: teacherSpec.version,
    teacherSpec,
    candidateIds: Object.freeze([spec.candidateA.id, spec.candidateB.id]),
    candidateRecordVersions: Object.freeze({A: spec.candidateA.recordVersion, B: spec.candidateB.recordVersion}),
    queryRecordVersion: spec.queryRecord.recordVersion,
    historySourceVersions: Object.freeze({query: spec.queryRecord.recordVersion, A: spec.candidateA.recordVersion, B: spec.candidateB.recordVersion}),
    candidateRecords: Object.freeze({A: frozenA, B: frozenB}),
    queryRecord: frozenQuery,
    incoming: Object.freeze([...spec.queryRecord.incoming]),
    muA: Object.freeze([...predictionA.mu]),
    muB: Object.freeze([...predictionB.mu]),
    incrementA: Object.freeze([...predictionA.increment]),
    incrementB: Object.freeze([...predictionB.increment]),
    margin: computed.value,
    marginGradient: Object.freeze([...computed.gradient]),
    eA: Object.freeze([...computed.eA]),
    eB: Object.freeze([...computed.eB]),
    associationSummary: Object.freeze({A: {value: computed.associationA.value, e: Object.freeze([...computed.associationA.e]), localTerms: Object.freeze(computed.associationA.localTerms.map(term => Object.freeze({...term})))}, B: {value: computed.associationB.value, e: Object.freeze([...computed.associationB.e]), localTerms: Object.freeze(computed.associationB.localTerms.map(term => Object.freeze({...term})))}}),
    options: Object.freeze({h: options.h, epsilon: options.epsilon, gamma: options.gamma, tau: options.tau, tolerance: options.tolerance, maxIterations: options.maxIterations, ...(typeof options.energyBias === "number" ? {energyBias: options.energyBias} : {})}),
    consumed: false,
  };
  return relationLockPendingSnapshot(pending);
}
function relationFrozenRecordCopy(record, includeY = true) {
  const copy = relationRecordCopy(record, includeY);
  Object.freeze(copy.x);
  Object.freeze(copy.incoming);
  if (copy.y !== undefined) Object.freeze(copy.y);
  return Object.freeze(copy);
}
function relationLockPendingSnapshot(pending) {
  for (const key of ["eventId", "thetaVersion", "thetaSnapshot", "action", "teacherVersion", "teacherSpec", "candidateIds", "candidateRecordVersions", "queryRecordVersion", "historySourceVersions", "candidateRecords", "queryRecord", "incoming", "muA", "muB", "incrementA", "incrementB", "margin", "marginGradient", "eA", "eB", "associationSummary", "options"]) {
    Object.defineProperty(pending, key, {value: pending[key], enumerable: true, writable: false, configurable: false});
  }
  return pending;
}
function relationCreateOnlineState({theta, thetaVersion = 0, liveState = relationZero()} = {}) {
  relationAssertTheta(theta);
  assert.ok(Number.isInteger(thetaVersion) && thetaVersion >= 0, "thetaVersion must be a nonnegative integer");
  relationAssertVector(liveState, 5, "liveState");
  return {theta: [...theta], thetaVersion, liveState: [...liveState], pending: null, consumedEventIds: []};
}
function relationBeginPending(session, spec) {
  assert.ok(session !== null && typeof session === "object", "online state is required");
  assert.equal(session.pending, null, "only one pending relation event is allowed");
  const pending = relationCreatePending({...spec, theta: session.theta, thetaVersion: session.thetaVersion});
  session.pending = pending;
  return pending;
}
function relationConsumePending(session, arrival, teacherSpec, eta) {
  assert.ok(session !== null && typeof session === "object", "online state is required");
  assert.ok(arrival !== null && typeof arrival === "object" && typeof arrival.eventId === "string", "arrival eventId is required");
  relationAssertVector(arrival.y, 2, "arrival.y");
  assert.ok(Number.isInteger(arrival.arrival) && arrival.arrival >= 0, "query arrival sequence must be explicit and nonnegative");
  assert.ok(typeof arrival.action === "string" && arrival.action.length > 0, "arrival action must be explicit and nonempty");
  assert.ok(typeof arrival.teacherVersion === "string" && arrival.teacherVersion.length > 0, "arrival teacher version must be explicit");
  if (session.consumedEventIds.includes(arrival.eventId)) throw new Error("duplicate arrival already consumed");
  assert.ok(session.pending !== null, "no pending relation event");
  const pending = session.pending;
  assert.equal(arrival.eventId, pending.eventId, "arrival eventId does not match the pending event");
  assert.equal(pending.thetaVersion, session.thetaVersion, "pending theta version mismatch");
  assert.ok(relationGap(session.theta, pending.thetaSnapshot) === 0, "pending theta value snapshot mismatch");
  assert.equal(arrival.action, pending.action, "arrival action does not match the pending action");
  assert.equal(arrival.teacherVersion, pending.teacherVersion, "arrival teacher version mismatch");
  const latestCandidateArrival = Math.max(pending.candidateRecords.A.arrival, pending.candidateRecords.B.arrival);
  assert.ok(arrival.arrival > latestCandidateArrival, "query arrival must follow completed candidate arrivals");
  assert.ok(arrival.arrival > pending.queryRecord.observedAt, "query arrival must follow query observedAt");
  if (teacherSpec !== undefined && teacherSpec !== null) {
    const suppliedTeacher = relationFreezeTeacherSpec(teacherSpec, true);
    assert.deepEqual(suppliedTeacher, pending.teacherSpec, "supplied teacher spec does not match the frozen pending teacher");
  }
  assert.ok(relationFinite(eta) && eta > 0, "eta must be positive and finite");
  const evidence = relationCoverageTeacher(pending, arrival, pending.teacherSpec), deltaTheta = pending.marginGradient.map(value => eta * evidence.delta * value);
  assert.ok(relationFiniteVector(deltaTheta, 4), "parameter write is not finite");
  const thetaBefore = [...session.theta], thetaAfter = relationAdd(session.theta, deltaTheta);
  assert.ok(relationFiniteVector(thetaAfter, 4), "updated theta is not finite");
  relationCoefficients(thetaAfter);
  pending.consumed = true;
  pending.arrivalRecord = {eventId: arrival.eventId, y: [...arrival.y], arrival: arrival.arrival, action: arrival.action, teacherVersion: arrival.teacherVersion};
  pending.evidence = evidence;
  pending.deltaTheta = deltaTheta;
  pending.thetaBefore = thetaBefore;
  pending.thetaAfter = thetaAfter;
  pending.thetaVersionAfter = session.thetaVersion + 1;
  session.theta = thetaAfter;
  session.thetaVersion += 1;
  session.consumedEventIds.push(arrival.eventId);
  session.pending = null;
  return {thetaBefore, thetaAfter, thetaVersionAfter: session.thetaVersion, evidence, deltaTheta, pending};
}

const freeToken = relationFreeToken;
const pairToken = relationPairToken;
const associationEnergy = relationAssociationEnergy;
const associationMargin = relationAssociationMargin;
const coverageTeacher = relationCoverageTeacher;
const createOnlineState = relationCreateOnlineState;
const beginRelationPending = relationBeginPending;
const consumeRelationPending = relationConsumePending;

function runRelationEnergyCheck() {
  const checks = [];
  const recordCheck = (name, passed, details = {}) => {
    const row = {name, passed: Boolean(passed), ...details};
    checks.push(row);
    return row.passed;
  };
  const tryError = operation => {
    try { operation(); return {threw: false, error: null}; }
    catch (error) { return {threw: true, error: String(error.message || error)}; }
  };
  const options = {h: 0.4, epsilon: 0, gamma: 1, tau: 1, tolerance: 1e-12, maxIterations: 80};
  const theta = [0.18, -0.11, 0.07, -0.05];
  const cold = relationZero();
  const warmHistory = (id, prefix) => {
    let incoming = [...cold];
    const steps = [];
    for (const x of prefix) {
      const record = {id: `${id}_${steps.length}`, x, incoming, recordVersion: `${id}-history-v1`};
      const free = freeToken(record, theta, options);
      steps.push({x: [...x], incoming: [...incoming], state: free.state, residual: free.residual});
      incoming = [...free.state];
    }
    return {id, prefix, incoming, steps};
  };
  const histories = {
    query: warmHistory("Q", [[1, 0], [0, 0.8], [-0.6, 0.2]]),
    A: warmHistory("A", [[0, 1], [0.7, -0.4], [0.2, 0.9]]),
    B: warmHistory("B", [[-0.8, 0.3], [0.4, 0.9], [-0.2, -0.7]]),
  };
  const query = {id: "query", x: [0.8, 0.2], action: "forward", arrival: null, completed: false, observedAt: 2, recordVersion: "query-live-v1", incoming: [...histories.query.incoming]};
  const candidateA = {id: "A", x: [0.3, 0.7], y: [0.4, 0.7], action: "forward", arrival: 1, completed: true, recordVersion: "candidate-A-v1", incoming: [...histories.A.incoming]};
  const candidateB = {id: "B", x: [0.7, -0.2], y: [0.6, -0.2], action: "forward", arrival: 2, completed: true, recordVersion: "candidate-B-v1", incoming: [...histories.B.incoming]};
  const freeQuery = freeToken(query, theta, options), expectedAlpha = theta.map(value => p.axial * Math.exp(value));
  recordCheck("Eq1 alpha is exactly k*exp(theta), finite and strictly positive, with no w field", relationGap(freeQuery.alpha, expectedAlpha) <= 1e-14 && freeQuery.alpha.every(value => value > 0 && relationFinite(value)) && !Object.prototype.hasOwnProperty.call(freeQuery, "weights"), {alpha: freeQuery.alpha, theta, k: p.axial});
  recordCheck("Eq2/Eq3 free conditional energy has a finite stationary solution", freeQuery.stationarityResidual <= 1e-9 && freeQuery.residual <= options.tolerance && relationFinite(freeQuery.energy), {energy: freeQuery.energy, stationarityResidual: freeQuery.stationarityResidual, residual: freeQuery.residual, histories});
  const associationA = associationEnergy(query, candidateA, theta, options), associationB = associationEnergy(query, candidateB, theta, options), margin = associationMargin(query, candidateA, candidateB, theta, options);
  const pair = associationA.pair;
  recordCheck("Eq8 A paired operator uses fixed S=[I4 0], two 5-state cells and a converged unique 10-state solve", pair.dimension.total === 10 && pair.dimension.query === 5 && pair.dimension.candidate === 5 && pair.unique && pair.residual <= options.tolerance && pair.portDifference.length === 4 && pair.couplingCost >= 0, {dimension: pair.dimension, port: "S=[I4 0]", gamma: options.gamma, pairResidual: pair.residual, pairGradientResidual: pair.gradientResidual, pairIterations: pair.iterations, portDifference: pair.portDifference, couplingCost: pair.couplingCost, stateErrorUpperBound: pair.stateErrorUpperBound, energyErrorUpperBound: pair.energyErrorUpperBound});
  recordCheck("Eq8 A net association energy is nonnegative for both candidates", associationA.value >= -1e-10 && associationB.value >= -1e-10 && associationA.pair.couplingCost >= 0 && associationB.pair.couplingCost >= 0, {A: associationA.value, B: associationB.value, couplingA: associationA.pair.couplingCost, couplingB: associationB.pair.couplingCost});
  recordCheck("Eq11 e contains all four paired/free query/candidate terms", associationA.localTerms.length === 4 && associationB.localTerms.length === 4 && associationA.localTerms.every(term => [term.pairedQuery, term.pairedCandidate, term.freeQuery, term.freeCandidate, term.value].every(relationFinite)) && relationFiniteVector(associationA.e, 4) && relationFiniteVector(associationB.e, 4), {eA: associationA.e, eB: associationB.e, localTermsA: associationA.localTerms, localTermsB: associationB.localTerms});
  recordCheck("Eq13 g uses f=(A_B-A_A)/tau and g=(e_B-e_A)/tau", Math.abs(margin.value - (associationB.value - associationA.value) / options.tau) <= 1e-14 && relationGap(margin.gradient, associationB.e.map((value, j) => (value - associationA.e[j]) / options.tau)) <= 1e-14, {f: margin.value, A_A: associationA.value, A_B: associationB.value, tau: options.tau, gradient: margin.gradient});
  const fdSteps = [1e-4, 3e-5, 1e-5], finiteDifferences = fdSteps.map(step => {
    const finite = theta.map((_, j) => {
      const plus = [...theta], minus = [...theta]; plus[j] += step; minus[j] -= step;
      return (associationMargin(query, candidateA, candidateB, plus, options).value - associationMargin(query, candidateA, candidateB, minus, options).value) / (2 * step);
    });
    return {step, finite, absoluteGap: relationGap(finite, margin.gradient)};
  });
  const analyticGradientNorm = relationEuclidean(margin.gradient);
  const options05 = {...options, epsilon: 0.05}, margin05 = associationMargin(query, candidateA, candidateB, theta, options05), finite05 = theta.map((_, j) => {
    const plus = [...theta], minus = [...theta]; plus[j] += 1e-5; minus[j] -= 1e-5;
    return (associationMargin(query, candidateA, candidateB, plus, options05).value - associationMargin(query, candidateA, candidateB, minus, options05).value) / 2e-5;
  });
  const finiteDifferencesWithRelative = finiteDifferences.map(row => ({...row, relativeGap: row.absoluteGap / Math.max(analyticGradientNorm, Number.MIN_VALUE)}));
  const finite05Gap = relationGap(finite05, margin05.gradient), finite05RelativeGap = finite05Gap / Math.max(relationEuclidean(margin05.gradient), Number.MIN_VALUE);
  recordCheck("Eq11 e/Eq13 g analytic gradient matches FD with incoming frozen, full re-solves, nonzero norm and epsilon=.05", analyticGradientNorm > 1e-6 && finiteDifferencesWithRelative.every(row => row.relativeGap <= 2e-7) && finite05RelativeGap <= 2e-7 && relationEuclidean(margin05.gradient) > 1e-6, {theta, frozenIncoming: true, analyticGradientNorm, finiteDifferences: finiteDifferencesWithRelative, analytic: margin.gradient, epsilon05: {analytic: margin05.gradient, gradientNorm: relationEuclidean(margin05.gradient), finite: finite05, absoluteGap: finite05Gap, relativeGap: finite05RelativeGap}});
  const swapped = associationEnergy(candidateA, query, theta, options);
  recordCheck("Eq8 A complete-record exchange symmetry holds", Math.abs(swapped.value - associationA.value) <= 2e-10, {forward: associationA.value, swapped: swapped.value, gap: Math.abs(swapped.value - associationA.value)});
  const biasOptions = {...options, energyBias: (record, th) => (record.id === "query" ? 3.1 : record.id === "A" ? -2.4 : 1.7) + relationDot(th, th)};
  const biasedMargin = associationMargin(query, candidateA, candidateB, theta, biasOptions);
  recordCheck("Eq8 A state-independent per-record energy biases cancel from A and its gradient", Math.abs(biasedMargin.value - margin.value) <= 2e-10 && relationGap(biasedMargin.gradient, margin.gradient) <= 2e-9, {unbiased: margin.value, biased: biasedMargin.value, gradientGap: relationGap(biasedMargin.gradient, margin.gradient)});
  const gammaZero = associationMargin(query, candidateA, candidateB, theta, {...options, gamma: 0});
  recordCheck("gamma=0 is a numerical habit control with actual port differences and energy residuals", Math.abs(gammaZero.value) <= 1e-12 && relationGap(gammaZero.associationA.pair.portDifference, relationSub(gammaZero.associationA.pair.queryState.slice(0, 4), gammaZero.associationA.pair.candidateState.slice(0, 4))) <= 1e-15 && relationFinite(gammaZero.associationA.pair.residual) && relationFinite(gammaZero.associationA.pair.gradientResidual) && gammaZero.associationA.pair.couplingCost === 0 && gammaZero.associationB.pair.couplingCost === 0, {margin: gammaZero.value, gradient: gammaZero.gradient, eA: gammaZero.eA, eB: gammaZero.eB, portDifferenceA: gammaZero.associationA.pair.portDifference, residualA: gammaZero.associationA.pair.residual, gradientResidualA: gammaZero.associationA.pair.gradientResidual});
  const pendingSpec = {eventId: "E0", queryRecord: query, candidateA, candidateB, teacherSpec: {version: "teacher-fixture-v1", pi0: 0.1, variance: 1, emptyMean: [2.5, -2.5], emptyVariance: 1}};
  const teacherFixture = {version: "teacher-fixture-v1", pi0: 0.1, variance: 1, emptyMean: [2.5, -2.5], emptyVariance: 1};
  const pending = relationCreatePending({...pendingSpec, theta, thetaVersion: 0, options});
  const arrival = {eventId: "E0", y: [0.9, 0.2], arrival: 3, action: "forward", teacherVersion: "teacher-fixture-v1"};
  const evidence = coverageTeacher(pending, arrival, teacherFixture);
  const teacherLossAtMargin = marginValue => coverageTeacher({...pending, margin: marginValue}, arrival, teacherFixture).loss;
  const teacherFd = [1e-4, 1e-5].map(step => {
    const finite = (teacherLossAtMargin(pending.margin + step) - teacherLossAtMargin(pending.margin - step)) / (2 * step), analytic = evidence.derivativeOfLossWrtMargin;
    return {step, finite, analytic, absoluteGap: Math.abs(finite - analytic), relativeGap: Math.abs(finite - analytic) / Math.max(Math.abs(analytic), Number.MIN_VALUE)};
  });
  const extremeTeacher = [40, 100, -1000].map((value, index) => coverageTeacher({...pending, margin: value}, {...arrival, eventId: `E_extreme_${index}`}, teacherFixture));
  const algebraicTeacher = {version: "teacher-log-algebra-v1", pi0: 0.1, variance: 1, logEmptyLikelihood: -3};
  const algebraicEvidence = coverageTeacher(pending, arrival, algebraicTeacher);
  const frozenTeacherMutation = tryError(() => { pending.teacherSpec.emptyMean[0] = 99; });
  recordCheck("Eq15 Gaussian outside is frozen at begin, dL/df=-delta matches loss FD and extreme margins stay finite", pending.teacherSpec.emptyMean[0] === 2.5 && Object.isFrozen(pending.teacherSpec) && Object.isFrozen(pending.teacherSpec.emptyMean) && frozenTeacherMutation.threw && [evidence.lambda, evidence.logMixture, evidence.logEvidence, evidence.omega, evidence.delta, evidence.derivativeOfLossWrtMargin].every(relationFinite) && Math.abs(evidence.derivativeOfLossWrtMargin + evidence.delta) <= 1e-15 && evidence.omega >= 0 && evidence.omega <= 1 && teacherFd.every(row => row.relativeGap <= 2e-7) && extremeTeacher.every(row => [row.logMixture, row.logEvidence, row.loss, row.delta].every(relationFinite)) && [algebraicEvidence.logEmptyLikelihood, algebraicEvidence.loss].every(relationFinite), {teacherFixture, evidence, teacherFd, extremeMargins: extremeTeacher.map(row => ({pBefore: row.pBefore, logMixture: row.logMixture, logEvidence: row.logEvidence, delta: row.delta})), algebraicLogDensity: algebraicEvidence.logEmptyLikelihood, frozenTeacherMutation});
  const equalQuery = {...query, id: "equal-query", x: [0.5, 0.25], recordVersion: "equal-query-v1"};
  const equalCandidateA = {...candidateA, id: "equal-A", x: [0.25, 0.25], y: [0.5, 0.25], recordVersion: "equal-A-v1"};
  const equalCandidateB = {...candidateB, id: "equal-B", x: [0.75, 0.25], y: [0.5, 0.25], recordVersion: "equal-B-v1"};
  const equalArrival = {eventId: "E_equal", y: [0.5, 0.25], arrival: 4, action: "forward", teacherVersion: "teacher-fixture-v1"};
  const equalSpec = {eventId: "E_equal", queryRecord: equalQuery, candidateA: equalCandidateA, candidateB: equalCandidateB, teacherSpec: {version: "teacher-fixture-v1", pi0: 0.1, variance: 1, emptyMean: [2.5, -2.5], emptyVariance: 1}};
  const equalPending = relationCreatePending({...equalSpec, theta, thetaVersion: 0, options});
  const equalEvidence = coverageTeacher(equalPending, equalArrival, teacherFixture);
  recordCheck("Eq15 equal likelihood produces exactly zero teacher delta", Math.abs(equalEvidence.lambda) <= 1e-14 && Math.abs(equalEvidence.delta) <= 1e-14, {equalEvidence});
  const outsideEvidence = coverageTeacher(pending, {eventId: "E_outside", y: [3, 3], arrival: 5}, {version: "teacher-outside-far-v1", pi0: 0.1, variance: 1, emptyMean: [3, 3], emptyVariance: 1});
  const noCoverageEvidence = coverageTeacher(pending, {eventId: "E_outside_no_coverage", y: [3, 3], arrival: 5}, {version: "teacher-outside-rare-v1", pi0: 1e-8, variance: 1, emptyMean: [3, 3], emptyVariance: 1});
  recordCheck("Eq15 candidate-set-outside Gaussian explanation suppresses rather than forces a large write", Math.abs(outsideEvidence.delta) < Math.abs(noCoverageEvidence.delta) && outsideEvidence.omega < noCoverageEvidence.omega, {outside: outsideEvidence, noCoverageReference: noCoverageEvidence});
  const session = createOnlineState({theta, thetaVersion: 0, liveState: [0.21, -0.17, 0.08, -0.03, 0.04]});
  const historyBefore = JSON.stringify({query, candidateA, candidateB}), liveBefore = [...session.liveState];
  const onlinePending = beginRelationPending(session, pendingSpec);
  const secondPending = tryError(() => beginRelationPending(session, {...pendingSpec, eventId: "E_second"}));
  const mismatchArrival = tryError(() => consumeRelationPending(session, {eventId: "wrong", y: arrival.y, arrival: 3, action: "forward", teacherVersion: "teacher-fixture-v1"}, teacherFixture, 0.02));
  const onlineUpdate = consumeRelationPending(session, arrival, teacherFixture, 0.02);
  const duplicateArrival = tryError(() => consumeRelationPending(session, arrival, teacherFixture, 0.02));
  const historyAfter = JSON.stringify({query, candidateA, candidateB});
  recordCheck("Eq17 single pending, event identity, action and teacher version are enforced", secondPending.threw && mismatchArrival.threw && onlineUpdate.thetaVersionAfter === 1 && session.pending === null && session.consumedEventIds.length === 1, {secondPending, mismatchArrival, thetaBefore: onlineUpdate.thetaBefore, thetaAfter: onlineUpdate.thetaAfter, version: session.thetaVersion});
  recordCheck("Eq17 duplicate arrival is rejected and history/live probes remain isolated", duplicateArrival.threw && historyBefore === historyAfter && relationGap(liveBefore, session.liveState) === 0 && relationFiniteVector(onlineUpdate.deltaTheta, 4), {duplicateArrival, historyUnchanged: historyBefore === historyAfter, liveUnchanged: relationGap(liveBefore, session.liveState) === 0, deltaTheta: onlineUpdate.deltaTheta});
  const staleSession = createOnlineState({theta, thetaVersion: 0, liveState: liveBefore});
  beginRelationPending(staleSession, {...pendingSpec, eventId: "E_stale"});
  staleSession.thetaVersion = 1;
  const staleVersion = tryError(() => consumeRelationPending(staleSession, {eventId: "E_stale", y: arrival.y, arrival: 3, action: "forward", teacherVersion: "teacher-fixture-v1"}, teacherFixture, 0.02));
  const staleValueSession = createOnlineState({theta, thetaVersion: 0, liveState: liveBefore});
  beginRelationPending(staleValueSession, {...pendingSpec, eventId: "E_stale_value"});
  staleValueSession.theta[0] += 1e-6;
  const staleValue = tryError(() => consumeRelationPending(staleValueSession, {eventId: "E_stale_value", y: arrival.y, arrival: 3, action: "forward", teacherVersion: "teacher-fixture-v1"}, teacherFixture, 0.02));
  const equalSession = createOnlineState({theta, thetaVersion: 0, liveState: liveBefore});
  beginRelationPending(equalSession, {...equalSpec, eventId: "E_equal_write"});
  const equalUpdate = consumeRelationPending(equalSession, {...equalArrival, eventId: "E_equal_write"}, teacherFixture, 0.02);
  const chronologySession = createOnlineState({theta, thetaVersion: 0, liveState: liveBefore});
  beginRelationPending(chronologySession, {...pendingSpec, eventId: "E_chronology"});
  const chronology = tryError(() => consumeRelationPending(chronologySession, {eventId: "E_chronology", y: arrival.y, arrival: 2, action: "forward", teacherVersion: "teacher-fixture-v1"}, teacherFixture, 0.02));
  recordCheck("Eq17 theta version/value snapshots, candidate chronology and equal-likelihood zero write are enforced", staleVersion.threw && staleValue.threw && chronology.threw && equalUpdate.deltaTheta.every(value => value === 0) && relationGap(equalUpdate.thetaAfter, theta) === 0, {staleVersion, staleValue, chronology, equalUpdate: {evidence: equalUpdate.evidence, deltaTheta: equalUpdate.deltaTheta}});
  const invalidDomain = {
    nanTheta: tryError(() => freeToken(query, [Number.NaN, 0, 0, 0], options)),
    overflowAlpha: tryError(() => freeToken(query, [1000, 0, 0, 0], options)),
    negativeLeak: tryError(() => freeToken({...query, x: [1, 0]}, theta, {...options, epsilon: 2})),
    negativeGamma: tryError(() => pairToken(query, candidateA, theta, {...options, gamma: -1})),
    zeroTau: tryError(() => associationMargin(query, candidateA, candidateB, theta, {...options, tau: 0})),
    nonconvergence: tryError(() => freeToken(query, theta, {...options, maxIterations: 0})),
    conflictingTeacher: tryError(() => relationCreatePending({...pendingSpec, eventId: "E_conflicting_teacher", theta, thetaVersion: 0, options, teacherSpec: {...teacherFixture, logEmptyLikelihood: -3}})),
    missingArrivalSequence: tryError(() => consumeRelationPending(createOnlineState({theta, thetaVersion: 0, liveState: liveBefore}), {eventId: "E_missing_sequence", y: arrival.y, action: "forward", teacherVersion: "teacher-fixture-v1"}, teacherFixture, 0.02)),
    futureCandidate: tryError(() => relationCreatePending({...pendingSpec, eventId: "E_future_candidate", candidateA: {...candidateA, arrival: 3, recordVersion: "candidate-A-future-v1"}, theta, thetaVersion: 0, options})),
  };
  recordCheck("domain, positive-leakage, alpha, tau/gamma and nonconvergence failures are hard rejects", Object.values(invalidDomain).every(row => row.threw), invalidDomain);
  const finiteEverything = [freeQuery.state, pair.state, margin.gradient, onlineUpdate.thetaAfter].flat().every(relationFinite);
  recordCheck("all primary path states, energies, gradients and writes are finite with no silent clipping", finiteEverything, {finiteEverything});
  const passed = checks.every(row => row.passed);
  return {
    scope: "Relation-energy T1 arithmetic only; no training, simulation, RGB/HM3D calibration, visual claim or novelty claim.",
    entryPoint: "--relation-energy-check",
    primary: true,
    legacyGate: {status: "historical regression only", currentRecommended: false, oldResultJsonWritten: false, oldResultJsonPaths: ["idea-stage/SELECTIVE_CORRECTION_CHECK_20260909.json", "idea-stage/STRUCTURE_COMPUTATION_CHECK_20260909.json"]},
    equations: {
      Eq1: "I_j=k d_j+alpha_j d_j^3; alpha_j=k exp(theta_j)>0; no w",
      Eq2: "Phi=1/2 sum L_j z_j^2+1/2 ell_s s^2+sum(k d_j^2/2+alpha_j d_j^4/4)",
      Eq3: "U=(2h)^(-1)||Z-Zminus||^2+Phi-v^T Z",
      Eq4: "free state is argmin U",
      Eq8: "A_C=E*_C-nu_q-nu_C, where E*_C=min_{Z,W}{U_q(Z)+U_C(W)+gamma/2||S Z-S W||^2}",
      Eq11: "e_C,j=B_j(Z_q|C)+B_j(Z_C|q)-B_j(Zbar_q)-B_j(Zbar_C), B_j=alpha_j d_j^4/4",
      Eq13: "g=(e_B-e_A)/tau and f=(A_B-A_A)/tau",
      Eq15: "M=pi P_A+(1-pi)P_B; Zev=(1-pi0)M+pi0 P_empty; dL/df=-delta in log-domain",
      Eq17: "Delta theta=eta*delta*g after one pending/action/version/chronology-checked arrival",
    },
    api: {
      freeToken: "freeToken(record, theta, options) -> {state, energy, alpha, stationarityResidual, ...}",
      pairToken: "pairToken(query, candidate, theta, options) -> {queryState, candidateState, energy, couplingCost, ...}",
      associationEnergy: "associationEnergy(query, candidate, theta, options) -> {value, e, localTerms, pair, ...}",
      associationMargin: "associationMargin(query, candidateA, candidateB, theta, options) -> {value, gradient, eA, eB, ...}",
      coverageTeacher: "coverageTeacher(pending, arrival, {version, pi0, variance, emptyMean, emptyVariance}) or finite logEmptyLikelihood for algebraic checks only -> {lambda, omega, delta, ...}",
      online: "createOnlineState -> beginRelationPending(session, spec with frozen Gaussian teacher) -> consumeRelationPending(session, arrival, optional matching teacher, eta); one pending and duplicate/version/value/chronology checks are enforced",
    },
    parameters: {...options, theta, k: p.axial, leak: p.leak, somaLeak: p.somaLeak, compartmentCount: 4, stateDimensionPerRecord: 5, pairStateDimension: 10, pi0: "explicit teacher fixture only: 0.1"},
    fixture: {query, candidateA, candidateB, histories, theta, evidenceArrival: arrival, equalArrival, teacherFixture, informationBoundary: {allowed: ["x", "actual action", "actual incoming", "completed candidate x/y records", "arrived y"], forbidden: ["pose", "depth", "overlap", "place identity", "reward", "success", "candidate truth"]}},
    complexity: {pairedSolve: "dense Newton on a 10x10 paired Hessian", perNewtonLinearSolve: "O((2(H+1))^3) time and O((2(H+1))^2) memory", nonlinearIterations: "reported per solve", schur: {implemented: false, claim: "no O(H) Schur claim"}},
    onlineSemantics: {kind: "conditional write transaction", continuousLiveDriverImplemented: false, callerMustCommitObservedFreeState: true, liveStateProbeWrites: false, chronology: "begin requires candidate arrivals <= query observedAt; consume requires query arrival > query observedAt and all completed candidates"},
    primaryObservations: {associationA: {value: associationA.value, e: associationA.e, pair: {energy: associationA.pair.energy, couplingCost: associationA.pair.couplingCost, residual: associationA.pair.residual, iterations: associationA.pair.iterations}}, associationB: {value: associationB.value, e: associationB.e, pair: {energy: associationB.pair.energy, couplingCost: associationB.pair.couplingCost, residual: associationB.pair.residual, iterations: associationB.pair.iterations}}, margin: {value: margin.value, gradient: margin.gradient}, evidence, online: {thetaBefore: onlineUpdate.thetaBefore, thetaAfter: onlineUpdate.thetaAfter, deltaTheta: onlineUpdate.deltaTheta, thetaVersionAfter: onlineUpdate.thetaVersionAfter}},
    checksPassed: checks.filter(row => row.passed).length,
    checksTotal: checks.length,
    checks,
    passed,
    verdict: passed ? "T1 relation-energy arithmetic and causal interface checks passed; this remains a conditional synthetic result." : "T1 failed; do not advance to T2 or visual claims.",
  };
}
function relationSelectivityClock() {
  return Number(process.hrtime.bigint()) / 1e6;
}
function relationSelectivityRelationSolveSummary(margin) {
  const queryFree = margin.queryFree, candidateAFree = margin.associationA.freeCandidate, candidateBFree = margin.associationB.freeCandidate, pairA = margin.associationA.pair, pairB = margin.associationB.pair;
  return {freeSolves: 3, pairSolves: 2, linearResponseSolves: 0, totalSolverCalls: 5, newtonIterations: queryFree.iterations + candidateAFree.iterations + candidateBFree.iterations + pairA.iterations + pairB.iterations};
}
function relationSelectivityReadValue(method, context, theta, options) {
  const started = relationSelectivityClock(), read = method === "relationEnergy"
    ? relationEnergyMarginValueOnly(context.queryRecord, context.candidateA, context.candidateB, theta, options)
    : relationDistanceMarginValueOnly(context.queryRecord, context.candidateA, context.candidateB, theta, options);
  const elapsedMs = relationSelectivityClock() - started;
  return {f: read.value, solveSummary: read.solveSummary, elapsedMs};
}
function relationSelectivityEvaluate(method, context, theta, options, y, teacherSpec) {
  const started = relationSelectivityClock(), read = method === "relationEnergy"
    ? relationAssociationMargin(context.queryRecord, context.candidateA, context.candidateB, theta, options)
    : relationDistanceMargin(context.queryRecord, context.candidateA, context.candidateB, theta, options);
  const predictions = context.predictions, teacherPending = {margin: read.value, muA: predictions.A.mu, muB: predictions.B.mu, queryRecord: {x: [...context.queryRecord.x]}}, teacher = relationCoverageTeacher(teacherPending, {y: [...y]}, teacherSpec), solveSummary = method === "relationEnergy" ? relationSelectivityRelationSolveSummary(read) : read.solveSummary, elapsedMs = relationSelectivityClock() - started;
  assert.ok(relationFinite(read.value) && relationFiniteVector(read.gradient, relationCompartmentCount) && relationFinite(teacher.delta), "selectivity evaluation is nonfinite");
  return {method, f: read.value, gradient: [...read.gradient], gradientNorm: relationEuclidean(read.gradient), delta: teacher.delta, sign: Math.sign(teacher.delta), teacher, solveSummary, elapsedMs};
}
function relationSelectivityMatchSelfCorrection({method, context, theta, options, y, teacherSpec, baseEvaluation, target}) {
  const targetTolerance = Math.max(1e-10, 1e-6 * target), gradientNorm = relationEuclidean(baseEvaluation.gradient), directionNorm = Math.abs(baseEvaluation.delta) * gradientNorm, lambdaLimit = 1, started = relationSelectivityClock(), trace = [];
  const common = {method, target, targetTolerance, lambdaLimit, directionNorm, trialUsesGradient: false, trialUsesTeacher: false, etaIsOfflineDiagnostic: true};
  const solveCount = () => trace.reduce((sum, row) => sum + (row.solveSummary == null ? 0 : row.solveSummary.totalSolverCalls), 0), responseSolveCount = () => trace.reduce((sum, row) => sum + (row.solveSummary == null ? 0 : row.solveSummary.linearResponseSolves), 0);
  const unmatched = (reason, lastTrial = null) => ({...common, matched: false, unmatchedReason: reason, eta: null, lambda: null, parameterStepL2: null, uMatched: null, matchError: null, solveCount: solveCount(), linearResponseSolves: responseSolveCount(), evaluations: trace.length, elapsedMs: relationSelectivityClock() - started, trace, lastTrial});
  if (baseEvaluation.delta === 0) return unmatched("delta_zero");
  if (gradientNorm === 0) return unmatched("gradient_zero");
  if (!relationFinite(directionNorm) || directionNorm <= 0) return unmatched("nonfinite_or_zero_search_direction");
  const evaluateLambda = lambda => {
    const eta = lambda / directionNorm, thetaTrial = theta.map((value, j) => value + eta * baseEvaluation.delta * baseEvaluation.gradient[j]);
    try {
      const read = relationSelectivityReadValue(method, context, thetaTrial, options), u = Math.sign(baseEvaluation.delta) * (read.f - baseEvaluation.f);
      if (!relationFinite(u) || !relationFinite(read.f)) return {ok: false, lambda, eta, theta: thetaTrial, error: "nonfinite self correction"};
      return {ok: true, lambda, eta, theta: thetaTrial, f: read.f, u, solveSummary: read.solveSummary, elapsedMs: read.elapsedMs};
    } catch (error) {
      return {ok: false, lambda, eta, theta: thetaTrial, error: String(error.message || error)};
    }
  };
  let lower = {ok: true, lambda: 0, eta: 0, theta: [...theta], f: baseEvaluation.f, u: 0, solveSummary: {totalSolverCalls: 0, linearResponseSolves: 0}, elapsedMs: 0}, upper = null, lambda = 1e-6;
  for (;;) {
    const trial = evaluateLambda(lambda), traceRow = {lambda: trial.lambda, eta: trial.eta, ok: trial.ok, f: trial.ok ? trial.f : null, u: trial.ok ? trial.u : null, solveSummary: trial.solveSummary || null, error: trial.error || null};
    trace.push(traceRow);
    if (!trial.ok) return unmatched("nonfinite_or_domain_during_search", trial);
    if (trial.u >= target) {
      upper = trial;
      break;
    }
    lower = trial;
    if (lambda >= lambdaLimit) return unmatched("search_budget_not_bracketed", trial);
    lambda = Math.min(lambdaLimit, lambda * 2);
  }
  let best = Math.abs(lower.u - target) <= Math.abs(upper.u - target) ? lower : upper;
  for (let iteration = 0; iteration < 80; iteration++) {
    if (Math.abs(best.u - target) <= targetTolerance) break;
    const midpoint = (lower.lambda + upper.lambda) / 2, trial = evaluateLambda(midpoint), traceRow = {lambda: trial.lambda, eta: trial.eta, ok: trial.ok, f: trial.ok ? trial.f : null, u: trial.ok ? trial.u : null, solveSummary: trial.solveSummary || null, error: trial.error || null};
    trace.push(traceRow);
    if (!trial.ok) return unmatched("nonfinite_or_domain_during_search", trial);
    if (Math.abs(trial.u - target) < Math.abs(best.u - target)) best = trial;
    if (trial.u >= target) upper = trial;
    else lower = trial;
  }
  const matchError = Math.abs(best.u - target), matched = matchError <= targetTolerance;
  return {...common, matched, unmatchedReason: matched ? null : "search_bracketed_but_bisection_tolerance_unmet", eta: matched ? best.eta : null, lambda: matched ? best.lambda : null, parameterStepL2: matched ? relationEuclidean(relationSub(best.theta, theta)) : null, uMatched: matched ? best.u : null, matchError: matched ? matchError : null, solveCount: solveCount(), linearResponseSolves: responseSolveCount(), evaluations: trace.length, elapsedMs: relationSelectivityClock() - started, thetaAfter: matched ? [...best.theta] : [...theta], fAfter: matched ? best.f : null, trace};
}
function relationSelectivityGram(gradients) {
  assert.ok(Array.isArray(gradients) && gradients.length === 2 && gradients.every(row => relationFiniteVector(row, relationCompartmentCount)), "selectivity Gram requires two finite length-4 gradients");
  const normG0 = relationEuclidean(gradients[0]), normG1 = relationEuclidean(gradients[1]), K = [[relationDot(gradients[0], gradients[0]), relationDot(gradients[0], gradients[1])], [relationDot(gradients[1], gradients[0]), relationDot(gradients[1], gradients[1])]], c = normG0 > 0 && normG1 > 0 ? K[0][1] / (normG0 * normG1) : null, absK01OverKii = {K00: Math.abs(K[0][1]) / Math.max(K[0][0], Number.MIN_VALUE), K11: Math.abs(K[0][1]) / Math.max(K[1][1], Number.MIN_VALUE)}, normalizedGramMinEigenvalue = c === null ? null : 1 - Math.abs(c), normalizedGramMaxEigenvalue = c === null ? null : 1 + Math.abs(c), Rmax = c === null ? null : Math.abs(c) * Math.max(normG0 / Math.max(normG1, Number.MIN_VALUE), normG1 / Math.max(normG0, Number.MIN_VALUE));
  return {G: gradients.map(row => [...row]), normG: [normG0, normG1], K, c, absK01OverKii, normalizedGramEigenvalues: c === null ? null : [normalizedGramMinEigenvalue, normalizedGramMaxEigenvalue], normalizedGramMinEigenvalue, Rmax, normalizedGramNote: "1-|c| is the normalized Gram minimum eigenvalue, not independent evidence"};
}
function relationSelectivityHistorySnapshot(id, prefix, theta, options) {
  let incoming = relationZero();
  const steps = [];
  for (let index = 0; index < prefix.length; index++) {
    const x = [...prefix[index]], record = {id: `${id}_${index}`, x, incoming: [...incoming], recordVersion: `${id}-theta-nonzero-v1`}, free = relationFreeToken(record, theta, options);
    steps.push({x, incoming: [...incoming], state: [...free.state], residual: free.residual, gradientResidual: free.gradientResidual});
    incoming = [...free.state];
  }
  return {id, prefix: prefix.map(x => [...x]), incoming: [...incoming], steps};
}
function relationSelectivityContext(basePointId, historyId, historySnapshots, candidateIncoming) {
  const queryRecord = {id: `query_${basePointId}_${historyId}`, x: [0.8, 0.2], action: "forward", arrival: null, completed: false, observedAt: 2, recordVersion: `query-${basePointId}-${historyId}-v1`, incoming: [...historySnapshots[historyId].incoming]};
  const candidateA = {id: "A", x: [0.3, 0.7], y: [0.4, 0.7], action: "forward", arrival: 1, completed: true, recordVersion: `candidate-A-${basePointId}-v1`, incoming: [...candidateIncoming.A]};
  const candidateB = {id: "B", x: [0.7, -0.2], y: [0.6, -0.2], action: "forward", arrival: 2, completed: true, recordVersion: `candidate-B-${basePointId}-v1`, incoming: [...candidateIncoming.B]};
  const predictionA = relationTransition(candidateA, queryRecord), predictionB = relationTransition(candidateB, queryRecord);
  return {basePointId, historyId, queryRecord, candidateA, candidateB, predictions: {A: predictionA, B: predictionB}};
}
function relationSelectivityProbabilitySummary(before, after) {
  return {priorBefore: before.teacher.pBefore, priorAfter: after.teacher.pBefore, priorChange: after.teacher.pBefore - before.teacher.pBefore, posteriorBefore: before.teacher.pPlus, posteriorAfter: after.teacher.pPlus, posteriorChange: after.teacher.pPlus - before.teacher.pPlus, logitBefore: before.f, logitAfter: after.f};
}
function relationSelectivityRunSequence({cellId, method, target, contexts, yByHistory, order, thetaInitial, options, teacherSpec}) {
  let theta = [...thetaInitial], evaluationCount = 0, totalSolverCalls = 0, totalLinearResponseSolves = 0, totalElapsedMs = 0;
  const evaluate = (historyId, thetaValue) => {
    const result = relationSelectivityEvaluate(method, contexts[historyId], thetaValue, options, yByHistory[historyId], teacherSpec);
    evaluationCount += 1;
    totalSolverCalls += result.solveSummary.totalSolverCalls;
    totalLinearResponseSolves += result.solveSummary.linearResponseSolves;
    totalElapsedMs += result.elapsedMs;
    return result;
  };
  const evaluateBoth = thetaValue => ({H0: evaluate("H0", thetaValue), H1: evaluate("H1", thetaValue)}), before = evaluateBoth(theta), writes = [], afterStates = [];
  let current = before;
  for (let stepIndex = 0; stepIndex < order.length; stepIndex++) {
    const historyId = order[stepIndex], startEvaluation = current[historyId], thetaBefore = [...theta], match = relationSelectivityMatchSelfCorrection({method, context: contexts[historyId], theta: thetaBefore, options, y: yByHistory[historyId], teacherSpec, baseEvaluation: startEvaluation, target}), applied = match.matched;
    const thetaAfter = applied ? [...match.thetaAfter] : [...thetaBefore], deltaTheta = applied ? relationSub(thetaAfter, thetaBefore) : [0, 0, 0, 0];
    if (applied) {
      const expectedDeltaTheta = startEvaluation.gradient.map(value => match.eta * startEvaluation.delta * value);
      assert.ok(relationGap(deltaTheta, expectedDeltaTheta) <= 2e-12, "matched selectivity write does not follow theta+eta*delta*g");
      relationAssertTheta(thetaAfter);
      relationCoefficients(thetaAfter);
    }
    theta = thetaAfter;
    current = evaluateBoth(theta);
    afterStates.push(current);
    const ownAfter = current[historyId], ownChange = Math.sign(startEvaluation.delta) === 0 ? null : Math.sign(startEvaluation.delta) * (ownAfter.f - startEvaluation.f);
    writes.push({step: stepIndex + 1, historyId, y: [...yByHistory[historyId]], thetaBefore, fBefore: startEvaluation.f, gradientBefore: [...startEvaluation.gradient], gradientNormBefore: startEvaluation.gradientNorm, delta: startEvaluation.delta, signFromActualDelta: Math.sign(startEvaluation.delta), teacherBefore: {lambda: startEvaluation.teacher.lambda, pBefore: startEvaluation.teacher.pBefore, pPlus: startEvaluation.teacher.pPlus, logit: startEvaluation.f}, match, applied, deltaTheta, parameterStepL2: relationEuclidean(deltaTheta), thetaAfter: [...thetaAfter], fSelfAfter: ownAfter.f, signedOwnCorrection: applied ? ownChange : null, teacherAfter: {lambda: ownAfter.teacher.lambda, pBefore: ownAfter.teacher.pBefore, pPlus: ownAfter.teacher.pPlus, logit: ownAfter.f}, probability: relationSelectivityProbabilitySummary(startEvaluation, ownAfter)});
  }
  const after1 = afterStates[0], after2 = afterStates[1], firstHistory = order[0], secondHistory = order[1], firstWrite = writes[0], secondWrite = writes[1], U0 = firstWrite.applied ? Math.sign(firstWrite.delta) * (after1[firstHistory].f - before[firstHistory].f) : null, U1 = secondWrite.applied ? Math.sign(secondWrite.delta) * (after2[secondHistory].f - after1[secondHistory].f) : null;
  // crossFirst is the first write's effect on the other association.  It is
  // reported separately.  The retained-first interference metric is the
  // second write's effect on the first association (crossSecond), signed by
  // the first write's actual teacher delta.
  const crossFirst = after1[secondHistory].f - before[secondHistory].f, signedCrossFirstByOtherDelta = Math.sign(secondWrite.delta) === 0 ? null : Math.sign(secondWrite.delta) * crossFirst, crossSecond = after2[firstHistory].f - after1[firstHistory].f, signedCrossSecondByFirstDelta = Math.sign(firstWrite.delta) === 0 ? null : Math.sign(firstWrite.delta) * crossSecond, H0 = U0 === null || signedCrossSecondByFirstDelta === null ? null : Math.max(0, -signedCrossSecondByFirstDelta), D0 = U0 === null || signedCrossSecondByFirstDelta === null ? null : U0 + signedCrossSecondByFirstDelta, HOverU0 = H0 === null || U0 <= 0 ? null : H0 / U0, DOverU0 = D0 === null || U0 <= 0 ? null : D0 / U0, absCrossOverU1 = U1 === null || U1 <= 0 ? null : Math.abs(crossSecond) / U1;
  const secondStepStartGram = relationSelectivityGram([after1.H0.gradient, after1.H1.gradient]);
  const scienceGate = {matchedBoth: firstWrite.applied && secondWrite.applied, U0, U1, crossFirst, signedCrossFirstByOtherDelta, crossSecond, signedCrossSecondByFirstDelta, H0, D0, HOverU0, DOverU0, absCrossOverU1, HOverUPassed: HOverU0 !== null && HOverU0 <= 0.25, DOverUPassed: DOverU0 !== null && DOverU0 >= 0.5, passed: firstWrite.applied && secondWrite.applied && HOverU0 !== null && HOverU0 <= 0.25 && DOverU0 !== null && DOverU0 >= 0.5};
  const result = {cellId, method, target, order: [...order], fTrajectory: {H0: {fBefore: before.H0.f, fAfter1: after1.H0.f, fAfter2: after2.H0.f}, H1: {fBefore: before.H1.f, fAfter1: after1.H1.f, fAfter2: after2.H1.f}}, initialTeacherDelta: {H0: before.H0.delta, H1: before.H1.delta}, writes, U0, U1, crossEffects: {firstWriteOnOther: {raw: crossFirst, signedByOtherDelta: signedCrossFirstByOtherDelta, interpretation: "reported separately; this is not the retained-first interference metric"}, retainedFirstAfterSecondWrite: {raw: crossSecond, signedByFirstDelta: signedCrossSecondByFirstDelta, H: H0, D: D0, HOverU: HOverU0, DOverU: DOverU0, absCrossOverUSecond: absCrossOverU1, interpretation: "second write's effect on the first association, signed by the first write's actual delta"}, secondWriteOnFirst: {raw: crossSecond, signedByFirstDelta: signedCrossSecondByFirstDelta, H: H0, D: D0}}, secondStepStart: {theta: [...writes[0].thetaAfter], f: {H0: after1.H0.f, H1: after1.H1.f}, gradients: {H0: [...after1.H0.gradient], H1: [...after1.H1.gradient]}, gram: secondStepStartGram, note: "actual post-first-write gradients/Gram; not the initial Gram finite-step prediction"}, solverCost: {fullEvaluationCount: evaluationCount, fullEvaluationSolverCalls: totalSolverCalls, fullEvaluationLinearResponseSolves: totalLinearResponseSolves, fullEvaluationElapsedMs: totalElapsedMs, matchingSolverCalls: writes.reduce((sum, write) => sum + write.match.solveCount, 0), matchingLinearResponseSolves: writes.reduce((sum, write) => sum + write.match.linearResponseSolves, 0), matchingEvaluations: writes.reduce((sum, write) => sum + write.match.evaluations, 0), matchingElapsedMs: writes.reduce((sum, write) => sum + write.match.elapsedMs, 0)}, scienceGate};
  result.trajectoryAudit = relationSelectivityTrajectoryAudit(result);
  return result;
}
function relationSelectivityRecomputeTrajectoryMetrics(sequence) {
  const firstHistory = sequence.order[0], secondHistory = sequence.order[1], firstWrite = sequence.writes[0], secondWrite = sequence.writes[1], firstTrajectory = sequence.fTrajectory[firstHistory], secondTrajectory = sequence.fTrajectory[secondHistory], signFirst = Math.sign(firstWrite.delta), signSecond = Math.sign(secondWrite.delta), U0 = firstWrite.applied ? signFirst * (firstTrajectory.fAfter1 - firstTrajectory.fBefore) : null, U1 = secondWrite.applied ? signSecond * (secondTrajectory.fAfter2 - secondTrajectory.fAfter1) : null, crossFirst = secondTrajectory.fAfter1 - secondTrajectory.fBefore, crossSecond = firstTrajectory.fAfter2 - firstTrajectory.fAfter1, signedCrossFirstByOtherDelta = signSecond === 0 ? null : signSecond * crossFirst, signedCrossSecondByFirstDelta = signFirst === 0 ? null : signFirst * crossSecond, H0 = U0 === null || signedCrossSecondByFirstDelta === null ? null : Math.max(0, -signedCrossSecondByFirstDelta), D0 = U0 === null || signedCrossSecondByFirstDelta === null ? null : U0 + signedCrossSecondByFirstDelta, HOverU0 = H0 === null || U0 <= 0 ? null : H0 / U0, DOverU0 = D0 === null || U0 <= 0 ? null : D0 / U0, absCrossOverU1 = U1 === null || U1 <= 0 ? null : Math.abs(crossSecond) / U1;
  return {firstHistory, secondHistory, signFirst, signSecond, U0, U1, crossFirst, signedCrossFirstByOtherDelta, crossSecond, signedCrossSecondByFirstDelta, H0, D0, HOverU0, DOverU0, absCrossOverU1};
}
function relationSelectivityScalarClose(left, right, tolerance = 2e-12) {
  if (left === null || right === null) return left === right;
  return relationFinite(left) && relationFinite(right) && Math.abs(left - right) <= tolerance * Math.max(1, Math.abs(left), Math.abs(right));
}
function relationSelectivityTrajectoryAudit(sequence) {
  const expected = relationSelectivityRecomputeTrajectoryMetrics(sequence), actual = sequence.scienceGate, fields = ["U0", "U1", "crossFirst", "signedCrossFirstByOtherDelta", "crossSecond", "signedCrossSecondByFirstDelta", "H0", "D0", "HOverU0", "DOverU0", "absCrossOverU1"], passed = fields.every(field => relationSelectivityScalarClose(expected[field], actual[field])) && relationSelectivityScalarClose(sequence.crossEffects.firstWriteOnOther.raw, expected.crossFirst) && relationSelectivityScalarClose(sequence.crossEffects.retainedFirstAfterSecondWrite.raw, expected.crossSecond) && relationSelectivityScalarClose(sequence.crossEffects.retainedFirstAfterSecondWrite.H, expected.H0) && relationSelectivityScalarClose(sequence.crossEffects.retainedFirstAfterSecondWrite.D, expected.D0);
  return {passed, fields, expected, actual: Object.fromEntries(fields.map(field => [field, actual[field]])), note: "independently recomputed from stored fBefore/fAfter1/fAfter2 and actual first/second teacher signs; crossFirst is not used for retained-first H/D"};
}
function relationSelectivitySyntheticTrajectoryUnitCheck() {
  const sequence = {order: ["H0", "H1"], fTrajectory: {H0: {fBefore: 0, fAfter1: 1, fAfter2: 0.8}, H1: {fBefore: 0, fAfter1: -0.3, fAfter2: -0.5}}, writes: [{delta: 1, applied: true}, {delta: -1, applied: true}]}, metrics = relationSelectivityRecomputeTrajectoryMetrics(sequence), expected = {U0: 1, U1: 0.2, crossFirst: -0.3, crossSecond: -0.2, signedCrossSecondByFirstDelta: -0.2, H0: 0.2, D0: 0.8, HOverU0: 0.2, DOverU0: 0.8, absCrossOverU1: 1};
  return {passed: Object.keys(expected).every(field => relationSelectivityScalarClose(metrics[field], expected[field])), sequence, metrics, expected, note: "opposite first/second delta signs and distinct fBefore/fAfter1/fAfter2 distinguish crossFirst from retained-first crossSecond"};
}
function relationSelectivityGradientCheck(method, context, theta, options, y, teacherSpec) {
  const analytic = relationSelectivityEvaluate(method, context, theta, options, y, teacherSpec), steps = [1e-3, 1e-4, 1e-5], finiteDifference = steps.map(step => {
    const gradient = theta.map((_, parameter) => {
      const plusTheta = theta.map((value, index) => value + (index === parameter ? step : 0)), minusTheta = theta.map((value, index) => value - (index === parameter ? step : 0));
      return (relationSelectivityReadValue(method, context, plusTheta, options).f - relationSelectivityReadValue(method, context, minusTheta, options).f) / (2 * step);
    });
    const absoluteError = relationMaxAbs(gradient.map((value, index) => value - analytic.gradient[index])), relativeError = absoluteError / Math.max(analytic.gradientNorm, Number.MIN_VALUE);
    return {step, finiteDifference: gradient, absoluteError, relativeError, passed: relativeError <= 5e-6};
  });
  return {method, analyticGradient: analytic.gradient, analyticGradientNorm: analytic.gradientNorm, minimumRequiredNorm: 1e-6, finiteDifference, passed: analytic.gradientNorm > 1e-6 && finiteDifference.every(row => row.passed)};
}
function relationSelectivityFiniteNumbers(value) {
  if (value === null || value === undefined || typeof value === "string" || typeof value === "boolean") return true;
  if (typeof value === "number") return Number.isFinite(value);
  if (Array.isArray(value)) return value.every(relationSelectivityFiniteNumbers);
  if (typeof value === "object") return Object.values(value).every(relationSelectivityFiniteNumbers);
  return true;
}
function runRelationEnergySelectivityCheck() {
  const thetaInitial = [0.18, -0.11, 0.07, -0.05], options = {h: 0.4, epsilon: 0, gamma: 1, tau: 1, tolerance: 1e-12, maxIterations: 80}, prefixes = {H0: [[1, 0], [0, 1]], H1: [[0, 1], [1, 0]]}, historySnapshots = {H0: relationSelectivityHistorySnapshot("H0", prefixes.H0, thetaInitial, options), H1: relationSelectivityHistorySnapshot("H1", prefixes.H1, thetaInitial, options)}, candidateIncoming = {cold: {A: relationZero(), B: relationZero()}, warm: {A: [...historySnapshots.H0.incoming], B: [...historySnapshots.H1.incoming]}, warmSwapped: {A: [...historySnapshots.H1.incoming], B: [...historySnapshots.H0.incoming]}}, basePointIds = ["cold", "warm", "warmSwapped"], historyIds = ["H0", "H1"], contexts = {}, teacherSpec = relationFreezeTeacherSpec({version: "selectivity-gaussian-empty-v1", pi0: 0.1, variance: 1, emptyMean: [2.5, -2.5], emptyVariance: 1}, true), yPlus = [0.9, 0.23], yMinus = [0.7, 0.23], evidenceDirections = [{id: "plus_plus", yByHistory: {H0: [...yPlus], H1: [...yPlus]}}, {id: "plus_minus", yByHistory: {H0: [...yPlus], H1: [...yMinus]}}, {id: "minus_plus", yByHistory: {H0: [...yMinus], H1: [...yPlus]}}, {id: "minus_minus", yByHistory: {H0: [...yMinus], H1: [...yMinus]}}], orders = [{id: "H0_then_H1", historyIds: ["H0", "H1"]}, {id: "H1_then_H0", historyIds: ["H1", "H0"]}], targets = [{id: "U_1e-4", value: 1e-4, role: "pressure_low"}, {id: "U_1e-3", value: 1e-3, role: "primary"}, {id: "U_1e-2", value: 1e-2, role: "pressure_high"}];
  for (const basePointId of basePointIds) {
    contexts[basePointId] = {};
    for (const historyId of historyIds) contexts[basePointId][historyId] = relationSelectivityContext(basePointId, historyId, historySnapshots, candidateIncoming[basePointId]);
  }
  const baseSensitivity = {}, gradientChecks = {};
  for (const basePointId of basePointIds) {
    baseSensitivity[basePointId] = {};
    for (const method of ["relationEnergy", "distance"]) {
      const evaluations = historyIds.map(historyId => relationSelectivityEvaluate(method, contexts[basePointId][historyId], thetaInitial, options, yPlus, teacherSpec));
      baseSensitivity[basePointId][method] = {...relationSelectivityGram(evaluations.map(evaluation => evaluation.gradient)), delta: evaluations.map(evaluation => evaluation.delta), f: evaluations.map(evaluation => evaluation.f), gradientNormMinimum: Math.min(...evaluations.map(evaluation => evaluation.gradientNorm))};
    }
  }
  gradientChecks.relationEnergy = relationSelectivityGradientCheck("relationEnergy", contexts.warm.H0, thetaInitial, options, yPlus, teacherSpec);
  gradientChecks.distance = relationSelectivityGradientCheck("distance", contexts.warm.H0, thetaInitial, options, yPlus, teacherSpec);
  const syntheticTrajectoryUnitTest = relationSelectivitySyntheticTrajectoryUnitCheck();
  const cells = [];
  for (const basePointId of basePointIds) for (const direction of evidenceDirections) for (const order of orders) {
    const cellId = basePointId + "__" + direction.id + "__" + order.id, methods = {};
    for (const method of ["relationEnergy", "distance"]) {
      methods[method] = {};
      for (const target of targets) methods[method][target.id] = relationSelectivityRunSequence({cellId, method, target: target.value, contexts: contexts[basePointId], yByHistory: direction.yByHistory, order: order.historyIds, thetaInitial, options, teacherSpec});
    }
    cells.push({cellId, basePoint: basePointId, evidenceDirection: direction.id, order: order.id, methods});
  }
  const allSequences = cells.flatMap(cell => Object.values(cell.methods).flatMap(methodResults => Object.values(methodResults))), trajectoryChecks = allSequences.map(sequence => ({cellId: sequence.cellId, method: sequence.method, target: sequence.target, ...relationSelectivityTrajectoryAudit(sequence)})), trajectoryChecksPassed = trajectoryChecks.every(check => check.passed), matchedWriteChecksPassed = allSequences.every(sequence => sequence.writes.every(write => !write.applied || (write.match.matched && write.match.lambda <= 1 + 1e-12 && write.match.parameterStepL2 <= 1 + 1e-12 && write.match.matchError <= write.match.targetTolerance))), completeCellsPassed = cells.length === 24 && cells.every(cell => Object.keys(cell.methods).length === 2 && ["relationEnergy", "distance"].every(method => Object.keys(cell.methods[method]).length === targets.length)), baseSensitivityPassed = Object.values(baseSensitivity).every(methods => Object.values(methods).every(summary => summary.gradientNormMinimum > 1e-6));
  const pairedComparisons = cells.flatMap(cell => targets.map(target => {
    const relationEnergy = cell.methods.relationEnergy[target.id], distance = cell.methods.distance[target.id], bothMatched = relationEnergy.scienceGate.matchedBoth && distance.scienceGate.matchedBoth;
    return {cellId: cell.cellId, basePoint: cell.basePoint, evidenceDirection: cell.evidenceDirection, order: cell.order, target: target.value, targetId: target.id, bothMatched, relationEnergy: {matched: relationEnergy.scienceGate.matchedBoth, HOverU: relationEnergy.scienceGate.HOverU0, DOverU: relationEnergy.scienceGate.DOverU0, absCrossOverUSecond: relationEnergy.scienceGate.absCrossOverU1}, distance: {matched: distance.scienceGate.matchedBoth, HOverU: distance.scienceGate.HOverU0, DOverU: distance.scienceGate.DOverU0, absCrossOverUSecond: distance.scienceGate.absCrossOverU1}, comparison: bothMatched ? {deltaAbsCrossOverUSecond: relationEnergy.scienceGate.absCrossOverU1 - distance.scienceGate.absCrossOverU1, noWinnerSelected: true} : {status: "not-compared-unmatched", noWinnerSelected: true}};
  }));
  const summary = {cellCount: cells.length, expectedCellCount: 24, methods: {}, pairedComparisonCount: pairedComparisons.length, pairedMatchedCount: pairedComparisons.filter(comparison => comparison.bothMatched).length, primaryTarget: "U_1e-3", targetRoles: Object.fromEntries(targets.map(target => [target.id, target.role]))};
  for (const method of ["relationEnergy", "distance"]) {
    summary.methods[method] = {};
    for (const target of targets) {
      const sequences = cells.map(cell => cell.methods[method][target.id]);
      summary.methods[method][target.id] = {total: sequences.length, firstMatched: sequences.filter(sequence => sequence.writes[0].applied).length, secondMatched: sequences.filter(sequence => sequence.writes[1].applied).length, bothMatched: sequences.filter(sequence => sequence.scienceGate.matchedBoth).length, thresholdPassed: sequences.filter(sequence => sequence.scienceGate.passed).length, unmatchedReasons: sequences.flatMap(sequence => sequence.writes.filter(write => !write.applied).map(write => write.match.unmatchedReason)).reduce((counts, reason) => ({...counts, [reason]: (counts[reason] || 0) + 1}), {})};
    }
  }
  const primaryTargetId = "U_1e-3", primaryRelationEnergySequences = cells.map(cell => cell.methods.relationEnergy[primaryTargetId]), scientificGatePassed = primaryRelationEnergySequences.every(sequence => sequence.scienceGate.matchedBoth && sequence.scienceGate.passed && sequence.scienceGate.absCrossOverU1 !== null && sequence.scienceGate.absCrossOverU1 <= 0.25), panelAllMethodsPrimaryThresholdsPassed = cells.every(cell => ["relationEnergy", "distance"].every(method => cell.methods[method][primaryTargetId].scienceGate.passed)), correctnessChecks = {gradientRelationEnergy: gradientChecks.relationEnergy.passed, gradientDistance: gradientChecks.distance.passed, nonzeroGradientNorms: baseSensitivityPassed, trajectoryMetricRecomputation: trajectoryChecksPassed, syntheticTrajectoryUnitTest: syntheticTrajectoryUnitTest.passed, matchedWriteBudgetAndTolerance: matchedWriteChecksPassed, complete24Cells: completeCellsPassed, allReportedNumbersFinite: relationSelectivityFiniteNumbers({baseSensitivity, gradientChecks, syntheticTrajectoryUnitTest, cells, pairedComparisons})}, correctnessPassed = Object.values(correctnessChecks).every(value => value === true);
  const fixture = {thetaInitial: [...thetaInitial], options, prefixes, historySnapshots, candidateIncoming, records: {query: {x: [0.8, 0.2], action: "forward", observedAt: 2, completed: false}, candidateA: {x: [0.3, 0.7], y: [0.4, 0.7], action: "forward", arrival: 1, completed: true}, candidateB: {x: [0.7, -0.2], y: [0.6, -0.2], action: "forward", arrival: 2, completed: true}}, evidenceDirections, orders, teacherSpec, candidateContextRule: "cold=(0,0) incoming; warm=(H0,H1); warmSwapped=(H1,H0)", learnerInputs: "only numeric y is passed to each evaluation; direction/order labels remain report metadata"};
  return {scope: "offline cross-association selectivity diagnostic after matching each method's own correction; not an online learning-rate result and not an automatic superiority proof", entryPoint: "--relation-energy-selectivity-check", diagnosticOnly: true, correctnessPassed, correctnessChecks, scientificGatePassed, panelAllMethodsPrimaryThresholdsPassed, equations: {relationEnergy: "existing associationMargin/e with full four-port paired energy", distance: "A_D(q,C)=gamma/2||S Z_q-S Z_C||^2; f_D=(A_D_B-A_D_A)/tau", update: "theta'=theta+eta*delta*g; lambda=eta*|delta|*||g|| is search coordinate only", interference: "retained-first cross= f_first(after2)-f_first(after1), signed by sign(first delta); first-write-on-other cross is separately reported"}, configuration: {theta: [...thetaInitial], h: options.h, epsilon: options.epsilon, gamma: options.gamma, tau: options.tau, targetTolerance: "max(1e-10,1e-6*U)", lambdaLimit: 1, absoluteCriterion: "|C|/U_second <= 0.25 for the primary E diagnostic gate", search: "nonnegative doubling then bisection; unmatched means no bracket within this budget, not global impossibility", teacherFreeze: "one cloned Gaussian outside model for all cells; no mutable callback; online teacher version is fixed", supplementalEpsilon005: {epsilon: 0.05, executed: false, reason: "main frozen diagnostic uses epsilon=0; no result selection by supplement"}}, fixture, gradientChecks, syntheticTrajectoryUnitTest, baseSensitivity, trajectoryAudit: {total: trajectoryChecks.length, passed: trajectoryChecks.filter(check => check.passed).length, failed: trajectoryChecks.filter(check => !check.passed)}, targets, summary, pairedComparisons, cells, interpretation: {scientificGate: "scientificGatePassed is only the current primary relation-energy E criterion: all primary E cells must be matched, pass H/U and D/U thresholds, and satisfy the absolute cross criterion; it is not a baseline or superiority gate", panelGate: "panelAllMethodsPrimaryThresholdsPassed reports the old all-method primary threshold panel separately and is not used to define scientificGatePassed", scientificThresholds: "per sequence H/U<=0.25 and D/U>=0.50 are retained in scienceGate; absCross/U_second is also reported", matching: "each lambda trial freezes the starting theta, actual delta and gradient and recomputes f only; next write recomputes gradients and teacher", probability: "probability changes and original f logits are retained to show the common energy-unit limitation", noAutomaticClaim: true}, limitations: {etaIsOfflineDiagnostic: true, continuousTrajectoryNotClaimed: true, noTrainingSimulationOrCalibration: true, densePairedSolve: "relation-energy uses the existing dense 10x10 Newton solve; distance uses free solves plus sharedResponse analytic chain; no O(H) Schur claim"}};
}
if (require.main === module) {
  if (relationEnergyOnly) {
    const relationEnergyQualification = runRelationEnergyCheck();
    const relationRendered = JSON.stringify(relationEnergyQualification, null, 2) + "\n";
    fs.writeFileSync("idea-stage/RELATION_ENERGY_CHECK_20260909.json", relationRendered, "utf8");
    if (!relationEnergyQualification.passed) process.exitCode = 1;
    process.stdout.write(relationRendered);
  } else if (relationEnergySelectivityOnly) {
    const selectivityQualification = runRelationEnergySelectivityCheck(), selectivityRendered = JSON.stringify(selectivityQualification, null, 2) + "\n";
    fs.writeFileSync("idea-stage/RELATION_ENERGY_SELECTIVITY_20260909.json", selectivityRendered, "utf8");
    if (!selectivityQualification.correctnessPassed) process.exitCode = 1;
    process.stdout.write(JSON.stringify({entryPoint: selectivityQualification.entryPoint, correctnessPassed: selectivityQualification.correctnessPassed, scientificGatePassed: selectivityQualification.scientificGatePassed, summary: selectivityQualification.summary, resultPath: "idea-stage/RELATION_ENERGY_SELECTIVITY_20260909.json"}, null, 2) + "\n");
  } else {
function voltage(segments, epsilon, step = 0.002, mismatch = 0, initial = [0, 0, 0, 0, 0]) {
  const bound = Math.max(...segments.flatMap(s => [Math.abs(s.a), Math.abs(s.b)]));
  // leak = fixed leak/2 + tonic shunt/2; even the modulated shunt stays positive.
  if (!(p.leak / 2 > Math.abs(epsilon) * bound)) throw new Error("modulated shunt is not strictly positive");
  return integrate(segments, initial, step, (v, a, b) => starFlow(v, a, b, epsilon, p.axial, mismatch));
}
function leadingOrder(segments, step = 0.002) {
  const lambda = p.leak + p.axial, somaLambda = p.somaLeak + 4 * p.axial;
  return integrate(segments, [0, 0, 0, 0], step, ([u, v, m1, s1], a, b) => [
    -lambda * u + a,
    -lambda * v + b,
    -lambda * m1 + 4 * p.axial * s1 - 2 * (b * u - a * v),
    -somaLambda * s1 + p.axial * m1,
  ]);
}
const base = voltage(path, 0.05), swapped = voltage(path.map(s => ({a: s.b, b: s.a, t: s.t})), 0.05);
check("port exchange reverses raw soma", Math.abs(base[4] + swapped[4]) < 1e-12);
const noModulation = voltage(path, 0);
check("without conductance modulation balanced soma is silent", Math.abs(noModulation[4]) < 1e-14 && norm(noModulation) > 0.01);
const sameStream = voltage(path.map(s => ({a: s.a, b: s.a, t: s.t})), 0.05);
check("identical streams have zero soma contrast", Math.abs(sameStream[4]) < 1e-14);
const negativeEpsilon = voltage(path, -0.05);
check("raw soma is odd in modulation strength", Math.abs(base[4] + negativeEpsilon[4]) < 1e-12);
const reference = -leadingOrder(path)[3] / 2;
const epsilonChecks = [0.2, 0.1, 0.05, 0.025].map(epsilon => {
  const state = voltage(path, epsilon), decoded = -state[4] / (2 * epsilon);
  const bound = epsilon ** 2 * 4 ** 4 * Math.exp(epsilon * 4) / 48;
  return {epsilon, rawSoma: state[4], decoded, reference, error: Math.abs(decoded - reference), errorBound: bound};
});
check("normalized small-modulation error is second order", epsilonChecks.slice(1).every((r, i) => epsilonChecks[i].error / r.error > 3.9 && epsilonChecks[i].error / r.error < 4.1));
check("finite-horizon conservative error bound holds on these inputs", epsilonChecks.every(r => r.error <= r.errorBound));
const refined = voltage(path, 0.05, 0.001), discretizationGap = gap(base, refined);
check("RK4 step refinement agrees", discretizationGap < 1e-10);
const stopped = voltage([{a: 0, b: 0, t: 40}], 0.05, 0.002, 0, base);
check("after input stops every activity state decays", norm(stopped) < 1e-12);
const slowerPath = path.map(s => ({a: s.a / 2, b: s.b / 2, t: s.t * 2}));
const slower = voltage(slowerPath, 0.05);
check("same integrated path at another speed does not preserve soma", Math.abs(base[4] - slower[4]) > 1e-5);
const mismatched = voltage(path, 0, 0.002, 0.01);
check("pair mismatch creates first-order soma leakage even at epsilon zero", Math.abs(mismatched[4]) > 1e-5);
assert.throws(() => voltage(path, 0.5), /shunt/);
checks.push("invalid conductance range fails explicitly");

// The attachment's compensated update is a contraction of an antisymmetric dyad.
const b = [1, 2, -1], e = [0.5, -1, 2], b2 = b.reduce((s, x) => s + x * x, 0);
const dot = (x, y) => x.reduce((s, v, i) => s + v * y[i], 0);
const omega = e.map((eh, h) => b.map((bk, k) => eh * bk - b[h] * e[k]));
const contraction = omega.map(row => dot(row, b) / b2), residual = e.map((x, i) => x - dot(b, e) / b2 * b[i]);
check("attachment wedge contraction equals compensated local residual", gap(contraction, residual) < 1e-14 && Math.abs(dot(b, contraction)) < 1e-14);
const ratio = (axis, q) => { const denominator = dot(q, q); return 1 - dot(axis, q) ** 2 / (dot(axis, axis) * denominator); };
const counterexamples = {
  hiddenButNoSafeDirection: {chi: ratio([1, 0], [0, 1]), kappa: ratio([0, 1], [0, 1])},
  visibleButSafeDirection: {chi: ratio([1, 0], [1, 0]), kappa: ratio([0, 1], [1, 0])},
};
check("chi and kappa do not imply one another", counterexamples.hiddenButNoSafeDirection.chi === 1 && counterexamples.hiddenButNoSafeDirection.kappa === 0 && counterexamples.visibleButSafeDirection.chi === 0 && counterexamples.visibleButSafeDirection.kappa === 1);
// Two identical soma histories, branch direction D=+-1 and legal future target Y=D.
const conditionalRisk = {somaOnlyMse: (1 ** 2 + (-1) ** 2) / 2, branchMse: 0, conditionalMeanGap: (1 ** 2 + (-1) ** 2) / 2};
check("finite conditional-prediction identity example", conditionalRisk.somaOnlyMse - conditionalRisk.branchMse === conditionalRisk.conditionalMeanGap);
const originalAxial = p.axial;
const couplingChecks = [0, 0.1, 0.4, 2, 10].map(axial => {
  p.axial = axial;
  return {axial, normalizedSoma: -voltage(path, 0.05, 0.001)[4] / 0.1};
});
p.axial = originalAxial;
check("finite coupling diagnostic: isolated zero, strong coupling attenuates this input", couplingChecks[0].normalizedSoma === 0 && Math.abs(couplingChecks[4].normalizedSoma) < Math.abs(couplingChecks[2].normalizedSoma));
function twoPulses(wait) {
  return [{a: 1, b: 0, t: 0.2}, {a: 0, b: 0, t: wait}, {a: 0, b: 1, t: 0.2}, {a: 0, b: 0, t: 0.4}];
}
const lagChecks = [0.1, 0.4, 2].map(axial => {
  p.axial = axial;
  const short = voltage(twoPulses(0.2), 0.05)[4], long = voltage(twoPulses(0.8), 0.05)[4];
  return {axial, short, long, ratio: long / short, expected: Math.exp(-(p.leak + axial) * 0.6)};
});
p.axial = originalAxial;
check("coupling changes lag selectivity, not just an overall gain", lagChecks.every(r => Math.abs(r.ratio - r.expected) < 1e-10) && Math.abs(lagChecks[0].ratio - lagChecks[2].ratio) > 0.1);

// J: conditional learning-interface qualification, not a visual benchmark.
// Same voltage() path as I. Hidden data-generation parameters never enter update().
const dim = 4;
const matrix = (f = () => 0) => Array.from({length: dim}, (_, i) => Array.from({length: dim}, (_, j) => f(i, j)));
const eye = matrix((i, j) => Number(i === j));
const tr = A => A[0].map((_, j) => A.map(row => row[j]));
const mv = (A, x) => A.map(row => dot(row, x));
const mm = (A, B) => A.map(row => B[0].map((_, j) => dot(row, B.map(r => r[j]))));
const mgap = (A, B) => gap(A.flat(), B.flat());
function solve(A, B) {
  const rows = A.map((row, i) => [...row, ...B[i]]);
  for (let j = 0; j < dim; j++) {
    let pivot = j;
    for (let i = j + 1; i < dim; i++) if (Math.abs(rows[i][j]) > Math.abs(rows[pivot][j])) pivot = i;
    if (Math.abs(rows[pivot][j]) < 1e-14) throw new Error("singular readout solve");
    [rows[j], rows[pivot]] = [rows[pivot], rows[j]];
    const scale = rows[j][j];
    rows[j] = rows[j].map(v => v / scale);
    for (let i = 0; i < dim; i++) if (i !== j) {
      const c = rows[i][j];
      rows[i] = rows[i].map((v, k) => v - c * rows[j][k]);
    }
  }
  return rows.map(row => row.slice(dim));
}
const cayley = (K, c = 1) => solve(matrix((i, j) => eye[i][j] - c * K[i][j] / 2), matrix((i, j) => eye[i][j] + c * K[i][j] / 2));
function observedResidual(K, x, y, c) {
  const m = x.map((v, i) => (v + y[i]) / 2), Km = mv(K, m);
  const r = x.map((v, i) => y[i] - v - c * Km[i]);
  return {m, r, loss: dot(r, r) / 2};
}
const phaseProtocol = (m0, m1, r0, r1) => [
  {a: m0, b: m1, t: 0.2}, {a: 0, b: 0, t: 0.2},
  {a: r0, b: r1, t: 0.2}, {a: 0, b: 0, t: 0.2},
];
const calibration = -leadingOrder(phaseProtocol(1, 0, 0, 1), 0.002)[3] / 2;
assert.ok(calibration > 0);
function update(K, x, y, c, kind, epsilon = 0.05, step = 0.01) {
  const {m, r} = observedResidual(K, x, y, c);
  if (kind === "ordinary") return matrix((i, j) => c * r[i] * m[j]);
  if (kind === "exact") return matrix((i, j) => c * (r[i] * m[j] - m[i] * r[j]) / 2);
  assert.equal(kind, "circuit");
  if (epsilon === 0) return matrix(); // structural ablation, never divide by zero
  const G = matrix();
  for (let i = 0; i < dim; i++) for (let j = i + 1; j < dim; j++) {
    const state = voltage(phaseProtocol(m[i], m[j], r[i], r[j]), epsilon, step);
    G[i][j] = c * state[4] / (4 * epsilon * calibration);
    G[j][i] = -G[i][j]; // explicit pair-tied plasticity assumption
  }
  return G;
}
const addUpdate = (K, G, eta) => matrix((i, j) => K[i][j] + eta * G[i][j]);
// Active research interface: read current teaching drive on the retained state.
// The reset/late-soma `update` above remains a frozen J regression reference.
const pairs = [];
for (let i = 0; i < dim; i++) for (let j = i + 1; j < dim; j++) pairs.push([i, j]);
const pairHistory = v => [(v[0] - v[1]) / 2, -(v[2] - v[3]) / 2];
function historyPrediction(K, states, c) {
  assert.equal(states.length, pairs.length);
  const prediction = Array(dim).fill(0);
  pairs.forEach(([i, j], edge) => {
    const [hi, hj] = pairHistory(states[edge]);
    prediction[i] += c * K[i][j] * hj;
    prediction[j] += c * K[j][i] * hi;
  });
  return prediction;
}
function historyLoss(K, states, target, c) {
  const prediction = historyPrediction(K, states, c);
  const r = target.map((value, i) => value - prediction[i]);
  return {prediction, r, loss: dot(r, r) / 2};
}
function plasticEvent(K, states, target, c, epsilon = 0.05, eta = 0.5) {
  assert.ok(epsilon >= 0 && Number.isFinite(epsilon) && eta > 0 && Number.isFinite(eta));
  assert.ok(c === 1 || c === -1);
  const before = historyLoss(K, states, target, c), delta = matrix();
  if (!(p.leak / 2 > epsilon * norm(before.r))) throw new Error("teaching shunt is not strictly positive");
  const currents = states.map((v, edge) => {
    const [i, j] = pairs[edge], local = teachingCurrent(v, before.r[i], before.r[j], epsilon);
    delta[i][j] = eta * c * local.reduce((sum, value) => sum + value, 0) / 4;
    delta[j][i] = -delta[i][j];
    return local;
  });
  const nextK = addUpdate(K, delta, 1);
  return {before, delta, nextK, currents, after: historyLoss(nextK, states, target, c)};
}
function observePair(K, states, x, y, c, epsilon = 0.05, eta = 0.5, step = 0.002) {
  const midpoint = x.map((v, i) => (v + y[i]) / 2), target = y.map((v, i) => v - x[i]);
  const atTeaching = states.map((v, edge) => {
    const [i, j] = pairs[edge];
    return voltage(phaseProtocol(midpoint[i], midpoint[j], 0, 0).slice(0, 2), epsilon, step, 0, v);
  });
  // Event reads pre-teaching voltages. The same sampled error then drives the
  // existing teaching/blank phases; it is not recomputed using newly written K.
  const event = plasticEvent(K, atTeaching, target, c, epsilon, eta);
  const nextStates = atTeaching.map((v, edge) => {
    const [i, j] = pairs[edge];
    return voltage(phaseProtocol(0, 0, event.before.r[i], event.before.r[j]).slice(2), epsilon, step, 0, v);
  });
  return {...event, atTeaching, nextStates};
}
function powers(T) {
  const inverse = solve(T, eye), result = [eye];
  let plus = eye, minus = eye;
  for (let n = 1; n <= 4; n++) {
    plus = mm(T, plus); minus = mm(inverse, minus);
    result.push(plus, minus);
  }
  return result;
}
const squaredDistance = (x, y) => x.reduce((s, v, i) => s + (v - y[i]) ** 2, 0);
const association = (transforms, x, y) => Math.min(...transforms.map(T => squaredDistance(y, mv(T, x))));
const exampleK = matrix((i, j) => (i - j) / 10);
const exampleX = [0.6, -0.2, 0.4, 0.1], exampleY = [0.3, 0.5, -0.1, 0.4];
const G = update(exampleK, exampleX, exampleY, 1, "exact");
const H = matrix((i, j) => Math.sin(i - j));
const finiteDiff = (observedResidual(addUpdate(exampleK, H, 1e-5), exampleX, exampleY, 1).loss - observedResidual(addUpdate(exampleK, H, -1e-5), exampleX, exampleY, 1).loss) / 2e-5;
check("J projected gradient matches finite difference with Frobenius convention", Math.abs(finiteDiff + dot(G.flat(), H.flat())) < 1e-9);
check("J reversing observations and signed action preserves exact update at fixed K", mgap(G, update(exampleK, exampleY, exampleX, -1, "exact")) < 1e-14);
check("J fixed step-size stress descends on the observed one-step objective", [0.25, 0.5, 0.75].every(eta => observedResidual(addUpdate(exampleK, G, eta), exampleX, exampleY, 1).loss < observedResidual(exampleK, exampleX, exampleY, 1).loss));
const exampleT = cayley(exampleK), reverseT = cayley(exampleK, -1);
check("J Cayley readout is orthogonal and signed action inverse", mgap(mm(tr(exampleT), exampleT), eye) < 1e-14 && mgap(mm(exampleT, reverseT), eye) < 1e-14);
const {r: exampleR} = observedResidual(exampleK, exampleX, exampleY, 1);
const explicitError = exampleY.map((y, i) => y - mv(exampleT, exampleX)[i]);
const residualTransport = mv(solve(matrix((i, j) => eye[i][j] - exampleK[i][j] / 2), eye), exampleR);
check("J implicit residual equals transported explicit prediction error", gap(explicitError, residualTransport) < 1e-14 && dot(explicitError, explicitError) <= dot(exampleR, exampleR));
check("J finite search association is symmetric for skew K", Math.abs(association(powers(exampleT), exampleX, exampleY) - association(powers(exampleT), exampleY, exampleX)) < 1e-14);
const gradientChecks = [0.1, 0.05, 0.025].map(epsilon => {
  const gc = update(exampleK, exampleX, exampleY, 1, "circuit", epsilon, 0.002);
  return {epsilon, maxError: mgap(G, gc), reverseMismatch: mgap(gc, update(exampleK, exampleY, exampleX, -1, "circuit", epsilon, 0.002))};
});
check("J finite-epsilon circuit gradient approaches exact local rule", gradientChecks.slice(1).every((r, i) => gradientChecks[i].maxError / r.maxError > 3.8 && gradientChecks[i].maxError / r.maxError < 4.2));
const zeroResidualY = mv(exampleT, exampleX);
const residualBias = norm(update(exampleK, exampleX, zeroResidualY, 1, "circuit").flat());
check("J finite-epsilon zero-residual bias is recorded, not called exact gradient", residualBias > 1e-10);
const circuitDiscretizationGap = mgap(update(exampleK, exampleX, exampleY, 1, "circuit", 0.05, 0.01), update(exampleK, exampleX, exampleY, 1, "circuit", 0.05, 0.005));
check("J circuit step refinement qualifies this update", circuitDiscretizationGap < 1e-7);
check("J no-modulation ablation has no persistent write", norm(update(exampleK, exampleX, exampleY, 1, "circuit", 0).flat()) === 0);

// Generator owns latent motion; learners only receive returned x/y/c triples.
const Q = [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]].map(row => row.map(x => x / 2));
const blockRotation = matrix((i, j) => {
  if (Math.floor(i / 2) !== Math.floor(j / 2)) return 0;
  const theta = i < 2 ? 0.35 : 0.7;
  return i === j ? Math.cos(theta) : (i % 2 ? 1 : -1) * Math.sin(theta);
});
const truthT = mm(mm(Q, blockRotation), tr(Q));
function rng(seed) { let state = seed >>> 0; return () => { state = (Math.imul(1664525, state) + 1013904223) >>> 0; return state / 4294967296; }; }
function randomUnit(random) { const x = Array.from({length: dim}, () => 2 * random() - 1), size = Math.sqrt(dot(x, x)); return x.map(v => v / size); }
function trainingStream(seed) {
  const random = rng(seed), rows = []; let x;
  for (let t = 0; t < 512; t++) {
    if (t % 32 === 0) x = randomUnit(random);
    const c = random() < 0.5 ? -1 : 1, y = mv(c === 1 ? truthT : tr(truthT), x);
    rows.push({x, y, c}); x = y;
  }
  return rows;
}
function evaluate(K, anchors) {
  const T = cayley(K), candidates = powers(T);
  const truthPowers = powers(truthT);
  let correct = 0, predictionError = 0;
  for (let id = 0; id < anchors.length; id++) {
    predictionError += squaredDistance(mv(truthT, anchors[id]), mv(T, anchors[id]));
    // powers() interleaves +n,-n; indices 5..8 are +/-3,+/-4.
    for (const trueTransform of truthPowers.slice(5, 9)) {
      const query = mv(trueTransform, anchors[id]);
      const scores = anchors.map(anchor => association(candidates, anchor, query));
      const predicted = scores.indexOf(Math.min(...scores));
      correct += Number(predicted === id);
    }
  }
  return {correct, total: anchors.length * 4, meanOneStepSquaredError: predictionError / anchors.length};
}
const learningRuns = creditKillOnly ? [] : [101, 202, 303].map(seed => {
  const rows = trainingStream(seed), random = rng(seed + 10000);
  const anchors = Array.from({length: 32}, () => randomUnit(random));
  const output = {seed, updates: rows.length, frozen: evaluate(matrix(), anchors)};
  const finalMatrices = {};
  for (const kind of ["exact", "circuit", "ordinary"]) {
    let K = matrix();
    for (const {x, y, c} of rows) K = addUpdate(K, update(K, x, y, c, kind), 0.5);
    finalMatrices[kind] = K;
    output[kind] = {...evaluate(K, anchors), persistentK: K};
    assert.ok(K.flat().every(Number.isFinite));
  }
  output.circuitVsExactMaxKGap = mgap(finalMatrices.circuit, finalMatrices.exact);
  return output;
});
if (!creditKillOnly) {
  check("J online updates are persistent and affect later held-out conditional association", learningRuns.every(run => norm(run.circuit.persistentK.flat()) > 0 && run.circuit.correct > run.frozen.correct));
  check("J ordinary exact learner remains the computational reference, not a defeated baseline", learningRuns.every(run => run.circuitVsExactMaxKGap < 0.01));
}
const aliasX = [0.5, 0.5, 0.5, 0.5], aliasY = mv(truthT, aliasX);
const aliasScore = association(powers(truthT), aliasX, aliasY);
check("J distinct identities on one observable orbit are not identifiable by this score", aliasScore < 1e-25);
const changedNorm = aliasX.map(v => v * 1.1);
const nonIsometricScore = association(powers(truthT), aliasX, changedNorm);
check("J same-identity non-isometric feature change need not receive zero score", nonIsometricScore >= 0.01 - 1e-14);
const learningQualification = {
  learningRunsSkipped: creditKillOnly,
  scope: "Constructed 4D orthogonal-action data only; not RGB, navigation performance, biology, or novelty evidence. Exact reference, phase buffers/reset, tied skew weights, and ordinary Cayley/search readout are explicit assumptions.",
  protocol: {seeds: [101, 202, 303], updates: 512, eta: 0.5, dimensions: dim, epsilon: 0.05, phase: 0.2, blank: 0.2, readDelay: 0.2, step: 0.01, gallery: 32, querySteps: [-4, -3, 3, 4], searchRange: [-4, 4], learnerInputs: ["x", "y", "executed c"], truthT, resetActivityEverySample: true, resetWeightsDuringStream: false},
  calibration, finiteDiff, gradientChecks, residualBias, circuitDiscretizationGap, learningRuns, aliasScore, nonIsometricScore,
};

// K: no-reset continuous state/sensitivity and topology-constrained time order.
// Matrix exponentials are reference arithmetic, not a second neural operator.
const zeros5 = [0, 0, 0, 0, 0];
const eye5 = zeros5.map((_, i) => zeros5.map((__, j) => Number(i === j)));
const scaleMatrix = (A, scale) => A.map(row => row.map(v => v * scale));
const addMatrix = (A, B, scale = 1) => A.map((row, i) => row.map((v, j) => v + scale * B[i][j]));
const commutator = (A, B) => addMatrix(mm(A, B), mm(B, A), -1);
function expMatrix(A) {
  const rowBound = Math.max(...A.map(row => row.reduce((s, v) => s + Math.abs(v), 0)));
  const squarings = Math.max(0, Math.ceil(Math.log2(Math.max(1, rowBound))));
  const scaled = scaleMatrix(A, 2 ** -squarings);
  let term = eye5, sum = eye5;
  for (let order = 1; order <= 40; order++) {
    term = scaleMatrix(mm(term, scaled), 1 / order);
    sum = addMatrix(sum, term);
  }
  for (let j = 0; j < squarings; j++) sum = mm(sum, sum);
  return sum;
}
function flowMatrix(a, b, epsilon = 0.2, axial = p.axial) {
  // Analytic Jacobian of the SAME starFlow, including both axial directions.
  const modulation = [b, -b, a, -a, 0];
  return eye5.map((row, i) => row.map((_, j) => {
    if (i === j) return i === 4 ? -(p.somaLeak + 4 * axial) : -(p.leak + axial + epsilon * modulation[i]);
    return (i === 4 || j === 4) ? axial : 0;
  }));
}
function propagator(segments, epsilon = 0.2, axial = p.axial) {
  return segments.reduce((U, seg) => mm(expMatrix(scaleMatrix(flowMatrix(seg.a, seg.b, epsilon, axial), seg.t)), U), eye5);
}
const phaseA = {a: 1, b: 0, t: 0.4}, phaseB = {a: 0, b: 1, t: 0.4};
const Uab = propagator([phaseA, phaseB]), Uba = propagator([phaseB, phaseA]);
const chronologicalGap = mgap(Uab, tr(Uab));
check("K instantaneous reciprocal matrices need not give reciprocal finite-time propagation", chronologicalGap > 1e-5);
check("K reversed gate chronology gives the transpose propagator", mgap(Uba, tr(Uab)) < 1e-14);
const stateExample = [0.3, -0.4, 0.2, 0.1, -0.2];
const analyticFlow = mv(flowMatrix(1, 0), stateExample);
check("K analysis Jacobian realizes the primary current RHS", gap(analyticFlow, starFlow(stateExample, 1, 0, 0.2).map((x, i) => x - [1, -1, 0, 0, 0][i])) < 1e-14);
const positive = voltage([phaseA, phaseB], 0.2, 0.002, 0, stateExample);
const zeroInitial = voltage([phaseA, phaseB], 0.2);
check("K matrix propagator matches no-reset physical state perturbation", gap(mv(Uab, stateExample), positive.map((v, i) => v - zeroInitial[i])) < 1e-10);
const FA = flowMatrix(1, 0), FB = flowMatrix(0, 1), bracket = commutator(FB, FA);
const deltaD = [0.2, -0.2, -0.2, 0.2, 0];
const starBracket = eye5.map((row, i) => row.map((_, j) => p.axial * (Number(i === 4) * deltaD[j] - deltaD[i] * Number(j === 4))));
check("K commutator is the claimed soma-axis outer-product difference", mgap(bracket, starBracket) < 1e-14);
check("K isolated compartments have no chronological matrix effect", mgap(propagator([phaseA, phaseB], 0.2, 0), propagator([phaseB, phaseA], 0.2, 0)) < 1e-14);
check("K constant conductance chronology remains reciprocal", mgap(propagator([phaseA, phaseA]), tr(propagator([phaseA, phaseA]))) < 1e-14);
function magnus(segments) {
  const T = segments.reduce((s, seg) => s + seg.t, 0);
  let omega1 = scaleMatrix(eye5, 0), omega2 = scaleMatrix(eye5, 0), time = 0;
  let temporalMoment = [...zeros5];
  for (let i = 0; i < segments.length; i++) {
    const seg = segments[i], F = flowMatrix(seg.a, seg.b);
    omega1 = addMatrix(omega1, F, seg.t);
    for (let j = 0; j < i; j++) omega2 = addMatrix(omega2, commutator(F, flowMatrix(segments[j].a, segments[j].b)), seg.t * segments[j].t / 2);
    const weightedDuration = (time + seg.t) ** 2 - time ** 2 - T * seg.t;
    const diag = [seg.b, -seg.b, seg.a, -seg.a, 0];
    temporalMoment = temporalMoment.map((v, j) => v + 0.2 * weightedDuration * diag[j]);
    time += seg.t;
  }
  const rankTwo = eye5.map((row, i) => row.map((_, j) => p.axial / 2 * (Number(i === 4) * temporalMoment[j] - temporalMoment[i] * Number(j === 4))));
  return {omega1, omega2, rankTwo, temporalMoment};
}
const threePhases = [{a: 1, b: 0, t: 1}, {a: 0, b: 1, t: 0.7}, {a: -0.6, b: 0.4, t: 0.9}];
const magnusCheck = magnus(threePhases);
check("K accumulated second Magnus term equals topology rank-two formula", mgap(magnusCheck.omega2, magnusCheck.rankTwo) < 1e-14);
const projectionInput = [0.4, -0.1, 0.8, -0.7, 0.3];
const rankTwoAction = mv(magnusCheck.omega2, projectionInput);
const expectedRankTwoAction = magnusCheck.temporalMoment.map((v, i) => p.axial / 2 * (Number(i === 4) * dot(magnusCheck.temporalMoment, projectionInput) - v * projectionInput[4]));
check("K rank-two term acts using two projections, not a scalar output gain", gap(rankTwoAction, expectedRankTwoAction) < 1e-14 && norm(rankTwoAction.slice(0, 4)) > 0);
const approximationChecks = [0.08, 0.04, 0.02].map(durationScale => {
  const segments = threePhases.map(seg => ({...seg, t: seg.t * durationScale}));
  const {omega1, omega2} = magnus(segments), U = propagator(segments);
  return {durationScale, unorderedError: mgap(U, expMatrix(omega1)), orderedError: mgap(U, expMatrix(addMatrix(omega1, omega2))), adjointError: mgap(tr(U), expMatrix(addMatrix(omega1, omega2, -1)))};
});
check("K second-order time term improves short-window propagation reference", approximationChecks.every(r => r.orderedError < r.unorderedError));
check("K short-window second-order Magnus error scales cubically", approximationChecks.slice(1).every((r, i) => approximationChecks[i].orderedError / r.orderedError > 6 && approximationChecks[i].orderedError / r.orderedError < 10));
check("K adjoint requires opposite skew correction", approximationChecks.every(r => Math.abs(r.adjointError - r.orderedError) < 1e-14));
function continuousSensitivity(segments, offset = 0, dropSomaFeedback = false) {
  let augmented = Array(10).fill(0);
  for (const seg of segments) {
    const axial = (seg.axial === undefined ? p.axial : seg.axial) + offset;
    if (axial < 0) throw new Error("negative axial conductance");
    const F = flowMatrix(seg.a, seg.b, 0.05, axial);
    augmented = integrate([seg], augmented, 0.002, values => {
      const z = values.slice(0, 5), P = values.slice(5);
      const source = z.slice(0, 4).map(v => z[4] - v);
      source.push(z.slice(0, 4).reduce((s, v) => s + v, 0) - 4 * z[4]);
      const dP = mv(F, P).map((v, i) => v + source[i] - (dropSomaFeedback && i < 4 ? axial * P[4] : 0));
      return [...starFlow(z, seg.a, seg.b, 0.05, axial), ...dP];
    });
  }
  return {state: augmented.slice(0, 5), sensitivity: augmented.slice(5)};
}
const fixedSensitivity = continuousSensitivity(path);
const differential = (plus, minus) => plus.map((v, i) => (v - minus[i]) / 2e-5);
const fixedFiniteDiff = differential(continuousSensitivity(path, 1e-5).state, continuousSensitivity(path, -1e-5).state);
check("K conductance sensitivity matches fixed-parameter finite differences", gap(fixedSensitivity.sensitivity, fixedFiniteDiff) < 1e-8);
check("K augmented sensitivity uses unchanged physical voltage path", gap(fixedSensitivity.state, base) < 1e-14);
const dropped = continuousSensitivity(path, 0, true);
const feedbackSensitivityGap = gap(dropped.sensitivity, fixedSensitivity.sensitivity);
check("K omitting soma feedback changes the conductance eligibility", feedbackSensitivityGap > 1e-7);
const schedule = path.map((seg, i) => ({...seg, axial: [0.2, 0.4, 0.7, 0.3][i]}));
const scheduleSensitivity = continuousSensitivity(schedule);
const scheduleFiniteDiff = differential(continuousSensitivity(schedule, 1e-5).state, continuousSensitivity(schedule, -1e-5).state);
check("K changing-parameter eligibility matches a common-offset history derivative", gap(scheduleSensitivity.sensitivity, scheduleFiniteDiff) < 1e-8);
const finalOnly = offset => continuousSensitivity(schedule.map((seg, i) => ({...seg, axial: seg.axial + (i === schedule.length - 1 ? offset : 0)}))).state;
const finalOnlyDerivative = differential(finalOnly(1e-5), finalOnly(-1e-5));
const replayDerivative = continuousSensitivity(path.map(seg => ({...seg, axial: 0.3}))).sensitivity;
check("K common-offset, last-transition and current-parameter replay derivatives differ", gap(scheduleSensitivity.sensitivity, finalOnlyDerivative) > 1e-4 && gap(scheduleSensitivity.sensitivity, replayDerivative) > 1e-4);
// Fast kill cases: same first two moments need not mean the same credit operator.
const zeroMomentBaseline = Array.from({length: 4}, () => ({a: 0, b: 0.2, t: 1}));
const zeroMomentCounterexample = [1, -1, -1, 1].map(sign => ({a: sign * 0.5, b: 0.2, t: 1}));
const baselineMoments = magnus(zeroMomentBaseline), counterMoments = magnus(zeroMomentCounterexample);
const equalMomentsPropagatorGap = mgap(propagator(zeroMomentBaseline), propagator(zeroMomentCounterexample));
check("K kill: identical first two Magnus terms do not determine the full propagator", mgap(baselineMoments.omega1, counterMoments.omega1) < 1e-14 && mgap(baselineMoments.omega2, counterMoments.omega2) < 1e-14 && equalMomentsPropagatorGap > 1e-6);
const arbitraryEdges = [0.2, 0.4, 0.6, 0.8];
const heteroF0 = eye5.map((row, i) => row.map((_, j) => {
  if (i === j) return i === 4 ? -p.somaLeak - arbitraryEdges.reduce((s, v) => s + v, 0) : -p.leak - arbitraryEdges[i];
  return i === 4 ? arbitraryEdges[j] : j === 4 ? arbitraryEdges[i] : 0;
}));
const heteroMoment = magnusCheck.temporalMoment.map((v, i) => i === 4 ? 0 : arbitraryEdges[i] * v);
const heteroCommutator = commutator(heteroF0, eye5.map((row, i) => row.map((_, j) => i === j ? magnusCheck.temporalMoment[i] : 0)));
const heteroRankTwo = eye5.map((row, i) => row.map((_, j) => Number(i === 4) * heteroMoment[j] - heteroMoment[i] * Number(j === 4)));
check("K kill: rank-two identity is generic fixed-hub algebra, not dependent on matched axial strengths", mgap(heteroCommutator, heteroRankTwo) < 1e-14);
const anti = addMatrix(propagator(threePhases), tr(propagator(threePhases)), -1);
const pfaffians = [];
for (let omit = 0; omit < 5; omit++) {
  const [i, j, k, l] = [0, 1, 2, 3, 4].filter(index => index !== omit);
  pfaffians.push(anti[i][j] * anti[k][l] - anti[i][k] * anti[j][l] + anti[i][l] * anti[j][k]);
}
const physicalZero = voltage(threePhases, 0.2);
const physicalU = tr(eye5.map(unit => voltage(threePhases, 0.2, 0.002, 0, unit).map((v, i) => v - physicalZero[i])));
const physicalPropagatorGap = mgap(physicalU, propagator(threePhases));
check("K kill: full temporal nonreciprocity is not rank two", norm(pfaffians) > 1e-10 && physicalPropagatorGap < 1e-10);
const killFindings = {zeroMomentBaseline, zeroMomentCounterexample, equalMomentsPropagatorGap, arbitraryEdges, heteroRankTwoGap: mgap(heteroCommutator, heteroRankTwo), fullAntisymmetricPrincipalPfaffians: pfaffians, physicalPropagatorGap};
const continuousCreditQualification = {
  scope: "Five-state conductance algebra and continuous sensitivity only; no HM3D/VPR run, no new-learning-rule or biological implementation claim.",
  parameters: {...p, shuntForPropagation: 0.2, shuntForSensitivity: 0.05}, phaseA, phaseB, threePhases,
  chronologicalGap, Uab, Uba, bracket, magnusCheck, approximationChecks, fixedSensitivity, fixedFiniteDiff,
  feedbackSensitivityGap, schedule, scheduleSensitivity, scheduleFiniteDiff, finalOnlyDerivative, replayDerivative,
  killFindings,
};
let continuityQualification;
if (continuityKillOnly) {
  // Audit J's existing per-pair readout on I's unchanged continuous voltage path.
  // No K update is committed, and no learning trajectory is run.
  const probeX = [0, 1, 0, 0], probeY = [...probeX];
  const probe = observedResidual(matrix(), probeX, probeY, 1);
  const segments = phaseProtocol(probe.m[0], probe.m[1], probe.r[0], probe.r[1]);
  const prefix = sign => [{a: sign, b: 0, t: 0.2}, {a: 0, b: 0, t: 0.2}];
  const cases = [0.1, 0.05, 0.025].map(epsilon => {
    const positiveInitial = voltage(prefix(1), epsilon);
    const negativeInitial = voltage(prefix(-1), epsilon);
    const reset = voltage(segments, epsilon);
    const positive = voltage(segments, epsilon, 0.002, 0, positiveInitial);
    const negative = voltage(segments, epsilon, 0.002, 0, negativeInitial);
    const joined = voltage([...prefix(1), ...segments], epsilon);
    const finerInitial = voltage(prefix(1), epsilon, 0.001);
    const finer = voltage(segments, epsilon, 0.001, 0, finerInitial);
    return {epsilon, positiveInitial, negativeInitial, reset, positive, negative,
      positiveWriteCoefficient: positive[4] / (4 * epsilon * calibration),
      negativeWriteCoefficient: negative[4] / (4 * epsilon * calibration),
      joinedGap: gap(joined, positive), stepGap: gap(finer, positive)};
  });
  const firstOrderWriteCoefficient = leadingOrder([...prefix(1), ...segments])[3] / (4 * calibration);
  check("L zero-residual current pair has zero exact J gradient", norm(probe.r) === 0 && norm(update(matrix(), probeX, probeY, 1, "exact").flat()) === 0);
  check("L reachable opposite branch histories have identical zero initial soma", cases.every(x => x.positiveInitial[4] === 0 && x.negativeInitial[4] === 0 && norm(x.positiveInitial) > 0.01 && gap(x.positiveInitial, x.negativeInitial.map(v => -v)) < 1e-14));
  check("L original reset interface has no write for this zero-residual pair", cases.every(x => norm(x.reset) > 0.01 && x.reset[4] === 0));
  check("L removing reset admits nonzero opposite write coefficients", cases.every(x => Math.abs(x.positiveWriteCoefficient) > 0.01 && Math.abs(x.positiveWriteCoefficient + x.negativeWriteCoefficient) < 1e-12));
  check("L continuation uses the original uninterrupted current dynamics", cases.every(x => x.joinedGap < 1e-14 && x.stepGap < 1e-10));
  const relativeErrors = cases.map(x => Math.abs(x.positiveWriteCoefficient - firstOrderWriteCoefficient));
  check("L normalized history contribution remains nonzero as modulation shrinks", Math.abs(firstOrderWriteCoefficient) > 0.01 && relativeErrors.slice(1).every((error, i) => relativeErrors[i] / error > 3.9 && relativeErrors[i] / error < 4.1));
  const unitInitial = cases[1].positiveInitial;
  const scaledResponses = [0, 0.5, 1].map(scale => voltage(segments, 0.05, 0.002, 0, unitInitial.map(v => scale * v))[4]);
  check("L dependence on incoming state is affine, not integration noise", Math.abs(scaledResponses[1] - (scaledResponses[0] + scaledResponses[2]) / 2) < 1e-13);
  continuityQualification = {scope: "Interface kill only: no committed synaptic updates, no training, no new mechanism. Refutes removal of J reset while retaining its current-pair gradient interpretation, not history-dependent plasticity in general.", probeX, probeY, probe, prefixPositive: prefix(1), segments, calibration, cases, firstOrderWriteCoefficient, relativeErrors, scaledResponses};
}
let repairQualification;
if (repairKillOnly) {
  const prefix = sign => [{a: sign, b: 0, t: 0.2}, {a: 0, b: 0, t: 0.2}];
  const incoming = sign => pairs.map((_, edge) => edge === 0 ? voltage(prefix(sign), 0.05) : [0, 0, 0, 0, 0]);
  const zeroX = [0, 1, 0, 0], zeroY = [...zeroX];
  const zeroCases = [1, -1].map(sign => observePair(matrix(), incoming(sign), zeroX, zeroY, 1));
  check("M repair: original zero-error histories cannot trigger any synaptic write", zeroCases.every(x => norm(x.before.r) === 0 && norm(x.delta.flat()) === 0));
  check("M repair: history and nonzero soma are retained while writing remains zero", zeroCases.every(x => norm(x.nextStates.flat()) > 0.01 && Math.abs(x.atTeaching[0][4]) > 1e-5));
  const fullPhysical = voltage([...prefix(1), ...phaseProtocol(0, 1, 0, 0)], 0.05);
  check("M repair: zero-error event does not reset or alter the original physical path", gap(zeroCases[0].nextStates[0], fullPhysical) < 1e-14);
  const states = pairs.map(([i, j]) => voltage([{a: (i + 1) / 4, b: -(j + 1) / 4, t: 0.3}, {a: -(j + 1) / 4, b: (i + 1) / 4, t: 0.2}], 0.05));
  const target = [0.2, -0.3, 0.1, 0.25], eta = 0.5, epsilon = 0.05;
  const priorK = exampleK.map(row => [...row]), priorStates = states.map(v => [...v]);
  const event = plasticEvent(priorK, states, target, 1, epsilon, eta);
  check("M repair: event leaves physical states and incoming weights untouched", mgap(states, priorStates) === 0 && mgap(priorK, exampleK) === 0);
  check("M repair: plastic current equals the teaching-modulated component of the same RHS", states.every((v, edge) => {
    const [i, j] = pairs[edge], a = event.before.r[i], b = event.before.r[j];
    const on = starFlow(v, a, b, epsilon), off = starFlow(v, a, b, 0);
    return gap(event.currents[edge], on.slice(0, 4).map((value, index) => value - off[index])) < 1e-14;
  }));
  const algebra = matrix();
  pairs.forEach(([i, j], edge) => {
    const [hi, hj] = pairHistory(states[edge]);
    algebra[i][j] = eta * epsilon * (event.before.r[i] * hj - event.before.r[j] * hi) / 2;
    algebra[j][i] = -algebra[i][j];
  });
  check("M repair: raw local currents match the history-based pair gradient", mgap(algebra, event.delta) < 1e-14);
  const derivative = (historyLoss(addUpdate(exampleK, H, 1e-5), states, target, 1).loss - historyLoss(addUpdate(exampleK, H, -1e-5), states, target, 1).loss) / 2e-5;
  const gradientGap = Math.abs(derivative + dot(event.delta.flat(), H.flat()) / (eta * epsilon));
  check("M repair: matched-history gradient passes finite difference", gradientGap < 1e-9);
  const stepChecks = [0.25, 0.5, 0.75].map(gain => ({eta: gain, before: event.before.loss, after: plasticEvent(exampleK, states, target, 1, epsilon, gain).after.loss}));
  check("M repair: fixed small one-step gains reduce the declared frozen-history loss", stepChecks.every(x => x.after < x.before));
  check("M repair: every committed pair remains antisymmetric", mgap(event.nextK, tr(event.nextK).map(row => row.map(v => -v))) < 1e-14);
  const exactTarget = historyPrediction(exampleK, states, 1);
  check("M repair: a nonzero matched prediction with zero residual also writes zero", norm(exactTarget) > 0 && norm(plasticEvent(exampleK, states, exactTarget, 1).delta.flat()) === 0);
  const balancedStates = pairs.map(() => [0.2, 0.2, -0.1, -0.1, 0.3]);
  const balancedEvent = plasticEvent(matrix(), balancedStates, target, 1);
  check("M repair: nonzero teaching with zero branch contrast has zero write", norm(balancedEvent.before.r) > 0 && norm(balancedEvent.currents.flat()) > 0 && norm(balancedEvent.delta.flat()) === 0);
  const shiftedStates = states.map(v => [v[0] + 0.2, v[1] + 0.2, v[2] - 0.1, v[3] - 0.1, v[4] + 0.5]);
  check("M repair: instantaneous write ignores soma stock and branch common mode", mgap(event.delta, plasticEvent(exampleK, shiftedStates, target, 1).delta) < 1e-14);
  const reverse = plasticEvent(exampleK, states, target.map(v => -v), -1);
  check("M repair: reversing target and action at fixed history preserves the write", mgap(event.delta, reverse.delta) < 1e-14);
  const zeroModulation = plasticEvent(exampleK, states, target, 1, 0);
  check("M repair: zero modulation means zero current and zero write without division", norm(zeroModulation.delta.flat()) === 0 && norm(zeroModulation.currents.flat()) === 0);
  assert.throws(() => plasticEvent(matrix(), states, [10, 0, 0, 0], 1, 0.05), /strictly positive/);
  checks.push("M repair: invalid teaching conductance fails explicitly");
  // Counterexample against changing only the write source while retaining J's
  // instantaneous-m residual: the two mismatched updates have opposite signs.
  const wrongStates = pairs.map((_, edge) => edge === 0 ? voltage([{a: -1, b: 0, t: 1}, {a: 1, b: 0, t: 0.2}], 0.05) : [0, 0, 0, 0, 0]);
  const wrongK = matrix(); wrongK[0][1] = 0.3; wrongK[1][0] = -0.3;
  const zeroTarget = [0, 0, 0, 0], instant = [1, 0, 0, 0];
  const matched = plasticEvent(wrongK, wrongStates, zeroTarget, 1);
  const staleError = observedResidual(wrongK, instant, instant, 1).r;
  const mismatchedDelta = matrix();
  pairs.forEach(([i, j], edge) => {
    mismatchedDelta[i][j] = eta * teachingCurrent(wrongStates[edge], staleError[i], staleError[j], epsilon).reduce((s, v) => s + v, 0) / 4;
    mismatchedDelta[j][i] = -mismatchedDelta[i][j];
  });
  const mismatchedAfter = historyLoss(addUpdate(wrongK, mismatchedDelta, 1), wrongStates, zeroTarget, 1).loss;
  check("M repair: keeping the old instantaneous residual can increase the intended history loss", matched.delta[0][1] * mismatchedDelta[0][1] < 0 && matched.after.loss < matched.before.loss && mismatchedAfter > matched.before.loss);
  const nonzeroEvents = [1, -1].map(sign => plasticEvent(matrix(), incoming(sign), [0, 0.2, 0, 0], 1));
  check("M repair: eligible history is used for nonzero teaching, not erased or ignored", nonzeroEvents[0].delta[0][1] * nonzeroEvents[1].delta[0][1] < 0);
  const first = observePair(exampleK, states, exampleX, exampleY, 1);
  const second = observePair(first.nextK, first.nextStates, exampleY, exampleX, -1);
  const coldSecond = observePair(first.nextK, pairs.map(() => [0, 0, 0, 0, 0]), exampleY, exampleX, -1);
  check("M repair: two-call continuation preserves a measurable history effect", norm(first.delta.flat()) > 0 && norm(second.delta.flat()) > 0 && mgap(second.atTeaching, coldSecond.atTeaching) > 0.001);
  const sharper = observePair(exampleK, states, exampleX, exampleY, 1, epsilon, eta, 0.001);
  const eventStepGap = mgap(first.delta, sharper.delta);
  check("M repair: caller write and voltage agree under step refinement", eventStepGap < 1e-10 && mgap(first.nextStates, sharper.nextStates) < 1e-10);
  repairQualification = {scope: "Research interface repair only. Events retain five physical states per pair and write on present teaching currents. No training sweep/VPR. The predictor now uses pair-specific history; J's instantaneous Cayley association guarantee is not inherited.", epsilon, eta, effectiveGradientStep: epsilon * eta, zeroCases, event, gradientGap, stepChecks, mismatch: {history: pairHistory(wrongStates[0]), before: matched.before.loss, matchedAfter: matched.after.loss, mismatchedAfter, matchedDelta: matched.delta[0][1], mismatchedDelta: mismatchedDelta[0][1]}, nonzeroHistoryWrites: nonzeroEvents.map(x => x.delta[0][1]), first, second, eventStepGap};
}
// N: fixed-parameter reduction tests on the SAME starFlow/voltage path.
// Mode coordinates and no-soma feedback are diagnostic controls, not operators.
let noveltyQualification;
if (noveltyKillOnly) {
  const lambda = p.leak + p.axial, somaLambda = p.somaLeak + 4 * p.axial;
  const modes = v => [...pairHistory(v), (v[0] + v[1]) / 2, (v[2] + v[3]) / 2, v[4]];
  const physical = ([H, V, P, Q, S]) => [P + H, P - H, Q - V, Q + V, S];
  const modeFlow = ([H, V, P, Q, S], a, b, epsilon, feedback = true) => [
    -lambda * H - epsilon * b * P + a,
    -lambda * V + epsilon * a * Q + b,
    -lambda * P - epsilon * b * H + (feedback ? p.axial * S : 0),
    -lambda * Q + epsilon * a * V + (feedback ? p.axial * S : 0),
    -somaLambda * S + 2 * p.axial * (P + Q),
  ];
  const modeRun = (segments, epsilon, feedback = true, step = 0.002) => integrate(segments, [0, 0, 0, 0, 0], step, (z, a, b) => modeFlow(z, a, b, epsilon, feedback));
  const expansion = segments => integrate(segments, [0, 0, 0, 0, 0, 0, 0], 0.002, ([H0, V0, P1, Q1, S1, H2, V2], a, b) => [
    -lambda * H0 + a, -lambda * V0 + b,
    -lambda * P1 - b * H0 + p.axial * S1,
    -lambda * Q1 + a * V0 + p.axial * S1,
    -somaLambda * S1 + 2 * p.axial * (P1 + Q1),
    -lambda * H2 - b * P1, -lambda * V2 + a * Q1,
  ]);
  const fixtures = [
    {name: "existing_short_phase", segments: phaseProtocol(0.7, -0.4, 0.2, 0.6)},
    {name: "existing_four_unit_loop", segments: path},
  ];
  const scaling = fixtures.map(({name, segments}) => {
    const series = expansion(segments), ema = series.slice(0, 2);
    const cases = [0.1, 0.05, 0.025].map(epsilon => {
      const raw = voltage(segments, epsilon), z = modes(raw), equivalent = modeRun(segments, epsilon);
      const noFeedback = modeRun(segments, epsilon, false), refined = modes(voltage(segments, epsilon, 0.001));
      const approximated = ema.map((value, i) => value + epsilon ** 2 * series[5 + i]);
      const featureGap = gap(z.slice(0, 2), ema), feedbackGap = gap(z.slice(0, 2), noFeedback.slice(0, 2));
      return {epsilon, physical: raw, modes: z, ema, featureGap, relativeFeatureGap: featureGap / norm(ema), feedbackGap, relativeFeedbackGap: feedbackGap / norm(ema), equivalentGap: gap(z, equivalent), expansionRemainder: gap(z.slice(0, 2), approximated), stepGap: gap(z, refined), negativeEpsilonModes: modes(voltage(segments, -epsilon))};
    });
    return {name, segments, series, cases};
  });
  check("N reduction: mode transform is invertible on actual reachable states", scaling.every(f => f.cases.every(r => gap(physical(r.modes), r.physical) < 1e-14)));
  check("N reduction: equal-budget five-state bilinear SSM exactly matches the circuit", scaling.every(f => f.cases.every(r => r.equivalentGap < 1e-12)));
  check("N symmetry: eligibility even and common modes odd in epsilon", scaling.every(f => f.cases.every(r => gap(r.modes, r.negativeEpsilonModes.map((x, i) => i < 2 ? x : -x)) < 1e-12)));
  check("N expansion: leading non-EMA eligibility scales quadratically", scaling.every(f => f.cases.slice(1).every((r, i) => f.cases[i].featureGap / r.featureGap > 3.9 && f.cases[i].featureGap / r.featureGap < 4.1)));
  check("N expansion: retaining second order leaves fourth-order error", scaling.every(f => f.cases.slice(1).every((r, i) => f.cases[i].expansionRemainder / r.expansionRemainder > 14 && f.cases[i].expansionRemainder / r.expansionRemainder < 18)));
  check("N integration: all fixed scale cases survive step refinement", scaling.every(f => f.cases.every(r => r.stepGap < 1e-10)));
  // Constructed collision, not data selection: solve the known EMA endpoint
  // analytically. Both streams equal => all antisymmetric-area filters are zero.
  const epsilon = 0.05, eta = 0.5, pulse = 0.2;
  const bipolar = amplitude => [{a: 1, b: 1, t: pulse}, {a: -amplitude, b: -amplitude, t: pulse}, {a: 0, b: 0, t: pulse}];
  const amplitude = Math.exp(-lambda * pulse), collisionSegments = bipolar(amplitude);
  const collisionRaw = voltage(collisionSegments, epsilon), collision = modes(collisionRaw), collisionTrace = leadingOrder(collisionSegments);
  const collisionNoFeedback = modeRun(collisionSegments, epsilon, false);
  check("N collision: analytic bipolar history cancels both EMA states and all leading-area states", norm(collisionTrace) < 1e-12);
  check("N collision: same EMA-area summary does not determine finite-shunt eligibility", norm(collision.slice(0, 2)) > 1e-9 && Math.abs(collision[4]) < 1e-14);
  check("N attribution kill: this EMA separation needs no soma feedback", gap(collision, collisionNoFeedback) < 1e-12);
  // A stronger reachable witness: find a pulse that also cancels the actual H,V.
  // Root finding constructs an equality constraint; it is not a trained parameter.
  let lo = 0, hi = 1;
  check("N witness: reachable zero-readout endpoint root is bracketed", pairHistory(voltage(bipolar(lo), epsilon))[0] > 0 && pairHistory(voltage(bipolar(hi), epsilon))[0] < 0);
  for (let iteration = 0; iteration < 55; iteration++) {
    const mid = (lo + hi) / 2;
    if (pairHistory(voltage(bipolar(mid), epsilon))[0] > 0) lo = mid; else hi = mid;
  }
  const rootAmplitude = (lo + hi) / 2, rootSegments = bipolar(rootAmplitude);
  const checkpointRaw = voltage(rootSegments, epsilon), checkpoint = modes(checkpointRaw);
  const probe = [{a: 1, b: 0, t: pulse}];
  const futureRaw = voltage(probe, epsilon, 0.002, 0, checkpointRaw), coldFutureRaw = voltage(probe, epsilon);
  const future = modes(futureRaw), coldFuture = modes(coldFutureRaw);
  const futureNoFeedback = modeRun([...rootSegments, ...probe], epsilon, false);
  const eventFor = raw => plasticEvent(matrix(), pairs.map((_, e) => e === 0 ? raw : [0, 0, 0, 0, 0]), [1, 0, 0, 0], 1, epsilon, eta);
  const futureWrite = eventFor(futureRaw).delta[0][1], coldFutureWrite = eventFor(coldFutureRaw).delta[0][1];
  check("N witness: identical readable history and soma can hide reachable common modes", norm([checkpoint[0], checkpoint[1], checkpoint[4]]) < 1e-12 && norm(checkpoint.slice(2, 4)) > 1e-5);
  check("N witness: same subsequent input reveals hidden common-mode state", gap(future.slice(0, 2), coldFuture.slice(0, 2)) > 1e-8);
  check("N witness: same K target and residual yield a different subsequent write", Math.abs(futureWrite - coldFutureWrite) > 1e-10);
  check("N witness: main separation also survives without soma feedback", gap(futureNoFeedback.slice(0, 2), coldFuture.slice(0, 2)) > 1e-8);
  check("N witness: refined reachable state and future agree", gap(future, modes(voltage([...rootSegments, ...probe], epsilon, 0.001))) < 1e-10);
  noveltyQualification = {
    scope: "Reduction and fixed arithmetic only, not proof of method novelty, biological plausibility or VPR utility. Same primary starFlow and plasticEvent. No training or new production path.",
    parameters: {epsilon, eta, pulse, lambda, somaLambda}, scaling,
    emaAreaCollision: {amplitude, segments: collisionSegments, baselineFourStates: collisionTrace, actualModes: collision, noSomaFeedbackModes: collisionNoFeedback, exactCurrentWrite: eventFor(collisionRaw).delta[0][1]},
    hiddenStateWitness: {rootAmplitude, rootSegments, checkpoint, probe, future, coldFuture, noSomaFeedbackFuture: futureNoFeedback, futureWrite, coldFutureWrite},
    verdict: "EMA is not exact, but this is insufficient to rescue novelty. Constructed memory witnesses do not require soma feedback; an equal-budget bilinear SSM is exactly equivalent. Optimizer remains LMS. No causal VPR/association value established.",
  };
}
// P: endpoint sensitivity of intrinsic tonic conductances, not a new operator.
// The 20 sensitivity scalars are a diagnostic oracle, NOT claimed local biology.
let visualCreditQualification, associationCreditQualification;
if (visualCreditKillOnly) {
  const theta = [0.15, -0.1, 0.08, -0.05], zeroTheta = [0, 0, 0, 0];
  const segments = [{a: 0.8, b: -0.3, t: 0.3}, {a: -0.2, b: 0.7, t: 0.4}];
  const epsilon = 0.05, initial = [0.1, -0.05, 0.03, -0.02, 0.01];
  const conductances = th => th.map(x => p.leak / 2 * Math.exp(x));
  function endpoint(th, segs = segments, z0 = initial, eps = epsilon, sensitivity = true, transport = true, step = 0.002) {
    const g = conductances(th);
    const state0 = sensitivity ? [...z0, ...Array(20).fill(0)] : z0;
    return integrate(segs, state0, step, (state, a, b) => {
      const z = state.slice(0, 5), dz = starFlow(z, a, b, eps, p.axial, 0, g);
      if (!sensitivity) return dz;
      const A = flowMatrix(a, b, eps).map((row, i) => row.map((v, j) => (transport || i === j ? v : 0) - (i === j && i < 4 ? g[i] - p.leak / 2 : 0)));
      const dE = Array.from({length: 20}, (_, k) => {
        const i = Math.floor(k / 4), j = k % 4;
        return A[i].reduce((s, v, n) => s + v * state[5 + 4 * n + j], 0) - (i === j ? g[j] * z[j] : 0);
      });
      return [...dz, ...dE];
    });
  }
  const baseG = conductances(zeroTheta);
  check("P continuity: baseline conductances recover the primary current exactly", gap(starFlow(initial, 0.8, -0.3, epsilon), starFlow(initial, 0.8, -0.3, epsilon, p.axial, 0, baseG)) === 0);
  const exact = endpoint(theta), z = exact.slice(0, 5);
  const E = Array.from({length: 5}, (_, i) => exact.slice(5 + 4 * i, 9 + 4 * i));
  const delta = 1e-5;
  const fdColumns = theta.map((_, j) => {
    const plus = [...theta], minus = [...theta]; plus[j] += delta; minus[j] -= delta;
    const zp = endpoint(plus, segments, initial, epsilon, false), zm = endpoint(minus, segments, initial, epsilon, false);
    return zp.map((v, i) => (v - zm[i]) / (2 * delta));
  });
  const derivativeGap = mgap(E, tr(fdColumns));
  check("P sensitivity: full transported eligibility matches finite differences", derivativeGap < 1e-9);
  const approximate = endpoint(theta, segments, initial, epsilon, true, false);
  const eligibilityGap = gap(exact.slice(5), approximate.slice(5));
  check("P locality kill: deleting cross-compartment sensitivity transport is not exact", eligibilityGap > 1e-5);
  const read = raw => pairHistory(raw.slice(0, 5));
  const readE = raw => [Array.from({length: 4}, (_, j) => (raw[5 + j] - raw[9 + j]) / 2), Array.from({length: 4}, (_, j) => -(raw[13 + j] - raw[17 + j]) / 2)];
  const target = [0.15, -0.1], prediction = read(exact), error = prediction.map((v, i) => v - target[i]);
  const gradient = tr(readE(exact)).map(col => dot(col, error));
  const loss = th => { const r = read(endpoint(th)).map((v, i) => v - target[i]); return dot(r, r) / 2; };
  const lossGradientFD = theta.map((_, j) => { const plus = [...theta], minus = [...theta]; plus[j] += delta; minus[j] -= delta; return (loss(plus) - loss(minus)) / (2 * delta); });
  check("P loss: fixed-start endpoint gradient matches finite differences", gap(gradient, lossGradientFD) < 1e-9);
  const eta = 0.1, nextTheta = theta.map((v, j) => v - eta * gradient[j]);
  check("P event: one diagnostic endpoint step decreases its own prediction loss", loss(nextTheta) < loss(theta));
  check("P teaching: zero current residual gives zero conductance write despite history", tr(readE(exact)).every(col => dot(col, [0, 0]) === 0) && norm(z) > 0);
  // Same parameter version and same cold read protocol for query AND reference.
  const cold = [0, 0, 0, 0, 0], q = [{a: 0.8, b: 0.2, t: 0.4}], ref = [{a: 0.3, b: 0.7, t: 0.4}];
  const cosine = (x, y) => { const nx = Math.sqrt(dot(x, x)), ny = Math.sqrt(dot(y, y)); assert.ok(nx > 0 && ny > 0, "undefined descriptor"); return dot(x, y) / (nx * ny); };
  function scoreAndGradient(th, query = q, reference = ref, eps = epsilon) {
    const u = endpoint(th, query, cold, eps), v = endpoint(th, reference, cold, eps), x = read(u), y = read(v), ex = readE(u), ey = readE(v);
    const nx = Math.sqrt(dot(x, x)), ny = Math.sqrt(dot(y, y)), score = cosine(x, y);
    const dx = y.map((val, i) => val / (nx * ny) - score * x[i] / (nx * nx));
    const dy = x.map((val, i) => val / (nx * ny) - score * y[i] / (ny * ny));
    return {score, gradient: theta.map((_, j) => dot(dx, ex.map(row => row[j])) + dot(dy, ey.map(row => row[j])))};
  }
  const sg = scoreAndGradient(theta);
  const scoreFD = theta.map((_, j) => { const plus = [...theta], minus = [...theta]; plus[j] += delta; minus[j] -= delta; return (scoreAndGradient(plus).score - scoreAndGradient(minus).score) / (2 * delta); });
  check("P readout: shared query-reference score has nonzero first-order parameter sensitivity", norm(sg.gradient) > 1e-4);
  check("P readout: normalized score includes BOTH query and reference derivatives", gap(sg.gradient, scoreFD) < 1e-9);
  const lin = (a, b) => read(endpoint(theta, [{a, b, t: 0.4}], cold, 0));
  check("P novelty kill: no input modulation reduces cold finite-time readout to a linear feature map", gap(lin(0.7, -0.2), lin(1, 0).map((v, i) => 0.7 * v - 0.2 * lin(0, 1)[i])) < 1e-12);
  check("P integration: endpoint and sensitivities survive step refinement", gap(exact, endpoint(theta, segments, initial, epsilon, true, true, 0.001)) < 1e-10);
  assert.throws(() => endpoint([-20, 0, 0, 0]), /conductance/);
  checks.push("P safety: invalid positive-but-too-small synaptic conductance fails explicitly");
  visualCreditQualification = {scope: "Intrinsic-conductance interface and conditional finite-horizon sensitivity checks only. Inputs are fixed numerical feature fixtures, not RGB/HM3D. No online learner, local-plasticity claim, retrieval-performance claim or novelty certification.", theta, epsilon, segments, initial, target, derivativeGap, eligibilityGap, gradient, lossBefore: loss(theta), lossAfter: loss(nextTheta), score: sg, scoreFD, extraSensitivityScalarsPerUnit: 20, physicalStatesPerUnit: 5, verdict: "Parameters now affect a shared normalized retrieval interface. Exact endpoint credit requires transported sensitivity; independent local traces are insufficient. Linear limit is an ordinary learned feature metric. Semantic association supervision and a novel implementable local credit rule remain open."};
  // Q: evidence-conditioned association credit. Likelihoods below are declared
  // arithmetic fixtures, NOT supplied by an implemented visual predictor.
  if (associationCreditKillOnly) {
    const sigmoid = f => 1 / (1 + Math.exp(-f));
    const alternate = [{a: 0.7, b: -0.2, t: 0.4}];
    function margin(th, query = q, first = ref, second = alternate, eps = epsilon) {
      const a = scoreAndGradient(th, query, first, eps), b = scoreAndGradient(th, query, second, eps);
      return {value: a.score - b.score, gradient: a.gradient.map((v, j) => v - b.gradient[j])};
    }
    const current = margin(theta), prior = sigmoid(current.value);
    const meanA = [0.15, -0.1], meanB = [-0.15, 0.1], observation = [...meanA], variance = 1;
    const squared = (x, y) => dot(x.map((v, i) => v - y[i]), x.map((v, i) => v - y[i]));
    const logDensity = (y, mu) => -squared(y, mu) / (2 * variance);
    const logA = logDensity(observation, meanA), logB = logDensity(observation, meanB);
    const logRatio = logA - logB, contrast = meanA.map((v, i) => v - meanB[i]);
    const midResidual = observation.map((v, i) => v - (meanA[i] + meanB[i]) / 2);
    check("Q evidence: Gaussian loss difference equals contrast-projected innovation", Math.abs(logRatio - dot(contrast, midResidual) / variance) < 1e-14);
    const posterior = sigmoid(current.value + logRatio), teaching = posterior - prior;
    check("Q evidence: identical predictive distributions produce no association write", sigmoid(current.value + logA - logA) === prior);
    const orthogonal = [contrast[1], -contrast[0]], shifted = observation.map((v, i) => v + orthogonal[i]);
    check("Q evidence: orthogonal innovation cancels, not arbitrary nuisance", Math.abs(logDensity(shifted, meanA) - logDensity(shifted, meanB) - logRatio) < 1e-14);
    const likelihoodLoss = th => {
      const prob = sigmoid(margin(th).value), pivot = Math.max(logA, logB);
      return -pivot - Math.log(prob * Math.exp(logA - pivot) + (1 - prob) * Math.exp(logB - pivot));
    };
    const evidenceGradient = current.gradient.map(v => -teaching * v);
    const evidenceFD = theta.map((_, j) => { const plus = [...theta], minus = [...theta]; plus[j] += delta; minus[j] -= delta; return (likelihoodLoss(plus) - likelihoodLoss(minus)) / (2 * delta); });
    const evidenceGap = gap(evidenceGradient, evidenceFD);
    check("Q conditional credit: cached-likelihood marginal gradient matches finite differences", evidenceGap < 1e-9);
    const updated = theta.map((v, j) => v + eta * teaching * current.gradient[j]);
    const ownChange = margin(updated).value - current.value;
    check("Q direction: supportive evidence increases this pair's association margin", teaching > 0 && ownChange > 0);
    check("Q objective: diagnostic step lowers cached-likelihood loss", likelihoodLoss(updated) < likelihoodLoss(theta));
    check("Q own response: first-order change uses the squared association sensitivity", Math.abs(ownChange - eta * teaching * dot(current.gradient, current.gradient)) < 1e-9);
    // Calibrated two-outcome law: exact finite sums, no sampled success labels.
    const densityA = [0.8, 0.2], densityB = [0.2, 0.8];
    const mixture = densityA.map((v, i) => prior * v + (1 - prior) * densityB[i]);
    const posteriors = densityA.map((v, i) => prior * v / mixture[i]);
    const shifts = posteriors.map(v => v - prior);
    const expectationA = dot(densityA, shifts), expectationB = dot(densityB, shifts);
    const posteriorVariance = dot(mixture, shifts.map(v => v * v));
    check("Q calibration: expected posterior teaching has correct class-conditional sign", expectationA > 0 && expectationB < 0);
    check("Q calibration: expectation equals posterior variance divided by prior mass", Math.abs(expectationA - posteriorVariance / prior) < 1e-14 && Math.abs(expectationB + posteriorVariance / (1 - prior)) < 1e-14);
    check("Q calibration: pooled expected teaching is zero, not evidence-free learning", Math.abs(dot(mixture, shifts)) < 1e-14);
    const swappedShifts = densityB.map((v, i) => sigmoid(current.value + Math.log(v / densityA[i])) - prior);
    check("Q teacher kill: swapped predictive models can reverse expected teaching", dot(densityA, swappedShifts) < 0);
    // Fixed contexts, declared before checking signs; no pose/place labels.
    const queryFixtures = [q, [{a: 0.2, b: 0.8, t: 0.4}], [{a: -0.4, b: 0.6, t: 0.4}], [{a: -0.7, b: -0.3, t: 0.4}], [{a: 0.6, b: -0.5, t: 0.4}], [{a: 0.1, b: -0.9, t: 0.4}]];
    const margins = queryFixtures.map(query => margin(theta, query));
    const responses = margins.map((m, i) => {
      const kernel = dot(m.gradient, current.gradient);
      const plus = theta.map((v, j) => v + delta * current.gradient[j]);
      const minus = theta.map((v, j) => v - delta * current.gradient[j]);
      const directionalFD = (margin(plus, queryFixtures[i]).value - margin(minus, queryFixtures[i]).value) / (2 * delta);
      return {query: queryFixtures[i], kernel, directionalFD, observedChange: margin(updated, queryFixtures[i]).value - m.value, predictedChange: eta * teaching * kernel};
    });
    check("Q transfer: cross-context sensitivity kernel matches directional finite differences", responses.every(r => Math.abs(r.kernel - r.directionalFD) < 1e-9));
    check("Q transfer: same update has the derived first-order cross-context response", responses.every(r => Math.abs(r.observedChange - r.predictedChange) < 1e-9));
    check("Q interference kill: fixed contexts include opposite signed score responses", responses.some(r => r.kernel > 1e-6) && responses.some(r => r.kernel < -1e-6));
    const normalized = margins.map(m => { const length = Math.sqrt(dot(m.gradient, m.gradient)); assert.ok(length > 0, "zero association sensitivity"); return m.gradient.map(v => v / length); });
    const gram = mm(normalized, tr(normalized));
    const count = gram.length, dimension = theta.length;
    const offDiagonalSquareSum = gram.reduce((sum, row, i) => sum + row.reduce((s, v, j) => s + (i === j ? 0 : v * v), 0), 0);
    const rankBound = count * (count - dimension) / dimension;
    check("Q capacity kill: normalized cross-talk obeys the finite-parameter Gram bound", offDiagonalSquareSum >= rankBound - 1e-10);
    check("Q capacity: each association sensitivity has exactly four parameter coordinates", normalized.every(v => v.length === 4) && gram.every((row, i) => Math.abs(row[i] - 1) < 1e-12));
    // Exact symmetric linear-limit kill: coupling remains ON, but all cold,
    // constant-input cosine-margin gradients lie on one gain-contrast axis.
    const axis = [1, 1, -1, -1];
    const linearMargins = queryFixtures.map(query => margin(zeroTheta, query, ref, alternate, 0));
    const linearAxisGap = Math.max(...linearMargins.map(m => gap(m.gradient, axis.map(v => v * dot(m.gradient, axis) / 4))));
    check("Q linear-limit kill: association sensitivities are rank one despite soma coupling", linearAxisGap < 1e-12 && linearMargins.every(m => norm(m.gradient) > 1e-6));
    const duration = 0.4, decay = p.leak + p.axial;
    const integral = (1 - (1 + decay * duration) * Math.exp(-decay * duration)) / (decay * decay);
    const gainDerivative = -p.leak / 4 * integral;
    const linearEA = readE(endpoint(zeroTheta, [{a: 1, b: 0, t: duration}], cold, 0));
    const linearEB = readE(endpoint(zeroTheta, [{a: 0, b: 1, t: duration}], cold, 0));
    const linearDerivativeGap = Math.max(mgap(linearEA, [[gainDerivative, gainDerivative, 0, 0], [0, 0, 0, 0]]), mgap(linearEB, [[0, 0, 0, 0], [0, 0, gainDerivative, gainDerivative]]));
    check("Q linear-limit derivation: finite-time gain derivatives match the closed form", linearDerivativeGap < 1e-11);
    const effectiveRank = count * count / (count + offDiagonalSquareSum);
    associationCreditQualification = {scope: "Conditional evidence-to-association algebra on fixed numerical contexts; no calibrated visual predictor, place labels, real retrieval experiment, new local rule or novelty claim.", likelihoodsFrozenDuringUpdate: true, meanA, meanB, observation, variance, current, prior, logRatio, posterior, teaching, evidenceGap, likelihoodLossBefore: likelihoodLoss(theta), likelihoodLossAfter: likelihoodLoss(updated), ownChange, calibration: {densityA, densityB, expectationA, expectationB, posteriorVariance, wrongModelExpectationA: dot(densityA, swappedShifts)}, crossContextResponses: responses, normalizedGram: gram, offDiagonalSquareSum, rankBound, effectiveRank, symmetricLinearLimit: {axis, linearAxisGap, gainDerivative, linearDerivativeGap, axial: p.axial}, verdict: "Evidence-conditioned credit can align the updated pair under frozen calibrated likelihood assumptions. Independent-context transfer still depends on signed sensitivity overlap; both signs occur. The symmetric linear limit is rank one even with soma coupling. Compartment coupling gives no automatic association-retention guarantee. The evidence predictor and a novel low-cost local realization remain unimplemented."};
  }
}
// R: constructive current redesign, on the SAME starFlow. The implicit token
// step is a changed discrete cell, not an exact replacement for the old RK4 ODE.
let constructiveCouplingQualification;
let structureComputationQualification;
let selectiveCorrectionQualification;
if (constructiveCouplingOnly) {
  const h = 0.4, epsilon = 0, theta = [0, 0, 0, 0], cold = [0, 0, 0, 0, 0];
  const connectionSigns = [1, -1, -1, 1];
  const read = state => pairHistory(state);
  const lift = b => [b[0] / 2, -b[0] / 2, -b[1] / 2, b[1] / 2, 0];
  // Reuse the shared five-state implicit dynamics and arrowhead response.  The
  // legacy closure remains only as names/parameters for byte-stable checks.
  const response = sharedResponse;
  const event = sharedImplicitCell;
  function localCredit(cell, stateError, details = false, audit = null) {
    const errorVoltage = response(cell.diagonal, cell.edge, stateError);
    if (audit !== null) {
      audit.responseCalls += 1;
      audit.creditResponseCalls += 1;
    }
    const gradient = cell.weights === undefined
      ? cell.coefficient.map((g, j) => cell.carrier === "axial"
        ? -h * g * cell.drop[j] ** 3 * (errorVoltage[j] - errorVoltage[4])
        : -h * g * cell.state[j] * errorVoltage[j])
      : cell.coefficient.map((g, j) => cell.carrier === "axial"
        ? -h * g * cell.weights[j] * cell.drop[j] ** 3 * (errorVoltage[j] - errorVoltage[4])
        : -h * g * cell.weights[j] * cell.state[j] * errorVoltage[j]);
    return details ? {gradient, errorVoltage} : gradient;
  }
  const delta = 1e-5, input = [0.8, -0.3], incoming = [0.1, -0.05, 0.03, -0.02, 0.01];
  const auditEpsilon = 0.05, auditTheta = [0.15, -0.1, 0.08, -0.05], cell = event(auditTheta, input, "axial", incoming, 1e-12, auditEpsilon);
  const stateError = [0.3, -0.2, 0.15, 0.1, -0.4];
  const gradient = localCredit(cell, stateError);
  const fd = auditTheta.map((_, j) => {
    const plus = [...auditTheta], minus = [...auditTheta]; plus[j] += delta; minus[j] -= delta;
    return (dot(stateError, event(plus, input, "axial", incoming, 1e-12, auditEpsilon).state) - dot(stateError, event(minus, input, "axial", incoming, 1e-12, auditEpsilon).state)) / (2 * delta);
  });
  const derivativeGap = gap(gradient, fd);
  check("R construction: one-soma local error response gives exact implicit-token gradient", derivativeGap < 1e-9);
  const zeroCubic = starFlow(incoming, ...input, epsilon, p.axial, 0, null, [0, 0, 0, 0]);
  check("R continuity: zero nonlinear current recovers the primary flow exactly", gap(zeroCubic, starFlow(incoming, ...input, epsilon)) === 0);
  const flow = starFlow(incoming, ...input, epsilon, p.axial, 0, null, cell.coefficient);
  const oldFlow = starFlow(incoming, ...input, epsilon), extra = flow.map((v, j) => v - oldFlow[j]);
  const expectedDissipation = -cell.coefficient.reduce((sum, g, j) => sum + g * (incoming[j] - incoming[4]) ** 4, 0);
  check("R conservation: nonlinear plastic current only transfers charge inside the cell", Math.abs(extra.reduce((sum, v) => sum + v, 0)) < 1e-14);
  check("R stability: added current has exactly nonpositive quartic energy contribution", Math.abs(dot(incoming, extra) - expectedDissipation) < 1e-14 && expectedDissipation < 0);
  check("R numerics: implicit state survives stricter convergence tolerance", gap(cell.state, event(auditTheta, input, "axial", incoming, 1e-14, auditEpsilon).state) < 1e-10);
  const shared = [0.3, 0.3, 0.3, 0.3, 0.3];
  check("R common mode: equal branch-soma voltages do not activate the new plastic current", gap(starFlow(shared, ...input, epsilon, p.axial, 0, null, cell.coefficient), starFlow(shared, ...input, epsilon)) === 0);
  function score(th, q, r, carrier, qStart = cold, rStart = cold, eps = epsilon) {
    const qc = event(th, q, carrier, qStart, 1e-12, eps), rc = event(th, r, carrier, rStart, 1e-12, eps);
    const x = read(qc.state), y = read(rc.state), nx = Math.sqrt(dot(x, x)), ny = Math.sqrt(dot(y, y));
    assert.ok(nx > 0 && ny > 0, "undefined normalized descriptor");
    const value = dot(x, y) / (nx * ny);
    const dx = y.map((v, i) => v / (nx * ny) - value * x[i] / (nx * nx));
    const dy = x.map((v, i) => v / (nx * ny) - value * y[i] / (ny * ny));
    const qCredit = localCredit(qc, lift(dx)), rCredit = localCredit(rc, lift(dy));
    return {value, gradient: qCredit.map((v, j) => v + rCredit[j])};
  }
  const referenceA = [0.3, 0.7], referenceB = [0.7, -0.2];
  const queries = [[0.8, 0.2], [0.2, 0.8], [-0.4, 0.6], [-0.7, -0.3], [0.6, -0.5], [0.1, -0.9]];
  function margin(th, q, carrier, start = cold, eps = epsilon) {
    const a = score(th, q, referenceA, carrier, start, cold, eps), b = score(th, q, referenceB, carrier, start, cold, eps);
    return {value: a.value - b.value, gradient: a.gradient.map((v, j) => v - b.gradient[j])};
  }
  function geometry(carrier) {
    const margins = queries.map(q => margin(theta, q, carrier));
    const directions = margins.map(m => { const length = Math.sqrt(dot(m.gradient, m.gradient)); assert.ok(length > 0); return m.gradient.map(v => v / length); });
    const gram = mm(directions, tr(directions)), squareSum = gram.flat().reduce((sum, v) => sum + v * v, 0);
    return {margins, gram, effectiveRank: queries.length ** 2 / squareSum};
  }
  const tonicGeometry = geometry("tonic"), axialGeometry = geometry("axial");
  const axialMargin = margin(auditTheta, queries[0], "axial");
  const marginFD = auditTheta.map((_, j) => {
    const plus = [...auditTheta], minus = [...auditTheta]; plus[j] += delta; minus[j] -= delta;
    return (margin(plus, queries[0], "axial").value - margin(minus, queries[0], "axial").value) / (2 * delta);
  });
  check("R association: exact local credit includes both normalized query and gallery paths", gap(axialMargin.gradient, marginFD) < 1e-9);
  check("R reference: the implicit tonic carrier retains the proved rank-one limit", Math.abs(tonicGeometry.effectiveRank - 1) < 1e-12);
  check("R constructive separation: nonlinear axial carrier breaks that rank-one limit", axialGeometry.effectiveRank > 1 + 1e-5);
  // Two previously fixed contexts receive contradictory evidence. Normalize
  // BOTH models' writes identically to isolate directions from gain magnitude.
  // This is an arithmetic learning assay, not a legal-visual-teacher claim.
  function twoWrites(carrier) {
    const learningRate = 0.02, evidenceMagnitude = 0.05, step = learningRate * evidenceMagnitude;
    const before = queries.slice(0, 2).map(q => margin(theta, q, carrier).value);
    const first = margin(theta, queries[0], carrier), firstLength = Math.sqrt(dot(first.gradient, first.gradient));
    const middle = theta.map((v, j) => v + step * first.gradient[j] / firstLength);
    const second = margin(middle, queries[1], carrier), secondLength = Math.sqrt(dot(second.gradient, second.gradient));
    const afterTheta = middle.map((v, j) => v - step * second.gradient[j] / secondLength);
    const after = queries.slice(0, 2).map(q => margin(afterTheta, q, carrier).value);
    return {learningRate, evidenceMagnitude, step, before, after, change: after.map((v, i) => v - before[i]), afterTheta};
  }
  const writeComparison = {tonic: twoWrites("tonic"), axial: twoWrites("axial")};
  check("R constructive learning: two opposite evidence writes separate the two fixed associations", writeComparison.axial.change[0] > 0 && writeComparison.axial.change[1] < 0);
  // Actual reachable prefixes under each carrier; same final observation.
  const prefixes = [[[1, 0], [0, 1]], [[0, 1], [1, 0]]];
  const historyCredit = prefixes.map(prefix => {
    let state = [...cold];
    for (const x of prefix) state = event(theta, x, "axial", state, 1e-12, auditEpsilon).state;
    return {prefix, incoming: state, margin: margin(theta, queries[0], "axial", state, auditEpsilon)};
  });
  check("R temporal interface: reachable different histories alter credit for the same current view", gap(historyCredit[0].margin.gradient, historyCredit[1].margin.gradient) > 1e-6);
  check("R temporal interface: original input modulation produces a nonzero soma context", historyCredit.every(item => Math.abs(item.incoming[4]) > 1e-7));
  const noSomaResponse = cell.coefficient.map((g, j) => -h * g * cell.drop[j] ** 3 * stateError[j] / cell.diagonal[j]);
  check("R local circuit: the computed soma error response is needed for this exact gradient", gap(noSomaResponse, fd) > 1e-6);
  check("R teaching: zero association evidence makes zero parameter write", gradient.every(v => 0 * v === 0));
  constructiveCouplingQualification = {scope: "Constructive nonlinear current and implicit-token local learning checks, not RGB navigation evidence or a novelty certificate.", h, epsilon, auditEpsilon, theta, input, incoming, derivativeGap, implicitIterations: cell.iterations, gradient, fd, noSomaResponse, expectedDissipation, referenceA, referenceB, queries, tonicGeometry, axialGeometry, writeComparison, historyCredit, physicalStates: 5, plasticParameters: 4, forwardSensitivityStatesForThisConditionalTokenGradient: 0, temporaryErrorCoordinatesPerLoss: 5, localDifferentialGainCoordinates: 4, inferenceCost: "Nonlinear implicit solve per token; linear-time star response per scalar loss. Different integrator/event definition from P, so no free speedup claim.", verdict: "The plastic carrier has been changed from tonic self-shunting to voltage-contrast cubic intercompartment transfer. Exact conditional token credit uses a one-soma Schur response, not 20 forward sensitivity states. Separation, teaching reliability and independent visual generalization must be assessed separately."};
  if (structureComputationOnly) {
    // S: the only new path is a fixed, synthetic computation fixture.  It
    // reuses R's event/response/localCredit machinery and never reads the old
    // JSON outputs or any visual/semantic teacher.
    const structureEpsilon = 0.05, structureTolerance = 1e-12, strictTolerance = 1e-14;
    const finiteDifferenceStep = 1e-5, finiteDifferenceStepFine = 5e-6;
    const mainEta = 0.02, rho = 0.001, rhoHalf = rho / 2;
    const theta0 = [0, 0, 0, 0];
    const validationChecks = [];
    const structureCheck = (name, passed, kind, details = {}) => {
      validationChecks.push({name, passed: Boolean(passed), kind, ...details});
      return Boolean(passed);
    };
    const vectorAdd = (a, b) => a.map((v, i) => v + b[i]);
    const vectorSub = (a, b) => a.map((v, i) => v - b[i]);
    const vectorScale = (a, c) => a.map(v => c * v);
    const vectorEqual = (a, b, tolerance = 0) => a.length === b.length && a.every((v, i) => Math.abs(v - b[i]) <= tolerance);
    const finiteVector = value => Array.isArray(value) && value.every(v => Number.isFinite(v));
    const sigmoid = value => 1 / (1 + Math.exp(-value));
    const logSumExp = values => {
      const maximum = Math.max(...values);
      return maximum + Math.log(values.reduce((sum, value) => sum + Math.exp(value - maximum), 0));
    };
    const numericScalarCount = value => {
      if (typeof value === "number") return 1;
      if (Array.isArray(value)) return value.reduce((sum, item) => sum + numericScalarCount(item), 0);
      if (value !== null && typeof value === "object") return Object.values(value).reduce((sum, item) => sum + numericScalarCount(item), 0);
      return 0;
    };
    const freshAudit = () => ({forwardCalls: 0, responseCalls: 0, newtonResponseCalls: 0, creditResponseCalls: 0, residualEvaluations: 0, newtonIterations: 0, backtrackSteps: 0, scoreCalls: 0});
    const mergeAudits = audits => audits.reduce((total, audit) => {
      for (const key of Object.keys(total)) total[key] += audit[key] || 0;
      return total;
    }, freshAudit());
    const describeCell = value => ({
      Zplus: value.state,
      drop: value.drop,
      alpha: value.coefficient,
      differentialConductance: value.edge.map(edge => edge / h),
      diagonal: value.diagonal,
      edge: value.edge,
      carrier: value.carrier,
      iterations: value.iterations,
      residual: value.residual,
    });
    const denseSolve = (matrixInput, rhsInput) => {
      const size = rhsInput.length, rows = matrixInput.map((row, i) => [...row, rhsInput[i]]);
      for (let column = 0; column < size; column++) {
        let pivot = column;
        for (let row = column + 1; row < size; row++) if (Math.abs(rows[row][column]) > Math.abs(rows[pivot][column])) pivot = row;
        if (Math.abs(rows[pivot][column]) < 1e-14) throw new Error("singular dense response reference");
        [rows[column], rows[pivot]] = [rows[pivot], rows[column]];
        const scale = rows[column][column];
        rows[column] = rows[column].map(value => value / scale);
        for (let row = 0; row < size; row++) if (row !== column) {
          const factor = rows[row][column];
          rows[row] = rows[row].map((value, index) => value - factor * rows[column][index]);
        }
      }
      return rows.map(row => row[size]);
    };
    const denseResponse = (cell, rhs) => {
      const matrixInput = Array.from({length: 5}, () => Array(5).fill(0));
      for (let j = 0; j < 4; j++) {
        matrixInput[j][j] = cell.diagonal[j];
        matrixInput[j][4] = -cell.edge[j];
        matrixInput[4][j] = -cell.edge[j];
      }
      matrixInput[4][4] = cell.diagonal[4];
      return denseSolve(matrixInput, rhs);
    };
    const detachedCredit = (cell, stateError) => {
      assert.equal(cell.carrier, "axial", "detached credit is only defined for the axial carrier");
      return cell.coefficient.map((alpha, j) => {
        const weight = cell.weights === undefined ? 1 : cell.weights[j];
        return -h * alpha * weight * cell.drop[j] ** 3 * stateError[j] / cell.diagonal[j];
      });
    };
    const scoreFromCells = (queryCell, referenceCell, creditMode, audit = null) => {
      if (audit !== null) audit.scoreCalls += 1;
      const queryDescriptor = read(queryCell.state), referenceDescriptor = read(referenceCell.state);
      const queryNorm = Math.sqrt(dot(queryDescriptor, queryDescriptor)), referenceNorm = Math.sqrt(dot(referenceDescriptor, referenceDescriptor));
      assert.ok(queryNorm > 0 && referenceNorm > 0, "undefined normalized descriptor");
      const value = dot(queryDescriptor, referenceDescriptor) / (queryNorm * referenceNorm);
      const queryDerivative = referenceDescriptor.map((v, i) => v / (queryNorm * referenceNorm) - value * queryDescriptor[i] / (queryNorm * queryNorm));
      const referenceDerivative = queryDescriptor.map((v, i) => v / (queryNorm * referenceNorm) - value * referenceDescriptor[i] / (referenceNorm * referenceNorm));
      const queryB = lift(queryDerivative), referenceB = lift(referenceDerivative);
      const queryFull = localCredit(queryCell, queryB, true, audit), referenceFull = localCredit(referenceCell, referenceB, true, audit);
      const queryDetached = queryCell.carrier === "axial" ? detachedCredit(queryCell, queryB) : queryFull.gradient;
      const referenceDetached = referenceCell.carrier === "axial" ? detachedCredit(referenceCell, referenceB) : referenceFull.gradient;
      return {
        value,
        gradient: creditMode === "detached" ? queryDetached.map((v, j) => v + referenceDetached[j]) : queryFull.gradient.map((v, j) => v + referenceFull.gradient[j]),
        fullGradient: queryFull.gradient.map((v, j) => v + referenceFull.gradient[j]),
        detachedGradient: queryDetached.map((v, j) => v + referenceDetached[j]),
        queryDescriptor,
        referenceDescriptor,
        queryB,
        referenceB,
        queryLambda: queryFull.errorVoltage,
        referenceLambda: referenceFull.errorVoltage,
        queryCredit: queryFull.gradient,
        referenceCredit: referenceFull.gradient,
        queryDetachedCredit: queryDetached,
        referenceDetachedCredit: referenceDetached,
      };
    };
    const marginFromCells = (queryCell, referenceACell, referenceBCell, creditMode, audit = null) => {
      const scoreA = scoreFromCells(queryCell, referenceACell, creditMode, audit);
      const scoreB = scoreFromCells(queryCell, referenceBCell, creditMode, audit);
      const makeRoute = (b, lambda, fullLocalCredit, detachedLocalCredit, cellKey) => ({
        cellKey,
        b,
        lambda,
        fullLocalCredit,
        detachedLocalCredit,
        localCredit: creditMode === "detached" ? detachedLocalCredit : fullLocalCredit,
      });
      const routes = {
        query: makeRoute(
          vectorSub(scoreA.queryB, scoreB.queryB),
          vectorSub(scoreA.queryLambda, scoreB.queryLambda),
          vectorSub(scoreA.queryCredit, scoreB.queryCredit),
          vectorSub(scoreA.queryDetachedCredit, scoreB.queryDetachedCredit),
          "query",
        ),
        referenceA: makeRoute(scoreA.referenceB, scoreA.referenceLambda, scoreA.referenceCredit, scoreA.referenceDetachedCredit, "referenceA"),
        referenceB: makeRoute(vectorScale(scoreB.referenceB, -1), vectorScale(scoreB.referenceLambda, -1), vectorScale(scoreB.referenceCredit, -1), vectorScale(scoreB.referenceDetachedCredit, -1), "referenceB"),
      };
      const fullGradient = vectorAdd(vectorAdd(routes.query.fullLocalCredit, routes.referenceA.fullLocalCredit), routes.referenceB.fullLocalCredit);
      const detachedGradient = vectorAdd(vectorAdd(routes.query.detachedLocalCredit, routes.referenceA.detachedLocalCredit), routes.referenceB.detachedLocalCredit);
      const gradient = creditMode === "detached" ? detachedGradient : fullGradient;
      return {
        value: scoreA.value - scoreB.value,
        gradient,
        fullGradient,
        detachedGradient,
        routes,
        rawCells: {query: queryCell, referenceA: referenceACell, referenceB: referenceBCell},
        cells: {query: describeCell(queryCell), referenceA: describeCell(referenceACell), referenceB: describeCell(referenceBCell)},
      };
    };
    const marginWithCells = (th, queryRecord, candidateA, candidateB, queryStart, referenceAStart, referenceBStart, carrier, creditMode, eps, audit = null, eventTolerance = structureTolerance) => {
      const queryCell = event(th, queryRecord.x, carrier, queryStart, eventTolerance, eps, audit);
      const referenceACell = event(th, candidateA.x, carrier, referenceAStart, eventTolerance, eps, audit);
      const referenceBCell = event(th, candidateB.x, carrier, referenceBStart, eventTolerance, eps, audit);
      return marginFromCells(queryCell, referenceACell, referenceBCell, creditMode, audit);
    };
    const transitionPrediction = (record, queryRecord) => {
      assert.equal(record.action, queryRecord.action, "candidate and query actions must be explicit and equal");
      const increment = vectorSub(record.y, record.x);
      return {candidateId: record.id, action: record.action, increment, mu: vectorAdd(queryRecord.x, increment)};
    };
    const relativeEvidence = (pending, y) => {
      const w = vectorSub(y, pending.queryRecord.x);
      const residualA = vectorSub(w, pending.incrementA), residualB = vectorSub(w, pending.incrementB);
      const squaredResidualA = dot(residualA, residualA), squaredResidualB = dot(residualB, residualB);
      const lambda = 0.5 * (squaredResidualB - squaredResidualA);
      const pBefore = sigmoid(pending.margin), pPlus = sigmoid(pending.margin + lambda), delta = pPlus - pBefore;
      const logEllA = -0.5 * (w.length * Math.log(2 * Math.PI) + squaredResidualA);
      const logEllB = -0.5 * (w.length * Math.log(2 * Math.PI) + squaredResidualB);
      const logLikelihood = logSumExp([Math.log(pBefore) + logEllA, Math.log1p(-pBefore) + logEllB]);
      return {y, w, residualA, residualB, squaredResidualA, squaredResidualB, logEllA, logEllB, ellA: Math.exp(logEllA), ellB: Math.exp(logEllB), lambda, pBefore, pPlus, delta, logLikelihood, loss: -logLikelihood};
    };
    const evidenceLossAtMargin = (marginValue, evidence) => -logSumExp([
      Math.log(sigmoid(marginValue)) + evidence.logEllA,
      Math.log1p(-sigmoid(marginValue)) + evidence.logEllB,
    ]);
    const buildPending = (spec, th, version, carrier = "axial", creditMode = "full", eventTolerance = structureTolerance) => {
      const audit = freshAudit();
      const predictionA = transitionPrediction(spec.candidateA, spec.queryRecord), predictionB = transitionPrediction(spec.candidateB, spec.queryRecord);
      const computed = marginWithCells(th, spec.queryRecord, spec.candidateA, spec.candidateB, spec.incoming, cold, cold, carrier, creditMode, structureEpsilon, audit, eventTolerance);
      const cacheScalarCountSubset = numericScalarCount({muA: predictionA.mu, muB: predictionB.mu, margin: computed.value, marginGradient: computed.gradient, marginGradientRoutes: computed.routes, cells: computed.cells});
      const cacheFields = {
        eventId: spec.eventId,
        thetaVersion: version,
        candidateIds: [spec.candidateA.id, spec.candidateB.id],
        candidateRecords: {A: spec.candidateA, B: spec.candidateB},
        queryRecord: {...spec.queryRecord, arrival: null},
        incoming: spec.incoming,
        historyId: spec.historyId,
        muA: predictionA.mu,
        muB: predictionB.mu,
        incrementA: predictionA.increment,
        incrementB: predictionB.increment,
        margin: computed.value,
        marginGradient: computed.gradient,
        fullMarginGradient: computed.fullGradient,
        detachedMarginGradient: computed.detachedGradient,
        marginGradientRoutes: computed.routes,
        cells: computed.cells,
        rawCells: computed.rawCells,
        cacheAudit: audit,
        carrier,
        creditMode,
        eventTolerance,
      };
      const cacheScalarCountTotal = numericScalarCount(cacheFields);
      const pending = {
        ...cacheFields,
        cacheScalarCount: cacheScalarCountTotal,
        cacheScalarCountSubset,
        cacheScalarCountDefinition: "cacheScalarCount counts all numeric leaves in cacheFields, including candidateRecords, queryRecord, incoming, rawCells, cells, route b/lambda/credits and audit counters; cacheScalarCountSubset is the earlier core-math subset.",
        consumed: false,
      };
      Object.defineProperty(pending, "_rawCells", {value: computed.rawCells, enumerable: false});
      return pending;
    };
    const consumePending = (pending, arrival, th, version, eta) => {
      if (pending.consumed) throw new Error("pending event already consumed");
      if (pending.thetaVersion !== version) throw new Error("pending theta version mismatch");
      const evidence = relativeEvidence(pending, arrival.y);
      const deltaTheta = pending.marginGradient.map(value => eta * evidence.delta * value);
      const thetaAfter = vectorAdd(th, deltaTheta);
      pending.consumed = true;
      pending.arrivalRecord = arrival;
      pending.evidence = evidence;
      pending.deltaTheta = deltaTheta;
      pending.thetaAfter = thetaAfter;
      pending.thetaVersionAfter = version + 1;
      return {thetaAfter, evidence, deltaTheta};
    };
    const consumeDiagnostic = (pending, arrival, th, version, radius) => {
      if (pending.consumed) throw new Error("pending event already consumed");
      if (pending.thetaVersion !== version) throw new Error("pending theta version mismatch");
      const evidence = relativeEvidence(pending, arrival.y);
      const gradientNorm = Math.sqrt(dot(pending.marginGradient, pending.marginGradient));
      const direction = gradientNorm === 0 ? pending.marginGradient.map(() => 0) : pending.marginGradient.map(value => Math.sign(evidence.delta) * value / gradientNorm);
      const deltaTheta = direction.map(value => radius * value), thetaAfter = vectorAdd(th, deltaTheta);
      pending.consumed = true;
      pending.arrivalRecord = arrival;
      pending.evidence = evidence;
      pending.deltaTheta = deltaTheta;
      pending.thetaAfter = thetaAfter;
      pending.thetaVersionAfter = version + 1;
      pending.diagnostic = {radius, gradientNorm, direction};
      return {thetaAfter, evidence, deltaTheta};
    };
    const consumeFrozen = (pending, arrival, th, version) => {
      if (pending.consumed) throw new Error("pending event already consumed");
      if (pending.thetaVersion !== version) throw new Error("pending theta version mismatch");
      const evidence = relativeEvidence(pending, arrival.y);
      const deltaTheta = pending.marginGradient.map(() => 0);
      pending.consumed = true;
      pending.arrivalRecord = arrival;
      pending.evidence = evidence;
      pending.deltaTheta = deltaTheta;
      pending.thetaAfter = [...th];
      pending.thetaVersionAfter = version + 1;
      pending.diagnostic = {radius: 0, gradientNorm: Math.sqrt(dot(pending.marginGradient, pending.marginGradient)), direction: deltaTheta, writeApplied: false, reason: "frozen_no_write_reference"};
      return {thetaAfter: [...th], evidence, deltaTheta};
    };
    const tryError = operation => {
      try { operation(); return {threw: false, error: null}; }
      catch (error) { return {threw: true, error: String(error.message || error)}; }
    };
    const pendingSignature = pending => ({
      thetaVersion: pending.thetaVersion,
      candidateIds: pending.candidateIds,
      queryRecord: pending.queryRecord,
      incoming: pending.incoming,
      historyId: pending.historyId,
      muA: pending.muA,
      muB: pending.muB,
      incrementA: pending.incrementA,
      incrementB: pending.incrementB,
      margin: pending.margin,
      marginGradient: pending.marginGradient,
      marginGradientRoutes: pending.marginGradientRoutes,
      cells: pending.cells,
    });
    const computePendingMargin = (pending, th, carrier = pending.carrier, creditMode = pending.creditMode, eps = structureEpsilon, audit = null, eventTolerance = pending.eventTolerance || structureTolerance) => marginWithCells(
      th,
      pending.queryRecord,
      pending.candidateRecords.A,
      pending.candidateRecords.B,
      pending.incoming,
      cold,
      cold,
      carrier,
      creditMode,
      eps,
      audit,
      eventTolerance,
    );
    const finiteDifferenceLoss = (pending, evidence, step, audit) => pending.marginGradient.map((_, j) => {
      const plus = [...theta0], minus = [...theta0];
      plus[j] += step; minus[j] -= step;
      return (evidenceLossAtMargin(computePendingMargin(pending, plus, "axial", "full", structureEpsilon, audit).value, evidence)
        - evidenceLossAtMargin(computePendingMargin(pending, minus, "axial", "full", structureEpsilon, audit).value, evidence)) / (2 * step);
    });
    const cosineValue = (queryCell, referenceCell) => {
      const x = read(queryCell.state), y = read(referenceCell.state), nx = Math.sqrt(dot(x, x)), ny = Math.sqrt(dot(y, y));
      assert.ok(nx > 0 && ny > 0, "undefined normalized descriptor in route finite difference");
      return dot(x, y) / (nx * ny);
    };
    const routeValueAtTheta = (pending, th, route, audit) => {
      const base = pending._rawCells;
      if (route === "query") {
        const queryCell = event(th, pending.queryRecord.x, "axial", pending.incoming, structureTolerance, structureEpsilon, audit);
        return cosineValue(queryCell, base.referenceA) - cosineValue(queryCell, base.referenceB);
      }
      if (route === "referenceA") {
        const referenceA = event(th, pending.candidateRecords.A.x, "axial", cold, structureTolerance, structureEpsilon, audit);
        return cosineValue(base.query, referenceA) - cosineValue(base.query, base.referenceB);
      }
      const referenceB = event(th, pending.candidateRecords.B.x, "axial", cold, structureTolerance, structureEpsilon, audit);
      return cosineValue(base.query, base.referenceA) - cosineValue(base.query, referenceB);
    };
    const finiteDifferenceRoute = (pending, route, step, audit) => pending.marginGradientRoutes[route].localCredit.map((_, j) => {
      const plus = [...theta0], minus = [...theta0]; plus[j] += step; minus[j] -= step;
      return (routeValueAtTheta(pending, plus, route, audit) - routeValueAtTheta(pending, minus, route, audit)) / (2 * step);
    });
    const s5Response = (cell, rhs) => {
      const chi = cell.diagonal.slice(0, 4).map((value, j) => value - cell.edge[j]);
      const chiSoma = cell.diagonal[4] - cell.edge.reduce((sum, value) => sum + value, 0);
      const denominator = chiSoma + cell.edge.reduce((sum, value, j) => sum + value * chi[j] / cell.diagonal[j], 0);
      const soma = (rhs[4] + cell.edge.reduce((sum, value, j) => sum + value * rhs[j] / cell.diagonal[j], 0)) / denominator;
      return {chi: [...chi, chiSoma], A: cell.diagonal.slice(0, 4), D: denominator, lambda: [...cell.edge.map((value, j) => (rhs[j] + value * soma) / cell.diagonal[j]), soma]};
    };
    const buildHistory = (historyId, prefix, carrier, eventTolerance = structureTolerance) => {
      let state = [...cold], steps = [];
      for (const inputValue of prefix) {
        const audit = freshAudit(), cellValue = event(theta0, inputValue, carrier, state, eventTolerance, structureEpsilon, audit);
        steps.push({input: inputValue, cell: describeCell(cellValue), audit});
        state = cellValue.state;
      }
      return {historyId, prefix, carrier, epsilon: structureEpsilon, eventTolerance, thetaVersion: 0, incoming: state, steps};
    };
    const probeSnapshot = (th, carrier, creditMode, historiesForCarrier, audit = null, eventTolerance = structureTolerance) => {
      const states = {};
      for (const [historyId, history] of Object.entries(historiesForCarrier)) {
        const rows = queries.map(queryValue => {
          const localAudit = audit || freshAudit();
          const computed = marginWithCells(th, {id: `probe_${historyId}`, x: queryValue, action: "forward", arrival: null}, candidateA, candidateB, history.incoming, cold, cold, carrier, creditMode, structureEpsilon, localAudit, eventTolerance);
          return {query: queryValue, margin: computed.value, gradient: computed.gradient};
        });
        const lengths = rows.map(row => Math.sqrt(dot(row.gradient, row.gradient)));
        const normalized = rows.map((row, i) => row.gradient.map(value => value / lengths[i]));
        states[historyId] = {historyId, rows, gradientGram: mm(normalized, tr(normalized))};
      }
      return {states, eventTolerance};
    };
    const routeCellGap = (a, b) => Math.max(gap(a.state, b.state), gap(a.drop, b.drop), gap(a.edge, b.edge), gap(a.diagonal, b.diagonal));
    const candidateA = {id: "A", x: [0.3, 0.7], y: [0.4, 0.7], action: "forward", arrival: 1, incoming: [...cold]};
    const candidateB = {id: "B", x: [0.7, -0.2], y: [0.6, -0.2], action: "forward", arrival: 2, incoming: [...cold]};
    const historyPrefixes = {H0: [[1, 0], [0, 1]], H1: [[0, 1], [1, 0]]};
    const axialHistories = {H0: buildHistory("H0", historyPrefixes.H0, "axial"), H1: buildHistory("H1", historyPrefixes.H1, "axial")};
    const queryRecord = {id: "query_current", x: [0.8, 0.2], action: "forward", arrival: null};
    const arrivalTemplates = [
      {eventId: "E0", queryRecord, y: [0.9, 0.23], arrival: 3, historyId: "H0"},
      {eventId: "E1", queryRecord, y: [0.7, 0.23], arrival: 4, historyId: "H1"},
      {eventId: "E0_prime", queryRecord, y: [0.9, 0.23], arrival: 5, historyId: "H0"},
      {eventId: "E1_prime", queryRecord, y: [0.7, 0.23], arrival: 6, historyId: "H1"},
    ];
    const materializeSpec = (template, historiesForCarrier, candidateBOverride = candidateB) => ({
      eventId: template.eventId,
      queryRecord: {...queryRecord},
      y: template.y,
      arrival: template.arrival,
      historyId: template.historyId,
      incoming: historiesForCarrier[template.historyId].incoming,
      candidateA,
      candidateB: candidateBOverride,
    });
    const axialSpecs = arrivalTemplates.map(template => materializeSpec(template, axialHistories));
    const blockOneStart = validationChecks.length;
    const pendingE0 = buildPending(axialSpecs[0], theta0, 0, "axial", "full");
    const pendingE1 = buildPending(axialSpecs[1], theta0, 0, "axial", "full");
    const evidenceE0 = relativeEvidence(pendingE0, axialSpecs[0].y), evidenceE1 = relativeEvidence(pendingE1, axialSpecs[1].y);
    const expectedMuA = vectorAdd(axialSpecs[0].queryRecord.x, vectorSub(candidateA.y, candidateA.x));
    const expectedMuB = vectorAdd(axialSpecs[0].queryRecord.x, vectorSub(candidateB.y, candidateB.x));
    structureCheck("S1 transition prediction uses completed-record increments", vectorEqual(pendingE0.muA, expectedMuA) && vectorEqual(pendingE0.muB, expectedMuB), "correctness", {muA: pendingE0.muA, muB: pendingE0.muB});
    structureCheck("S2 E0 relative evidence matches the fixed density identity", Math.abs(evidenceE0.lambda - 0.02) <= 1e-10, "correctness", {lambda: evidenceE0.lambda, expected: 0.02});
    structureCheck("S2 E1 relative evidence reverses on the fixed successor", Math.abs(evidenceE1.lambda + 0.02) <= 1e-10, "correctness", {lambda: evidenceE1.lambda, expected: -0.02});
    structureCheck("S3 evidence changes only after the pending successor arrives", pendingE0.consumed === false && pendingE0.queryRecord.arrival === null && Number.isFinite(pendingE0.margin), "correctness", {consumed: pendingE0.consumed, queryArrival: pendingE0.queryRecord.arrival});
    const changedSuccessorPending = buildPending({...axialSpecs[0], y: [0.9, 0.30]}, theta0, 0, "axial", "full");
    structureCheck("S3 changing unrevealed y leaves frozen pending cache unchanged", JSON.stringify(pendingSignature(pendingE0)) === JSON.stringify(pendingSignature(changedSuccessorPending)), "correctness", {originalEventId: pendingE0.eventId, changedSuccessorIgnoredUntilArrival: true});
    const orthogonalEvidence = relativeEvidence(pendingE0, [0.9, 0.30]);
    structureCheck("S2 orthogonal successor perturbation leaves Lambda unchanged", Math.abs(orthogonalEvidence.lambda - evidenceE0.lambda) <= 1e-10, "correctness", {baseLambda: evidenceE0.lambda, perturbedLambda: orthogonalEvidence.lambda, perturbation: [0, 0.07]});
    const equalIncrement = vectorSub(candidateA.y, candidateA.x);
    const equalCandidateB = {...candidateB, id: "B_equal_increment", y: vectorAdd(candidateB.x, equalIncrement)}, equalSpec = materializeSpec(axialSpecs[0], axialHistories, equalCandidateB);
    const equalPending = buildPending(equalSpec, theta0, 0, "axial", "full"), equalEvidence = relativeEvidence(equalPending, equalSpec.y);
    const equalConsumed = consumePending(equalPending, equalSpec, theta0, 0, mainEta);
    const equalQualificationNorm = Math.sqrt(dot(equalPending.marginGradient, equalPending.marginGradient));
    structureCheck("S2 equal candidate increments preserve distinct B input and nonzero association qualification", gap(equalCandidateB.x, candidateA.x) > 0 && gap(equalIncrement, vectorSub(equalCandidateB.y, equalCandidateB.x)) <= 1e-14 && equalQualificationNorm > 1e-12, "correctness", {candidateAX: candidateA.x, candidateBX: equalCandidateB.x, incrementA: equalIncrement, incrementB: vectorSub(equalCandidateB.y, equalCandidateB.x), qualificationMargin: equalPending.margin, qualificationGradient: equalPending.marginGradient, qualificationNorm: equalQualificationNorm});
    structureCheck("S2 equal candidate increments produce zero relative evidence", Math.abs(equalEvidence.lambda) <= 1e-14 && Math.abs(equalEvidence.delta) <= 1e-14, "correctness", {lambda: equalEvidence.lambda, delta: equalEvidence.delta, w: equalEvidence.w, residualA: equalEvidence.residualA, residualB: equalEvidence.residualB});
    structureCheck("S4 zero evidence leaves every parameter unchanged", vectorEqual(equalConsumed.thetaAfter, theta0), "correctness", {thetaBefore: theta0, thetaAfter: equalConsumed.thetaAfter});
    const coarseLossAudit = freshAudit(), fineLossAudit = freshAudit();
    const analyticLossGradient = vectorScale(pendingE0.marginGradient, -evidenceE0.delta);
    const coarseLossGradient = finiteDifferenceLoss(pendingE0, evidenceE0, finiteDifferenceStep, coarseLossAudit);
    const fineLossGradient = finiteDifferenceLoss(pendingE0, evidenceE0, finiteDifferenceStepFine, fineLossAudit);
    const coarseLossGap = gap(analyticLossGradient, coarseLossGradient), fineLossGap = gap(analyticLossGradient, fineLossGradient);
    const lossRelativeScale = Math.max(norm(analyticLossGradient), norm(fineLossGradient), 1e-300);
    structureCheck("S4 full edge likelihood gradient matches central finite differences", coarseLossGap <= 1e-9 && fineLossGap <= 1e-9, "correctness", {analyticLossGradient, coarseLossGradient, fineLossGradient, coarseAbsoluteGap: coarseLossGap, fineAbsoluteGap: fineLossGap, coarseRelativeGap: coarseLossGap / lossRelativeScale, fineRelativeGap: fineLossGap / lossRelativeScale});
    structureCheck("S4 both finite-difference spacings agree on the gradient", gap(coarseLossGradient, fineLossGradient) <= 1e-9, "correctness", {spacingGap: gap(coarseLossGradient, fineLossGradient), finiteDifferenceStep, finiteDifferenceStepFine});
    const routeFiniteDifference = {}, routeFiniteDifferenceAudit = freshAudit();
    for (const route of ["query", "referenceA", "referenceB"]) {
      const coarse = finiteDifferenceRoute(pendingE0, route, finiteDifferenceStep, routeFiniteDifferenceAudit);
      const fine = finiteDifferenceRoute(pendingE0, route, finiteDifferenceStepFine, routeFiniteDifferenceAudit);
      const analytic = pendingE0.marginGradientRoutes[route].localCredit;
      routeFiniteDifference[route] = {analytic, coarse, fine, coarseAbsoluteGap: gap(analytic, coarse), fineAbsoluteGap: gap(analytic, fine), spacingGap: gap(coarse, fine)};
    }
    structureCheck("S4 query route margin derivative matches finite differences", routeFiniteDifference.query.coarseAbsoluteGap <= 1e-9 && routeFiniteDifference.query.fineAbsoluteGap <= 1e-9, "correctness", routeFiniteDifference.query);
    structureCheck("S4 candidate-A route margin derivative matches finite differences", routeFiniteDifference.referenceA.coarseAbsoluteGap <= 1e-9 && routeFiniteDifference.referenceA.fineAbsoluteGap <= 1e-9, "correctness", routeFiniteDifference.referenceA);
    structureCheck("S4 candidate-B route margin derivative matches finite differences", routeFiniteDifference.referenceB.coarseAbsoluteGap <= 1e-9 && routeFiniteDifference.referenceB.fineAbsoluteGap <= 1e-9, "correctness", routeFiniteDifference.referenceB);
    const routeSum = vectorAdd(vectorAdd(pendingE0.marginGradientRoutes.query.localCredit, pendingE0.marginGradientRoutes.referenceA.localCredit), pendingE0.marginGradientRoutes.referenceB.localCredit);
    structureCheck("S4 three route credits sum to the cached margin gradient", gap(routeSum, pendingE0.marginGradient) <= 1e-12, "correctness", {routeSum, marginGradient: pendingE0.marginGradient, gap: gap(routeSum, pendingE0.marginGradient)});
    const denseReference = {}, s5Rows = {};
    for (const route of ["query", "referenceA", "referenceB"]) {
      const routeData = pendingE0.marginGradientRoutes[route], cellValue = pendingE0._rawCells[routeData.cellKey];
      const arrowhead = response(cellValue.diagonal, cellValue.edge, routeData.b), dense = denseResponse(cellValue, routeData.b), schur = s5Response(cellValue, routeData.b);
      denseReference[route] = {arrowhead, dense, gap: gap(arrowhead, dense)};
      s5Rows[route] = {lambda: routeData.lambda, lambdaS5: schur.lambda, gap: gap(routeData.lambda, schur.lambda), chi: schur.chi, A: schur.A, D: schur.D};
    }
    structureCheck("S5 arrowhead response agrees with dense 5x5 reference", Object.values(denseReference).every(row => row.gap <= 1e-12), "correctness", denseReference);
    structureCheck("S5 closed Schur response agrees with the local response on all routes", Object.values(s5Rows).every(row => row.gap <= 1e-10), "correctness", s5Rows);
    const s6Route = "query", s6Cell = pendingE0._rawCells.query, s6B = pendingE0.marginGradientRoutes[s6Route].b;
    const s6Chi = s6Cell.diagonal.slice(0, 4).map((value, j) => value - s6Cell.edge[j]), s6ChiSoma = s6Cell.diagonal[4] - s6Cell.edge.reduce((sum, value) => sum + value, 0);
    const s6D = s6ChiSoma + s6Cell.edge.reduce((sum, value, j) => sum + value * s6Chi[j] / s6Cell.diagonal[j], 0);
    const s6Analytic = Array.from({length: 4}, () => Array(4).fill(0));
    for (let j = 0; j < 4; j++) for (let k = 0; k < 4; k++) if (j !== k) s6Analytic[j][k] = h * s6Cell.coefficient[j] * s6Cell.drop[j] ** 3 * s6Chi[j] * s6Cell.edge[k] / (s6Cell.diagonal[j] * s6Cell.diagonal[k] * s6D);
    const s6FiniteDifference = Array.from({length: 4}, () => Array(4).fill(0));
    for (let k = 0; k < 4; k++) {
      const plus = [...s6B], minus = [...s6B]; plus[k] += finiteDifferenceStep; minus[k] -= finiteDifferenceStep;
      const gPlus = localCredit(s6Cell, plus), gMinus = localCredit(s6Cell, minus);
      for (let j = 0; j < 4; j++) s6FiniteDifference[j][k] = (gPlus[j] - gMinus[j]) / (2 * finiteDifferenceStep);
    }
    const s6CrossFiniteDifference = s6FiniteDifference.map((row, j) => row.map((value, k) => j === k ? 0 : value));
    const s6CrossGap = mgap(s6Analytic, s6CrossFiniteDifference);
    structureCheck("S6 single-coordinate cross-branch credit transport matches finite differences", s6CrossGap <= 1e-9, "correctness", {route: s6Route, analyticCrossBranch: s6Analytic, finiteDifferenceFull: s6FiniteDifference, finiteDifferenceCrossBranch: s6CrossFiniteDifference, diagonalFiniteDifference: s6FiniteDifference.map((row, j) => row[j]), absoluteCrossBranchGap: s6CrossGap, step: finiteDifferenceStep});
    const commonScale = 0.7, commonCell = pendingE0._rawCells.query, commonChi = [...s5Rows.query.chi], commonB = commonChi.map(value => commonScale * value), commonCredit = localCredit(commonCell, commonB);
    structureCheck("S5 common weighted error has zero axial local response", norm(commonCredit) <= 1e-12, "correctness", {commonScale, commonB, localCredit: commonCredit});
    const fullDetachedNumerator = norm(vectorSub(pendingE0.fullMarginGradient, pendingE0.detachedMarginGradient)), fullGradientNorm = norm(pendingE0.fullMarginGradient);
    const fullDetachedFraction = fullGradientNorm === 0 ? null : fullDetachedNumerator / fullGradientNorm;
    structureCheck("S4 full and detached gradients are separately recorded on the actual margin", fullDetachedFraction !== null && Number.isFinite(fullDetachedFraction), "correctness", {fullGradient: pendingE0.fullMarginGradient, detachedGradient: pendingE0.detachedMarginGradient, numerator: fullDetachedNumerator, denominator: fullGradientNorm, fraction: fullDetachedFraction});
    const paperA = pendingE0.marginGradient, paperB = pendingE1.marginGradient;
    const paperANorm = Math.sqrt(dot(paperA, paperA)), paperBNorm = Math.sqrt(dot(paperB, paperB));
    const paperD0 = evidenceE0.delta, paperD1 = -evidenceE1.delta;
    const paperC = dot(paperA, paperB) / (paperANorm * paperBNorm), paperR = paperD0 * paperANorm / (paperD1 * paperBNorm);
    const paperPrediction = {a: paperA, b: paperB, d0: paperD0, d1: paperD1, c: paperC, r: paperR, lowerMargin: paperR - paperC, upperMargin: 1 / paperC - paperR, condition: paperC < paperR && paperR < 1 / paperC, source: "theta0 H0/H1 margin gradients and arrived E0/E1 relative evidence; this is the registered first-order prediction, not a visual result."};
    const strictE0Audit = freshAudit();
    const strictE0 = marginWithCells(theta0, axialSpecs[0].queryRecord, candidateA, candidateB, axialSpecs[0].incoming, cold, cold, "axial", "full", structureEpsilon, strictE0Audit, strictTolerance);
    const strictE0Gap = Math.max(Math.abs(strictE0.value - pendingE0.margin), gap(strictE0.gradient, pendingE0.marginGradient));
    structureCheck("S4 strict 1e-14 event solve agrees with the 1e-12 pending solve", strictE0Gap <= 1e-12, "correctness", {strictTolerance, structureTolerance, valueGap: Math.abs(strictE0.value - pendingE0.margin), gradientGap: gap(strictE0.gradient, pendingE0.marginGradient), audit: strictE0Audit});
    const stalePending = buildPending(axialSpecs[0], theta0, 0, "axial", "full");
    const staleVersion = tryError(() => consumePending(stalePending, {eventId: "E0", y: axialSpecs[0].y, arrival: axialSpecs[0].arrival}, theta0, 1, mainEta));
    structureCheck("S3 pending consumption rejects a stale theta version", staleVersion.threw, "correctness", staleVersion);
    const duplicatePending = buildPending(axialSpecs[0], theta0, 0, "axial", "full");
    consumePending(duplicatePending, {eventId: "E0", y: axialSpecs[0].y, arrival: axialSpecs[0].arrival}, theta0, 0, mainEta);
    const duplicateConsumption = tryError(() => consumePending(duplicatePending, {eventId: "E0", y: axialSpecs[0].y, arrival: axialSpecs[0].arrival}, vectorAdd(theta0, duplicatePending.deltaTheta), 1, mainEta));
    structureCheck("S3 pending consumption rejects a duplicate event", duplicateConsumption.threw, "correctness", duplicateConsumption);
    const pendingRecord = pending => ({...pending, rawCells: pending._rawCells});
    const blockOneChecks = validationChecks.slice(blockOneStart);
    const blockOne = {
      id: "validation_block_one_source_evidence_and_exact_credit",
      tolerance: structureTolerance,
      strictTolerance,
      finiteDifferenceStep,
      finiteDifferenceStepFine,
      inputs: {candidateA, candidateB, queryRecord, histories: historyPrefixes, epsilon: structureEpsilon, theta0, eta: mainEta},
      historyLibrary: axialHistories,
      events: {
        E0: {spec: axialSpecs[0], pending: pendingRecord(pendingE0), evidence: evidenceE0},
        E1: {spec: axialSpecs[1], pending: pendingRecord(pendingE1), evidence: evidenceE1},
      },
      boundaryCases: {
        changedSuccessorPending: pendingRecord(changedSuccessorPending),
        orthogonalEvidence,
        equalCandidate: {candidateB: equalCandidateB, pending: pendingRecord(equalPending), evidence: equalEvidence, consumption: equalConsumed},
        staleVersion,
        duplicateConsumption,
      },
      finiteDifference: {
        analyticLossGradient,
        coarseLossGradient,
        fineLossGradient,
        coarseAbsoluteGap: coarseLossGap,
        fineAbsoluteGap: fineLossGap,
        coarseRelativeGap: coarseLossGap / lossRelativeScale,
        fineRelativeGap: fineLossGap / lossRelativeScale,
        spacingGap: gap(coarseLossGradient, fineLossGradient),
        coarseAudit: coarseLossAudit,
        fineAudit: fineLossAudit,
      },
      routeFiniteDifference,
      denseReference,
      s5Rows,
      s6: {route: s6Route, analyticCrossBranch: s6Analytic, finiteDifferenceFull: s6FiniteDifference, finiteDifferenceCrossBranch: s6CrossFiniteDifference, absoluteCrossBranchGap: s6CrossGap, step: finiteDifferenceStep},
      commonWeightedError: {commonScale, commonB, localCredit: commonCredit, reason: "b is proportional to the full Schur chi vector; this is an algebra audit, not a teacher."},
      fullDetached: {fullGradient: pendingE0.fullMarginGradient, detachedGradient: pendingE0.detachedMarginGradient, numerator: fullDetachedNumerator, denominator: fullGradientNorm, fraction: fullDetachedFraction},
      paperPrediction,
      strictAudit: {E0: {margin: strictE0.value, gradient: strictE0.gradient, audit: strictE0Audit, gap: strictE0Gap}},
      checks: blockOneChecks,
      correctnessPassed: blockOneChecks.every(item => item.passed),
    };
    const associationSnapshot = (th, carrier, creditMode, historiesForCarrier, audit = null, eventTolerance = structureTolerance) => {
      const localAudit = audit || freshAudit(), rows = {};
      for (const [historyId, history] of Object.entries(historiesForCarrier)) {
        const computed = marginWithCells(th, queryRecord, candidateA, candidateB, history.incoming, cold, cold, carrier, creditMode, structureEpsilon, localAudit, eventTolerance);
        rows[historyId] = {
          historyId,
          query: queryRecord.x,
          margin: computed.value,
          gradient: computed.gradient,
          fullGradient: computed.fullGradient,
          detachedGradient: computed.detachedGradient,
          routes: computed.routes,
          cells: computed.cells,
          rawCells: computed.rawCells,
        };
      }
      return {rows, audit: localAudit, eventTolerance};
    };
    const forwardEquivalence = (historiesForCarrier, eventTolerance = structureTolerance) => {
      const rows = {};
      for (const [historyId, history] of Object.entries(historiesForCarrier)) {
        const fullAudit = freshAudit(), detachedAudit = freshAudit();
        const full = marginWithCells(theta0, queryRecord, candidateA, candidateB, history.incoming, cold, cold, "axial", "full", structureEpsilon, fullAudit, eventTolerance);
        const detached = marginWithCells(theta0, queryRecord, candidateA, candidateB, history.incoming, cold, cold, "axial", "detached", structureEpsilon, detachedAudit, eventTolerance);
        const stateGap = Math.max(...["query", "referenceA", "referenceB"].map(route => routeCellGap(full.rawCells[route], detached.rawCells[route])));
        const valueGap = Math.abs(full.value - detached.value);
        rows[historyId] = {historyId, valueGap, stateGap, full: {value: full.value, cells: full.cells, audit: fullAudit}, detached: {value: detached.value, cells: detached.cells, audit: detachedAudit}};
      }
      return {eventTolerance, rows};
    };
    const euclideanNorm = values => Math.sqrt(dot(values, values));
    const fTable = snapshot => Object.fromEntries(Object.entries(snapshot.rows).map(([historyId, row]) => [historyId, {margin: row.margin, gradient: row.gradient, fullGradient: row.fullGradient, detachedGradient: row.detachedGradient}]));
    const runSequence = ({label, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates, mode, radius, eventTolerance = structureTolerance}) => {
      let currentTheta = [...theta0], currentVersion = 0;
      const initialAssociationAudit = freshAudit(), initialAssociation = associationSnapshot(currentTheta, carrier, creditMode, historiesForCarrier, initialAssociationAudit, eventTolerance);
      const initialProbeAudit = freshAudit(), initialProbe = probeSnapshot(currentTheta, carrier, creditMode, probeHistoriesForCarrier, initialProbeAudit, eventTolerance);
      const eventRows = [];
      for (const template of templates) {
        const spec = materializeSpec(template, historiesForCarrier);
        const pending = buildPending(spec, currentTheta, currentVersion, carrier, creditMode, eventTolerance);
        const beforeAssociationAudit = freshAudit(), beforeAssociation = associationSnapshot(currentTheta, carrier, creditMode, historiesForCarrier, beforeAssociationAudit, eventTolerance);
        const beforeProbeAudit = freshAudit(), beforeProbe = probeSnapshot(currentTheta, carrier, creditMode, probeHistoriesForCarrier, beforeProbeAudit, eventTolerance);
        const arrivalRecord = {eventId: spec.eventId, y: spec.y, arrival: spec.arrival, queryRecord: spec.queryRecord, historyId: spec.historyId};
        const update = mode === "main"
          ? consumePending(pending, arrivalRecord, currentTheta, currentVersion, mainEta)
          : mode === "diagnostic"
            ? consumeDiagnostic(pending, arrivalRecord, currentTheta, currentVersion, radius)
            : consumeFrozen(pending, arrivalRecord, currentTheta, currentVersion);
        const afterAssociationAudit = freshAudit(), afterAssociation = associationSnapshot(update.thetaAfter, carrier, creditMode, historiesForCarrier, afterAssociationAudit, eventTolerance);
        const afterProbeAudit = freshAudit(), afterProbe = probeSnapshot(update.thetaAfter, carrier, creditMode, probeHistoriesForCarrier, afterProbeAudit, eventTolerance);
        const targetEffects = {};
        for (const historyId of Object.keys(historiesForCarrier)) {
          const beforeRow = beforeAssociation.rows[historyId], afterRow = afterAssociation.rows[historyId];
          const actualChange = afterRow.margin - beforeRow.margin;
          const readoutGradient = beforeRow.fullGradient;
          const updateGradient = pending.marginGradient;
          const gradientDot = dot(readoutGradient, update.deltaTheta);
          const kernel = dot(readoutGradient, updateGradient);
          targetEffects[historyId] = {
            historyId,
            fBefore: beforeRow.margin,
            fAfter: afterRow.margin,
            actualChange,
            firstOrderFromDeltaTheta: gradientDot,
            firstOrderEtaDeltaKernel: mode === "main" ? mainEta * update.evidence.delta * kernel : null,
            firstOrderNormalizedKernel: mode === "diagnostic" ? radius * Math.sign(update.evidence.delta) * kernel / Math.max(euclideanNorm(updateGradient), 1e-300) : null,
            readoutGradientFull: readoutGradient,
            updateGradientSelected: updateGradient,
            updateGradientFull: pending.fullMarginGradient,
            updateGradientDetached: pending.detachedMarginGradient,
          };
        }
        const row = {
          eventId: spec.eventId,
          arrival: {arrival: spec.arrival, y: spec.y, action: spec.queryRecord.action, historyId: spec.historyId},
          carrier,
          creditMode,
          mode,
          eventTolerance,
          thetaVersionBefore: currentVersion,
          thetaBefore: currentTheta,
          fBefore: fTable(beforeAssociation),
          gBefore: Object.fromEntries(Object.entries(beforeAssociation.rows).map(([historyId, value]) => [historyId, value.gradient])),
          pending: pendingRecord(pending),
          evidence: update.evidence,
          deltaTheta: update.deltaTheta,
          parameterWriteNorm: euclideanNorm(update.deltaTheta),
          parameterWriteMaxAbs: norm(update.deltaTheta),
          thetaAfter: update.thetaAfter,
          thetaVersionAfter: currentVersion + 1,
          fAfter: fTable(afterAssociation),
          gAfter: Object.fromEntries(Object.entries(afterAssociation.rows).map(([historyId, value]) => [historyId, value.gradient])),
          targetEffects,
          diagnostics: mode === "main" ? {eta: mainEta, delta: update.evidence.delta, gradient: pending.marginGradient} : pending.diagnostic,
          audits: {cache: pending.cacheAudit, beforeAssociation: beforeAssociation.audit, beforeProbe: beforeProbe.audit, afterAssociation: afterAssociation.audit, afterProbe: afterProbe.audit},
          beforeAssociation,
          afterAssociation,
          beforeProbe,
          afterProbe,
        };
        eventRows.push(row);
        currentTheta = update.thetaAfter;
        currentVersion += 1;
      }
      const finalAssociationAudit = freshAudit(), finalAssociation = associationSnapshot(currentTheta, carrier, creditMode, historiesForCarrier, finalAssociationAudit, eventTolerance);
      const finalProbeAudit = freshAudit(), finalProbe = probeSnapshot(currentTheta, carrier, creditMode, probeHistoriesForCarrier, finalProbeAudit, eventTolerance);
      const firstLambda = {};
      for (const row of eventRows) if (!(row.arrival.historyId in firstLambda)) firstLambda[row.arrival.historyId] = row.evidence.lambda;
      const retention = {};
      for (const historyId of Object.keys(historiesForCarrier)) {
        const lambda = firstLambda[historyId] ?? null, sign = lambda === null || lambda === 0 ? null : Math.sign(lambda);
        const initialValue = initialAssociation.rows[historyId].margin, finalValue = finalAssociation.rows[historyId].margin;
        const rawFinalChange = finalValue - initialValue, D = sign === null ? null : sign * rawFinalChange;
        const own = eventRows.find(row => row.arrival.historyId === historyId);
        const ownRawChange = own ? own.targetEffects[historyId].actualChange : null, B = sign === null || own === undefined ? null : sign * ownRawChange;
        let keepRatio = null, keepReason = null;
        if (B === null) keepReason = sign === null ? "zero_or_missing_first_lambda" : "missing_first_event";
        else if (B <= 1e-12) keepReason = "B_not_above_1e-12";
        else keepRatio = D / B;
        retention[historyId] = {historyId, firstLambda: lambda, sign, fInitial: initialValue, fFinal: finalValue, rawFinalChange, D, ownEventId: own ? own.eventId : null, ownRawChange, B, keepRatio, keepReason};
      }
      const targetEffects = eventRows.flatMap(row => Object.values(row.targetEffects).map(effect => ({eventId: row.eventId, eventHistoryId: row.arrival.historyId, ...effect, otherEvent: row.arrival.historyId !== effect.historyId})));
      const retractions = targetEffects.filter(effect => effect.otherEvent && effect.firstOrderFromDeltaTheta !== null && ((retention[effect.historyId].sign !== null && retention[effect.historyId].sign * effect.actualChange < -1e-12))).map(effect => ({eventId: effect.eventId, historyId: effect.historyId, directionalChange: retention[effect.historyId].sign * effect.actualChange, actualChange: effect.actualChange, firstOrderFromDeltaTheta: effect.firstOrderFromDeltaTheta}));
      const probeDrift = {};
      for (const historyId of Object.keys(probeHistoriesForCarrier)) {
        const initialRows = initialProbe.states[historyId].rows, finalRows = finalProbe.states[historyId].rows;
        probeDrift[historyId] = initialRows.map((initialRow, index) => ({
          queryIndex: index,
          query: initialRow.query,
          initialMargin: initialRow.margin,
          finalMargin: finalRows[index].margin,
          marginChange: finalRows[index].margin - initialRow.margin,
          initialGradient: initialRow.gradient,
          finalGradient: finalRows[index].gradient,
          gradientChange: vectorSub(finalRows[index].gradient, initialRow.gradient),
          gradientGap: gap(finalRows[index].gradient, initialRow.gradient),
        }));
      }
      return {
        label,
        carrier,
        creditMode,
        mode,
        radius: mode === "main" ? mainEta : radius,
        eventTolerance,
        thetaInitial: theta0,
        thetaFinal: currentTheta,
        parameterPathLength: eventRows.reduce((sum, row) => sum + row.parameterWriteNorm, 0),
        parameterWriteNorms: eventRows.map(row => row.parameterWriteNorm),
        initialAssociation,
        finalAssociation,
        initialProbe,
        finalProbe,
        probeDrift,
        fInitial: fTable(initialAssociation),
        fFinal: fTable(finalAssociation),
        events: eventRows,
        retention,
        retractions,
        firstOrderPredictions: eventRows.map(row => ({eventId: row.eventId, targetEffects: row.targetEffects})),
        finalVersion: currentVersion,
      };
    };
    const makeColdProbeHistory = carrier => ({historyId: "cold", prefix: [], carrier, epsilon: structureEpsilon, eventTolerance: structureTolerance, thetaVersion: 0, incoming: [...cold], steps: []});
    const runWarmGroup = ({label, carrier, creditMode, historiesForCarrier, strictHistoriesForCarrier = null}) => {
      const probeHistoriesForCarrier = {cold: makeColdProbeHistory(carrier), ...historiesForCarrier};
      return {
        label,
        carrier,
        creditMode,
        historyLibrary: historiesForCarrier,
        probeHistoryLibrary: probeHistoriesForCarrier,
        main: runSequence({label: `${label}_main`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates.slice(0, 2), mode: "main", radius: mainEta}),
        rho: runSequence({label: `${label}_rho`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates.slice(0, 2), mode: "diagnostic", radius: rho}),
        rhoHalf: runSequence({label: `${label}_rho_half`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates.slice(0, 2), mode: "diagnostic", radius: rhoHalf}),
        frozen: runSequence({label: `${label}_frozen`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates.slice(0, 2), mode: "frozen", radius: 0}),
        interleavedMain: runSequence({label: `${label}_interleaved_main`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates, mode: "main", radius: mainEta}),
        interleavedRho: runSequence({label: `${label}_interleaved_rho`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates, mode: "diagnostic", radius: rho}),
        interleavedFrozen: runSequence({label: `${label}_interleaved_frozen`, carrier, creditMode, historiesForCarrier, probeHistoriesForCarrier, templates: arrivalTemplates, mode: "frozen", radius: 0}),
        strictMain: strictHistoriesForCarrier === null ? null : runSequence({label: `${label}_strict_main`, carrier, creditMode, historiesForCarrier: strictHistoriesForCarrier, probeHistoriesForCarrier: {cold: makeColdProbeHistory(carrier), ...strictHistoriesForCarrier}, templates: arrivalTemplates.slice(0, 2), mode: "main", radius: mainEta, eventTolerance: strictTolerance}),
        strictRho: strictHistoriesForCarrier === null ? null : runSequence({label: `${label}_strict_rho`, carrier, creditMode, historiesForCarrier: strictHistoriesForCarrier, probeHistoriesForCarrier: {cold: makeColdProbeHistory(carrier), ...strictHistoriesForCarrier}, templates: arrivalTemplates.slice(0, 2), mode: "diagnostic", radius: rho, eventTolerance: strictTolerance}),
      };
    };
    const tonicHistories = {H0: buildHistory("H0", historyPrefixes.H0, "tonic"), H1: buildHistory("H1", historyPrefixes.H1, "tonic")};
    const strictAxialHistories = {H0: buildHistory("H0", historyPrefixes.H0, "axial", strictTolerance), H1: buildHistory("H1", historyPrefixes.H1, "axial", strictTolerance)};
    const blockTwoStart = validationChecks.length;
    const forwardEquivalenceData = forwardEquivalence(axialHistories);
    for (const [historyId, row] of Object.entries(forwardEquivalenceData.rows)) {
      structureCheck(`block2 ${historyId} detached credit keeps the full axial forward cells`, row.valueGap <= 1e-12 && row.stateGap <= 1e-12, "correctness", row);
    }
    const axialFullGroup = runWarmGroup({label: "axial_full", carrier: "axial", creditMode: "full", historiesForCarrier: axialHistories, strictHistoriesForCarrier: strictAxialHistories});
    const axialDetachedGroup = runWarmGroup({label: "axial_detached", carrier: "axial", creditMode: "detached", historiesForCarrier: axialHistories});
    const tonicGroup = runWarmGroup({label: "tonic", carrier: "tonic", creditMode: "full", historiesForCarrier: tonicHistories});
    for (const group of [axialFullGroup, axialDetachedGroup, tonicGroup]) {
      for (const frozenSequence of [group.frozen, group.interleavedFrozen]) {
        structureCheck(`block2 ${frozenSequence.label} is a no-write reference`, frozenSequence.events.every(row => row.deltaTheta.every(value => value === 0)) && vectorEqual(frozenSequence.thetaFinal, theta0) && frozenSequence.parameterPathLength === 0 && frozenSequence.events.every(row => row.pending.diagnostic.writeApplied === false), "correctness", {label: frozenSequence.label, thetaFinal: frozenSequence.thetaFinal, parameterPathLength: frozenSequence.parameterPathLength, eventCount: frozenSequence.events.length});
      }
    }
    const signStatus = value => value === null || Math.abs(value) <= 1e-12 ? "unresolved" : Math.sign(value);
    const compareDSigns = (left, right) => Object.keys(left.retention).every(historyId => signStatus(left.retention[historyId].D) === signStatus(right.retention[historyId].D));
    structureCheck("block2 1e-12 and 1e-14 axial main signs agree", compareDSigns(axialFullGroup.main, axialFullGroup.strictMain), "correctness", {coarse: axialFullGroup.main.retention, strict: axialFullGroup.strictMain.retention});
    structureCheck("block2 1e-12 and 1e-14 axial rho signs agree", compareDSigns(axialFullGroup.rho, axialFullGroup.strictRho), "correctness", {coarse: axialFullGroup.rho.retention, strict: axialFullGroup.strictRho.retention});
    const directionCosine = (left, right, index) => {
      const a = left.events[index].deltaTheta, bValue = right.events[index].deltaTheta;
      return dot(a, bValue) / Math.max(euclideanNorm(a) * euclideanNorm(bValue), 1e-300);
    };
    const rhoHalfDirection = [0, 1].map(index => ({eventId: axialFullGroup.rhoHalf.events[index].eventId, cosineWithRho: directionCosine(axialFullGroup.rhoHalf, axialFullGroup.rho, index), cosineWithMain: directionCosine(axialFullGroup.rhoHalf, axialFullGroup.main, index)}));
    const rhoHalfSignStability = ["H0", "H1"].map(historyId => {
      const rhoHalfD = axialFullGroup.rhoHalf.retention[historyId].D, rhoD = axialFullGroup.rho.retention[historyId].D, mainD = axialFullGroup.main.retention[historyId].D;
      return {historyId, rhoHalfD, rhoD, mainD, rhoHalfStatus: signStatus(rhoHalfD), rhoStatus: signStatus(rhoD), mainStatus: signStatus(mainD), stableWithRho: signStatus(rhoHalfD) === signStatus(rhoD), stableWithMain: signStatus(rhoHalfD) === signStatus(mainD)};
    });
    structureCheck("block2 rho/2 first two direction diagnostics are finite and recorded", rhoHalfDirection.length === 2 && rhoHalfDirection.every(row => Number.isFinite(row.cosineWithRho) && Number.isFinite(row.cosineWithMain)), "correctness", {rhoHalfDirection, rho, rhoHalf});
    structureCheck("block2 rho/2 final D signs are recorded for both evidence directions", rhoHalfSignStability.length === 2 && rhoHalfSignStability.every(row => row.rhoHalfStatus !== undefined), "correctness", {rhoHalfSignStability});
    const blockTwoChecks = validationChecks.slice(blockTwoStart);
    const blockTwo = {
      id: "validation_block_two_warm_opposite_corrections_and_interleaved_retention",
      tolerance: structureTolerance,
      strictTolerance,
      parameters: {eta: mainEta, rho, rhoHalf, epsilon: structureEpsilon, theta0, arrivalOrder: arrivalTemplates.map(item => ({eventId: item.eventId, arrival: item.arrival}))},
      fixedContexts: {axial: axialHistories, tonic: tonicHistories, strictAxial: strictAxialHistories},
      forwardEquivalence: forwardEquivalenceData,
      groups: {axialFull: axialFullGroup, axialDetached: axialDetachedGroup, tonic: tonicGroup},
      strictCrossCheck: {axialMain: axialFullGroup.strictMain, axialRho: axialFullGroup.strictRho, rhoHalfDirection, rhoHalfSignStability},
      checks: blockTwoChecks,
      correctnessPassed: blockTwoChecks.every(item => item.passed),
    };
    const scientificChecks = [];
    const scientificCheck = (name, passed, details = {}) => {
      scientificChecks.push({name, passed: Boolean(passed), kind: "scientific", ...details});
      return Boolean(passed);
    };
    const retentionDPositive = sequence => Object.values(sequence.retention).every(item => item.D !== null && item.D > 1e-12);
    scientificCheck("legal interface: block one correctness layer", blockOne.correctnessPassed, {reason: "reported separately from scientific outcomes"});
    scientificCheck("constructive slice: axial full main D0 and D1 exceed 1e-12", retentionDPositive(axialFullGroup.main), {retention: axialFullGroup.main.retention});
    scientificCheck("constructive slice: axial full equal-rho D0 and D1 exceed 1e-12", retentionDPositive(axialFullGroup.rho), {retention: axialFullGroup.rho.retention});
    scientificCheck("constructive slice: rho/2 final D0/D1 evidence-direction signs are stable", rhoHalfSignStability.every(row => row.rhoHalfStatus !== "unresolved" && row.stableWithRho && row.stableWithMain), {rhoHalfSignStability, cosineDescription: rhoHalfDirection});
    scientificCheck("interleaved retention: axial full main D0 and D1 exceed 1e-12", retentionDPositive(axialFullGroup.interleavedMain), {retention: axialFullGroup.interleavedMain.retention});
    const keepSummary = sequence => {
      const rows = Object.values(sequence.retention), invalid = rows.filter(item => item.keepRatio === null || !Number.isFinite(item.keepRatio));
      return {valid: invalid.length === 0 && rows.length === 2, min: invalid.length === 0 && rows.length === 2 ? Math.min(...rows.map(item => item.keepRatio)) : null, ratios: Object.fromEntries(rows.map(item => [item.historyId, item.keepRatio])), invalid: invalid.map(item => ({historyId: item.historyId, reason: item.keepReason, B: item.B, D: item.D}))};
    };
    const fullMainKeep = keepSummary(axialFullGroup.interleavedMain), detachedMainKeep = keepSummary(axialDetachedGroup.interleavedMain), tonicMainKeep = keepSummary(tonicGroup.interleavedMain);
    const fullRhoKeep = keepSummary(axialFullGroup.interleavedRho), detachedRhoKeep = keepSummary(axialDetachedGroup.interleavedRho), tonicRhoKeep = keepSummary(tonicGroup.interleavedRho);
    const fullDetachedAttribution = blockOne.fullDetached.fraction !== null && blockOne.fullDetached.fraction > 1e-6;
    scientificCheck("structural attribution: full axial and detached credit differ on the actual margin", fullDetachedAttribution, {fullDetached: blockOne.fullDetached});
    scientificCheck("structural attribution: equal-rho full axial minimum interleaved keep exceeds detached", fullRhoKeep.valid && detachedRhoKeep.valid && fullRhoKeep.min > detachedRhoKeep.min, {fullRhoKeep, detachedRhoKeep, mainComparison: {fullMainKeep, detachedMainKeep}});
    scientificCheck("structural attribution: equal-rho full axial minimum interleaved keep exceeds tonic", fullRhoKeep.valid && tonicRhoKeep.valid && fullRhoKeep.min > tonicRhoKeep.min, {fullRhoKeep, tonicRhoKeep, mainComparison: {fullMainKeep, tonicMainKeep}});
    const scientificFailures = scientificChecks.filter(item => !item.passed);
    const correctnessFailures = validationChecks.filter(item => !item.passed && item.kind === "correctness");
    const metrics = {
      axialFull: {main: axialFullGroup.main.retention, rho: axialFullGroup.rho.retention, interleavedMain: axialFullGroup.interleavedMain.retention, interleavedRho: axialFullGroup.interleavedRho.retention, frozen: axialFullGroup.frozen.retention, interleavedFrozen: axialFullGroup.interleavedFrozen.retention, mainKeep: fullMainKeep, rhoKeep: fullRhoKeep},
      axialDetached: {main: axialDetachedGroup.main.retention, rho: axialDetachedGroup.rho.retention, interleavedMain: axialDetachedGroup.interleavedMain.retention, interleavedRho: axialDetachedGroup.interleavedRho.retention, frozen: axialDetachedGroup.frozen.retention, interleavedFrozen: axialDetachedGroup.interleavedFrozen.retention, mainKeep: detachedMainKeep, rhoKeep: detachedRhoKeep},
      tonic: {main: tonicGroup.main.retention, rho: tonicGroup.rho.retention, interleavedMain: tonicGroup.interleavedMain.retention, interleavedRho: tonicGroup.interleavedRho.retention, frozen: tonicGroup.frozen.retention, interleavedFrozen: tonicGroup.interleavedFrozen.retention, mainKeep: tonicMainKeep, rhoKeep: tonicRhoKeep},
      rhoHalfDirection,
      rhoHalfSignStability,
      forwardEquivalence: forwardEquivalenceData,
    };
    structureComputationQualification = {
      scope: "Fixed synthetic conductance-circuit and record-derived attachment arithmetic only; no RGB/HM3D navigation evidence, no visual teacher, no novelty claim.",
      parameters: {h, leak: p.leak, somaLeak: p.somaLeak, axial: p.axial, theta0, epsilon: structureEpsilon, tolerance: structureTolerance, strictTolerance, newtonMaxIterations: 40, eta: mainEta, rho, rhoHalf, finiteDifferenceStep, finiteDifferenceStepFine},
      fixture: {candidateA, candidateB, historyPrefixes, queryRecord, arrivals: arrivalTemplates, incomingContexts: {axial: axialHistories, tonic: tonicHistories}, teacherSource: "completed-record transition increments plus arrived successor observations only; no target sign or candidate identity is stored.", informationBoundary: {allowed: ["candidate x/y/action/arrival records", "query x/action record", "fixed incoming context", "theta and thetaVersion"], forbidden: ["pose", "depth", "overlap", "place/instance/scene identity", "reward", "success", "target sign", "candidate identity"], auditMeaning: "These declarations document the access contract; none is read as a fixture value or used to choose an update."}},
      validationBlocks: [blockOne, blockTwo],
      validationChecks,
      correctnessPassed: correctnessFailures.length === 0,
      correctnessFailures,
      scientificChecks,
      scientificFailures,
      metrics,
      verdict: "Layer 1 is a legal computation-interface result only when all correctness checks pass. Layers 2–4 report this fixed synthetic slice, interleaved retention, and separate structural attribution; none is a navigation, visual-generalization, or novelty result.",
    };
    if (selectiveCorrectionStageA) {
      const stageACstarExpected = 0.019658382399647234;
      const stageAFiniteDifferenceSteps = [1e-5, 3e-6];
      const stageAMainEta = 0.02, stageARho = 0.001;
      const stageGate = {maxSCross: 0.25, maxAbsC: 0.25, minNormalizedGramEigenvalue: 0.75, maxRelativeFiniteDifferenceError: 1e-4, maxAbsoluteFiniteDifferenceError: 1e-9, maxHOverU: 0.25, minDOverU: 0.5, minRhoBOverOriginalR: 0.25};
      const stageChecks = [], stageCorrectnessChecks = [], stageScientificChecks = [];
      const stageCheck = (name, passed, kind, details = {}) => {
        const row = {name, passed: Boolean(passed), kind, ...details};
        stageChecks.push(row);
        if (kind === "scientific") stageScientificChecks.push(row);
        else stageCorrectnessChecks.push(row);
        return Boolean(passed);
      };
      const stageQ = incoming => incoming.slice(0, 4).map(value => (value - incoming[4]) ** 2);
      const stageP = qValues => qValues.reduce((sum, value) => sum + value, 0) / qValues.length;
      const stageCstarQ = stageQ(axialHistories.H0.incoming), stageCstarP = stageP(stageCstarQ);
      const stageACstarComputed = Math.sqrt(stageCstarQ.reduce((sum, value) => sum + (value - stageCstarP) ** 2, 0) / stageCstarQ.length);
      stageCheck("Stage A frozen C* comes from the theta=0 H0 prefix before E0", Math.abs(stageACstarComputed - stageACstarExpected) <= 1e-12, "correctness", {expected: stageACstarExpected, computed: stageACstarComputed, q: stageCstarQ, P: stageCstarP, source: "H0 prefix after its first two observation tokens; before E0"});
      const stageModes = {
        full: {weightMode: "full", creditMode: "full", description: "alpha=k exp(theta*(1+p*e)); full conditional three-route credit"},
        differenceOnly: {weightMode: "difference", creditMode: "full", description: "alpha=k exp(theta*(p*e)); full conditional three-route credit"},
        detachedError: {weightMode: "full", creditMode: "detached", description: "full forward formula; detached branch-only error response diagnostic"},
        commonOnlyOriginalR: {weightMode: "common", creditMode: "full", description: "w=1 original-R coefficient path; full conditional three-route credit"},
      };
      const stageBasePoints = [
        {id: "cold_candidates", candidateA, candidateB, candidateStarts: {A: [...cold], B: [...cold]}},
        {id: "warm_candidates_H0_H1", candidateA, candidateB, candidateStarts: {A: [...axialHistories.H0.incoming], B: [...axialHistories.H1.incoming]}},
        {id: "warm_candidates_H1_H0", candidateA, candidateB, candidateStarts: {A: [...axialHistories.H1.incoming], B: [...axialHistories.H0.incoming]}},
      ];
      const stageCellDetails = cell => ({
        state: cell.state,
        drop: cell.drop,
        q: cell.weightData.q,
        P: cell.weightData.P,
        Cstar: cell.weightData.Cstar,
        e: cell.weightData.e,
        p: cell.weightData.connectionSigns,
        w: cell.weightData.weights,
        alpha: cell.coefficient,
        differentialConductance: cell.edge.map(edge => edge / h),
        diagonal: cell.diagonal,
        edge: cell.edge,
        carrier: cell.carrier,
        iterations: cell.iterations,
        residual: cell.residual,
      });
      const stageRouteDetails = routes => Object.fromEntries(Object.entries(routes).map(([route, value]) => [route, {
        cellKey: value.cellKey,
        b: value.b,
        lambda: value.lambda,
        fullLocalCredit: value.fullLocalCredit,
        detachedLocalCredit: value.detachedLocalCredit,
        localCredit: value.localCredit,
      }]));
      const stageMargin = (th, historyId, basePoint, weightMode, creditMode, audit = null, eventTolerance = structureTolerance) => {
        const context = {Cstar: stageACstarExpected, mode: weightMode};
        const queryCell = event(th, queryRecord.x, "axial", axialHistories[historyId].incoming, eventTolerance, structureEpsilon, audit, context);
        const referenceACell = event(th, basePoint.candidateA.x, "axial", basePoint.candidateStarts.A, eventTolerance, structureEpsilon, audit, context);
        const referenceBCell = event(th, basePoint.candidateB.x, "axial", basePoint.candidateStarts.B, eventTolerance, structureEpsilon, audit, context);
        return {...marginFromCells(queryCell, referenceACell, referenceBCell, creditMode, audit), queryCell, referenceACell, referenceBCell};
      };
      const stageSelectedGradient = (computed, creditMode) => creditMode === "detached" ? computed.detachedGradient : computed.fullGradient;
      const stageGradientErrorRows = (analytic, finite) => analytic.map((value, index) => {
        const absoluteError = Math.abs(value - finite[index]), nearZero = Math.abs(value) <= stageGate.maxAbsoluteFiniteDifferenceError;
        const relativeError = nearZero ? null : absoluteError / Math.abs(value);
        return {analytic: value, finite: finite[index], absoluteError, relativeError, nearZero, passed: absoluteError <= stageGate.maxAbsoluteFiniteDifferenceError && (nearZero || relativeError <= stageGate.maxRelativeFiniteDifferenceError)};
      });
      const stageFiniteDifference = (basePoint, historyId, weightMode, creditMode, step) => {
        const plus = [...theta0], minus = [...theta0];
        for (let j = 0; j < plus.length; j++) { plus[j] += step; minus[j] -= step; }
        const plusValue = stageMargin(plus, historyId, basePoint, weightMode, "full").value;
        const minusValue = stageMargin(minus, historyId, basePoint, weightMode, "full").value;
        const finite = theta0.map((_, j) => {
          const plusTheta = [...theta0], minusTheta = [...theta0];
          plusTheta[j] += step; minusTheta[j] -= step;
          return (stageMargin(plusTheta, historyId, basePoint, weightMode, "full").value - stageMargin(minusTheta, historyId, basePoint, weightMode, "full").value) / (2 * step);
        });
        return {step, plusValue, minusValue, finite, note: creditMode === "detached" ? "finite difference is compared with full conditional gradient; detached is retained as a diagnostic" : "central finite difference of the same fixed-incoming forward path"};
      };
      const stageSpectrum = gradients => {
        const K = gradients.map(left => gradients.map(right => dot(left, right)));
        const diagonalPositive = K.every((row, index) => Number.isFinite(row[index]) && row[index] > 0);
        if (!diagonalPositive) return {K, normalizedGram: null, eigenvalues: null, conditionNumber: null, effectiveRank: null, c: null, minEigenvalue: null, diagonalPositive};
        const normalizedGram = [[1, K[0][1] / Math.sqrt(K[0][0] * K[1][1])], [K[1][0] / Math.sqrt(K[0][0] * K[1][1]), 1]];
        const c = normalizedGram[0][1], minEigenvalue = 1 - Math.abs(c), maxEigenvalue = 1 + Math.abs(c);
        return {K, normalizedGram, eigenvalues: [minEigenvalue, maxEigenvalue], conditionNumber: maxEigenvalue / minEigenvalue, effectiveRank: 2 / (1 + c * c), c, minEigenvalue, diagonalPositive};
      };
      const stagePanels = {};
      for (const basePoint of stageBasePoints) {
        stagePanels[basePoint.id] = {};
        for (const [modeId, mode] of Object.entries(stageModes)) {
          const rows = {};
          for (const historyId of ["H0", "H1"]) {
            const audit = freshAudit(), computed = stageMargin(theta0, historyId, basePoint, mode.weightMode, mode.creditMode, audit);
            const selectedGradient = stageSelectedGradient(computed, mode.creditMode), finiteDifferences = stageAFiniteDifferenceSteps.map(step => {
              const fd = stageFiniteDifference(basePoint, historyId, mode.weightMode, mode.creditMode, step);
              const target = computed.fullGradient, errors = stageGradientErrorRows(target, fd.finite);
              return {...fd, target, errors, passed: errors.every(row => row.passed)};
            });
            rows[historyId] = {
              historyId,
              margin: computed.value,
              gradient: selectedGradient,
              fullConditionalGradient: computed.fullGradient,
              detachedGradient: computed.detachedGradient,
              routeCredits: stageRouteDetails(computed.routes),
              cells: {query: stageCellDetails(computed.queryCell), referenceA: stageCellDetails(computed.referenceACell), referenceB: stageCellDetails(computed.referenceBCell)},
              finiteDifferences,
              finiteDifferencePassed: finiteDifferences.every(item => item.passed),
              selectedGradientFiniteDifferenceGap: gap(selectedGradient, finiteDifferences[0].finite),
              audit,
            };
          }
          const spectrum = stageSpectrum([rows.H0.gradient, rows.H1.gradient]);
          stagePanels[basePoint.id][modeId] = {basePoint: basePoint.id, mode: modeId, weightMode: mode.weightMode, creditMode: mode.creditMode, description: mode.description, rows, ...spectrum, finiteDifferencePassed: Object.values(rows).every(row => row.finiteDifferencePassed)};
        }
      }
      const stageBasePointSummary = Object.fromEntries(stageBasePoints.map(basePoint => {
        const panel = stagePanels[basePoint.id].full;
        const passed = panel.diagonalPositive && panel.minEigenvalue >= stageGate.minNormalizedGramEigenvalue && panel.c !== null && Math.abs(panel.c) <= stageGate.maxAbsC && Math.max(Math.abs(panel.K[0][1]) / panel.K[0][0], Math.abs(panel.K[0][1]) / panel.K[1][1]) <= stageGate.maxSCross;
        stageCheck(`${basePoint.id} full panel passes S_cross, |c| and normalized Gram spectrum`, passed, "scientific", {basePoint: basePoint.id, S_cross: Math.max(Math.abs(panel.K[0][1]) / panel.K[0][0], Math.abs(panel.K[0][1]) / panel.K[1][1]), c: panel.c, minNormalizedGramEigenvalue: panel.minEigenvalue, thresholds: {S_cross: stageGate.maxSCross, absC: stageGate.maxAbsC, minEigenvalue: stageGate.minNormalizedGramEigenvalue}});
        stageCheck(`${basePoint.id} full panel finite differences satisfy absolute/relative tolerances`, panel.finiteDifferencePassed, "correctness", {basePoint: basePoint.id, finiteDifferencePassed: panel.finiteDifferencePassed});
        return [basePoint.id, {S_cross: Math.max(Math.abs(panel.K[0][1]) / panel.K[0][0], Math.abs(panel.K[0][1]) / panel.K[1][1]), c: panel.c, minNormalizedGramEigenvalue: panel.minEigenvalue, eigenvalues: panel.eigenvalues, conditionNumber: panel.conditionNumber, effectiveRank: panel.effectiveRank, passed}];
      }));
      const stageForwardChecks = [];
      for (const basePoint of stageBasePoints) {
        for (const historyId of ["H0", "H1"]) {
          const newComputed = stageMargin(theta0, historyId, basePoint, "full", "full");
          const oldCells = {
            query: event(theta0, queryRecord.x, "axial", axialHistories[historyId].incoming, structureTolerance, structureEpsilon),
            referenceA: event(theta0, basePoint.candidateA.x, "axial", basePoint.candidateStarts.A, structureTolerance, structureEpsilon),
            referenceB: event(theta0, basePoint.candidateB.x, "axial", basePoint.candidateStarts.B, structureTolerance, structureEpsilon),
          };
          const forwardGap = Math.max(...["query", "referenceA", "referenceB"].map(route => routeCellGap(oldCells[route], newComputed[`${route}Cell`])));
          stageForwardChecks.push({basePoint: basePoint.id, historyId, gap: forwardGap});
        }
      }
      stageCheck("theta=0 selective forward cells equal the original R cells at all three basepoints", stageForwardChecks.every(item => item.gap <= 1e-12), "correctness", {checks: stageForwardChecks});
      const coldPanel = stagePanels.cold_candidates.full;
      const coldWeights = [coldPanel.rows.H0.cells.referenceA.w, coldPanel.rows.H0.cells.referenceB.w, coldPanel.rows.H1.cells.referenceA.w, coldPanel.rows.H1.cells.referenceB.w];
      const coldCandidateGradients = [coldPanel.rows.H0.routeCredits.referenceA.fullLocalCredit, coldPanel.rows.H0.routeCredits.referenceB.fullLocalCredit, coldPanel.rows.H1.routeCredits.referenceA.fullLocalCredit, coldPanel.rows.H1.routeCredits.referenceB.fullLocalCredit];
      stageCheck("cold candidates retain w=1 and nonzero A/B conditional derivatives", coldWeights.every(values => values.every(value => value === 1)) && coldCandidateGradients.every(values => norm(values) > 1e-12), "correctness", {weights: coldWeights, candidateGradientNorms: coldCandidateGradients.map(norm)});
      const stagePending = (th, historyId, basePoint, weightMode, creditMode, version = 0, candidateBOverride = null) => {
        const effectiveBasePoint = candidateBOverride === null ? basePoint : {...basePoint, candidateB: candidateBOverride};
        const computed = stageMargin(th, historyId, effectiveBasePoint, weightMode, creditMode);
        const effectiveCandidateB = effectiveBasePoint.candidateB;
        return {
          ...computed,
          thetaVersion: version,
          candidateRecords: {A: effectiveBasePoint.candidateA, B: effectiveCandidateB},
          queryRecord: {...queryRecord, arrival: null},
          incoming: [...axialHistories[historyId].incoming],
          historyId,
          incrementA: vectorSub(effectiveBasePoint.candidateA.y, effectiveBasePoint.candidateA.x),
          incrementB: vectorSub(effectiveCandidateB.y, effectiveCandidateB.x),
          margin: computed.value,
          marginGradient: stageSelectedGradient(computed, creditMode),
          consumed: false,
        };
      };
      const stageConsume = (pending, y, th, version, budget) => {
        if (pending.consumed) throw new Error("stage A pending event already consumed");
        if (pending.thetaVersion !== version) throw new Error("stage A pending theta version mismatch");
        const evidence = relativeEvidence(pending, y), gradient = pending.marginGradient, gradientNorm = Math.sqrt(dot(gradient, gradient));
        let deltaTheta;
        if (budget.kind === "eta") deltaTheta = gradient.map(value => budget.value * evidence.delta * value);
        else {
          const direction = gradientNorm === 0 || evidence.delta === 0 ? gradient.map(() => 0) : gradient.map(value => Math.sign(evidence.delta) * value / gradientNorm);
          deltaTheta = direction.map(value => budget.value * value);
        }
        if (!Number.isFinite(evidence.delta) || !gradient.every(Number.isFinite) || !deltaTheta.every(Number.isFinite)) throw new Error(`stage A update out of domain: budget=${JSON.stringify(budget)} evidence=${JSON.stringify({margin: pending.margin, lambda: evidence.lambda, pBefore: evidence.pBefore, pPlus: evidence.pPlus, delta: evidence.delta})} gradient=${JSON.stringify(gradient)} deltaTheta=${JSON.stringify(deltaTheta)}`);
        pending.consumed = true;
        return {evidence, gradient, gradientNorm, deltaTheta, thetaAfter: vectorAdd(th, deltaTheta), thetaVersionAfter: version + 1};
      };
      const tryStageError = operation => { try { operation(); return {threw: false, error: null}; } catch (error) { return {threw: true, error: String(error.message || error)}; } };
      const equalCandidateB = {...candidateB, id: "B_equal_increment_stage_a", y: vectorAdd(candidateB.x, vectorSub(candidateA.y, candidateA.x))};
      const equalPending = stagePending(theta0, "H0", stageBasePoints[0], "full", "full", 0, equalCandidateB), equalEvidence = relativeEvidence(equalPending, equalCandidateB.y), equalConsumed = stageConsume(equalPending, equalCandidateB.y, theta0, 0, {kind: "eta", value: stageAMainEta});
      stageCheck("zero evidence leaves the selective parameter write exactly zero", Math.abs(equalEvidence.lambda) <= 1e-14 && Math.abs(equalEvidence.delta) <= 1e-14 && equalConsumed.deltaTheta.every(value => value === 0), "correctness", {lambda: equalEvidence.lambda, delta: equalEvidence.delta, deltaTheta: equalConsumed.deltaTheta, candidateBInputRemainsDistinct: gap(equalCandidateB.x, candidateA.x) > 0});
      const stalePending = stagePending(theta0, "H0", stageBasePoints[0], "full", "full", 0);
      const staleVersion = tryStageError(() => stageConsume(stalePending, arrivalTemplates[0].y, theta0, 1, {kind: "eta", value: stageAMainEta}));
      stageCheck("selective pending consumption rejects a stale theta version", staleVersion.threw, "correctness", staleVersion);
      const duplicatePending = stagePending(theta0, "H0", stageBasePoints[0], "full", "full", 0);
      stageConsume(duplicatePending, arrivalTemplates[0].y, theta0, 0, {kind: "eta", value: stageAMainEta});
      const duplicateConsumption = tryStageError(() => stageConsume(duplicatePending, arrivalTemplates[0].y, vectorAdd(theta0, duplicatePending.marginGradient), 1, {kind: "eta", value: stageAMainEta}));
      stageCheck("selective pending consumption rejects duplicate consumption", duplicateConsumption.threw, "correctness", duplicateConsumption);
      stageCheck("selective pending cache has no unrevealed successor y", !Object.prototype.hasOwnProperty.call(stalePending, "y"), "correctness", {pendingKeys: Object.keys(stalePending).filter(key => key === "y")});
      const stageSnapshot = (th, basePoint, weightMode, creditMode) => {
        const rows = {};
        for (const historyId of ["H0", "H1"]) {
          const computed = stageMargin(th, historyId, basePoint, weightMode, creditMode);
          rows[historyId] = {margin: computed.value, fullGradient: computed.fullGradient, detachedGradient: computed.detachedGradient, selectedGradient: stageSelectedGradient(computed, creditMode)};
        }
        return {rows};
      };
      const compactStageEvidence = evidence => ({y: evidence.y, w: evidence.w, residualA: evidence.residualA, residualB: evidence.residualB, squaredResidualA: evidence.squaredResidualA, squaredResidualB: evidence.squaredResidualB, lambda: evidence.lambda, pBefore: evidence.pBefore, pPlus: evidence.pPlus, delta: evidence.delta});
      const stageRunTwoEvent = ({basePoint, relation, order, weightMode, creditMode, budget}) => {
        let currentTheta = [...theta0], currentVersion = 0;
        const initial = stageSnapshot(currentTheta, basePoint, weightMode, creditMode), eventRows = [];
        for (let step = 0; step < order.length; step++) {
          const historyId = order[step], pending = stagePending(currentTheta, historyId, basePoint, weightMode, creditMode, currentVersion);
          const before = stageSnapshot(currentTheta, basePoint, weightMode, creditMode), y = relation.yByHistory[historyId];
          const update = stageConsume(pending, y, currentTheta, currentVersion, budget), after = stageSnapshot(update.thetaAfter, basePoint, weightMode, creditMode);
          const targetEffects = {};
          for (const targetHistoryId of ["H0", "H1"]) {
            const actualChange = after.rows[targetHistoryId].margin - before.rows[targetHistoryId].margin;
            targetEffects[targetHistoryId] = {historyId: targetHistoryId, fBefore: before.rows[targetHistoryId].margin, fAfter: after.rows[targetHistoryId].margin, actualChange, firstOrderFromActualDeltaTheta: dot(before.rows[targetHistoryId].fullGradient, update.deltaTheta), ownEvent: targetHistoryId === historyId};
          }
          const logAlphaDelta = Object.fromEntries(Object.entries(pending.rawCells).map(([route, cell]) => [route, cell.weights.map((weight, j) => weight * update.deltaTheta[j])]));
          eventRows.push({
            eventId: `${relation.id}_${order.join("_")}_${historyId}_${step}`,
            historyId,
            arrival: 3 + step,
            action: queryRecord.action,
            y,
            thetaVersionBefore: currentVersion,
            thetaBefore: currentTheta,
            fBefore: Object.fromEntries(Object.entries(before.rows).map(([id, row]) => [id, row.margin])),
            gReadoutFull: before.rows[historyId].fullGradient,
            gUpdateFull: pending.fullGradient,
            gUpdateDetached: pending.detachedGradient,
            gUpdateSelected: update.gradient,
            routeGradients: stageRouteDetails(pending.routes),
            cells: {query: stageCellDetails(pending.queryCell), referenceA: stageCellDetails(pending.referenceACell), referenceB: stageCellDetails(pending.referenceBCell)},
            evidence: compactStageEvidence(update.evidence),
            deltaTheta: update.deltaTheta,
            deltaLogAlpha: logAlphaDelta,
            parameterWriteNorm: Math.sqrt(dot(update.deltaTheta, update.deltaTheta)),
            thetaAfter: update.thetaAfter,
            thetaVersionAfter: update.thetaVersionAfter,
            fAfter: Object.fromEntries(Object.entries(after.rows).map(([id, row]) => [id, row.margin])),
            targetEffects,
            diagnostic: budget.kind === "rho" ? {rho: budget.value, gradientNorm: update.gradientNorm} : {eta: budget.value},
          });
          currentTheta = update.thetaAfter;
          currentVersion += 1;
        }
        const final = stageSnapshot(currentTheta, basePoint, weightMode, creditMode), signs = Object.fromEntries(["H0", "H1"].map(historyId => {
          const own = eventRows.find(row => row.historyId === historyId);
          return [historyId, own === undefined || own.evidence.lambda === 0 ? null : Math.sign(own.evidence.lambda)];
        }));
        const retention = {};
        for (const historyId of ["H0", "H1"]) {
          const sign = signs[historyId], ownEvents = eventRows.filter(row => row.historyId === historyId), otherEvents = eventRows.filter(row => row.historyId !== historyId);
          const U = sign === null ? null : ownEvents.reduce((sum, row) => sum + sign * row.targetEffects[historyId].actualChange, 0);
          const H = sign === null ? null : otherEvents.reduce((sum, row) => sum + Math.max(0, -sign * row.targetEffects[historyId].actualChange), 0);
          const D = sign === null ? null : sign * (final.rows[historyId].margin - initial.rows[historyId].margin);
          const B = U;
          retention[historyId] = {historyId, sign, fInitial: initial.rows[historyId].margin, fFinal: final.rows[historyId].margin, U, H, D, B, HOverU: U === null || U <= 0 ? null : H / U, DOverU: U === null || U <= 0 ? null : D / U};
        }
        const directionalEffects = eventRows.flatMap(row => ["H0", "H1"].map(historyId => ({eventId: row.eventId, eventHistoryId: row.historyId, historyId, ownEvent: row.historyId === historyId, actualChange: row.targetEffects[historyId].actualChange, directionalChange: signs[historyId] === null ? null : signs[historyId] * row.targetEffects[historyId].actualChange})));
        const gateRows = Object.values(retention).map(row => ({historyId: row.historyId, UPositive: row.U !== null && row.U > 0, DPositive: row.D !== null && row.D > 0, HOverUPassed: row.HOverU !== null && row.HOverU <= stageGate.maxHOverU, DOverUPassed: row.DOverU !== null && row.DOverU >= stageGate.minDOverU}));
        return {budget, weightMode, creditMode, thetaInitial: theta0, thetaFinal: currentTheta, initialMargins: Object.fromEntries(Object.entries(initial.rows).map(([id, row]) => [id, row.margin])), finalMargins: Object.fromEntries(Object.entries(final.rows).map(([id, row]) => [id, row.margin])), events: eventRows, directionalEffects, retention, gateRows, gatePassed: gateRows.every(row => row.UPositive && row.DPositive && row.HOverUPassed && row.DOverUPassed), finalVersion: currentVersion};
      };
      const compactReferenceSequence = sequence => ({retention: sequence.retention, gateRows: sequence.gateRows, gatePassed: sequence.gatePassed, firstEventB: Object.fromEntries(Object.entries(sequence.retention).map(([historyId, row]) => [historyId, row.B]))});
      const yPlus = [0.9, 0.23], yMinus = [0.7, 0.23];
      const stageRelations = [
        {id: "yplus_yplus", yByHistory: {H0: [...yPlus], H1: [...yPlus]}},
        {id: "yminus_yminus", yByHistory: {H0: [...yMinus], H1: [...yMinus]}},
        {id: "yplus_yminus", yByHistory: {H0: [...yPlus], H1: [...yMinus]}},
        {id: "yminus_yplus", yByHistory: {H0: [...yMinus], H1: [...yPlus]}},
      ];
      const stageOrders = [{id: "H0_then_H1", historyIds: ["H0", "H1"]}, {id: "H1_then_H0", historyIds: ["H1", "H0"]}];
      const stageUnits = [];
      for (const basePoint of stageBasePoints) for (const relation of stageRelations) for (const order of stageOrders) {
        const etaSequence = stageRunTwoEvent({basePoint, relation, order: order.historyIds, weightMode: "full", creditMode: "full", budget: {kind: "eta", value: stageAMainEta}});
        const rhoSequence = stageRunTwoEvent({basePoint, relation, order: order.historyIds, weightMode: "full", creditMode: "full", budget: {kind: "rho", value: stageARho}});
        const originalRho = stageRunTwoEvent({basePoint, relation, order: order.historyIds, weightMode: "common", creditMode: "full", budget: {kind: "rho", value: stageARho}});
        const rhoBComparison = Object.fromEntries(["H0", "H1"].map(historyId => {
          const scB = rhoSequence.retention[historyId].B, originalB = originalRho.retention[historyId].B;
          return [historyId, {scB, originalRhoB: originalB, ratio: originalB === null || originalB <= 0 ? null : scB / originalB, passed: originalB !== null && originalB > 0 && scB >= stageGate.minRhoBOverOriginalR * originalB}];
        }));
        const gatePassed = etaSequence.gatePassed && rhoSequence.gatePassed && Object.values(rhoBComparison).every(row => row.passed);
        const unit = {
          id: `${basePoint.id}__${relation.id}__${order.id}`,
          basePoint: basePoint.id,
          relation: relation.id,
          yByHistory: relation.yByHistory,
          order: order.historyIds,
          eta: etaSequence,
          rho: rhoSequence,
          originalRhoReference: compactReferenceSequence(originalRho),
          rhoBComparison,
          gatePassed,
        };
        stageUnits.push(unit);
        stageCheck(`${unit.id} satisfies eta/rho U,D,H/U,D/U and rho B reference gates`, gatePassed, "scientific", {unit: unit.id, eta: {retention: etaSequence.retention, gateRows: etaSequence.gateRows}, rho: {retention: rhoSequence.retention, gateRows: rhoSequence.gateRows}, rhoBComparison});
      }
      const correctnessPassed = stageCorrectnessChecks.every(item => item.passed), scientificGatePassed = stageScientificChecks.every(item => item.passed), stageAGatePassed = correctnessPassed && scientificGatePassed;
      selectiveCorrectionQualification = {
        scope: "Stage A fixed synthetic selective-correction arithmetic only; no RGB/HM3D evidence, no visual teacher, no Stage B implementation.",
        formula: {alpha: "k*exp(theta_j*(1+p_j*e_j))", current: "I_j=k*(z_j-s)+alpha_j*(z_j-s)^3", connectionSigns: [...connectionSigns], q: "(z_j_minus-s_minus)^2", P: "mean(q)", e: "(q-P)/Cstar", Cstar: stageACstarExpected, incomingFixedDuringNewton: true, wComputedOnceAtEventEntrance: true, theta0RestoresOriginalR: true},
        parameters: {h, epsilon: structureEpsilon, theta0, eta: stageAMainEta, rho: stageARho, finiteDifferenceSteps: stageAFiniteDifferenceSteps, tolerances: {event: structureTolerance, strict: strictTolerance}, gate: stageGate},
        fixture: {candidateA: {id: candidateA.id, x: candidateA.x, y: candidateA.y, action: candidateA.action, arrival: candidateA.arrival}, candidateB: {id: candidateB.id, x: candidateB.x, y: candidateB.y, action: candidateB.action, arrival: candidateB.arrival}, queryRecord, historyPrefixes, historyIncoming: {H0: axialHistories.H0.incoming, H1: axialHistories.H1.incoming}, basePoints: stageBasePoints.map(basePoint => ({id: basePoint.id, candidateStarts: basePoint.candidateStarts})), evidenceRelations: stageRelations, orders: stageOrders, CstarSource: {expected: stageACstarExpected, computed: stageACstarComputed, q: stageCstarQ, P: stageCstarP, source: "H0 prefix after two theta=0 observation tokens, before E0"}},
        modes: stageModes,
        basePointSummary: stageBasePointSummary,
        basePointPanels: stagePanels,
        units: stageUnits,
        counts: {basePoints: stageBasePoints.length, relations: stageRelations.length, orders: stageOrders.length, units: stageUnits.length, budgetsPerUnit: 2, expectedUnits: 24},
        checks: stageChecks,
        correctnessChecks: stageCorrectnessChecks,
        scientificChecks: stageScientificChecks,
        correctnessPassed,
        scientificGatePassed,
        coreGatePassed: stageAGatePassed,
        stageB: {implemented: false, eligibleByCoreGate: stageAGatePassed, requiresExplicitRelease: true, stoppedOnFailure: !stageAGatePassed},
        legacyRegressionReference: {oldConstructiveChecks: 132, oldAssociationChecks: 118, oldStructureOutputSha256: "CE124DC26B855A5E96EB91BB0F812D89BC3C5BDA8839C4CF4B2F8C20691E6E56"},
        verdict: stageAGatePassed ? "Stage A core arithmetic and frozen scientific gates passed; this is not Stage B, visual evidence, navigation evidence, or Goal completion." : "Stage A retained all raw basepoint/unit quantities and at least one frozen correctness/scientific gate failed; stop before Stage B and return to Astra.",
      };
    }
  }
}
const output = {
  scope: "Fixed nondimensional conductance-circuit and attachment arithmetic checks only; no visual learning, biological measurement, or novelty proof.",
  parameters: p, path, checksPassed: checks.length, checks,
  results: {base, swapped, noModulation, sameStream, negativeEpsilon, epsilonChecks, discretizationGap, stopped, slower, mismatched, omega, contraction, residual, counterexamples, conditionalRisk, couplingChecks, lagChecks},
  learningQualification,
  ...(continuityKillOnly ? {continuityQualification} : {}),
  ...(repairKillOnly ? {repairQualification} : {}),
  ...(noveltyKillOnly ? {noveltyQualification} : {}),
  ...(visualCreditKillOnly ? {visualCreditQualification} : {}),
  ...(associationCreditKillOnly ? {associationCreditQualification} : {}),
  ...(constructiveCouplingOnly ? {constructiveCouplingQualification} : {}),
  ...(structureComputationOnly ? {structureComputationQualification} : {}),
  continuousCreditQualification,
};
const rendered = JSON.stringify(selectiveCorrectionStageA ? {
  scope: "Stage A selective-correction raw core only; old arithmetic outputs are referenced, not embedded.",
  selectiveCorrectionQualification,
} : output, null, 2) + "\n";
if (selectiveCorrectionStageA) {
  if (!noWriteLegacy) fs.writeFileSync("idea-stage/SELECTIVE_CORRECTION_CHECK_20260909.json", rendered, "utf8");
  if (!selectiveCorrectionQualification.correctnessPassed || !selectiveCorrectionQualification.scientificGatePassed) process.exitCode = 2;
} else if (structureComputationOnly) {
  const structureRendered = JSON.stringify(output, null, 2) + "\n";
  if (!noWriteLegacy) fs.writeFileSync("idea-stage/STRUCTURE_COMPUTATION_CHECK_20260909.json", structureRendered, "utf8");
  if (!structureComputationQualification.correctnessPassed) process.exitCode = 1;
}
process.stdout.write(rendered);
  }
}

module.exports = {
  freeToken,
  pairToken,
  associationEnergy,
  associationMargin,
  coverageTeacher,
  createOnlineState,
  beginRelationPending,
  consumeRelationPending,
  distanceMargin: relationDistanceMargin,
  runRelationEnergyCheck,
  runRelationEnergySelectivityCheck,
  sharedImplicitCell,
  sharedResponse,
};
