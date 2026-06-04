/**
 * RMEE 相対軌道可視化 — MEE / RMEE 定義・変換・軌道可視化
 */

const CONFIG = {
  MU_EARTH_KM3_S2: 398600.4418,
  R_EARTH_KM: 6378.137,
  NUM_ORBIT_POINTS: 400,
  ECI_CHIEF_ORBIT_POINTS: 360,
  EARTH_EQUATOR_POINTS: 128,
  DEFAULT_NUM_DRIFT_ORBITS: 5,
  DEFAULT_ECI_EXAGGERATION_LOG: 3,
  ECI_ZOOM_PADDING_RATIO: 0.2,
  AXIS_PADDING: 0.15,
  DEFAULT_CHIEF_COE: {
    a_km: 42164,
    e: 0,
    i_deg: 0,
    raan_deg: 0,
    argp_deg: 0,
    M_deg: 0,
  },
  DEFAULT_DEPUTY_COE: {
    a_km: 42165,
    e: 0.000237,
    i_deg: 0.0407,
    raan_deg: 0,
    argp_deg: 0,
    M_deg: 0.1359,
  },
};

const TWO_PI = 2 * Math.PI;
const MEE_SPECS = [
  { key: "p_km", suffix: "p", label: "p", unit: "km", step: "1" },
  { key: "L_deg", suffix: "L", label: "L", unit: "°", step: "0.001" },
  { key: "f", suffix: "f", label: "f", unit: "—", step: "0.000001" },
  { key: "g", suffix: "g", label: "g", unit: "—", step: "0.000001" },
  { key: "h", suffix: "h", label: "h", unit: "—", step: "0.000001" },
  { key: "k", suffix: "k", label: "k", unit: "—", step: "0.000001" },
];
const RMEE_SPECS = [
  { key: "dp_over_pc", suffix: "dp-pc", label: "δp/p_c", unit: "—", step: "0.00000001" },
  { key: "dL_deg", suffix: "dL", label: "δL", unit: "°", step: "0.001" },
  { key: "df", suffix: "df", label: "δf", unit: "—", step: "0.000001" },
  { key: "dg", suffix: "dg", label: "δg", unit: "—", step: "0.000001" },
  { key: "dh", suffix: "dh", label: "δh", unit: "—", step: "0.000001" },
  { key: "dk", suffix: "dk", label: "δk", unit: "—", step: "0.000001" },
];
const DRIFT_LAMBDA_FACTOR = 1.5;
const RTN_PLOT_IDS = ["plot-3d", "plot-tr", "plot-tn", "plot-rn"];
const ECI_PLOT_IDS = ["plot-eci-3d", "plot-eci-xy", "plot-eci-xz", "plot-eci-yz"];
const PLOT_DIV_IDS = [...RTN_PLOT_IDS, ...ECI_PLOT_IDS];
const VIEW_PLOT_IDS = { rtn: RTN_PLOT_IDS, eci: ECI_PLOT_IDS };
const PLOT_LAYOUT_BASE = {
  paper_bgcolor: "#1c2333",
  plot_bgcolor: "#1c2333",
  font: { color: "#e6edf3", family: "Segoe UI, sans-serif", size: 12 },
  margin: { l: 50, r: 20, t: 30, b: 45 },
};

let activeMode = "coe-mee";
/** @type {"rtn" | "eci"} */
let activeView = "eci";
/** RMEE 系タブで最後に選んだ RTN / ECI */
let relativePlotView = "rtn";
let plotsInitialized = false;
let lastPlotStructureKey = null;

function degToRad(d) {
  return (d * Math.PI) / 180;
}
function radToDeg(r) {
  return (r * 180) / Math.PI;
}
function normalizeAngleRad(r) {
  return ((r % TWO_PI) + TWO_PI) % TWO_PI;
}
function normalizeAngleSigned(r) {
  let x = normalizeAngleRad(r);
  if (x > Math.PI) x -= TWO_PI;
  return x;
}
function formatNum(v, digits = 6) {
  if (!Number.isFinite(v)) return "—";
  if (Math.abs(v) < 1e-12) return "0";
  if (Math.abs(v) >= 1e4 || Math.abs(v) < 1e-4) return v.toExponential(4);
  return v.toFixed(digits);
}

function formatKm(value) {
  if (!Number.isFinite(value)) return "—";
  if (Math.abs(value) < 1e-9) return "0";
  const abs = Math.abs(value);
  if (abs >= 100 || abs < 0.01) return value.toExponential(2);
  return value.toFixed(2);
}

function coeKmDegToRad(c) {
  return {
    a_km: c.a_km,
    e: c.e,
    i_rad: degToRad(c.i_deg),
    raan_rad: degToRad(c.raan_deg),
    argp_rad: degToRad(c.argp_deg),
    M_rad: degToRad(c.M_deg),
  };
}
function coeRadToKmDeg(c) {
  return {
    a_km: c.a_km,
    e: c.e,
    i_deg: radToDeg(c.i_rad),
    raan_deg: radToDeg(c.raan_rad),
    argp_deg: radToDeg(c.argp_rad),
    M_deg: radToDeg(c.M_rad),
  };
}

function solveKeplerE(M_rad, e, tol = 1e-12) {
  const M = normalizeAngleRad(M_rad);
  if (e < 1e-12) return M;
  let E = e < 0.8 ? M : Math.PI;
  for (let i = 0; i < 60; i++) {
    const f = E - e * Math.sin(E) - M;
    const fp = 1 - e * Math.cos(E);
    E -= f / fp;
    if (Math.abs(f / fp) < tol) break;
  }
  return E;
}

function trueAnomalyFromMeanAnomaly(M_rad, e) {
  if (e < 1e-12) return normalizeAngleRad(M_rad);
  const E = solveKeplerE(M_rad, e);
  const beta = Math.sqrt(1 - e * e);
  return Math.atan2(beta * Math.sin(E), Math.cos(E) - e);
}

function meanAnomalyFromTrueAnomaly(nu_rad, e) {
  if (e < 1e-12) return normalizeAngleRad(nu_rad);
  const E = Math.atan2(
    Math.sqrt(1 - e * e) * Math.sin(nu_rad),
    e + Math.cos(nu_rad)
  );
  return normalizeAngleRad(E - e * Math.sin(E));
}

/** COE → Modified Equinoctial Elements α = [p, L, f, g, h, k] */
function coeToMee(coe) {
  const { a_km, e, i_rad, raan_rad, argp_rad, M_rad } = coe;
  const Omw = raan_rad + argp_rad;
  return {
    p_km: a_km * (1 - e * e),
    f: e * Math.cos(Omw),
    g: e * Math.sin(Omw),
    h: Math.tan(i_rad / 2) * Math.cos(raan_rad),
    k: Math.tan(i_rad / 2) * Math.sin(raan_rad),
    L_rad: normalizeAngleRad(raan_rad + argp_rad + M_rad),
  };
}

function meeToKmDeg(mee) {
  return {
    p_km: mee.p_km,
    f: mee.f,
    g: mee.g,
    h: mee.h,
    k: mee.k,
    L_deg: radToDeg(mee.L_rad),
  };
}

/** MEE → COE（L = Ω + ω + M から M を復元） */
function meeToCoe(mee) {
  const { p_km, f, g, h, k, L_rad } = mee;
  const e = Math.hypot(f, g);
  const denom = Math.max(1 - e * e, 1e-12);
  const a_km = p_km / denom;
  const i_rad = 2 * Math.atan(Math.hypot(h, k));
  const raan_rad = Math.atan2(k, h);
  const argp_rad = normalizeAngleRad(Math.atan2(g, f) - raan_rad);
  const M_rad = normalizeAngleRad(L_rad - raan_rad - argp_rad);
  return { a_km, e, i_rad, raan_rad, argp_rad, M_rad };
}

function meeKmDegToRad(m) {
  return {
    p_km: m.p_km,
    f: m.f,
    g: m.g,
    h: m.h,
    k: m.k,
    L_rad: degToRad(m.L_deg),
  };
}

/**
 * RMEE: δα = [δp/p_c, δL, δf, δg, δh, δk]^T
 * （第1成分のみ chief の p で正規化）
 */
function computeRmee(chiefMee, deputyMee) {
  const pc = chiefMee.p_km;
  return {
    dp_over_pc: pc !== 0 ? (deputyMee.p_km - pc) / pc : 0,
    df: deputyMee.f - chiefMee.f,
    dg: deputyMee.g - chiefMee.g,
    dh: deputyMee.h - chiefMee.h,
    dk: deputyMee.k - chiefMee.k,
    dL_rad: normalizeAngleSigned(deputyMee.L_rad - chiefMee.L_rad),
  };
}

/** RMEE 入力から Deputy MEE を復元（δp/p_c 正規化成分に対応） */
function deputyMeeFromRmee(chiefMee, rmee) {
  return {
    p_km: chiefMee.p_km * (1 + rmee.dp_over_pc),
    f: chiefMee.f + rmee.df,
    g: chiefMee.g + rmee.dg,
    h: chiefMee.h + rmee.dh,
    k: chiefMee.k + rmee.dk,
    L_rad: normalizeAngleRad(chiefMee.L_rad + rmee.dL_rad),
  };
}

function coeResultRows(coeKm) {
  return [
    { label: "<var>a</var>", value: formatNum(coeKm.a_km, 2), unit: "km" },
    { label: "<var>e</var>", value: formatNum(coeKm.e, 6) },
    { label: "<var>i</var>", value: formatNum(coeKm.i_deg, 4), unit: "°" },
    { label: "<var>Ω</var>", value: formatNum(coeKm.raan_deg, 2), unit: "°" },
    { label: "<var>ω</var>", value: formatNum(coeKm.argp_deg, 2), unit: "°" },
    { label: "<var>M</var>", value: formatNum(coeKm.M_deg, 4), unit: "°" },
  ];
}

function meeResultRows(mee) {
  const m = meeToKmDeg(mee);
  return [
    { label: "<var>p</var>", value: formatNum(m.p_km, 2), unit: "km" },
    { label: "<var>L</var>", value: formatNum(m.L_deg, 4), unit: "°" },
    { label: "<var>f</var>", value: formatNum(m.f) },
    { label: "<var>g</var>", value: formatNum(m.g) },
    { label: "<var>h</var>", value: formatNum(m.h) },
    { label: "<var>k</var>", value: formatNum(m.k) },
  ];
}

/** Gim–Alfriend 準非特異 ROE（無次元） */
function coeToRoe(chief, deputy) {
  const { a_km: a, e, i_rad: i, raan_rad: Om, argp_rad: w, M_rad: M } = chief;
  const {
    a_km: ad,
    e: ed,
    i_rad: id,
    raan_rad: Omd,
    argp_rad: wd,
    M_rad: Md,
  } = deputy;
  return {
    delta_a: (ad - a) / a,
    delta_lambda: Md - M + (Omd - Om) * Math.cos(i),
    delta_ex: ed * Math.cos(wd) - e * Math.cos(w),
    delta_ey: ed * Math.sin(wd) - e * Math.sin(w),
    delta_ix: id - i,
    delta_iy: (Omd - Om) * Math.sin(i),
  };
}

function roeKmFromDimensionless(roe, a_km) {
  return {
    delta_a: roe.delta_a * a_km,
    delta_lambda: roe.delta_lambda * a_km,
    delta_ex: roe.delta_ex * a_km,
    delta_ey: roe.delta_ey * a_km,
    delta_ix: roe.delta_ix * a_km,
    delta_iy: roe.delta_iy * a_km,
  };
}

function roeResultRows(roeKm) {
  return [
    { label: "a·δa", value: formatKm(roeKm.delta_a), unit: "km" },
    { label: "a·δλ", value: formatKm(roeKm.delta_lambda), unit: "km" },
    { label: "a·δe_x", value: formatKm(roeKm.delta_ex), unit: "km" },
    { label: "a·δe_y", value: formatKm(roeKm.delta_ey), unit: "km" },
    { label: "a·δi_x", value: formatKm(roeKm.delta_ix), unit: "km" },
    { label: "a·δi_y", value: formatKm(roeKm.delta_iy), unit: "km" },
  ];
}

function rmeeResultRows(rmee) {
  return [
    { label: "δ<var>p</var>/<var>p</var><sub>c</sub>", value: formatNum(rmee.dp_over_pc, 8) },
    { label: "δ<var>L</var>", value: formatNum(radToDeg(rmee.dL_rad), 4), unit: "°" },
    { label: "δ<var>f</var>", value: formatNum(rmee.df) },
    { label: "δ<var>g</var>", value: formatNum(rmee.dg) },
    { label: "δ<var>h</var>", value: formatNum(rmee.dh) },
    { label: "δ<var>k</var>", value: formatNum(rmee.dk) },
  ];
}

function resultListHtml(rows) {
  return rows
    .map(({ label, value, unit }) => {
      const unitHtml = unit
        ? `<span class="result-unit">${unit}</span>`
        : "";
      return (
        `<li>` +
        `<span class="result-label">${label}</span>` +
        `<output class="result-value">${value}</output>${unitHtml}` +
        `</li>`
      );
    })
    .join("");
}

function renderResultList(listEl, rows) {
  if (!listEl) return;
  listEl.innerHTML = resultListHtml(rows);
}

function renderResultStack(containerEl, blocks) {
  if (!containerEl) return;
  containerEl.innerHTML = blocks
    .map(
      (block) =>
        `<div class="result-block">` +
        `<p class="result-block__title">${block.title}</p>` +
        `<ul class="result-list">${resultListHtml(block.rows)}</ul>` +
        `</div>`
    )
    .join("");
}

function trueAnomalyFromU(u_rad, argp_rad) {
  return u_rad - argp_rad;
}
function radiusFromTrueAnomaly(a_km, e, nu) {
  return (a_km * (1 - e * e)) / (1 + e * Math.cos(nu));
}
function positionEciKm(coe, u_rad) {
  const { a_km, e, i_rad, raan_rad, argp_rad } = coe;
  const nu = trueAnomalyFromU(u_rad, argp_rad);
  const r = radiusFromTrueAnomaly(a_km, e, nu);
  const cu = Math.cos(u_rad);
  const su = Math.sin(u_rad);
  const cO = Math.cos(raan_rad);
  const sO = Math.sin(raan_rad);
  const ci = Math.cos(i_rad);
  const si = Math.sin(i_rad);
  return [
    r * (cO * cu - sO * su * ci),
    r * (sO * cu + cO * su * ci),
    r * su * si,
  ];
}
function initialUFromCoe(coe) {
  const nu = trueAnomalyFromMeanAnomaly(coe.M_rad, coe.e);
  return normalizeAngleRad(coe.argp_rad + nu);
}

function buildOrbitEci(coe, numPoints) {
  const u0 = initialUFromCoe(coe);
  const x_km = [];
  const y_km = [];
  const z_km = [];
  for (let i = 0; i < numPoints; i++) {
    const u = u0 + (TWO_PI * i) / (numPoints - 1);
    const [x, y, z] = positionEciKm(coe, u);
    x_km.push(x);
    y_km.push(y);
    z_km.push(z);
  }
  return { x_km, y_km, z_km };
}

function chiefMeanMotionRadS(a_km) {
  return Math.sqrt(CONFIG.MU_EARTH_KM3_S2 / (a_km * a_km * a_km));
}

function deltaLambdaAtTimeSec(roe, n_rad_s, t_sec) {
  return roe.delta_lambda - DRIFT_LAMBDA_FACTOR * n_rad_s * roe.delta_a * t_sec;
}

function rtnAtLatitude(a_km, roe, u_rad, deltaLambdaOverride = null) {
  const delta_lambda =
    deltaLambdaOverride !== null ? deltaLambdaOverride : roe.delta_lambda;
  const { delta_a, delta_ex, delta_ey, delta_ix, delta_iy } = roe;
  const cu = Math.cos(u_rad);
  const su = Math.sin(u_rad);
  const R_km = a_km * (delta_a - delta_ex * cu - delta_ey * su);
  const T_km = a_km * (delta_lambda + 2 * delta_ex * su - 2 * delta_ey * cu);
  const N_km = a_km * (delta_ix * su - delta_iy * cu);
  return { R_km, T_km, N_km };
}

function buildSingleOrbitSamples(a_km, roe, numPoints) {
  const samples = [];
  for (let i = 0; i < numPoints; i++) {
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    samples.push({ u_rad, ...rtnAtLatitude(a_km, roe, u_rad) });
  }
  return samples;
}

function buildDriftOrbitSamples(a_km, roe, numPointsPerOrbit, numOrbits) {
  const n_rad_s = chiefMeanMotionRadS(a_km);
  const orbitPeriodSec = TWO_PI / n_rad_s;
  const samples = [];
  for (let k = 0; k < numOrbits; k++) {
    for (let i = 0; i < numPointsPerOrbit; i++) {
      if (k > 0 && i === 0) continue;
      const u_rad = (TWO_PI * i) / (numPointsPerOrbit - 1);
      const t_sec = k * orbitPeriodSec + u_rad / n_rad_s;
      const delta_lambda_t = deltaLambdaAtTimeSec(roe, n_rad_s, t_sec);
      samples.push({ u_rad, t_sec, ...rtnAtLatitude(a_km, roe, u_rad, delta_lambda_t) });
    }
  }
  return { samples, n_rad_s, orbitPeriodSec };
}

function samplesToRtnOrbit(samples) {
  return {
    R_km: samples.map((s) => s.R_km),
    T_km: samples.map((s) => s.T_km),
    N_km: samples.map((s) => s.N_km),
  };
}

function chiefVelocityEciKmS(coe, nu_rad, r_km) {
  const { a_km, e, i_rad, raan_rad, argp_rad } = coe;
  const mu = CONFIG.MU_EARTH_KM3_S2;
  const h = Math.sqrt(mu * a_km * (1 - e * e));
  const p = a_km * (1 - e * e);
  const vr = (mu / h) * e * Math.sin(nu_rad);
  const vt = (mu / h) * (p / r_km);
  const cnu = Math.cos(nu_rad);
  const snu = Math.sin(nu_rad);
  const v_pf_x = vr * cnu - vt * snu;
  const v_pf_y = vr * snu + vt * cnu;
  const cO = Math.cos(raan_rad);
  const sO = Math.sin(raan_rad);
  const ci = Math.cos(i_rad);
  const si = Math.sin(i_rad);
  const cw = Math.cos(argp_rad);
  const sw = Math.sin(argp_rad);
  const R11 = cO * cw - sO * sw * ci;
  const R12 = -cO * sw - sO * cw * ci;
  const R21 = sO * cw + cO * sw * ci;
  const R22 = -sO * sw + cO * cw * ci;
  const R31 = sw * si;
  const R32 = cw * si;
  return [
    R11 * v_pf_x + R12 * v_pf_y,
    R21 * v_pf_x + R22 * v_pf_y,
    R31 * v_pf_x + R32 * v_pf_y,
  ];
}

function normalizeVec3(v) {
  const n = Math.hypot(v[0], v[1], v[2]);
  if (n < 1e-12) return [1, 0, 0];
  return [v[0] / n, v[1] / n, v[2] / n];
}

function crossVec3(a, b) {
  return [
    a[1] * b[2] - a[2] * b[1],
    a[2] * b[0] - a[0] * b[2],
    a[0] * b[1] - a[1] * b[0],
  ];
}

function rtnBasisEci(coe, u_rad) {
  const { argp_rad } = coe;
  const nu = trueAnomalyFromU(u_rad, argp_rad);
  const r_km = radiusFromTrueAnomaly(coe.a_km, coe.e, nu);
  const r_vec = positionEciKm(coe, u_rad);
  const v_vec = chiefVelocityEciKmS(coe, nu, r_km);
  const R_hat = normalizeVec3(r_vec);
  const N_hat = normalizeVec3(crossVec3(r_vec, v_vec));
  const T_hat = normalizeVec3(crossVec3(N_hat, R_hat));
  return { R_hat, T_hat, N_hat };
}

function relativeRtnToEciKm(R_km, T_km, N_km, coe, u_rad) {
  const { R_hat, T_hat, N_hat } = rtnBasisEci(coe, u_rad);
  return [
    R_km * R_hat[0] + T_km * T_hat[0] + N_km * N_hat[0],
    R_km * R_hat[1] + T_km * T_hat[1] + N_km * N_hat[1],
    R_km * R_hat[2] + T_km * T_hat[2] + N_km * N_hat[2],
  ];
}

function buildChiefOrbitEci(coe, numPoints) {
  const x_km = [];
  const y_km = [];
  const z_km = [];
  for (let i = 0; i < numPoints; i++) {
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    const [x, y, z] = positionEciKm(coe, u_rad);
    x_km.push(x);
    y_km.push(y);
    z_km.push(z);
  }
  return { x_km, y_km, z_km };
}

function buildEarthEquatorCrossSectionEci() {
  const r = CONFIG.R_EARTH_KM;
  const n = CONFIG.EARTH_EQUATOR_POINTS;
  const x_km = [];
  const y_km = [];
  const z_km = [];
  for (let i = 0; i <= n; i++) {
    const t = (TWO_PI * i) / n;
    x_km.push(r * Math.cos(t));
    y_km.push(r * Math.sin(t));
    z_km.push(0);
  }
  return { x_km, y_km, z_km };
}
function safeAxisRange(min, max, fallbackSpan = 10) {
  if (!Number.isFinite(min) || !Number.isFinite(max)) {
    return [-fallbackSpan / 2, fallbackSpan / 2];
  }
  if (min === max) {
    const half = Math.max(Math.abs(min) * 0.01, fallbackSpan / 2);
    return [min - half, max + half];
  }
  return [min, max];
}

function samplesToEciExaggeratedOrbit(coe, samples, exaggerationK) {
  const x_km = [];
  const y_km = [];
  const z_km = [];
  for (const sample of samples) {
    const chief = positionEciKm(coe, sample.u_rad);
    const delta = relativeRtnToEciKm(
      sample.R_km,
      sample.T_km,
      sample.N_km,
      coe,
      sample.u_rad
    );
    x_km.push(chief[0] + exaggerationK * delta[0]);
    y_km.push(chief[1] + exaggerationK * delta[1]);
    z_km.push(chief[2] + exaggerationK * delta[2]);
  }
  return { x_km, y_km, z_km };
}

function readEciExaggerationFactor() {
  const el = document.getElementById("eci-exaggeration");
  const logExp = el ? parseInt(el.value, 10) : CONFIG.DEFAULT_ECI_EXAGGERATION_LOG;
  return Math.pow(10, Number.isFinite(logExp) ? logExp : CONFIG.DEFAULT_ECI_EXAGGERATION_LOG);
}

function buildEciExaggeratedOrbitData(coe, primarySamples, referenceSamples) {
  const k = readEciExaggerationFactor();
  return {
    exaggerationK: k,
    primaryEci: samplesToEciExaggeratedOrbit(coe, primarySamples, k),
    referenceEci: referenceSamples
      ? samplesToEciExaggeratedOrbit(coe, referenceSamples, k)
      : null,
    chiefOrbitEci: buildChiefOrbitEci(coe, CONFIG.ECI_CHIEF_ORBIT_POINTS),
  };
}

function computeEciAxisRanges(orbits) {
  const coords = { x: [], y: [], z: [] };
  for (const orbit of orbits) {
    if (!orbit) continue;
    coords.x.push(...orbit.x_km);
    coords.y.push(...orbit.y_km);
    coords.z.push(...orbit.z_km);
  }
  const padRatio = CONFIG.ECI_ZOOM_PADDING_RATIO;
  const axisRange = (arr) => {
    if (arr.length === 0) return safeAxisRange(-1, 1);
    const min = Math.min(...arr);
    const max = Math.max(...arr);
    const span = Math.max(max - min, 0.5);
    const pad = span * padRatio;
    return safeAxisRange(min - pad, max + pad);
  };
  return { x: axisRange(coords.x), y: axisRange(coords.y), z: axisRange(coords.z) };
}

/** 実スケール 3D 用：軌道中心の立方体レンジ（Z が潰れないよう等尺） */
function computeEciAxisRangesCube(orbits, padRatio = CONFIG.AXIS_PADDING) {
  const xs = [];
  const ys = [];
  const zs = [];
  for (const orbit of orbits) {
    if (!orbit) continue;
    xs.push(...orbit.x_km);
    ys.push(...orbit.y_km);
    zs.push(...orbit.z_km);
  }
  if (xs.length === 0) {
    return { x: safeAxisRange(-1, 1), y: safeAxisRange(-1, 1), z: safeAxisRange(-1, 1) };
  }
  const minX = Math.min(...xs);
  const maxX = Math.max(...xs);
  const minY = Math.min(...ys);
  const maxY = Math.max(...ys);
  const minZ = Math.min(...zs);
  const maxZ = Math.max(...zs);
  const cx = (minX + maxX) / 2;
  const cy = (minY + maxY) / 2;
  const cz = (minZ + maxZ) / 2;
  const maxSpan = Math.max(
    maxX - minX,
    maxY - minY,
    maxZ - minZ,
    CONFIG.R_EARTH_KM * 2,
    1
  );
  const half = (maxSpan * (1 + padRatio)) / 2;
  return {
    x: [cx - half, cx + half],
    y: [cy - half, cy + half],
    z: [cz - half, cz + half],
  };
}

function computeEciAxisRangesForPlane(orbits, axisKey) {
  const values = [];
  for (const orbit of orbits) {
    if (!orbit) continue;
    values.push(...orbit[axisKey]);
  }
  if (values.length === 0) return safeAxisRange(-1, 1);
  const min = Math.min(...values);
  const max = Math.max(...values);
  const span = Math.max(max - min, 0.5);
  const pad = span * CONFIG.ECI_ZOOM_PADDING_RATIO;
  return safeAxisRange(min - pad, max + pad);
}

function readNumDriftOrbits() {
  const el = document.getElementById("num-drift-orbits");
  const v = parseInt(el?.value, 10);
  return Number.isFinite(v) && v >= 1 ? v : CONFIG.DEFAULT_NUM_DRIFT_ORBITS;
}

function computeOrbitsForDisplay(a_km, roe) {
  const numPoints = CONFIG.NUM_ORBIT_POINTS;
  const numOrbits = readNumDriftOrbits();
  const { samples: primarySamples } = buildDriftOrbitSamples(
    a_km,
    roe,
    numPoints,
    numOrbits
  );
  const showRef = document.getElementById("show-single-orbit-ref")?.checked ?? true;
  const referenceSamples = showRef
    ? buildSingleOrbitSamples(a_km, roe, numPoints)
    : null;
  return {
    primary: samplesToRtnOrbit(primarySamples),
    reference: referenceSamples ? samplesToRtnOrbit(referenceSamples) : null,
    primarySamples,
    referenceSamples,
  };
}

function updateDriftRateDisplay(a_km, roeKm) {
  const el = document.getElementById("drift-rate-display");
  if (!el) return;
  const deltaAKm = roeKm.delta_a;
  if (Math.abs(deltaAKm) < 1e-9) {
    el.textContent =
      "Δa ≈ 0 km → along-track ドリフトなし。\n" +
      "複数周表示は同じ閉曲線に重なります（周回数を 1 にすると見やすいです）。";
    el.classList.add("drift-rate--warn");
    return;
  }
  el.classList.remove("drift-rate--warn");
  const roeDim = {
    delta_a: roeKm.delta_a / a_km,
    delta_lambda: roeKm.delta_lambda / a_km,
    delta_ex: roeKm.delta_ex / a_km,
    delta_ey: roeKm.delta_ey / a_km,
    delta_ix: roeKm.delta_ix / a_km,
    delta_iy: roeKm.delta_iy / a_km,
  };
  const n = chiefMeanMotionRadS(a_km);
  const dDeltaLambdaDt = -DRIFT_LAMBDA_FACTOR * n * roeDim.delta_a;
  const driftPerOrbitKm = -DRIFT_LAMBDA_FACTOR * TWO_PI * deltaAKm;
  const periodHr = TWO_PI / n / 3600;
  el.textContent =
    `δλ̇ ≈ ${formatKm(dDeltaLambdaDt * a_km)}/s（along-track スケール）\n` +
    `1 周あたり ΔT ≈ ${formatKm(driftPerOrbitKm)} km\n` +
    `T_orbit ≈ ${periodHr.toFixed(2)} h`;
}

const ORBIT_ENDPOINT_MARKER = {
  start: { color: "#58a6ff", size: 9, symbol: "circle", line: { color: "#e6edf3", width: 1 } },
  end: { color: "#f85149", size: 9, symbol: "square", line: { color: "#e6edf3", width: 1 } },
};
const RTN_AXIS_LABELS_3D = { x: "T", y: "R", z: "N" };
const ECI_AXIS_LABELS_3D = { x: "X", y: "Y", z: "Z" };

function pushOrbitStartMarker3D(traces, coords, labelPrefix, axisLabels) {
  const { x, y, z } = coords;
  if (!x || x.length < 1) return;
  const axes = axisLabels || ECI_AXIS_LABELS_3D;
  traces.push({
    type: "scatter3d",
    mode: "markers",
    name: `${labelPrefix} 始点`,
    x: [x[0]],
    y: [y[0]],
    z: [z[0]],
    marker: { ...ORBIT_ENDPOINT_MARKER.start },
    hovertemplate:
      `${labelPrefix} 始点 (M から)<br>` +
      `${axes.x}: %{x:.3f} ${axes.y}: %{y:.3f} ${axes.z}: %{z:.3f} km<extra></extra>`,
  });
}

function pushOrbitEndpointMarkers3D(traces, coords, labelPrefix, axisLabels) {
  const { x, y, z } = coords;
  if (!x || x.length < 1) return;
  const axes = axisLabels || ECI_AXIS_LABELS_3D;
  const last = x.length - 1;
  traces.push({
    type: "scatter3d",
    mode: "markers",
    name: `${labelPrefix} 始点`,
    x: [x[0]],
    y: [y[0]],
    z: [z[0]],
    marker: { ...ORBIT_ENDPOINT_MARKER.start },
    hovertemplate:
      `${labelPrefix} 始点 (u ≈ 0)<br>` +
      `${axes.x}: %{x:.3f} ${axes.y}: %{y:.3f} ${axes.z}: %{z:.3f} km<extra></extra>`,
  });
  traces.push({
    type: "scatter3d",
    mode: "markers",
    name: `${labelPrefix} 終点`,
    x: [x[last]],
    y: [y[last]],
    z: [z[last]],
    marker: { ...ORBIT_ENDPOINT_MARKER.end },
    hovertemplate:
      `${labelPrefix} 終点<br>` +
      `${axes.x}: %{x:.3f} ${axes.y}: %{y:.3f} ${axes.z}: %{z:.3f} km<extra></extra>`,
  });
}

function pushOrbitEndpointMarkers2D(traces, x, y, labelPrefix) {
  if (!x || x.length < 1) return;
  const last = x.length - 1;
  traces.push({
    type: "scatter",
    mode: "markers",
    name: `${labelPrefix} 始点`,
    x: [x[0]],
    y: [y[0]],
    marker: { ...ORBIT_ENDPOINT_MARKER.start },
    hovertemplate: `${labelPrefix} 始点: %{x:.3f}, %{y:.3f} km<extra></extra>`,
  });
  traces.push({
    type: "scatter",
    mode: "markers",
    name: `${labelPrefix} 終点`,
    x: [x[last]],
    y: [y[last]],
    marker: { ...ORBIT_ENDPOINT_MARKER.end },
    hovertemplate: `${labelPrefix} 終点: %{x:.3f}, %{y:.3f} km<extra></extra>`,
  });
}

function buildChiefTrace3D() {
  return {
    type: "scatter3d",
    mode: "markers",
    name: "Chief（原点）",
    x: [0],
    y: [0],
    z: [0],
    marker: { color: "#f0c14b", size: 6, symbol: "diamond" },
    hovertemplate: "Chief @ 原点<extra></extra>",
  };
}

function buildDeputyTrace3D(orbit, name, color, width) {
  return {
    type: "scatter3d",
    mode: "lines",
    name,
    x: orbit.T_km,
    y: orbit.R_km,
    z: orbit.N_km,
    line: { color, width },
    hovertemplate: "T: %{x:.3f} km<br>R: %{y:.3f} km<br>N: %{z:.3f} km<extra></extra>",
  };
}

function buildOrbitTraces3D(primaryOrbit, referenceOrbit) {
  const traces = [];
  if (referenceOrbit) {
    traces.push(
      buildDeputyTrace3D(referenceOrbit, "参考: 1 周（δλ 固定）", "rgba(139, 148, 158, 0.55)", 2)
    );
  }
  traces.push(buildDeputyTrace3D(primaryOrbit, "Deputy（相対軌道）", "#3fb950", 4));
  pushOrbitEndpointMarkers3D(
    traces,
    { x: primaryOrbit.T_km, y: primaryOrbit.R_km, z: primaryOrbit.N_km },
    "Deputy",
    RTN_AXIS_LABELS_3D
  );
  traces.push(buildChiefTrace3D());
  return traces;
}

function layout3D() {
  return {
    ...PLOT_LAYOUT_BASE,
    height: 420,
    showlegend: true,
    legend: { x: 0, y: 1, bgcolor: "rgba(0,0,0,0.3)" },
    scene: {
      xaxis: { title: "T [km]", gridcolor: "#30363d", zerolinecolor: "#484f58" },
      yaxis: { title: "R [km]", gridcolor: "#30363d", zerolinecolor: "#484f58" },
      zaxis: { title: "N [km]", gridcolor: "#30363d", zerolinecolor: "#484f58" },
      bgcolor: "#1c2333",
      aspectmode: "data",
      camera: { eye: { x: 1.4, y: 1.2, z: 0.9 } },
    },
  };
}

function build2DTrace(x, y, name, color, width = 2, dash = null) {
  const line = { color, width };
  if (dash) line.dash = dash;
  return {
    type: "scatter",
    mode: "lines",
    name,
    x,
    y,
    line,
    hovertemplate: "%{fullData.name}: %{x:.3f}, %{y:.3f} km<extra></extra>",
  };
}

function buildChief2D() {
  return {
    type: "scatter",
    mode: "markers",
    name: "Chief",
    x: [0],
    y: [0],
    marker: { color: "#f0c14b", size: 8, symbol: "diamond" },
    showlegend: false,
    hovertemplate: "Chief @ 原点<extra></extra>",
  };
}

const PLANE_TR = {
  axisX: "T_km",
  axisY: "R_km",
  xTitle: "T [km]",
  yTitle: "R [km]",
  color: "#3fb950",
  id: "tr",
};
const PLANE_TN = {
  axisX: "T_km",
  axisY: "N_km",
  xTitle: "T [km]",
  yTitle: "N [km]",
  color: "#58a6ff",
  id: "tn",
};
const PLANE_RN = {
  axisX: "R_km",
  axisY: "N_km",
  xTitle: "R [km]",
  yTitle: "N [km]",
  color: "#a371f7",
  id: "rn",
};

function layout2D(xTitle, yTitle, height, planeId) {
  const useEqualScale = !(planeId === "tr" || planeId === "tn");
  const xaxis = {
    title: xTitle,
    gridcolor: "#30363d",
    zerolinecolor: "#484f58",
  };
  if (useEqualScale) {
    xaxis.scaleanchor = "y";
    xaxis.scaleratio = 1;
  }
  return {
    ...PLOT_LAYOUT_BASE,
    height,
    xaxis,
    yaxis: { title: yTitle, gridcolor: "#30363d", zerolinecolor: "#484f58" },
    showlegend: true,
    legend: { font: { size: 10 } },
  };
}

function build2DTraces(orbit, plane, referenceOrbit) {
  const traces = [];
  const { axisX, axisY, xTitle, yTitle, color, id } = plane;
  if (referenceOrbit) {
    traces.push(
      build2DTrace(
        referenceOrbit[axisX],
        referenceOrbit[axisY],
        "参考: 1 周",
        "rgba(139, 148, 158, 0.7)",
        1.5,
        "dot"
      )
    );
  }
  traces.push(build2DTrace(orbit[axisX], orbit[axisY], "相対軌道", color, 2));
  pushOrbitEndpointMarkers2D(traces, orbit[axisX], orbit[axisY], "Deputy");
  traces.push(buildChief2D());
  return { traces, layout: layout2D(xTitle, yTitle, 260, id) };
}

function buildEciLineTrace3D(orbit, name, color, width, dash) {
  const line = { color, width };
  if (dash) line.dash = dash;
  return {
    type: "scatter3d",
    mode: "lines",
    name,
    x: orbit.x_km,
    y: orbit.y_km,
    z: orbit.z_km,
    line,
    hovertemplate: "X: %{x:.3f}<br>Y: %{y:.3f}<br>Z: %{z:.3f} km<extra></extra>",
  };
}

function buildEarthEquatorTrace3D() {
  const earth = buildEarthEquatorCrossSectionEci();
  return buildEciLineTrace3D(
    earth,
    `地球（赤道断面, R⊕=${CONFIG.R_EARTH_KM.toFixed(1)} km）`,
    "rgba(110, 159, 239, 0.9)",
    2,
    "dot"
  );
}

function buildEciTraces3D(primaryEci, referenceEci, chiefOrbitEci, exaggerationK) {
  const earthTrace = buildEarthEquatorTrace3D();
  earthTrace.hovertemplate =
    `R⊕ = ${CONFIG.R_EARTH_KM.toFixed(1)} km（赤道断面）<extra></extra>`;
  const traces = [
    earthTrace,
    buildEciLineTrace3D(chiefOrbitEci, "Chief 軌道", "rgba(240, 193, 75, 0.85)", 2, null),
  ];
  if (referenceEci) {
    traces.push(
      buildEciLineTrace3D(
        referenceEci,
        "参考: 1 周（δλ 固定）",
        "rgba(139, 148, 158, 0.55)",
        2,
        "dot"
      )
    );
  }
  traces.push(
    buildEciLineTrace3D(primaryEci, `Deputy（k=${exaggerationK}）`, "#3fb950", 4, null)
  );
  pushOrbitEndpointMarkers3D(
    traces,
    { x: primaryEci.x_km, y: primaryEci.y_km, z: primaryEci.z_km },
    "Deputy",
    ECI_AXIS_LABELS_3D
  );
  return traces;
}

function layoutEci3D(axisRanges, { aspectCube = false } = {}) {
  const axisStyle = { gridcolor: "#30363d", zerolinecolor: "#484f58" };
  return {
    ...PLOT_LAYOUT_BASE,
    height: 420,
    showlegend: true,
    legend: { x: 0, y: 1, bgcolor: "rgba(0,0,0,0.3)" },
    scene: {
      xaxis: { title: "X [km]", range: axisRanges.x, ...axisStyle },
      yaxis: { title: "Y [km]", range: axisRanges.y, ...axisStyle },
      zaxis: { title: "Z [km]", range: axisRanges.z, ...axisStyle },
      bgcolor: "#1c2333",
      aspectmode: aspectCube ? "cube" : "data",
      camera: { eye: { x: 1.35, y: 1.25, z: 0.85 }, center: { x: 0, y: 0, z: 0 } },
    },
  };
}

function buildEci2DTrace(orbit, axisX, axisY, name, color, width, dash) {
  const line = { color, width };
  if (dash) line.dash = dash;
  return {
    type: "scatter",
    mode: "lines",
    name,
    x: orbit[axisX],
    y: orbit[axisY],
    line,
    hovertemplate: "%{x:.3f}, %{y:.3f} km<extra></extra>",
  };
}

const PLANE_ECI_XY = {
  axisX: "x_km",
  axisY: "y_km",
  xTitle: "X [km]",
  yTitle: "Y [km]",
  color: "#3fb950",
};
const PLANE_ECI_XZ = {
  axisX: "x_km",
  axisY: "z_km",
  xTitle: "X [km]",
  yTitle: "Z [km]",
  color: "#58a6ff",
};
const PLANE_ECI_YZ = {
  axisX: "y_km",
  axisY: "z_km",
  xTitle: "Y [km]",
  yTitle: "Z [km]",
  color: "#a371f7",
};

function buildEci2DTraces(primaryEci, referenceEci, chiefOrbitEci, plane, exaggerationK) {
  const { axisX, axisY, xTitle, yTitle } = plane;
  const earthEquator = buildEarthEquatorCrossSectionEci();
  const traces = [
    buildEci2DTrace(earthEquator, axisX, axisY, "地球（R⊕ 断面）", "rgba(110, 159, 239, 0.9)", 1.5, "dot"),
    buildEci2DTrace(chiefOrbitEci, axisX, axisY, "Chief 軌道", "rgba(240, 193, 75, 0.85)", 1.5, "dot"),
  ];
  if (referenceEci) {
    traces.push(
      buildEci2DTrace(
        referenceEci,
        axisX,
        axisY,
        "参考: 1 周",
        "rgba(139, 148, 158, 0.7)",
        1.5,
        "dot"
      )
    );
  }
  traces.push(
    buildEci2DTrace(primaryEci, axisX, axisY, `Deputy（k=${exaggerationK}）`, plane.color, 2, null)
  );
  pushOrbitEndpointMarkers2D(traces, primaryEci[axisX], primaryEci[axisY], "Deputy");
  const allOrbits = [earthEquator, chiefOrbitEci, primaryEci, referenceEci];
  const axisRanges = {
    x: computeEciAxisRangesForPlane(allOrbits, axisX),
    y: computeEciAxisRangesForPlane(allOrbits, axisY),
  };
  return {
    traces,
    layout: {
      ...PLOT_LAYOUT_BASE,
      height: 260,
      xaxis: {
        title: xTitle,
        range: axisRanges.x,
        gridcolor: "#30363d",
        zerolinecolor: "#484f58",
      },
      yaxis: {
        title: yTitle,
        range: axisRanges.y,
        gridcolor: "#30363d",
        zerolinecolor: "#484f58",
      },
      showlegend: true,
      legend: { font: { size: 10 } },
    },
  };
}

function updateEciCaption(exaggerationK) {
  const el = document.getElementById("plot-eci-caption");
  if (!el) return;
  if (exaggerationK <= 1) {
    el.textContent =
      "実スケール（k=1）。deputy は chief 軌道にほぼ重なって見えます。";
    return;
  }
  el.textContent =
    `表示専用: r = r_chief + k·(r_deputy−r_chief)、k=${exaggerationK}。` +
    "実際の相対距離は km オーダ（RTN プロット参照）。";
}

function plotDivExists(id) {
  return Boolean(document.getElementById(id));
}

function drawPlot(divId, traces, layout, useNewPlot, plotOpts) {
  if (!plotDivExists(divId)) return Promise.resolve();
  if (useNewPlot) return Plotly.newPlot(divId, traces, layout, plotOpts);
  return Plotly.react(divId, traces, layout);
}

function getPlotIdsForView(view = activeView) {
  return VIEW_PLOT_IDS[view] ?? RTN_PLOT_IDS;
}

function schedulePlotResize(view = activeView) {
  requestAnimationFrame(() => {
    for (const id of getPlotIdsForView(view)) {
      if (plotDivExists(id)) Plotly.Plots.resize(id);
    }
  });
}

function getChiefPlotKey(coe) {
  if (!coe) return "";
  return [
    coe.a_km.toFixed(3),
    coe.e.toFixed(6),
    radToDeg(coe.i_rad).toFixed(3),
    radToDeg(coe.raan_rad).toFixed(2),
    radToDeg(coe.argp_rad).toFixed(2),
    radToDeg(coe.M_rad).toFixed(3),
  ].join("|");
}

function getPlotStructureKey(chiefPlotKey = "") {
  const showRef = document.getElementById("show-single-orbit-ref")?.checked ?? true;
  return [
    `mode=${activeMode}`,
    `ref=${showRef ? 1 : 0}`,
    `orb=${readNumDriftOrbits()}`,
    `k=${activeMode === "coe-mee" ? 1 : readEciExaggerationFactor()}`,
    chiefPlotKey,
  ].join("|");
}

function buildPlotBundleFromCoePair(chiefCoe, deputyCoe) {
  const roe = coeToRoe(chiefCoe, deputyCoe);
  const a_km = chiefCoe.a_km;
  const { primary, reference, primarySamples, referenceSamples } =
    computeOrbitsForDisplay(a_km, roe);
  const {
    primaryEci,
    referenceEci,
    chiefOrbitEci,
    exaggerationK,
  } = buildEciExaggeratedOrbitData(chiefCoe, primarySamples, referenceSamples);
  const earthEquator = buildEarthEquatorCrossSectionEci();
  const axisRanges3d = computeEciAxisRanges([
    earthEquator,
    chiefOrbitEci,
    primaryEci,
    referenceEci,
  ]);
  return {
    traces3d: buildOrbitTraces3D(primary, reference),
    layout3d: layout3D(),
    tr: build2DTraces(primary, PLANE_TR, reference),
    tn: build2DTraces(primary, PLANE_TN, reference),
    rn: build2DTraces(primary, PLANE_RN, reference),
    tracesEci3d: buildEciTraces3D(primaryEci, referenceEci, chiefOrbitEci, exaggerationK),
    layoutEci3d: layoutEci3D(axisRanges3d),
    eciXy: buildEci2DTraces(primaryEci, referenceEci, chiefOrbitEci, PLANE_ECI_XY, exaggerationK),
    eciXz: buildEci2DTraces(primaryEci, referenceEci, chiefOrbitEci, PLANE_ECI_XZ, exaggerationK),
    eciYz: buildEci2DTraces(primaryEci, referenceEci, chiefOrbitEci, PLANE_ECI_YZ, exaggerationK),
    roeKm: roeKmFromDimensionless(roe, a_km),
  };
}

function drawPlots(plotBundle, { eciOnly = false, chiefPlotKey = "", forceNewPlot = false } = {}) {
  const plotOpts = { responsive: true, displayModeBar: true };
  const structureKey = getPlotStructureKey(chiefPlotKey);
  const useNewPlot =
    forceNewPlot || !plotsInitialized || structureKey !== lastPlotStructureKey;
  lastPlotStructureKey = structureKey;

  if (useNewPlot && plotsInitialized) {
    const idsToPurge = eciOnly ? PLOT_DIV_IDS : PLOT_DIV_IDS;
    for (const id of idsToPurge) {
      if (plotDivExists(id)) Plotly.purge(id);
    }
  }

  const { traces3d, layout3d, tr, tn, rn, tracesEci3d, layoutEci3d, eciXy, eciXz, eciYz } =
    plotBundle;

  const drawTasks = eciOnly
    ? [
        drawPlot("plot-eci-3d", tracesEci3d, layoutEci3d, useNewPlot, plotOpts),
        drawPlot("plot-eci-xy", eciXy.traces, eciXy.layout, useNewPlot, plotOpts),
        drawPlot("plot-eci-xz", eciXz.traces, eciXz.layout, useNewPlot, plotOpts),
        drawPlot("plot-eci-yz", eciYz.traces, eciYz.layout, useNewPlot, plotOpts),
      ]
    : [
        drawPlot("plot-3d", traces3d, layout3d, useNewPlot, plotOpts),
        drawPlot("plot-tr", tr.traces, tr.layout, useNewPlot, plotOpts),
        drawPlot("plot-tn", tn.traces, tn.layout, useNewPlot, plotOpts),
        drawPlot("plot-rn", rn.traces, rn.layout, useNewPlot, plotOpts),
        drawPlot("plot-eci-3d", tracesEci3d, layoutEci3d, useNewPlot, plotOpts),
        drawPlot("plot-eci-xy", eciXy.traces, eciXy.layout, useNewPlot, plotOpts),
        drawPlot("plot-eci-xz", eciXz.traces, eciXz.layout, useNewPlot, plotOpts),
        drawPlot("plot-eci-yz", eciYz.traces, eciYz.layout, useNewPlot, plotOpts),
      ];

  Promise.all(drawTasks)
    .then(() => {
      schedulePlotResize(eciOnly ? "eci" : activeView);
      if (eciOnly) {
        requestAnimationFrame(() => {
          requestAnimationFrame(() => {
            for (const id of ECI_PLOT_IDS) {
              if (plotDivExists(id)) Plotly.Plots.resize(id);
            }
          });
        });
      }
    })
    .catch((err) => console.error("Plotly draw failed:", err));
  plotsInitialized = true;
}

/** @deprecated use drawPlots */
function drawAllPlots(plotBundle, forceNewPlot = false) {
  drawPlots(plotBundle, { eciOnly: false, forceNewPlot });
}

function buildChiefOnlyEci2D(chiefOrbitEci, plane) {
  const { axisX, axisY, xTitle, yTitle } = plane;
  const earth = buildEarthEquatorCrossSectionEci();
  const traces = [
    buildEci2DTrace(earth, axisX, axisY, "地球（R⊕ 断面）", "rgba(110, 159, 239, 0.9)", 1.5, "dot"),
    buildEci2DTrace(chiefOrbitEci, axisX, axisY, "Chief 軌道", "rgba(240, 193, 75, 0.85)", 1.5, "dot"),
  ];
  const allOrbits = [earth, chiefOrbitEci];
  const axisRanges = {
    x: computeEciAxisRangesForPlane(allOrbits, axisX),
    y: computeEciAxisRangesForPlane(allOrbits, axisY),
  };
  return {
    traces,
    layout: {
      ...PLOT_LAYOUT_BASE,
      height: 260,
      xaxis: {
        title: xTitle,
        range: axisRanges.x,
        gridcolor: "#30363d",
        zerolinecolor: "#484f58",
      },
      yaxis: {
        title: yTitle,
        range: axisRanges.y,
        gridcolor: "#30363d",
        zerolinecolor: "#484f58",
      },
      showlegend: true,
      legend: { font: { size: 10 } },
    },
  };
}

function updatePlotsFromCoePair(chiefCoe, deputyCoe) {
  const cap3d = document.getElementById("plot-3d-caption");
  if (cap3d) {
    cap3d.textContent =
      "横軸: T、縦軸: R、奥行き: N [km]（chief は原点）。● 始点（u≈0）、■ 終点。";
  }
  const bundle = buildPlotBundleFromCoePair(chiefCoe, deputyCoe);
  updateDriftRateDisplay(chiefCoe.a_km, bundle.roeKm);
  updateEciCaption(readEciExaggerationFactor());
  drawAllPlots(bundle);
}

function plotChiefOnlyEci(chiefCoe) {
  const chiefOrbitEci = buildOrbitEci(chiefCoe, CONFIG.ECI_CHIEF_ORBIT_POINTS);
  const axisRanges3d = computeEciAxisRangesCube([chiefOrbitEci]);
  const cap = document.getElementById("plot-eci-caption");
  if (cap) {
    cap.textContent =
      "COE↔MEE タブ — Chief 絶対軌道のみ（実スケール k=1、Deputy / RTN なし）。";
  }
  const tracesEci3d = [
    buildEarthEquatorTrace3D(),
    buildEciLineTrace3D(chiefOrbitEci, "Chief 軌道", "rgba(240, 193, 75, 0.85)", 2, null),
  ];
  pushOrbitStartMarker3D(
    tracesEci3d,
    { x: chiefOrbitEci.x_km, y: chiefOrbitEci.y_km, z: chiefOrbitEci.z_km },
    "Chief",
    ECI_AXIS_LABELS_3D
  );
  drawPlots(
    {
      tracesEci3d,
      layoutEci3d: layoutEci3D(axisRanges3d, { aspectCube: true }),
      eciXy: buildChiefOnlyEci2D(chiefOrbitEci, PLANE_ECI_XY),
      eciXz: buildChiefOnlyEci2D(chiefOrbitEci, PLANE_ECI_XZ),
      eciYz: buildChiefOnlyEci2D(chiefOrbitEci, PLANE_ECI_YZ),
    },
    { eciOnly: true, chiefPlotKey: getChiefPlotKey(chiefCoe), forceNewPlot: !plotsInitialized }
  );
}

function syncOrbitDisplaySection() {
  const section = document.getElementById("orbit-display-section");
  if (section) section.hidden = activeMode === "coe-mee";
}

function syncPlotUiForMode() {
  const isCoeMee = activeMode === "coe-mee";
  const viewTabs = document.getElementById("view-tabs");
  const tabRtn = document.getElementById("view-tab-rtn");
  const exaggerationSection = document.getElementById("eci-exaggeration-section");
  const eciTitle = document.getElementById("plot-eci-title");
  const eciProjTitle = document.getElementById("plot-eci-proj-title");

  if (viewTabs) viewTabs.hidden = isCoeMee;
  if (tabRtn) tabRtn.hidden = isCoeMee;
  if (exaggerationSection) exaggerationSection.hidden = isCoeMee;

  if (eciTitle) {
    eciTitle.textContent = isCoeMee
      ? "3D — ECI（慣性系・実スケール）"
      : "3D — ECI（慣性系・誇張表示）";
  }
  if (eciProjTitle) {
    eciProjTitle.textContent = isCoeMee
      ? "ECI — 2D 投影"
      : "ECI — 2D 投影（誇張表示）";
  }

  if (isCoeMee) {
    setActiveView("eci");
  } else {
    setActiveView(relativePlotView);
    updateEciCaption(readEciExaggerationFactor());
  }
}

function setActiveView(view) {
  if (view !== "rtn" && view !== "eci") return;
  if (activeMode !== "coe-mee") {
    relativePlotView = view;
  }
  activeView = view;

  const tabRtn = document.getElementById("view-tab-rtn");
  const tabEci = document.getElementById("view-tab-eci");
  const panelRtn = document.getElementById("view-panel-rtn");
  const panelEci = document.getElementById("view-panel-eci");
  if (!tabRtn || !tabEci || !panelRtn || !panelEci) return;

  const isRtn = view === "rtn" && activeMode !== "coe-mee";
  const isEci = !isRtn;

  tabRtn.classList.toggle("is-active", isRtn);
  tabEci.classList.toggle("is-active", isEci);
  tabRtn.setAttribute("aria-selected", String(isRtn));
  tabEci.setAttribute("aria-selected", String(isEci));
  tabRtn.tabIndex = isRtn ? 0 : -1;
  tabEci.tabIndex = isEci ? 0 : -1;

  panelRtn.classList.toggle("is-active", isRtn);
  panelEci.classList.toggle("is-active", isEci);
  panelRtn.hidden = !isRtn;
  panelEci.hidden = !isEci;

  schedulePlotResize(isEci ? "eci" : "rtn");
}

function initViewTabs() {
  document.querySelectorAll(".view-tab[data-view]").forEach((tab) => {
    tab.addEventListener("click", () => {
      const view = tab.getAttribute("data-view");
      if (view && view !== activeView) setActiveView(view);
    });
  });
  setActiveView(activeView);
}

function initOrbitControls() {
  const orbitSlider = document.getElementById("num-drift-orbits");
  const orbitOutput = document.getElementById("val-num-drift-orbits");
  if (!orbitSlider || !orbitOutput) return;
  orbitSlider.value = String(CONFIG.DEFAULT_NUM_DRIFT_ORBITS);
  orbitOutput.textContent = String(CONFIG.DEFAULT_NUM_DRIFT_ORBITS);
  orbitSlider.addEventListener("input", () => {
    orbitOutput.textContent = orbitSlider.value;
    updateAll();
  });
  document.getElementById("show-single-orbit-ref")?.addEventListener("change", () => updateAll());
}

function initEciExaggerationControl() {
  const slider = document.getElementById("eci-exaggeration");
  const output = document.getElementById("val-eci-exaggeration");
  if (!slider || !output) return;
  slider.value = String(CONFIG.DEFAULT_ECI_EXAGGERATION_LOG);
  output.textContent = String(readEciExaggerationFactor());
  slider.addEventListener("input", () => {
    output.textContent = String(readEciExaggerationFactor());
    updateEciCaption(readEciExaggerationFactor());
    updateAll();
  });
}

function appendRoeResultBlock(blocks, chiefCoe, deputyCoe) {
  const roeKm = roeKmFromDimensionless(coeToRoe(chiefCoe, deputyCoe), chiefCoe.a_km);
  blocks.push({
    title: "ROE（Gim–Alfriend, a·δ）[km]",
    rows: roeResultRows(roeKm),
  });
}

function readNumInput(id, fallback) {
  const el = document.getElementById(id);
  const v = parseFloat(el?.value, 10);
  return Number.isFinite(v) ? v : fallback;
}

function readChiefCoeFromInputs() {
  const defaults = CONFIG.DEFAULT_CHIEF_COE;
  return coeKmDegToRad({
    a_km: readNumInput("chief-a", defaults.a_km),
    e: readNumInput("chief-e", defaults.e),
    i_deg: readNumInput("chief-i", defaults.i_deg),
    raan_deg: readNumInput("chief-raan", defaults.raan_deg),
    argp_deg: readNumInput("chief-argp", defaults.argp_deg),
    M_deg: readNumInput("chief-M", defaults.M_deg),
  });
}

function readRmeeChiefCoeFromInputs() {
  const defaults = CONFIG.DEFAULT_CHIEF_COE;
  return coeKmDegToRad({
    a_km: readNumInput("rmee-chief-a", defaults.a_km),
    e: readNumInput("rmee-chief-e", defaults.e),
    i_deg: readNumInput("rmee-chief-i", defaults.i_deg),
    raan_deg: readNumInput("rmee-chief-raan", defaults.raan_deg),
    argp_deg: readNumInput("rmee-chief-argp", defaults.argp_deg),
    M_deg: readNumInput("rmee-chief-M", defaults.M_deg),
  });
}

function readDeputyCoeFromInputs() {
  const defaults = CONFIG.DEFAULT_DEPUTY_COE;
  return coeKmDegToRad({
    a_km: readNumInput("deputy-a", defaults.a_km),
    e: readNumInput("deputy-e", defaults.e),
    i_deg: readNumInput("deputy-i", defaults.i_deg),
    raan_deg: readNumInput("deputy-raan", defaults.raan_deg),
    argp_deg: readNumInput("deputy-argp", defaults.argp_deg),
    M_deg: readNumInput("deputy-M", defaults.M_deg),
  });
}

function buildMeeInputs(container, prefix, defaultsKmDeg) {
  container.innerHTML = "";
  for (const spec of MEE_SPECS) {
    const id = `${prefix}-${spec.suffix}`;
    const label = document.createElement("label");
    label.className = "field-label";
    label.htmlFor = id;
    label.innerHTML = `<var>${spec.label}</var> [${spec.unit}]`;
    const input = document.createElement("input");
    input.type = "number";
    input.id = id;
    input.className = "number-input mee-field";
    input.dataset.prefix = prefix;
    input.dataset.key = spec.key;
    input.step = spec.step;
    input.value = String(defaultsKmDeg[spec.key]);
    container.appendChild(label);
    container.appendChild(input);
  }
}

function buildRmeeInputs(container, prefix, defaultsDisplay) {
  container.innerHTML = "";
  for (const spec of RMEE_SPECS) {
    const id = `${prefix}-${spec.suffix}`;
    const label = document.createElement("label");
    label.className = "field-label";
    label.htmlFor = id;
    label.innerHTML = `<var>${spec.label}</var> [${spec.unit}]`;
    const input = document.createElement("input");
    input.type = "number";
    input.id = id;
    input.className = "number-input rmee-field";
    input.dataset.prefix = prefix;
    input.dataset.key = spec.key;
    input.step = spec.step;
    input.value = String(defaultsDisplay[spec.key]);
    container.appendChild(label);
    container.appendChild(input);
  }
}

function readMeeGroup(prefix, defaultsKmDeg) {
  const out = { ...defaultsKmDeg };
  document.querySelectorAll(`.mee-field[data-prefix="${prefix}"]`).forEach((el) => {
    const v = parseFloat(el.value, 10);
    if (Number.isFinite(v)) out[el.dataset.key] = v;
  });
  return meeKmDegToRad(out);
}

function getDefaultChiefMeeKmDeg() {
  return meeToKmDeg(coeToMee(readChiefCoeFromInputs()));
}

function getDefaultDeputyMeeKmDeg() {
  return meeToKmDeg(coeToMee(coeKmDegToRad(CONFIG.DEFAULT_DEPUTY_COE)));
}

function getDefaultRmeeChiefMeeKmDeg() {
  return meeToKmDeg(coeToMee(readRmeeChiefCoeFromInputs()));
}

function getDefaultDeputyRmeeDisplay() {
  const chiefMee = coeToMee(coeKmDegToRad(CONFIG.DEFAULT_CHIEF_COE));
  const deputyMee = coeToMee(coeKmDegToRad(CONFIG.DEFAULT_DEPUTY_COE));
  const rmee = computeRmee(chiefMee, deputyMee);
  return {
    dp_over_pc: rmee.dp_over_pc,
    dL_deg: radToDeg(rmee.dL_rad),
    df: rmee.df,
    dg: rmee.dg,
    dh: rmee.dh,
    dk: rmee.dk,
  };
}

function readRmeeGroup(prefix, defaultsDisplay) {
  const out = { ...defaultsDisplay };
  document.querySelectorAll(`.rmee-field[data-prefix="${prefix}"]`).forEach((el) => {
    const v = parseFloat(el.value, 10);
    if (Number.isFinite(v)) out[el.dataset.key] = v;
  });
  return {
    dp_over_pc: out.dp_over_pc,
    df: out.df,
    dg: out.dg,
    dh: out.dh,
    dk: out.dk,
    dL_rad: degToRad(out.dL_deg),
  };
}

function getChiefDefMode() {
  const el = document.querySelector('input[name="chief-def-mode"]:checked');
  return el?.value === "mee" ? "mee" : "coe";
}

function getDeputyDefMode() {
  const el = document.querySelector('input[name="deputy-def-mode"]:checked');
  if (el?.value === "mee") return "mee";
  if (el?.value === "rmee") return "rmee";
  return "coe";
}

function readChiefMeeForRmeeTab() {
  if (getChiefDefMode() === "mee") {
    return readMeeGroup("chief-mee", getDefaultRmeeChiefMeeKmDeg());
  }
  return coeToMee(readRmeeChiefCoeFromInputs());
}

function readDeputyMeeForRmeeTab(chiefMee) {
  const mode = getDeputyDefMode();
  if (mode === "mee") {
    return readMeeGroup("deputy-mee", getDefaultDeputyMeeKmDeg());
  }
  if (mode === "rmee") {
    return deputyMeeFromRmee(
      chiefMee,
      readRmeeGroup("deputy-rmee", getDefaultDeputyRmeeDisplay())
    );
  }
  return coeToMee(readDeputyCoeFromInputs());
}

function getCoeMeeInputMode() {
  const el = document.querySelector('input[name="coe-mee-input-mode"]:checked');
  return el?.value === "mee-to-coe" ? "mee-to-coe" : "coe-to-mee";
}

function syncCoeMeeInputPanels() {
  const mode = getCoeMeeInputMode();
  document.getElementById("coe-to-mee-panel").hidden = mode !== "coe-to-mee";
  document.getElementById("mee-to-coe-panel").hidden = mode !== "mee-to-coe";
  document.getElementById("coe-mee-result-mee").hidden = mode !== "coe-to-mee";
  document.getElementById("coe-mee-result-coe").hidden = mode !== "mee-to-coe";
}

function syncRmeeInputPanels() {
  if (activeMode !== "rmee") return;
  const chiefMode = getChiefDefMode();
  const deputyMode = getDeputyDefMode();
  document.getElementById("chief-def-coe-panel").hidden = chiefMode !== "coe";
  document.getElementById("chief-def-mee-panel").hidden = chiefMode !== "mee";
  document.getElementById("deputy-def-coe-panel").hidden = deputyMode !== "coe";
  document.getElementById("deputy-def-mee-panel").hidden = deputyMode !== "mee";
  document.getElementById("deputy-def-rmee-panel").hidden = deputyMode !== "rmee";
}

function updateCoeMeeMode() {
  syncCoeMeeInputPanels();
  const inputMode = getCoeMeeInputMode();

  if (inputMode === "coe-to-mee") {
    const coe = readChiefCoeFromInputs();
    const mee = coeToMee(coe);
    renderResultList(document.getElementById("list-coe-mee-mee"), meeResultRows(mee));
    plotChiefOnlyEci(coe);
    return;
  }

  const mee = readMeeGroup("conv-mee", getDefaultChiefMeeKmDeg());
  const coe = meeToCoe(mee);
  renderResultList(
    document.getElementById("list-coe-mee-coe"),
    coeResultRows(coeRadToKmDeg(coe))
  );
  plotChiefOnlyEci(coe);
}

function updateRmeeDefMode() {
  const chief = coeKmDegToRad(CONFIG.DEFAULT_CHIEF_COE);
  const deputy = coeKmDegToRad(CONFIG.DEFAULT_DEPUTY_COE);
  const chiefMee = coeToMee(chief);
  const deputyMee = coeToMee(deputy);
  const rmee = computeRmee(chiefMee, deputyMee);
  const roeKm = roeKmFromDimensionless(coeToRoe(chief, deputy), chief.a_km);
  renderResultStack(document.getElementById("result-rmee-def"), [
    { title: "Chief COE（例）", rows: coeResultRows(CONFIG.DEFAULT_CHIEF_COE) },
    { title: "Deputy COE（例）", rows: coeResultRows(CONFIG.DEFAULT_DEPUTY_COE) },
    { title: "Chief MEE α", rows: meeResultRows(chiefMee) },
    { title: "Deputy MEE α", rows: meeResultRows(deputyMee) },
    { title: "ROE（Gim–Alfriend, a·δ）[km]", rows: roeResultRows(roeKm) },
    {
      title: "RMEE δα = [ δp/p<sub>c</sub>, δL, δf, δg, δh, δk ]",
      rows: rmeeResultRows(rmee),
    },
  ]);
  updatePlotsFromCoePair(chief, deputy);
}

function buildRmeeResultStack(chiefMode, deputyMode, chiefCoe, deputyCoe, chiefMee, deputyMee, rmee) {
  const blocks = [];
  if (chiefMode === "coe") {
    blocks.push({ title: "Chief COE", rows: coeResultRows(coeRadToKmDeg(chiefCoe)) });
  }
  blocks.push({ title: "Chief MEE α", rows: meeResultRows(chiefMee) });
  if (deputyMode === "coe") {
    blocks.push({ title: "Deputy COE", rows: coeResultRows(coeRadToKmDeg(deputyCoe)) });
  }
  blocks.push({
    title: deputyMode === "mee" ? "Deputy MEE α（入力）" : "Deputy MEE α（算出）",
    rows: meeResultRows(deputyMee),
  });
  appendRoeResultBlock(blocks, chiefCoe, deputyCoe);
  blocks.push({
    title: "RMEE δα = [ δp/p<sub>c</sub>, δL, δf, δg, δh, δk ]",
    rows: rmeeResultRows(rmee),
  });
  return blocks;
}

function updateRmeeMode() {
  syncRmeeInputPanels();
  const chiefMode = getChiefDefMode();
  const deputyMode = getDeputyDefMode();
  const chiefMee = readChiefMeeForRmeeTab();
  const deputyMee = readDeputyMeeForRmeeTab(chiefMee);
  const chiefCoe = meeToCoe(chiefMee);
  const deputyCoe = meeToCoe(deputyMee);
  const rmee = computeRmee(chiefMee, deputyMee);

  renderResultStack(
    document.getElementById("result-rmee"),
    buildRmeeResultStack(
      chiefMode,
      deputyMode,
      chiefCoe,
      deputyCoe,
      chiefMee,
      deputyMee,
      rmee
    )
  );

  updatePlotsFromCoePair(chiefCoe, deputyCoe);
}

function updateAll() {
  syncOrbitDisplaySection();
  if (activeMode === "coe-mee") updateCoeMeeMode();
  else if (activeMode === "rmee-def") updateRmeeDefMode();
  else updateRmeeMode();
}

function setActiveMode(mode) {
  activeMode = mode;
  document.querySelectorAll(".mode-tab").forEach((tab) => {
    const on = tab.dataset.mode === mode;
    tab.classList.toggle("is-active", on);
  });
  document.querySelectorAll(".mode-panel").forEach((panel) => {
    panel.hidden = panel.id !== `panel-${mode}`;
  });
  syncOrbitDisplaySection();
  syncPlotUiForMode();
  updateAll();
}

function initModeTabs() {
  document.querySelectorAll(".mode-tab").forEach((tab) => {
    tab.addEventListener("click", () => setActiveMode(tab.dataset.mode));
  });
}

function initRmeeDefModes() {
  document.querySelectorAll('input[name="chief-def-mode"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      syncRmeeInputPanels();
      updateAll();
    });
  });
  document.querySelectorAll('input[name="deputy-def-mode"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      syncRmeeInputPanels();
      updateAll();
    });
  });
}

function initCoeMeeInputMode() {
  document.querySelectorAll('input[name="coe-mee-input-mode"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      syncCoeMeeInputPanels();
      updateAll();
    });
  });
}

function bindInputHandlers() {
  document.querySelectorAll(".number-input").forEach((input) => {
    input.addEventListener("input", updateAll);
    input.addEventListener("change", updateAll);
  });
}

function init() {
  buildMeeInputs(
    document.querySelector('[data-mee-group="conv-mee"]'),
    "conv-mee",
    getDefaultChiefMeeKmDeg()
  );
  buildMeeInputs(
    document.querySelector('[data-mee-group="chief-mee"]'),
    "chief-mee",
    getDefaultRmeeChiefMeeKmDeg()
  );
  buildMeeInputs(
    document.querySelector('[data-mee-group="deputy-mee"]'),
    "deputy-mee",
    getDefaultDeputyMeeKmDeg()
  );
  buildRmeeInputs(
    document.querySelector('[data-rmee-group="deputy-rmee"]'),
    "deputy-rmee",
    getDefaultDeputyRmeeDisplay()
  );
  initModeTabs();
  initCoeMeeInputMode();
  initRmeeDefModes();
  initViewTabs();
  initOrbitControls();
  initEciExaggerationControl();
  bindInputHandlers();
  syncOrbitDisplaySection();
  syncPlotUiForMode();
  setActiveMode("coe-mee");
  window.addEventListener("resize", () => {
    if (plotsInitialized) {
      for (const id of getPlotIdsForView()) Plotly.Plots.resize(id);
    }
  });
}

document.addEventListener("DOMContentLoaded", init);
