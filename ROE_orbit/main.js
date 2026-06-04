/**
 * ROE 相対軌道可視化
 * 近円軌道・小さい相対運動の線形近似（ブラウザ内計算のみ）
 */

// ---------------------------------------------------------------------------
// 設定定数（研究用途で後から変更しやすいよう集約）
// ---------------------------------------------------------------------------

const CONFIG = {
  /** ROE 数値入力の推奨 step [km] */
  ROE_INPUT_STEP_KM: 0.5,

  /** 1 周あたりのサンプル点数 */
  NUM_ORBIT_POINTS: 400,

  /** 時間ドリフト表示のデフォルト周回数 */
  DEFAULT_NUM_DRIFT_ORBITS: 5,

  /** ECI 表示用 chief 軌道のサンプル点数 */
  ECI_CHIEF_ORBIT_POINTS: 360,

  /** ECI 誇張表示のデフォルト k = 10^3 */
  DEFAULT_ECI_EXAGGERATION_LOG: 3,

  /** 軸余白比率 */
  ECI_ZOOM_PADDING_RATIO: 0.2,

  /** 地球赤道半径 [km]（WGS84 赤道半径） */
  R_EARTH_KM: 6378.137,

  /** 地球赤道断面円のサンプル点数 */
  EARTH_EQUATOR_POINTS: 128,

  /** 地球重力定数 [km³/s²]（平均運動 n の計算用） */
  MU_EARTH_KM3_S2: 398600.4418,

  /** Chief COE 初期値（GEO 近円赤道） */
  DEFAULT_CHIEF_COE: {
    a_km: 42164,
    e: 0,
    i_deg: 0,
    raan_deg: 0,
    argp_deg: 0,
    M_deg: 0,
  },

  /**
   * ROE 初期値 [km]（a = DEFAULT_CHIEF_COE.a_km 時の従来無次元値 × a）
   * 内部計算では δ = (km 値) / a に変換
   */
  ROE_DEFAULTS_KM: {
    delta_a: 1,
    delta_lambda: 100,
    delta_ex: 10,
    delta_ey: 0,
    delta_ix: 30,
    delta_iy: 0,
  },
};

/** ROE 数値入力（UI は km、内部は無次元 δ） */
const ROE_FIELDS = [
  { key: "delta_a", inputId: "roe-delta-a" },
  { key: "delta_lambda", inputId: "roe-delta-lambda" },
  { key: "delta_ex", inputId: "roe-delta-ex" },
  { key: "delta_ey", inputId: "roe-delta-ey" },
  { key: "delta_ix", inputId: "roe-delta-ix" },
  { key: "delta_iy", inputId: "roe-delta-iy" },
];

/** COE → ROE 変換後の読み取り専用表示（Deputy COE モード） */
const COE_DERIVED_ROE_OUTPUTS = [
  { key: "delta_a", outputId: "coe-val-delta-a" },
  { key: "delta_lambda", outputId: "coe-val-delta-lambda" },
  { key: "delta_ex", outputId: "coe-val-delta-ex" },
  { key: "delta_ey", outputId: "coe-val-delta-ey" },
  { key: "delta_ix", outputId: "coe-val-delta-ix" },
  { key: "delta_iy", outputId: "coe-val-delta-iy" },
];

/** @type {"roe" | "coe"} */
let inputMode = "roe";

const TWO_PI = 2 * Math.PI;
const DRIFT_LAMBDA_FACTOR = 1.5; // δλ̇ = -(3/2) n δa

/** Plotly が描画する div id 一覧 */
const RTN_PLOT_IDS = ["plot-3d", "plot-tr", "plot-tn", "plot-rn"];
const ECI_PLOT_IDS = [
  "plot-eci-3d",
  "plot-eci-xy",
  "plot-eci-xz",
  "plot-eci-yz",
];
const PLOT_DIV_IDS = [...RTN_PLOT_IDS, ...ECI_PLOT_IDS];

const VIEW_PLOT_IDS = {
  rtn: RTN_PLOT_IDS,
  eci: ECI_PLOT_IDS,
};

/** @type {"rtn" | "eci"} */
let activeView = "rtn";

// Plotly 用の共通スタイル（ダークテーマに合わせる）
const PLOT_LAYOUT_BASE = {
  paper_bgcolor: "#1c2333",
  plot_bgcolor: "#1c2333",
  font: { color: "#e6edf3", family: "Segoe UI, sans-serif", size: 12 },
  margin: { l: 50, r: 20, t: 30, b: 45 },
};

// ---------------------------------------------------------------------------
// 軌道力学
// ---------------------------------------------------------------------------

/** 近円 chief 軌道の平均運動 [rad/s] */
function chiefMeanMotionRadS(a_km) {
  return Math.sqrt(CONFIG.MU_EARTH_KM3_S2 / (a_km * a_km * a_km));
}

/**
 * δa による along-track セクラー変化（近円・線形 ROE）
 * @returns {number} δλ(t) [無次元]
 */
function deltaLambdaAtTimeSec(roe, n_rad_s, t_sec) {
  return roe.delta_lambda - DRIFT_LAMBDA_FACTOR * n_rad_s * roe.delta_a * t_sec;
}

/**
 * 引数緯度 u [rad] における RTN 相対位置 [km]
 * @param {number|null} deltaLambdaOverride - 指定時は roe.delta_lambda の代わりに使用
 */
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

// ---------------------------------------------------------------------------
// 軌道サンプル（RTN / ECI 共通）
// ---------------------------------------------------------------------------

/** u = 0 … 2π、δλ 固定の 1 周サンプル */
function buildSingleOrbitSamples(a_km, roe, numPoints) {
  const samples = [];
  for (let i = 0; i < numPoints; i++) {
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    const pos = rtnAtLatitude(a_km, roe, u_rad);
    samples.push({ u_rad, ...pos });
  }
  return samples;
}

/** 時間ドリフト付き複数周サンプル */
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
      const pos = rtnAtLatitude(a_km, roe, u_rad, delta_lambda_t);
      samples.push({ u_rad, t_sec, ...pos });
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

function computeRelativeOrbitRTN(a_km, roe, numPoints) {
  return samplesToRtnOrbit(buildSingleOrbitSamples(a_km, roe, numPoints));
}

function computeDriftTrajectoryRTN(a_km, roe, numPointsPerOrbit, numOrbits) {
  const { samples, n_rad_s, orbitPeriodSec } = buildDriftOrbitSamples(
    a_km,
    roe,
    numPointsPerOrbit,
    numOrbits
  );
  return { ...samplesToRtnOrbit(samples), n_rad_s, orbitPeriodSec };
}

// ---------------------------------------------------------------------------
// ECI 変換（一般 COE の楕円・傾斜 chief 軌道）
// ---------------------------------------------------------------------------

function normalizeAngleRad(rad) {
  return ((rad % TWO_PI) + TWO_PI) % TWO_PI;
}

/** ケプラー方程式 E − e sin E = M を Newton 法で解く */
function solveKeplerE(M_rad, e, tol = 1e-12) {
  const M = normalizeAngleRad(M_rad);
  if (e < 1e-12) return M;

  let E = e < 0.8 ? M : Math.PI;
  for (let k = 0; k < 60; k++) {
    const f = E - e * Math.sin(E) - M;
    const fp = 1 - e * Math.cos(E);
    const dE = f / fp;
    E -= dE;
    if (Math.abs(dE) < tol) break;
  }
  return E;
}

function trueAnomalyFromArgumentOfLatitude(u_rad, argp_rad) {
  return u_rad - argp_rad;
}

function radiusFromTrueAnomaly(a_km, e, nu_rad) {
  return (a_km * (1 - e * e)) / (1 + e * Math.cos(nu_rad));
}

/** 引数緯度 u における chief 位置 [km]（ECI） */
function chiefPositionEciKm(coe, u_rad) {
  const { a_km, e, i_rad, raan_rad, argp_rad } = coe;
  const nu = trueAnomalyFromArgumentOfLatitude(u_rad, argp_rad);
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

/** 真近点角 ν における chief 速度 [km/s]（ECI, PQW→ECI 回転） */
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

/**
 * 引数緯度 u における RTN 基底単位ベクトル（ECI 成分）
 * R̂: 半径方向, T̂: 沿軌道方向, N̂ = R̂ × T̂
 */
function rtnBasisEci(coe, u_rad) {
  const { argp_rad } = coe;
  const nu = trueAnomalyFromArgumentOfLatitude(u_rad, argp_rad);
  const r_km = radiusFromTrueAnomaly(coe.a_km, coe.e, nu);
  const r_vec = chiefPositionEciKm(coe, u_rad);
  const v_vec = chiefVelocityEciKmS(coe, nu, r_km);

  const R_hat = normalizeVec3(r_vec);
  const h_vec = crossVec3(r_vec, v_vec);
  const N_hat = normalizeVec3(h_vec);
  const T_hat = normalizeVec3(crossVec3(N_hat, R_hat));

  return { R_hat, T_hat, N_hat };
}

/** RTN 相対変位を ECI ベクトル [km] に変換 */
function relativeRtnToEciKm(R_km, T_km, N_km, coe, u_rad) {
  const { R_hat, T_hat, N_hat } = rtnBasisEci(coe, u_rad);
  return [
    R_km * R_hat[0] + T_km * T_hat[0] + N_km * N_hat[0],
    R_km * R_hat[1] + T_km * T_hat[1] + N_km * N_hat[1],
    R_km * R_hat[2] + T_km * T_hat[2] + N_km * N_hat[2],
  ];
}

/**
 * 地球赤道断面（ECI 赤道面 z=0 上の円, 半径 R⊕）
 * ワイヤーフレーム球の赤道断面として表示
 */
function buildEarthEquatorCrossSectionEci() {
  const r = CONFIG.R_EARTH_KM;
  const n = CONFIG.EARTH_EQUATOR_POINTS;
  const x_km = [];
  const y_km = [];
  const z_km = [];

  for (let i = 0; i <= n; i++) {
    const theta = (TWO_PI * i) / n;
    x_km.push(r * Math.cos(theta));
    y_km.push(r * Math.sin(theta));
    z_km.push(0);
  }

  return { x_km, y_km, z_km };
}

/** chief の 1 周軌道（ECI, 引数緯度 u でサンプル） */
function buildChiefOrbitEci(coe, numPoints) {
  const x_km = [];
  const y_km = [];
  const z_km = [];

  for (let i = 0; i < numPoints; i++) {
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    const [x, y, z] = chiefPositionEciKm(coe, u_rad);
    x_km.push(x);
    y_km.push(y);
    z_km.push(z);
  }

  return { x_km, y_km, z_km };
}

/**
 * ECI 誇張表示: r_display = r_chief + k · (r_deputy − r_chief)
 * @param {number} exaggerationK - 表示専用の倍率（1 = 実スケール）
 */
function samplesToEciExaggeratedOrbit(coe, samples, exaggerationK) {
  const x_km = [];
  const y_km = [];
  const z_km = [];

  for (const sample of samples) {
    const chief = chiefPositionEciKm(coe, sample.u_rad);
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

/** 有効な Plotly 軸範囲 [min, max] */
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

/** ECI プロットの軸範囲 [km]（chief 軌道 + 誇張 deputy を含める） */
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

  return {
    x: axisRange(coords.x),
    y: axisRange(coords.y),
    z: axisRange(coords.z),
  };
}

// ---------------------------------------------------------------------------
// COE ↔ ROE
// ---------------------------------------------------------------------------

function degToRad(deg) {
  return (deg * Math.PI) / 180;
}

function radToDeg(rad) {
  return (rad * 180) / Math.PI;
}

function defaultChiefCoeKmDeg() {
  return { ...CONFIG.DEFAULT_CHIEF_COE };
}

function chiefCoeKmDegToRad(coeKmDeg) {
  return {
    a_km: coeKmDeg.a_km,
    e: coeKmDeg.e,
    i_rad: degToRad(coeKmDeg.i_deg),
    raan_rad: degToRad(coeKmDeg.raan_deg),
    argp_rad: degToRad(coeKmDeg.argp_deg),
    M_rad: degToRad(coeKmDeg.M_deg),
  };
}

function chiefCoeRadToKmDeg(coe) {
  return {
    a_km: coe.a_km,
    e: coe.e,
    i_deg: radToDeg(coe.i_rad),
    raan_deg: radToDeg(coe.raan_rad),
    argp_deg: radToDeg(coe.argp_rad),
    M_deg: radToDeg(coe.M_rad),
  };
}

function readChiefCoeFromInputs() {
  const defaults = defaultChiefCoeKmDeg();

  const readNum = (id, fallback) => {
    const el = document.getElementById(id);
    const v = parseFloat(el?.value, 10);
    return Number.isFinite(v) ? v : fallback;
  };

  return chiefCoeKmDegToRad({
    a_km: readNum("chief-a", defaults.a_km),
    e: readNum("chief-e", defaults.e),
    i_deg: readNum("chief-i", defaults.i_deg),
    raan_deg: readNum("chief-raan", defaults.raan_deg),
    argp_deg: readNum("chief-argp", defaults.argp_deg),
    M_deg: readNum("chief-M", defaults.M_deg),
  });
}

function getChiefCoe() {
  return readChiefCoeFromInputs();
}

function applyChiefCoeToInputs(coeKmDeg) {
  const fields = [
    ["chief-a", coeKmDeg.a_km, (v) => String(Math.round(v))],
    ["chief-e", coeKmDeg.e, (v) => v.toFixed(6)],
    ["chief-i", coeKmDeg.i_deg, (v) => v.toFixed(5)],
    ["chief-raan", coeKmDeg.raan_deg, (v) => v.toFixed(2)],
    ["chief-argp", coeKmDeg.argp_deg, (v) => v.toFixed(2)],
    ["chief-M", coeKmDeg.M_deg, (v) => v.toFixed(5)],
  ];

  for (const [id, value, fmt] of fields) {
    const el = document.getElementById(id);
    if (el) el.value = fmt(value);
  }
}

/**
 * Gim–Alfriend 準非特異 ROE（無次元）
 * @see D'Amico, Gim & Alfriend (2003)
 */
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

/** 現在の ROE 数値入力 [km] から Deputy COE を逆算 */
function deputyCoeFromRoeKm(roeKm, chiefCoe) {
  const a_km = chiefCoe.a_km;
  const roe = roeDimensionlessFromKm(roeKm, a_km);
  const { e, i_rad: i, raan_rad: Om, argp_rad: w, M_rad: M } = chiefCoe;

  const ex = roe.delta_ex + e * Math.cos(w);
  const ey = roe.delta_ey + e * Math.sin(w);

  let Omd = Om;
  if (Math.abs(Math.sin(i)) > 1e-10) {
    Omd = Om + roe.delta_iy / Math.sin(i);
  }

  const Md = roe.delta_lambda + M - (Omd - Om) * Math.cos(i);

  return {
    a_km: a_km * (1 + roe.delta_a),
    e: Math.hypot(ex, ey),
    i_deg: radToDeg(i + roe.delta_ix),
    raan_deg: radToDeg(normalizeAngleRad(Omd)),
    argp_deg: radToDeg(normalizeAngleRad(Math.atan2(ey, ex))),
    M_deg: radToDeg(normalizeAngleRad(Md)),
  };
}

function getDefaultDeputyCoe() {
  return deputyCoeFromRoeKm(CONFIG.ROE_DEFAULTS_KM, getChiefCoe());
}

function readInputMode() {
  const el = document.querySelector('input[name="input-mode"]:checked');
  return el && el.value === "coe" ? "coe" : "roe";
}

function readDeputyCoeFromInputs() {
  const defaults = getDefaultDeputyCoe();

  const readNum = (id, fallback) => {
    const el = document.getElementById(id);
    const v = parseFloat(el?.value, 10);
    return Number.isFinite(v) ? v : fallback;
  };

  return {
    a_km: readNum("deputy-a", defaults.a_km),
    e: readNum("deputy-e", defaults.e),
    i_rad: degToRad(readNum("deputy-i", defaults.i_deg)),
    raan_rad: degToRad(readNum("deputy-raan", defaults.raan_deg)),
    argp_rad: degToRad(readNum("deputy-argp", defaults.argp_deg)),
    M_rad: degToRad(readNum("deputy-M", defaults.M_deg)),
  };
}

function roeKmFromDimensionless(roe, a_km) {
  const roeKm = {};
  for (const { key } of ROE_FIELDS) {
    roeKm[key] = roe[key] * a_km;
  }
  return roeKm;
}

function readRoeForPlot() {
  const chief = getChiefCoe();
  if (readInputMode() === "coe") {
    return coeToRoe(chief, readDeputyCoeFromInputs());
  }
  return readRoeFromInputs();
}

function readRoeKmForDisplay() {
  return roeKmFromDimensionless(readRoeForPlot(), getChiefCoe().a_km);
}

function updateCoeDerivedDisplay() {
  const roeKm = readRoeKmForDisplay();
  for (const { key, outputId } of COE_DERIVED_ROE_OUTPUTS) {
    const output = document.getElementById(outputId);
    if (output) output.textContent = formatKm(roeKm[key]);
  }
}

function applyRoeKmToInputs(roeKm) {
  for (const { key, inputId } of ROE_FIELDS) {
    const input = document.getElementById(inputId);
    if (!input) continue;
    input.value = String(roeKm[key]);
  }
}

function applyDeputyCoeToInputs(deputyCoe) {
  const fields = [
    ["deputy-a", deputyCoe.a_km, (v) => String(Math.round(v))],
    ["deputy-e", deputyCoe.e, (v) => v.toFixed(6)],
    ["deputy-i", deputyCoe.i_deg, (v) => v.toFixed(5)],
    ["deputy-raan", deputyCoe.raan_deg, (v) => v.toFixed(2)],
    ["deputy-argp", deputyCoe.argp_deg, (v) => v.toFixed(2)],
    ["deputy-M", deputyCoe.M_deg, (v) => v.toFixed(5)],
  ];

  for (const [id, value, fmt] of fields) {
    const el = document.getElementById(id);
    if (el) el.value = fmt(value);
  }
}

function syncDeputyCoeFromRoeInputs() {
  applyDeputyCoeToInputs(
    deputyCoeFromRoeKm(readRoeKmFromInputs(), getChiefCoe())
  );
}

function syncRoeInputsFromDeputyCoe() {
  applyRoeKmToInputs(readRoeKmForDisplay());
}

function setInputMode(mode) {
  if (mode !== "roe" && mode !== "coe") return;

  inputMode = mode;

  const roePanel = document.getElementById("roe-input-panel");
  const coePanel = document.getElementById("coe-input-panel");
  if (!roePanel || !coePanel) return;

  const isRoe = mode === "roe";
  roePanel.hidden = !isRoe;
  coePanel.hidden = isRoe;

  const radio = document.querySelector(`input[name="input-mode"][value="${mode}"]`);
  if (radio) radio.checked = true;

  if (isRoe) {
    syncRoeInputsFromDeputyCoe();
  } else {
    syncDeputyCoeFromRoeInputs();
    updateCoeDerivedDisplay();
  }
}

function initInputModeControl() {
  const radios = document.querySelectorAll('input[name="input-mode"]');
  if (radios.length === 0) return;

  radios.forEach((radio) => {
    radio.addEventListener("change", () => {
      if (!radio.checked) return;
      const nextMode = radio.value === "coe" ? "coe" : "roe";
      if (nextMode === inputMode) return;
      setInputMode(nextMode);
      updateDriftRateDisplay(readSemiMajorAxisKm(), readRoeKmForDisplay());
      updateAllPlots();
    });
  });

  setInputMode(inputMode);
}

function initDeputyCoeInputs() {
  applyDeputyCoeToInputs(getDefaultDeputyCoe());

  document.querySelectorAll(".coe-input").forEach((input) => {
    const handler = () => {
      if (readInputMode() !== "coe") return;
      updateCoeDerivedDisplay();
      updateDriftRateDisplay(getChiefCoe().a_km, readRoeKmForDisplay());
      updateAllPlots();
    };
    input.addEventListener("input", handler);
    input.addEventListener("change", handler);
  });
}

function initChiefCoeInputs() {
  applyChiefCoeToInputs(defaultChiefCoeKmDeg());

  document.querySelectorAll(".chief-input").forEach((input) => {
    const handler = () => {
      refreshRelativeInputDisplays();
      updateAllPlots();
    };
    input.addEventListener("input", handler);
    input.addEventListener("change", handler);
  });
}

function refreshRelativeInputDisplays() {
  if (readInputMode() === "coe") {
    updateCoeDerivedDisplay();
  }
  updateDriftRateDisplay(getChiefCoe().a_km, readRoeKmForDisplay());
}

function getChiefStructureKey() {
  const c = getChiefCoe();
  return [
    c.a_km.toFixed(3),
    c.e.toFixed(6),
    radToDeg(c.i_rad).toFixed(3),
    radToDeg(c.raan_rad).toFixed(2),
    radToDeg(c.argp_rad).toFixed(2),
    radToDeg(c.M_rad).toFixed(3),
  ].join("|");
}

// ---------------------------------------------------------------------------
// UI ヘルパ
// ---------------------------------------------------------------------------

function formatKm(value) {
  if (!Number.isFinite(value)) return "—";
  if (Math.abs(value) < 1e-9) return "0";
  const abs = Math.abs(value);
  if (abs >= 100 || abs < 0.01) return value.toExponential(2);
  return value.toFixed(2);
}

/** 数値入力 [km] → 無次元 ROE（δ ≈ 線形スケール値 / a） */
function roeDimensionlessFromKm(roeKm, a_km) {
  const roe = {};
  for (const { key } of ROE_FIELDS) {
    roe[key] = roeKm[key] / a_km;
  }
  return roe;
}

function readSemiMajorAxisKm() {
  return readChiefCoeFromInputs().a_km;
}

function readRoeKmFromInputs() {
  const defaults = CONFIG.ROE_DEFAULTS_KM;
  const roeKm = {};
  for (const { key, inputId } of ROE_FIELDS) {
    const el = document.getElementById(inputId);
    const v = parseFloat(el?.value, 10);
    roeKm[key] = Number.isFinite(v) ? v : defaults[key];
  }
  return roeKm;
}

function readRoeFromInputs() {
  const a_km = readSemiMajorAxisKm();
  return roeDimensionlessFromKm(readRoeKmFromInputs(), a_km);
}

function readNumDriftOrbits() {
  const el = document.getElementById("num-drift-orbits");
  const v = parseInt(el.value, 10);
  return Number.isFinite(v) && v >= 1 ? v : CONFIG.DEFAULT_NUM_DRIFT_ORBITS;
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

  const roe = roeDimensionlessFromKm(roeKm, a_km);
  const n = chiefMeanMotionRadS(a_km);
  const dDeltaLambdaDt = -DRIFT_LAMBDA_FACTOR * n * roe.delta_a;
  const driftPerOrbitKm = -DRIFT_LAMBDA_FACTOR * TWO_PI * deltaAKm;
  const periodHr = TWO_PI / n / 3600;

  el.textContent =
    `δλ̇ ≈ ${formatKm(dDeltaLambdaDt * a_km)}/s（along-track スケール）\n` +
    `1 周あたり ΔT ≈ ${formatKm(driftPerOrbitKm)} km\n` +
    `T_orbit ≈ ${periodHr.toFixed(2)} h`;
}

function initRoeInputs() {
  applyRoeKmToInputs(CONFIG.ROE_DEFAULTS_KM);

  document.querySelectorAll(".roe-input").forEach((input) => {
    const handler = () => {
      if (readInputMode() !== "roe") return;
      updateDriftRateDisplay(readSemiMajorAxisKm(), readRoeKmFromInputs());
      updateAllPlots();
    };
    input.addEventListener("input", handler);
    input.addEventListener("change", handler);
  });
}

/** Plotly.react はトレース数変更に弱い → 構造が変わったら newPlot */
function getPlotStructureKey() {
  const showRef = document.getElementById("show-single-orbit-ref").checked;
  return (
    `mode=${readInputMode()}|chief=${getChiefStructureKey()}|ref=${showRef ? 1 : 0}|orb=${readNumDriftOrbits()}|k=${readEciExaggerationFactor()}`
  );
}

function initOrbitControls() {
  const orbitSlider = document.getElementById("num-drift-orbits");
  const orbitOutput = document.getElementById("val-num-drift-orbits");
  orbitSlider.value = String(CONFIG.DEFAULT_NUM_DRIFT_ORBITS);
  orbitOutput.textContent = String(CONFIG.DEFAULT_NUM_DRIFT_ORBITS);

  orbitSlider.addEventListener("input", () => {
    orbitOutput.textContent = orbitSlider.value;
    updateAllPlots();
  });

  document
    .getElementById("show-single-orbit-ref")
    .addEventListener("change", () => updateAllPlots(true));

  initEciExaggerationControl();
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
    updateAllPlots();
  });
}

// ---------------------------------------------------------------------------
// Plotly 描画
// ---------------------------------------------------------------------------

let plotsInitialized = false;
let lastPlotStructureKey = null;

/** 軌道線の始点・終点マーカー（全プロット共通） */
const ORBIT_ENDPOINT_MARKER = {
  start: { color: "#58a6ff", size: 9, symbol: "circle", line: { color: "#e6edf3", width: 1 } },
  end: { color: "#f85149", size: 9, symbol: "square", line: { color: "#e6edf3", width: 1 } },
};

function pushOrbitEndpointMarkers3D(traces, coords, labelPrefix, axisLabels) {
  const { x, y, z } = coords;
  if (!x || x.length < 1) return;

  const axes = axisLabels || { x: "X", y: "Y", z: "Z" };
  const last = x.length - 1;
  const hoverStart =
    `${labelPrefix} 始点 (u ≈ 0)<br>` +
    `${axes.x}: %{x:.3f} ${axes.y}: %{y:.3f} ${axes.z}: %{z:.3f} km<extra></extra>`;
  const hoverEnd =
    `${labelPrefix} 終点<br>` +
    `${axes.x}: %{x:.3f} ${axes.y}: %{y:.3f} ${axes.z}: %{z:.3f} km<extra></extra>`;

  traces.push({
    type: "scatter3d",
    mode: "markers",
    name: `${labelPrefix} 始点`,
    x: [x[0]],
    y: [y[0]],
    z: [z[0]],
    marker: { ...ORBIT_ENDPOINT_MARKER.start },
    hovertemplate: hoverStart,
  });
  traces.push({
    type: "scatter3d",
    mode: "markers",
    name: `${labelPrefix} 終点`,
    x: [x[last]],
    y: [y[last]],
    z: [z[last]],
    marker: { ...ORBIT_ENDPOINT_MARKER.end },
    hovertemplate: hoverEnd,
  });
}

const RTN_AXIS_LABELS_3D = { x: "T", y: "R", z: "N" };
const ECI_AXIS_LABELS_3D = { x: "X", y: "Y", z: "Z" };

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
    hovertemplate:
      "T: %{x:.3f} km<br>R: %{y:.3f} km<br>N: %{z:.3f} km<extra></extra>",
  };
}

function buildOrbitTraces3D(primaryOrbit, referenceOrbit) {
  const traces = [];

  if (referenceOrbit) {
    traces.push(
      buildDeputyTrace3D(
        referenceOrbit,
        "参考: 1 周（δλ 固定）",
        "rgba(139, 148, 158, 0.55)",
        2
      )
    );
  }

  traces.push(
    buildDeputyTrace3D(
      primaryOrbit,
      "Deputy（相対軌道）",
      "#3fb950",
      4
    )
  );
  pushOrbitEndpointMarkers3D(
    traces,
    {
      x: primaryOrbit.T_km,
      y: primaryOrbit.R_km,
      z: primaryOrbit.N_km,
    },
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

function build2DTraces(orbit, plane, referenceOrbit) {
  const traces = [];
  const { axisX, axisY, xTitle, yTitle, color } = plane;

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

  traces.push(
    build2DTrace(
      orbit[axisX],
      orbit[axisY],
      "相対軌道",
      color,
      2
    )
  );
  pushOrbitEndpointMarkers2D(traces, orbit[axisX], orbit[axisY], "Deputy");
  traces.push(buildChief2D());

  const planeId =
    plane === PLANE_TR ? "tr" : plane === PLANE_TN ? "tn" : "rn";
  return { traces, layout: layout2D(xTitle, yTitle, 260, planeId) };
}

const PLANE_TR = {
  axisX: "T_km",
  axisY: "R_km",
  xTitle: "T [km]",
  yTitle: "R [km]",
  color: "#3fb950",
};
const PLANE_TN = {
  axisX: "T_km",
  axisY: "N_km",
  xTitle: "T [km]",
  yTitle: "N [km]",
  color: "#58a6ff",
};
const PLANE_RN = {
  axisX: "R_km",
  axisY: "N_km",
  xTitle: "R [km]",
  yTitle: "N [km]",
  color: "#a371f7",
};

function layout2D(xTitle, yTitle, height, planeId) {
  // T 方向ドリフトで範囲が伸びやすい → T 軸を含む面は等尺を外す
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
    yaxis: {
      title: yTitle,
      gridcolor: "#30363d",
      zerolinecolor: "#484f58",
    },
    showlegend: true,
    legend: { font: { size: 10 } },
  };
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

  const showRef = document.getElementById("show-single-orbit-ref").checked;
  const referenceSamples = showRef
    ? buildSingleOrbitSamples(a_km, roe, numPoints)
    : null;

  return {
    primary: samplesToRtnOrbit(primarySamples),
    reference: referenceSamples
      ? samplesToRtnOrbit(referenceSamples)
      : null,
    primarySamples,
    referenceSamples,
  };
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
    hovertemplate:
      "X: %{x:.3f}<br>Y: %{y:.3f}<br>Z: %{z:.3f} km<extra></extra>",
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
  const deputyLabel = `Deputy（k=${exaggerationK}）`;

  const earthTrace = buildEarthEquatorTrace3D();
  earthTrace.hovertemplate =
    `R⊕ = ${CONFIG.R_EARTH_KM.toFixed(1)} km（赤道断面）<extra></extra>`;

  const traces = [
    earthTrace,
    buildEciLineTrace3D(
      chiefOrbitEci,
      "Chief 軌道",
      "rgba(240, 193, 75, 0.85)",
      2,
      null
    ),
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
    buildEciLineTrace3D(primaryEci, deputyLabel, "#3fb950", 4, null)
  );
  pushOrbitEndpointMarkers3D(
    traces,
    {
      x: primaryEci.x_km,
      y: primaryEci.y_km,
      z: primaryEci.z_km,
    },
    "Deputy",
    ECI_AXIS_LABELS_3D
  );

  return traces;
}

function layoutEci3D(axisRanges) {
  const axisStyle = {
    gridcolor: "#30363d",
    zerolinecolor: "#484f58",
  };

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
      aspectmode: "data",
      camera: { eye: { x: 1.35, y: 1.25, z: 0.85 } },
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

function layoutEci2D(xTitle, yTitle, height, axisRanges) {
  return {
    ...PLOT_LAYOUT_BASE,
    height,
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
  };
}

function buildEci2DTraces(
  primaryEci,
  referenceEci,
  chiefOrbitEci,
  plane,
  exaggerationK
) {
  const { axisX, axisY, xTitle, yTitle } = plane;
  const earthEquator = buildEarthEquatorCrossSectionEci();
  const traces = [
    buildEci2DTrace(
      earthEquator,
      axisX,
      axisY,
      `地球（R⊕ 断面）`,
      "rgba(110, 159, 239, 0.9)",
      1.5,
      "dot"
    ),
    buildEci2DTrace(
      chiefOrbitEci,
      axisX,
      axisY,
      "Chief 軌道",
      "rgba(240, 193, 75, 0.85)",
      1.5,
      "dot"
    ),
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

  const deputyLabel = `Deputy（k=${exaggerationK}）`;

  traces.push(
    buildEci2DTrace(primaryEci, axisX, axisY, deputyLabel, plane.color, 2, null)
  );
  pushOrbitEndpointMarkers2D(
    traces,
    primaryEci[axisX],
    primaryEci[axisY],
    "Deputy"
  );

  const allOrbits = [
    buildEarthEquatorCrossSectionEci(),
    chiefOrbitEci,
    primaryEci,
    referenceEci,
  ];
  const axisRanges = {
    x: computeEciAxisRangesForPlane(allOrbits, axisX),
    y: computeEciAxisRangesForPlane(allOrbits, axisY),
  };

  return { traces, layout: layoutEci2D(xTitle, yTitle, 260, axisRanges) };
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
  if (!plotDivExists(divId)) {
    console.warn(`Plot container #${divId} not found`);
    return Promise.resolve();
  }

  if (useNewPlot) {
    return Plotly.newPlot(divId, traces, layout, plotOpts);
  }
  return Plotly.react(divId, traces, layout);
}

function getPlotIdsForView(view = activeView) {
  return VIEW_PLOT_IDS[view] ?? RTN_PLOT_IDS;
}

function resizePlots(plotIds) {
  for (const id of plotIds) {
    if (plotDivExists(id)) {
      Plotly.Plots.resize(id);
    }
  }
}

function schedulePlotResize(view = activeView) {
  requestAnimationFrame(() => {
    resizePlots(getPlotIdsForView(view));
  });
}

function setActiveView(view) {
  if (view !== "rtn" && view !== "eci") return;

  activeView = view;

  const tabRtn = document.getElementById("view-tab-rtn");
  const tabEci = document.getElementById("view-tab-eci");
  const panelRtn = document.getElementById("view-panel-rtn");
  const panelEci = document.getElementById("view-panel-eci");
  if (!tabRtn || !tabEci || !panelRtn || !panelEci) return;

  const isRtn = view === "rtn";

  tabRtn.classList.toggle("is-active", isRtn);
  tabEci.classList.toggle("is-active", !isRtn);
  tabRtn.setAttribute("aria-selected", String(isRtn));
  tabEci.setAttribute("aria-selected", String(!isRtn));
  tabRtn.tabIndex = isRtn ? 0 : -1;
  tabEci.tabIndex = isRtn ? -1 : 0;

  panelRtn.classList.toggle("is-active", isRtn);
  panelEci.classList.toggle("is-active", !isRtn);
  panelRtn.hidden = !isRtn;
  panelEci.hidden = isRtn;

  schedulePlotResize(view);
}

function initViewTabs() {
  const tabs = document.querySelectorAll(".view-tab[data-view]");
  if (tabs.length === 0) return;

  tabs.forEach((tab) => {
    tab.addEventListener("click", () => {
      const view = tab.getAttribute("data-view");
      if (view === activeView) return;
      setActiveView(view);
    });
  });

  setActiveView(activeView);
}

function drawAllPlots(plotBundle, forceNewPlot) {
  const plotOpts = { responsive: true, displayModeBar: true };
  const structureKey = getPlotStructureKey();
  const useNewPlot =
    forceNewPlot ||
    !plotsInitialized ||
    structureKey !== lastPlotStructureKey;

  lastPlotStructureKey = structureKey;

  if (useNewPlot && plotsInitialized) {
    for (const id of PLOT_DIV_IDS) {
      if (plotDivExists(id)) {
        Plotly.purge(id);
      }
    }
  }

  const {
    traces3d,
    layout3d,
    tr,
    tn,
    rn,
    tracesEci3d,
    layoutEci3d,
    eciXy,
    eciXz,
    eciYz,
  } = plotBundle;

  Promise.all([
    drawPlot("plot-3d", traces3d, layout3d, useNewPlot, plotOpts),
    drawPlot("plot-tr", tr.traces, tr.layout, useNewPlot, plotOpts),
    drawPlot("plot-tn", tn.traces, tn.layout, useNewPlot, plotOpts),
    drawPlot("plot-rn", rn.traces, rn.layout, useNewPlot, plotOpts),
    drawPlot("plot-eci-3d", tracesEci3d, layoutEci3d, useNewPlot, plotOpts),
    drawPlot("plot-eci-xy", eciXy.traces, eciXy.layout, useNewPlot, plotOpts),
    drawPlot("plot-eci-xz", eciXz.traces, eciXz.layout, useNewPlot, plotOpts),
    drawPlot("plot-eci-yz", eciYz.traces, eciYz.layout, useNewPlot, plotOpts),
  ])
    .then(schedulePlotResize)
    .catch((err) => {
      console.error("Plotly draw failed:", err);
      const cap = document.getElementById("plot-eci-caption");
      if (cap) {
        cap.textContent = `プロット描画エラー: ${err.message}`;
      }
    });

  plotsInitialized = true;
}

function updateAllPlots(forceNewPlot = false) {
  const chief = getChiefCoe();
  const a_km = chief.a_km;
  const roe = readRoeForPlot();
  const { primary, reference, primarySamples, referenceSamples } =
    computeOrbitsForDisplay(a_km, roe);

  refreshRelativeInputDisplays();

  const traces3d = buildOrbitTraces3D(primary, reference);
  const layout3d = layout3D();
  const tr = build2DTraces(primary, PLANE_TR, reference);
  const tn = build2DTraces(primary, PLANE_TN, reference);
  const rn = build2DTraces(primary, PLANE_RN, reference);

  const {
    primaryEci,
    referenceEci,
    chiefOrbitEci,
    exaggerationK,
  } = buildEciExaggeratedOrbitData(chief, primarySamples, referenceSamples);
  updateEciCaption(exaggerationK);

  const earthEquator = buildEarthEquatorCrossSectionEci();
  const axisRanges3d = computeEciAxisRanges([
    earthEquator,
    chiefOrbitEci,
    primaryEci,
    referenceEci,
  ]);

  const plotBundle = {
    traces3d,
    layout3d,
    tr,
    tn,
    rn,
    tracesEci3d: buildEciTraces3D(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      exaggerationK
    ),
    layoutEci3d: layoutEci3D(axisRanges3d),
    eciXy: buildEci2DTraces(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      PLANE_ECI_XY,
      exaggerationK
    ),
    eciXz: buildEci2DTraces(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      PLANE_ECI_XZ,
      exaggerationK
    ),
    eciYz: buildEci2DTraces(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      PLANE_ECI_YZ,
      exaggerationK
    ),
  };

  drawAllPlots(plotBundle, forceNewPlot);
}

// ---------------------------------------------------------------------------
// 初期化
// ---------------------------------------------------------------------------

function init() {
  initRoeInputs();
  initChiefCoeInputs();
  initDeputyCoeInputs();
  initInputModeControl();
  initOrbitControls();

  initViewTabs();

  window.addEventListener("resize", () => {
    if (plotsInitialized) {
      resizePlots(getPlotIdsForView());
    }
  });

  updateEciCaption(readEciExaggerationFactor());
  refreshRelativeInputDisplays();
  updateAllPlots();
}

document.addEventListener("DOMContentLoaded", init);
