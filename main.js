/**
 * ROE → RTN 相対軌道可視化
 * 近円軌道・小さい相対運動の線形近似（ブラウザ内計算のみ）
 */

// ---------------------------------------------------------------------------
// 設定定数（研究用途で後から変更しやすいよう集約）
// ---------------------------------------------------------------------------

const CONFIG = {
  /** スライダー範囲 [km]（線形近似: 無次元 ROE ≈ 値/a） */
  SLIDER_MIN_KM: -50,
  SLIDER_MAX_KM: 50,
  SLIDER_STEP_KM: 0.5,

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

  /** chief 半長軸の初期値 [km]（GEO 想定） */
  DEFAULT_SEMI_MAJOR_AXIS_KM: 42164,

  /**
   * スライダー初期値 [km]（a = DEFAULT_SEMI_MAJOR_AXIS_KM 時の従来無次元値 × a）
   * 内部計算では δ = (km 値) / a に変換
   */
  ROE_DEFAULTS_KM: {
    delta_a: 0,
    delta_lambda: 4.2164,
    delta_ex: 4.2164,
    delta_ey: 0,
    delta_ix: 4.2164,
    delta_iy: 0,
  },
};

/** ROE スライダー（UI は km、内部は無次元 δ） */
const ROE_FIELDS = [
  { key: "delta_a", sliderId: "roe-delta-a", outputId: "val-delta-a" },
  { key: "delta_lambda", sliderId: "roe-delta-lambda", outputId: "val-delta-lambda" },
  { key: "delta_ex", sliderId: "roe-delta-ex", outputId: "val-delta-ex" },
  { key: "delta_ey", sliderId: "roe-delta-ey", outputId: "val-delta-ey" },
  { key: "delta_ix", sliderId: "roe-delta-ix", outputId: "val-delta-ix" },
  { key: "delta_iy", sliderId: "roe-delta-iy", outputId: "val-delta-iy" },
];

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
// ECI 変換（赤道面内近円 chief 軌道を仮定）
// ---------------------------------------------------------------------------

/**
 * 引数緯度 u における RTN 基底単位ベクトル（ECI 成分）
 * R̂: 地心→chief 方向, T̂: 速度方向, N̂ = R̂×T̂
 */
function rtnBasisEci(u_rad) {
  const cu = Math.cos(u_rad);
  const su = Math.sin(u_rad);
  return {
    R_hat: [cu, su, 0],
    T_hat: [-su, cu, 0],
    N_hat: [0, 0, 1],
  };
}

/** chief 位置 [km]（ECI） */
function chiefPositionEciKm(a_km, u_rad) {
  const { R_hat } = rtnBasisEci(u_rad);
  return [a_km * R_hat[0], a_km * R_hat[1], a_km * R_hat[2]];
}

/** RTN 相対変位を ECI ベクトル [km] に変換 */
function relativeRtnToEciKm(R_km, T_km, N_km, u_rad) {
  const { R_hat, T_hat, N_hat } = rtnBasisEci(u_rad);
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

/** chief の 1 周軌道（ECI） */
function buildChiefOrbitEci(a_km, numPoints) {
  const x_km = [];
  const y_km = [];
  const z_km = [];

  for (let i = 0; i < numPoints; i++) {
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    const [x, y, z] = chiefPositionEciKm(a_km, u_rad);
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
function samplesToEciExaggeratedOrbit(a_km, samples, exaggerationK) {
  const x_km = [];
  const y_km = [];
  const z_km = [];

  for (const sample of samples) {
    const chief = chiefPositionEciKm(a_km, sample.u_rad);
    const delta = relativeRtnToEciKm(
      sample.R_km,
      sample.T_km,
      sample.N_km,
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

function buildEciExaggeratedOrbitData(a_km, primarySamples, referenceSamples) {
  const k = readEciExaggerationFactor();
  return {
    exaggerationK: k,
    primaryEci: samplesToEciExaggeratedOrbit(a_km, primarySamples, k),
    referenceEci: referenceSamples
      ? samplesToEciExaggeratedOrbit(a_km, referenceSamples, k)
      : null,
    chiefOrbitEci: buildChiefOrbitEci(a_km, CONFIG.ECI_CHIEF_ORBIT_POINTS),
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
// UI ヘルパ
// ---------------------------------------------------------------------------

function formatKm(value) {
  if (!Number.isFinite(value)) return "—";
  if (Math.abs(value) < 1e-9) return "0";
  const abs = Math.abs(value);
  if (abs >= 100 || abs < 0.01) return value.toExponential(2);
  return value.toFixed(2);
}

/** スライダー [km] → 無次元 ROE（δ ≈ 線形スケール値 / a） */
function roeDimensionlessFromKm(roeKm, a_km) {
  const roe = {};
  for (const { key } of ROE_FIELDS) {
    roe[key] = roeKm[key] / a_km;
  }
  return roe;
}

function readSemiMajorAxisKm() {
  const el = document.getElementById("semi-major-axis");
  const v = parseFloat(el.value, 10);
  return Number.isFinite(v) && v > 0 ? v : CONFIG.DEFAULT_SEMI_MAJOR_AXIS_KM;
}

function readRoeKmFromSliders() {
  const roeKm = {};
  for (const { key, sliderId } of ROE_FIELDS) {
    const slider = document.getElementById(sliderId);
    roeKm[key] = parseFloat(slider.value, 10);
  }
  return roeKm;
}

function readRoeFromSliders() {
  const a_km = readSemiMajorAxisKm();
  return roeDimensionlessFromKm(readRoeKmFromSliders(), a_km);
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

function initSliders() {
  const { SLIDER_MIN_KM, SLIDER_MAX_KM, SLIDER_STEP_KM, ROE_DEFAULTS_KM } = CONFIG;

  for (const { key, sliderId, outputId } of ROE_FIELDS) {
    const slider = document.getElementById(sliderId);
    const output = document.getElementById(outputId);

    slider.min = String(SLIDER_MIN_KM);
    slider.max = String(SLIDER_MAX_KM);
    slider.step = String(SLIDER_STEP_KM);
    slider.value = String(ROE_DEFAULTS_KM[key]);
    output.textContent = formatKm(ROE_DEFAULTS_KM[key]);

    slider.addEventListener("input", () => {
      const km = parseFloat(slider.value, 10);
      output.textContent = formatKm(km);
      updateDriftRateDisplay(readSemiMajorAxisKm(), readRoeKmFromSliders());
      updateAllPlots();
    });
  }
}

/** Plotly.react はトレース数変更に弱い → 構造が変わったら newPlot */
function getPlotStructureKey() {
  const showRef = document.getElementById("show-single-orbit-ref").checked;
  return (
    `ref=${showRef ? 1 : 0}|orb=${readNumDriftOrbits()}|k=${readEciExaggerationFactor()}`
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
      "実スケール（k=1）。deputy は chief 軌道にほぼ重なって見えます。" +
      " 軸 [km] は表示座標＝実 ECI 座標です。";
    return;
  }

  el.textContent =
    `表示式: r = r_chief + k·(r_deputy−r_chief)、k=${exaggerationK}。` +
    " 軸 [km] は表示座標（chief・地球は実スケール、deputy の相対ずれのみ ×k）。" +
    " 実相対距離は下のパネルを参照。";
}

/** 軌道サンプルから実相対距離 [km] の統計を算出 */
function computeRelativeDistanceStats(samples) {
  if (!samples || samples.length === 0) return null;

  let maxAbsR = 0;
  let maxAbsT = 0;
  let maxAbsN = 0;
  let maxNorm = 0;
  let minR = Infinity;
  let maxR = -Infinity;
  let minT = Infinity;
  let maxT = -Infinity;
  let minN = Infinity;
  let maxN = -Infinity;

  for (const sample of samples) {
    const { R_km, T_km, N_km } = sample;
    maxAbsR = Math.max(maxAbsR, Math.abs(R_km));
    maxAbsT = Math.max(maxAbsT, Math.abs(T_km));
    maxAbsN = Math.max(maxAbsN, Math.abs(N_km));
    maxNorm = Math.max(maxNorm, Math.hypot(R_km, T_km, N_km));
    minR = Math.min(minR, R_km);
    maxR = Math.max(maxR, R_km);
    minT = Math.min(minT, T_km);
    maxT = Math.max(maxT, T_km);
    minN = Math.min(minN, N_km);
    maxN = Math.max(maxN, N_km);
  }

  return {
    maxAbsR,
    maxAbsT,
    maxAbsN,
    maxNorm,
    spanR: maxR - minR,
    spanT: maxT - minT,
    spanN: maxN - minN,
  };
}

function updateEciTrueDistanceDisplay(samples, exaggerationK) {
  const el = document.getElementById("eci-true-distance");
  if (!el) return;

  const stats = computeRelativeDistanceStats(samples);
  if (!stats) {
    el.textContent = "";
    return;
  }

  const k = exaggerationK;
  const lines = [
    "実相対距離（ROE 線形近似・RTN、表示軌道全体）",
    `  |Δr|max = ${formatKm(stats.maxNorm)} km`,
    `  R: 最大 |R| = ${formatKm(stats.maxAbsR)} km、範囲 ${formatKm(stats.spanR)} km`,
    `  T: 最大 |T| = ${formatKm(stats.maxAbsT)} km、範囲 ${formatKm(stats.spanT)} km`,
    `  N: 最大 |N| = ${formatKm(stats.maxAbsN)} km、範囲 ${formatKm(stats.spanN)} km`,
  ];

  if (k > 1) {
    lines.push(
      "",
      `表示上の deputy 相対オフセット（目安）≈ 実距離 × k = ${formatKm(stats.maxNorm * k)} km`,
      "（chief 半径・地球サイズは誇張されません）"
    );
  } else {
    lines.push("", "k=1: 表示座標＝実 ECI 座標（相対ずれは chief に重なって見えにくい場合があります）");
  }

  el.textContent = lines.join("\n");
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
  const a_km = readSemiMajorAxisKm();
  const roe = readRoeFromSliders();
  const { primary, reference, primarySamples, referenceSamples } =
    computeOrbitsForDisplay(a_km, roe);

  updateDriftRateDisplay(a_km, readRoeKmFromSliders());

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
  } = buildEciExaggeratedOrbitData(a_km, primarySamples, referenceSamples);
  updateEciCaption(exaggerationK);
  updateEciTrueDistanceDisplay(primarySamples, exaggerationK);

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
  initSliders();
  initOrbitControls();

  const semiMajorInput = document.getElementById("semi-major-axis");
  semiMajorInput.addEventListener("input", () => {
    updateDriftRateDisplay(readSemiMajorAxisKm(), readRoeKmFromSliders());
    updateAllPlots();
  });
  semiMajorInput.addEventListener("change", () => {
    updateDriftRateDisplay(readSemiMajorAxisKm(), readRoeKmFromSliders());
    updateAllPlots();
  });

  initViewTabs();

  window.addEventListener("resize", () => {
    if (plotsInitialized) {
      resizePlots(getPlotIdsForView());
    }
  });

  updateDriftRateDisplay(
    readSemiMajorAxisKm(),
    readRoeKmFromSliders()
  );
  updateEciCaption(readEciExaggerationFactor());
  updateAllPlots();
}

document.addEventListener("DOMContentLoaded", init);
