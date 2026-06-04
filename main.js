/**
 * ROE → RTN 相対軌道可視化
 * 近円軌道・小さい相対運動の線形近似（ブラウザ内計算のみ）
 */

// ---------------------------------------------------------------------------
// 設定定数（研究用途で後から変更しやすいよう集約）
// ---------------------------------------------------------------------------

const CONFIG = {
  /** スライダー範囲 [無次元] */
  SLIDER_MIN: -1.0e-3,
  SLIDER_MAX: 1.0e-3,
  SLIDER_STEP: 1.0e-5,

  /** 1 周あたりのサンプル点数 */
  NUM_ORBIT_POINTS: 400,

  /** 時間ドリフト表示のデフォルト周回数 */
  DEFAULT_NUM_DRIFT_ORBITS: 5,

  /** ECI 表示用 chief 軌道のサンプル点数 */
  ECI_CHIEF_ORBIT_POINTS: 360,

  /** deputy 付近ズーム時の軸余白比率 */
  ECI_ZOOM_PADDING_RATIO: 0.2,

  /** 地球重力定数 [km³/s²]（平均運動 n の計算用） */
  MU_EARTH_KM3_S2: 398600.4418,

  /** chief 半長軸の初期値 [km]（GEO 想定） */
  DEFAULT_SEMI_MAJOR_AXIS_KM: 42164,

  /** ROE 初期値 [無次元] */
  ROE_DEFAULTS: {
    delta_a: 0,
    delta_lambda: 1.0e-4,
    delta_ex: 1.0e-4,
    delta_ey: 0,
    delta_ix: 1.0e-4,
    delta_iy: 0,
  },
};

/** ROE キーと DOM id の対応 */
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
const PLOT_DIV_IDS = [
  "plot-3d",
  "plot-tr",
  "plot-tn",
  "plot-rn",
  "plot-eci-3d",
  "plot-eci-xy",
  "plot-eci-xz",
  "plot-eci-yz",
];

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

/** deputy 絶対位置 [km]（ECI） */
function deputyPositionEciKm(a_km, sample) {
  const chief = chiefPositionEciKm(a_km, sample.u_rad);
  const delta = relativeRtnToEciKm(
    sample.R_km,
    sample.T_km,
    sample.N_km,
    sample.u_rad
  );
  return [
    chief[0] + delta[0],
    chief[1] + delta[1],
    chief[2] + delta[2],
  ];
}

function samplesToEciOrbit(a_km, samples) {
  const x_km = [];
  const y_km = [];
  const z_km = [];

  for (const sample of samples) {
    const [x, y, z] = deputyPositionEciKm(a_km, sample);
    x_km.push(x);
    y_km.push(y);
    z_km.push(z);
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

function isEciFullScale() {
  const el = document.getElementById("eci-full-scale");
  return el ? el.checked : false;
}

/** ECI プロットの軸範囲 [km] */
function computeEciAxisRanges(a_km, eciOrbits, fullScale) {
  if (fullScale) {
    const pad = a_km * 0.08;
    const lim = a_km + pad;
    return {
      x: [-lim, lim],
      y: [-lim, lim],
      z: [-lim, lim],
    };
  }

  const coords = { x: [], y: [], z: [] };
  for (const orbit of eciOrbits) {
    if (!orbit) continue;
    coords.x.push(...orbit.x_km);
    coords.y.push(...orbit.y_km);
    coords.z.push(...orbit.z_km);
  }

  const padRatio = CONFIG.ECI_ZOOM_PADDING_RATIO;
  const axisRange = (arr) => {
    const min = Math.min(...arr);
    const max = Math.max(...arr);
    const span = Math.max(max - min, 0.5);
    const pad = span * padRatio;
    return [min - pad, max + pad];
  };

  return {
    x: axisRange(coords.x),
    y: axisRange(coords.y),
    z: axisRange(coords.z),
  };
}

/** 1 軌道周期あたりの along-track ドリフト量 [km] */
function alongTrackDriftPerOrbitKm(a_km, delta_a) {
  return -DRIFT_LAMBDA_FACTOR * TWO_PI * a_km * delta_a;
}

// ---------------------------------------------------------------------------
// UI ヘルパ
// ---------------------------------------------------------------------------

function formatRoeValue(value) {
  const abs = Math.abs(value);
  if (abs === 0) return "0";
  if (abs >= 0.01 || abs < 1e-6) return value.toExponential(2);
  return value.toExponential(4);
}

function readSemiMajorAxisKm() {
  const el = document.getElementById("semi-major-axis");
  const v = parseFloat(el.value, 10);
  return Number.isFinite(v) && v > 0 ? v : CONFIG.DEFAULT_SEMI_MAJOR_AXIS_KM;
}

function readRoeFromSliders() {
  const roe = {};
  for (const { key, sliderId } of ROE_FIELDS) {
    const slider = document.getElementById(sliderId);
    roe[key] = parseFloat(slider.value, 10);
  }
  return roe;
}

function isDriftMode() {
  const selected = document.querySelector('input[name="display-mode"]:checked');
  return selected && selected.value === "drift";
}

function readNumDriftOrbits() {
  const el = document.getElementById("num-drift-orbits");
  const v = parseInt(el.value, 10);
  return Number.isFinite(v) && v >= 1 ? v : CONFIG.DEFAULT_NUM_DRIFT_ORBITS;
}

function updateDriftRateDisplay(a_km, roe) {
  const el = document.getElementById("drift-rate-display");
  if (!el) return;

  if (Math.abs(roe.delta_a) < 1e-12) {
    el.textContent =
      "δa ≈ 0 のため along-track ドリフトはありません。\n" +
      "δa スライダーを 0 以外に動かしてください。";
    el.classList.add("drift-rate--warn");
    return;
  }

  el.classList.remove("drift-rate--warn");

  const n = chiefMeanMotionRadS(a_km);
  const dDeltaLambdaDt = -DRIFT_LAMBDA_FACTOR * n * roe.delta_a;
  const driftPerOrbitKm = alongTrackDriftPerOrbitKm(a_km, roe.delta_a);
  const periodHr = TWO_PI / n / 3600;

  el.textContent =
    `δλ̇ = ${formatRoeValue(dDeltaLambdaDt)} /s\n` +
    `1 周あたり ΔT ≈ ${driftPerOrbitKm.toFixed(3)} km\n` +
    `T_orbit ≈ ${periodHr.toFixed(2)} h`;
}

function updateModeUI() {
  const driftMode = isDriftMode();
  const driftControls = document.getElementById("drift-controls");
  const modeHint = document.getElementById("mode-hint");
  const caption = document.getElementById("plot-3d-caption");

  driftControls.hidden = !driftMode;

  if (driftMode) {
    modeHint.textContent =
      "δλ(t) = δλ₀ − (3/2)nδa·t を用い、複数周の開いた軌跡を表示";
    caption.textContent =
      "時間ドリフト: δa ≠ 0 で T 方向にセクラー変化（複数周分）";
    updateDriftRateDisplay(readSemiMajorAxisKm(), readRoeFromSliders());
  } else {
    modeHint.textContent =
      "固定 ROE で u を 1 周 → 閉じた相対軌道";
    caption.textContent =
      "横軸: T、縦軸: R、奥行き: N（chief は原点）";
  }
}

function initSliders() {
  const { SLIDER_MIN, SLIDER_MAX, SLIDER_STEP, ROE_DEFAULTS } = CONFIG;

  for (const { key, sliderId, outputId } of ROE_FIELDS) {
    const slider = document.getElementById(sliderId);
    const output = document.getElementById(outputId);

    slider.min = String(SLIDER_MIN);
    slider.max = String(SLIDER_MAX);
    slider.step = String(SLIDER_STEP);
    slider.value = String(ROE_DEFAULTS[key]);
    output.textContent = formatRoeValue(ROE_DEFAULTS[key]);

    slider.addEventListener("input", () => {
      output.textContent = formatRoeValue(parseFloat(slider.value, 10));
      updateModeUI();
      updateAllPlots();
    });
  }
}

/** Plotly.react はトレース数変更に弱い → 構造が変わったら newPlot */
function getPlotStructureKey() {
  const drift = isDriftMode();
  const showRef =
    drift && document.getElementById("show-single-orbit-ref").checked;
  return (
    `${drift ? "drift" : "single"}|ref=${showRef ? 1 : 0}|` +
    `orb=${drift ? readNumDriftOrbits() : 0}|eci=${isEciFullScale() ? 1 : 0}`
  );
}

function initModeControls() {
  document.querySelectorAll('input[name="display-mode"]').forEach((radio) => {
    const onModeChange = () => {
      updateModeUI();
      updateAllPlots(true);
    };
    radio.addEventListener("change", onModeChange);
    radio.addEventListener("click", onModeChange);
  });

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

  const eciScaleCheckbox = document.getElementById("eci-full-scale");
  if (eciScaleCheckbox) {
    eciScaleCheckbox.addEventListener("change", () => updateAllPlots(true));
  }

  updateModeUI();
}

// ---------------------------------------------------------------------------
// Plotly 描画
// ---------------------------------------------------------------------------

let plotsInitialized = false;
let lastPlotStructureKey = null;

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
      isDriftMode() ? "Deputy（ドリフト軌跡）" : "Deputy（相対軌道）",
      "#3fb950",
      4
    ),
    buildChiefTrace3D()
  );

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
      isDriftMode() ? "ドリフト軌跡" : "相対軌道",
      color,
      2
    ),
    buildChief2D()
  );

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
  const drift = isDriftMode();
  // ドリフト時は T 範囲が大きくなりやすい → T 軸を含む面は等尺を外してドリフトを見やすく
  const useEqualScale = !(drift && (planeId === "tr" || planeId === "tn"));

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
    showlegend: drift,
    legend: { font: { size: 10 } },
  };
}

function computeOrbitsForDisplay(a_km, roe) {
  const numPoints = CONFIG.NUM_ORBIT_POINTS;

  if (!isDriftMode()) {
    const primarySamples = buildSingleOrbitSamples(a_km, roe, numPoints);
    return {
      primary: samplesToRtnOrbit(primarySamples),
      reference: null,
      primarySamples,
      referenceSamples: null,
    };
  }

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

function buildEciTraces3D(primaryEci, referenceEci, chiefOrbitEci) {
  const traces = [
    {
      type: "scatter3d",
      mode: "markers",
      name: "地球（原点）",
      x: [0],
      y: [0],
      z: [0],
      marker: { color: "#6e9fef", size: 5, symbol: "circle" },
      hovertemplate: "地心 ECI 原点<extra></extra>",
    },
    buildEciLineTrace3D(
      chiefOrbitEci,
      "Chief 軌道",
      "rgba(240, 193, 75, 0.65)",
      2,
      "dot"
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
    buildEciLineTrace3D(
      primaryEci,
      isDriftMode() ? "Deputy（ドリフト）" : "Deputy",
      "#3fb950",
      4,
      null
    )
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
      aspectmode: "cube",
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

function buildEci2DTraces(primaryEci, referenceEci, chiefOrbitEci, plane) {
  const { axisX, axisY, xTitle, yTitle } = plane;
  const traces = [
    buildEci2DTrace(
      chiefOrbitEci,
      axisX,
      axisY,
      "Chief 軌道",
      "rgba(240, 193, 75, 0.65)",
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

  traces.push(
    buildEci2DTrace(
      primaryEci,
      axisX,
      axisY,
      isDriftMode() ? "Deputy（ドリフト）" : "Deputy",
      plane.color,
      2,
      null
    )
  );

  const axisRanges = {
    x: computeEciAxisRangesForPlane(
      [primaryEci, referenceEci, chiefOrbitEci],
      axisX,
      isEciFullScale()
    ),
    y: computeEciAxisRangesForPlane(
      [primaryEci, referenceEci, chiefOrbitEci],
      axisY,
      isEciFullScale()
    ),
  };

  return { traces, layout: layoutEci2D(xTitle, yTitle, 260, axisRanges) };
}

function computeEciAxisRangesForPlane(orbits, axisKey, fullScale) {
  const a_km = readSemiMajorAxisKm();
  if (fullScale) {
    const lim = a_km * 1.08;
    return [-lim, lim];
  }

  const values = [];
  for (const orbit of orbits) {
    if (!orbit) continue;
    values.push(...orbit[axisKey]);
  }

  const min = Math.min(...values);
  const max = Math.max(...values);
  const span = Math.max(max - min, 0.5);
  const pad = span * CONFIG.ECI_ZOOM_PADDING_RATIO;
  return [min - pad, max + pad];
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

function updateEciCaption(fullScale) {
  const el = document.getElementById("plot-eci-caption");
  if (!el) return;

  el.textContent = fullScale
    ? "chief 全軌道スケール。deputy の相対運動は円に対して非常に小さく見えます。"
    : "deputy 軌跡付近に自動ズーム。δi 成分は Z 方向の振動として現れます。";
}

function drawAllPlots(plotBundle, forceNewPlot) {
  const plotOpts = { responsive: true, displayModeBar: true };
  const structureKey = getPlotStructureKey();
  const useNewPlot =
    forceNewPlot ||
    !plotsInitialized ||
    structureKey !== lastPlotStructureKey;

  lastPlotStructureKey = structureKey;

  const draw = useNewPlot ? Plotly.newPlot : Plotly.react;

  if (useNewPlot && plotsInitialized) {
    for (const id of PLOT_DIV_IDS) {
      Plotly.purge(id);
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

  draw("plot-3d", traces3d, layout3d, plotOpts);
  draw("plot-tr", tr.traces, tr.layout, plotOpts);
  draw("plot-tn", tn.traces, tn.layout, plotOpts);
  draw("plot-rn", rn.traces, rn.layout, plotOpts);
  draw("plot-eci-3d", tracesEci3d, layoutEci3d, plotOpts);
  draw("plot-eci-xy", eciXy.traces, eciXy.layout, plotOpts);
  draw("plot-eci-xz", eciXz.traces, eciXz.layout, plotOpts);
  draw("plot-eci-yz", eciYz.traces, eciYz.layout, plotOpts);
  plotsInitialized = true;
}

function updateAllPlots(forceNewPlot = false) {
  const a_km = readSemiMajorAxisKm();
  const roe = readRoeFromSliders();
  const { primary, reference, primarySamples, referenceSamples } =
    computeOrbitsForDisplay(a_km, roe);

  if (isDriftMode()) {
    updateDriftRateDisplay(a_km, roe);
  }

  const traces3d = buildOrbitTraces3D(primary, reference);
  const layout3d = layout3D();
  const tr = build2DTraces(primary, PLANE_TR, reference);
  const tn = build2DTraces(primary, PLANE_TN, reference);
  const rn = build2DTraces(primary, PLANE_RN, reference);

  const primaryEci = samplesToEciOrbit(a_km, primarySamples);
  const referenceEci = referenceSamples
    ? samplesToEciOrbit(a_km, referenceSamples)
    : null;
  const chiefOrbitEci = buildChiefOrbitEci(
    a_km,
    CONFIG.ECI_CHIEF_ORBIT_POINTS
  );

  const eciFullScale = isEciFullScale();
  updateEciCaption(eciFullScale);

  const axisRanges3d = computeEciAxisRanges(
    a_km,
    [primaryEci, referenceEci, chiefOrbitEci],
    eciFullScale
  );

  const plotBundle = {
    traces3d,
    layout3d,
    tr,
    tn,
    rn,
    tracesEci3d: buildEciTraces3D(primaryEci, referenceEci, chiefOrbitEci),
    layoutEci3d: layoutEci3D(axisRanges3d),
    eciXy: buildEci2DTraces(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      PLANE_ECI_XY
    ),
    eciXz: buildEci2DTraces(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      PLANE_ECI_XZ
    ),
    eciYz: buildEci2DTraces(
      primaryEci,
      referenceEci,
      chiefOrbitEci,
      PLANE_ECI_YZ
    ),
  };

  drawAllPlots(plotBundle, forceNewPlot);
}

// ---------------------------------------------------------------------------
// 初期化
// ---------------------------------------------------------------------------

function init() {
  initSliders();
  initModeControls();

  const semiMajorInput = document.getElementById("semi-major-axis");
  semiMajorInput.addEventListener("input", () => {
    updateModeUI();
    updateAllPlots();
  });
  semiMajorInput.addEventListener("change", () => {
    updateModeUI();
    updateAllPlots();
  });

  window.addEventListener("resize", () => {
    if (plotsInitialized) {
      for (const id of PLOT_DIV_IDS) {
        Plotly.Plots.resize(id);
      }
    }
  });

  updateAllPlots();
}

document.addEventListener("DOMContentLoaded", init);
