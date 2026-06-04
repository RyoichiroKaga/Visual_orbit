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

/**
 * u = 0 … 2π で 1 周分（δλ 固定の閉曲線）
 */
function computeRelativeOrbitRTN(a_km, roe, numPoints) {
  const R_km = [];
  const T_km = [];
  const N_km = [];

  for (let i = 0; i < numPoints; i++) {
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    const pos = rtnAtLatitude(a_km, roe, u_rad);
    R_km.push(pos.R_km);
    T_km.push(pos.T_km);
    N_km.push(pos.N_km);
  }

  return { R_km, T_km, N_km };
}

/**
 * 複数周回の時間ドリフト軌跡
 * 各点で t = k·(2π/n) + u/n、δλ(t) を用いて T 方向ドリフトを表現
 */
function computeDriftTrajectoryRTN(a_km, roe, numPointsPerOrbit, numOrbits) {
  const n_rad_s = chiefMeanMotionRadS(a_km);
  const orbitPeriodSec = TWO_PI / n_rad_s;
  const R_km = [];
  const T_km = [];
  const N_km = [];

  for (let k = 0; k < numOrbits; k++) {
    for (let i = 0; i < numPointsPerOrbit; i++) {
      // 最後の周の終端以外で点を重複させない
      if (k > 0 && i === 0) continue;

      const u_rad = (TWO_PI * i) / (numPointsPerOrbit - 1);
      const t_sec = k * orbitPeriodSec + u_rad / n_rad_s;
      const delta_lambda_t = deltaLambdaAtTimeSec(roe, n_rad_s, t_sec);
      const pos = rtnAtLatitude(a_km, roe, u_rad, delta_lambda_t);
      R_km.push(pos.R_km);
      T_km.push(pos.T_km);
      N_km.push(pos.N_km);
    }
  }

  return { R_km, T_km, N_km, n_rad_s, orbitPeriodSec };
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

  const n = chiefMeanMotionRadS(a_km);
  const dDeltaLambdaDt = -DRIFT_LAMBDA_FACTOR * n * roe.delta_a;
  const driftPerOrbitKm = alongTrackDriftPerOrbitKm(a_km, roe.delta_a);
  const periodHr = (TWO_PI / n) / 3600;

  el.textContent =
    `δλ̇ = ${formatRoeValue(dDeltaLambdaDt)} /s\n` +
    `1 周あたり ΔT ≈ ${driftPerOrbitKm.toFixed(3)} km\n` +
    `T_orbit ≈ ${(periodHr).toFixed(2)} h`;
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

function initModeControls() {
  document.querySelectorAll('input[name="display-mode"]').forEach((radio) => {
    radio.addEventListener("change", () => {
      updateModeUI();
      updateAllPlots();
    });
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
    .addEventListener("change", updateAllPlots);

  updateModeUI();
}

// ---------------------------------------------------------------------------
// Plotly 描画
// ---------------------------------------------------------------------------

let plotsInitialized = false;

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

  return { traces, layout: layout2D(xTitle, yTitle, 260) };
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

function layout2D(xTitle, yTitle, height) {
  return {
    ...PLOT_LAYOUT_BASE,
    height,
    xaxis: {
      title: xTitle,
      gridcolor: "#30363d",
      zerolinecolor: "#484f58",
      scaleanchor: "y",
      scaleratio: 1,
    },
    yaxis: {
      title: yTitle,
      gridcolor: "#30363d",
      zerolinecolor: "#484f58",
    },
    showlegend: isDriftMode(),
    legend: { font: { size: 10 } },
  };
}

function computeOrbitsForDisplay(a_km, roe) {
  const numPoints = CONFIG.NUM_ORBIT_POINTS;

  if (!isDriftMode()) {
    const orbit = computeRelativeOrbitRTN(a_km, roe, numPoints);
    return { primary: orbit, reference: null };
  }

  const numOrbits = readNumDriftOrbits();
  const primary = computeDriftTrajectoryRTN(
    a_km,
    roe,
    numPoints,
    numOrbits
  );

  const showRef = document.getElementById("show-single-orbit-ref").checked;
  const reference = showRef
    ? computeRelativeOrbitRTN(a_km, roe, numPoints)
    : null;

  return { primary, reference };
}

function updateAllPlots() {
  const a_km = readSemiMajorAxisKm();
  const roe = readRoeFromSliders();
  const { primary, reference } = computeOrbitsForDisplay(a_km, roe);

  if (isDriftMode()) {
    updateDriftRateDisplay(a_km, roe);
  }

  const traces3d = buildOrbitTraces3D(primary, reference);
  const layout3d = layout3D();
  const tr = build2DTraces(primary, PLANE_TR, reference);
  const tn = build2DTraces(primary, PLANE_TN, reference);
  const rn = build2DTraces(primary, PLANE_RN, reference);
  const plotOpts = { responsive: true, displayModeBar: true };

  if (!plotsInitialized) {
    Plotly.newPlot("plot-3d", traces3d, layout3d, plotOpts);
    Plotly.newPlot("plot-tr", tr.traces, tr.layout, plotOpts);
    Plotly.newPlot("plot-tn", tn.traces, tn.layout, plotOpts);
    Plotly.newPlot("plot-rn", rn.traces, rn.layout, plotOpts);
    plotsInitialized = true;
  } else {
    Plotly.react("plot-3d", traces3d, layout3d);
    Plotly.react("plot-tr", tr.traces, tr.layout);
    Plotly.react("plot-tn", tn.traces, tn.layout);
    Plotly.react("plot-rn", rn.traces, rn.layout);
  }
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
      Plotly.Plots.resize("plot-3d");
      Plotly.Plots.resize("plot-tr");
      Plotly.Plots.resize("plot-tn");
      Plotly.Plots.resize("plot-rn");
    }
  });

  updateAllPlots();
}

document.addEventListener("DOMContentLoaded", init);
