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

  /** 1 周分のサンプル点数 */
  NUM_ORBIT_POINTS: 400,

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

// Plotly 用の共通スタイル（ダークテーマに合わせる）
const PLOT_LAYOUT_BASE = {
  paper_bgcolor: "#1c2333",
  plot_bgcolor: "#1c2333",
  font: { color: "#e6edf3", family: "Segoe UI, sans-serif", size: 12 },
  margin: { l: 50, r: 20, t: 30, b: 45 },
};

// ---------------------------------------------------------------------------
// ROE → RTN 変換
// ---------------------------------------------------------------------------

/**
 * 引数緯度 u [rad] における RTN 相対位置 [km] を計算
 * @param {number} a_km - chief 半長軸 [km]
 * @param {object} roe - 無次元 ROE
 * @param {number} u_rad - argument of latitude [rad]
 */
function rtnAtLatitude(a_km, roe, u_rad) {
  const { delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy } = roe;
  const cu = Math.cos(u_rad);
  const su = Math.sin(u_rad);

  const R_km =
    a_km * (delta_a - delta_ex * cu - delta_ey * su);
  const T_km =
    a_km * (delta_lambda + 2 * delta_ex * su - 2 * delta_ey * cu);
  const N_km = a_km * (delta_ix * su - delta_iy * cu);

  return { R_km, T_km, N_km };
}

/**
 * u = 0 … 2π で 1 周分の相対軌道をサンプリング
 * @returns {{ R_km: number[], T_km: number[], N_km: number[] }}
 */
function computeRelativeOrbitRTN(a_km, roe, numPoints) {
  const R_km = [];
  const T_km = [];
  const N_km = [];

  for (let i = 0; i < numPoints; i++) {
    // 端点を含めて閉じた軌道にする（i=0 → u=0, i=n-1 → u=2π）
    const u_rad = (TWO_PI * i) / (numPoints - 1);
    const pos = rtnAtLatitude(a_km, roe, u_rad);
    R_km.push(pos.R_km);
    T_km.push(pos.T_km);
    N_km.push(pos.N_km);
  }

  return { R_km, T_km, N_km };
}

// ---------------------------------------------------------------------------
// UI ヘルパ
// ---------------------------------------------------------------------------

/** 科学記数法風の表示（小さい ROE 向け） */
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
      updateAllPlots();
    });
  }
}

// ---------------------------------------------------------------------------
// Plotly 描画
// ---------------------------------------------------------------------------

let plotsInitialized = false;

function buildOrbitTraces(orbit) {
  const { T_km, R_km, N_km } = orbit;

  const deputyTrace = {
    type: "scatter3d",
    mode: "lines",
    name: "Deputy（相対軌道）",
    x: T_km,
    y: R_km,
    z: N_km,
    line: { color: "#3fb950", width: 4 },
    hovertemplate:
      "T: %{x:.3f} km<br>R: %{y:.3f} km<br>N: %{z:.3f} km<extra></extra>",
  };

  const chiefTrace = {
    type: "scatter3d",
    mode: "markers",
    name: "Chief（原点）",
    x: [0],
    y: [0],
    z: [0],
    marker: { color: "#f0c14b", size: 6, symbol: "diamond" },
    hovertemplate: "Chief @ 原点<extra></extra>",
  };

  return [deputyTrace, chiefTrace];
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
      camera: {
        eye: { x: 1.4, y: 1.2, z: 0.9 },
      },
    },
  };
}

function build2DTrace(x, y, name, color) {
  return {
    type: "scatter",
    mode: "lines",
    name,
    x,
    y,
    line: { color, width: 2 },
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
    showlegend: false,
  };
}

function updateAllPlots() {
  const a_km = readSemiMajorAxisKm();
  const roe = readRoeFromSliders();
  const orbit = computeRelativeOrbitRTN(
    a_km,
    roe,
    CONFIG.NUM_ORBIT_POINTS
  );

  const traces3d = buildOrbitTraces(orbit);
  const layout3d = layout3D();

  const plotOpts = { responsive: true, displayModeBar: true };

  if (!plotsInitialized) {
    Plotly.newPlot("plot-3d", traces3d, layout3d, plotOpts);

    Plotly.newPlot(
      "plot-tr",
      [
        build2DTrace(orbit.T_km, orbit.R_km, "相対軌道", "#3fb950"),
        buildChief2D(),
      ],
      layout2D("T [km]", "R [km]", 260),
      plotOpts
    );

    Plotly.newPlot(
      "plot-tn",
      [
        build2DTrace(orbit.T_km, orbit.N_km, "相対軌道", "#58a6ff"),
        buildChief2D(),
      ],
      layout2D("T [km]", "N [km]", 260),
      plotOpts
    );

    Plotly.newPlot(
      "plot-rn",
      [
        build2DTrace(orbit.R_km, orbit.N_km, "相対軌道", "#a371f7"),
        buildChief2D(),
      ],
      layout2D("R [km]", "N [km]", 260),
      plotOpts
    );

    plotsInitialized = true;
  } else {
    Plotly.react("plot-3d", traces3d, layout3d);
    Plotly.react(
      "plot-tr",
      [
        build2DTrace(orbit.T_km, orbit.R_km, "相対軌道", "#3fb950"),
        buildChief2D(),
      ],
      layout2D("T [km]", "R [km]", 260)
    );
    Plotly.react(
      "plot-tn",
      [
        build2DTrace(orbit.T_km, orbit.N_km, "相対軌道", "#58a6ff"),
        buildChief2D(),
      ],
      layout2D("T [km]", "N [km]", 260)
    );
    Plotly.react(
      "plot-rn",
      [
        build2DTrace(orbit.R_km, orbit.N_km, "相対軌道", "#a371f7"),
        buildChief2D(),
      ],
      layout2D("R [km]", "N [km]", 260)
    );
  }
}

// ---------------------------------------------------------------------------
// 初期化
// ---------------------------------------------------------------------------

function init() {
  initSliders();

  const semiMajorInput = document.getElementById("semi-major-axis");
  semiMajorInput.addEventListener("input", updateAllPlots);
  semiMajorInput.addEventListener("change", updateAllPlots);

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
