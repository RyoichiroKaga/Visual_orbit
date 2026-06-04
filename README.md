# ROE → RTN 相対軌道可視化

Relative Orbital Elements（ROE）の各パラメータをスライダーで調整しながら、chief 衛星に対する deputy 衛星の相対軌道を **RTN 座標系**でリアルタイムに確認するための、教育・研究補助用 Web デモです。

厳密な軌道伝播ツールではなく、「ROE の各成分が相対軌道の形状にどう効くか」を直感的に理解するための可視化に特化しています。

## この Web アプリの目的

- 近円軌道・小さい相対運動を仮定した **線形近似**により、ROE から RTN 相対位置を即座に計算する
- スライダー操作と連動した **RTN** および **ECI（慣性系）** の 3D / 2D プロットで相対・絶対軌道を確認する
- 研究室のデモや発表資料向けに、ブラウザだけで動く **静的サイト**として配布する

## ROE（Relative Orbital Elements）とは

ROE は、chief 衛星の軌道に対する deputy 衛星の相対運動を、6 つの無次元パラメータで表現する記述です。本アプリのスライダーは **km スケール**（\(a \cdot \delta\) に相当）で操作し、内部で \(\delta = \text{値}/a\) に変換します。

| UI [km] | 対応する δ | 主な RTN への効き |
|---------|-----------|------------------|
| \(a \cdot \delta a\) | 半長軸差 | R オフセット、along-track ドリフト |
| \(a \cdot \delta\lambda\) | 平均経度差 | T オフセット |
| \(a \cdot \delta e_x,\ \delta e_y\) | 離心率ベクトル | R–T 平面の楕円 |
| \(a \cdot \delta i_x,\ \delta i_y\) | 傾斜角ベクトル | N 方向の振動 |

## RTN 座標系

chief 衛星を基準とした局所直交座標系です。

- **R（Radial）** — chief から地心方向（半径方向）
- **T（Along-track）** — 速度方向に近い軌道面内の接線方向
- **N（Cross-track）** — 軌道角運動量に沿った法線方向（軌道面の外側）

本アプリの RTN 3D プロットでは、横軸 **T**、縦軸 **R**、奥行き **N** [km] とし、chief を原点に表示します。

### ECI（Earth-Centered Inertial）— 誇張表示

GEO スケールでは実際の相対変位（km オーダ）が chief 軌道（約 4 万 km）に対して小さすぎるため、**表示専用の誇張倍率** \(k\) を用います。

\[
\mathbf{r}_{\mathrm{display}}^{\mathrm{ECI}}
= \mathbf{r}_{\mathrm{chief}}^{\mathrm{ECI}}
+ k\left(\mathbf{r}_{\mathrm{deputy}}^{\mathrm{ECI}} - \mathbf{r}_{\mathrm{chief}}^{\mathrm{ECI}}\right)
\]

- \(k = 1\) … 実スケール（deputy は chief にほぼ重なる）
- \(k = 10^3\) … デフォルト（デモ向け）
- 相対形状の厳密な km 値は **RTN プロット**を参照
- **地球** … 赤道面（\(z=0\)）上に \(R_\oplus \approx 6378\) km の円（ワイヤーフレーム球の赤道断面、実スケール）

## 使用している近似式

近円軌道・小さい相対運動を仮定した線形モデルです。引数緯度 \(u\) [rad] に対し：

\[
\begin{aligned}
R &= a\,(\delta a - \delta e_x \cos u - \delta e_y \sin u) \\
T &= a\,(\delta\lambda + 2\delta e_x \sin u - 2\delta e_y \cos u) \\
N &= a\,(\delta i_x \sin u - \delta i_y \cos u)
\end{aligned}
\]

- \(a\) … chief 軌道の半長軸 [km]
- \(u\) … 各周で 0 から \(2\pi\)、複数周分を連続表示

### 時間ドリフト（常時適用）

\[
\delta\lambda(t) = \delta\lambda_0 - \frac{3}{2}\, n\, \delta a\, t, \quad
t = k\frac{2\pi}{n} + \frac{u}{n}
\]

- \(n=\sqrt{\mu/a^3}\) [rad/s]
- **Δa = 0 km** のときは各周が同じ閉曲線に重なる（周回数 1 推奨）
- **Δa ≠ 0** のとき along-track 方向へ周回ごとにずれる

## ファイル構成

```
.
├── index.html   # ページ構造・UI
├── style.css    # レイアウト・ダークテーマ
├── main.js      # ROE→RTN 計算・Plotly 描画
└── README.md
```

## ローカルでの表示

リポジトリのルートで簡易 HTTP サーバを起動し、ブラウザで `index.html` を開きます。

```bash
# Python 3
python3 -m http.server 8080

# または Node.js (npx)
npx --yes serve -l 8080
```

ブラウザで `http://localhost:8080` にアクセスしてください。

`file://` で直接開くことも可能ですが、CDN（Plotly.js）の読み込み環境によっては HTTP サーバの利用を推奨します。

## GitHub Pages での公開方法

1. このプロジェクトを GitHub リポジトリに push する
2. リポジトリの **Settings → Pages**
3. **Build and deployment → Source** で **Deploy from a branch** を選択
4. **Branch** を `main`（または使用ブランチ）、フォルダを **`/ (root)`** に設定して Save
5. 数分後、`https://<ユーザー名>.github.io/<リポジトリ名>/` で公開される

ルート以外のサブディレクトリに置く場合は、そのフォルダを Pages のソースに指定するか、リポジトリ構成に合わせてパスを調整してください。

## 注意事項

- **近円軌道**を仮定した式であり、大きな離心率や強い相対運動には適用できません
- **線形化**されたモデルのため、長期伝播や非線形効果は表現しません（ドリフトモードも \(\dot{\delta\lambda}\propto\delta a\) の一次近似のみ）
- **摂動**（J2、第三体、SRP など）は考慮していません
- 単位は UI 上で明示していますが、数値はデモ用のスケール感です。実ミッション設計には専用の伝播・制御ツールを使用してください

## 技術スタック

- HTML / CSS / JavaScript（フレームワークなし）
- [Plotly.js](https://plotly.com/javascript/)（CDN）
- 計算はすべてブラウザ内で完結（サーバー不要）

## ライセンス

研究・教育用途での利用を想定しています。必要に応じてリポジトリにライセンスを追加してください。
