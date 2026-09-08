const { useState, useRef, useEffect, useCallback, useMemo } = React;

/* ============================================================
   跟驰稳定性实验台 — 环道微观仿真 × 解析判据
   场景：S0 基准 IDM / S1 前车加速度前馈 / S5 时间延迟
   ============================================================ */

/* ---------- 复数工具 ---------- */
const C = (re, im = 0) => ({ re, im });
const cAdd = (a, b) => C(a.re + b.re, a.im + b.im);
const cSub = (a, b) => C(a.re - b.re, a.im - b.im);
const cMul = (a, b) => C(a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re);
const cScale = (a, k) => C(a.re * k, a.im * k);
const cDiv = (a, b) => {
  const d = b.re * b.re + b.im * b.im;
  return C((a.re * b.re + a.im * b.im) / d, (a.im * b.re - a.re * b.im) / d);
};
const cExp = (a) => {
  const e = Math.exp(a.re);
  return C(e * Math.cos(a.im), e * Math.sin(a.im));
};
const cAbs = (a) => Math.hypot(a.re, a.im);
const cSqrt = (a) => {
  const r = cAbs(a);
  const re = Math.sqrt(Math.max(0, (r + a.re) / 2));
  let im = Math.sqrt(Math.max(0, (r - a.re) / 2));
  if (a.im < 0) im = -im;
  return C(re, im);
};

/* ---------- IDM 平衡态与偏导 ---------- */
function idmEquilibriumSpeed(se, P) {
  // 解 se = (s0 + v T)/sqrt(1-(v/v0)^delta)
  const f = (v) => {
    const z2 = 1 - Math.pow(v / P.v0, P.delta);
    if (z2 <= 1e-12) return 1e9;
    return (P.s0 + v * P.T) / Math.sqrt(z2) - se;
  };
  if (f(0) > 0) return 0; // 间距比最小停车间距还小
  let lo = 0, hi = P.v0 * 0.999999;
  for (let i = 0; i < 80; i++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) < 0) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

function idmPartials(ve, P) {
  const z2 = 1 - Math.pow(ve / P.v0, P.delta);
  const z = Math.sqrt(Math.max(z2, 1e-12));
  const se = (P.s0 + ve * P.T) / z;
  const A = ve > 0 ? (P.a * P.delta / P.v0) * Math.pow(ve / P.v0, P.delta - 1) : 0;
  const B = (2 * P.a * z * P.T) / se;
  const Cc = Math.sqrt(P.a / P.b) * (z * ve) / se;
  return {
    fs: (2 * P.a * z * z) / se,
    fv: -(A + B) - Cc,
    fvl: Cc,
    se, z, A, B, C: Cc,
  };
}

/* ---------- 解析判据 ---------- */
// Psi(omega) >= 0  <=>  |G(i omega)| <= 1   （单向 + 前馈 kappa + 延迟）
function psiOmega(w, p) {
  const { fs, fv, fvl } = p.part;
  const kappa = p.fa ?? p.kappa ?? 0;
  const { tau0, tauA } = p;
  const ph = w * (tauA - tau0);
  return (
    (1 - kappa * kappa) * w * w +
    fv * fv - fvl * fvl -
    2 * fs * Math.cos(w * tau0) +
    2 * fv * w * Math.sin(w * tau0) +
    2 * kappa * fs * Math.cos(ph) -
    2 * kappa * fvl * w * Math.sin(ph)
  );
}

// 任意有限加速度核：前方用 e^{-ijk}，后方用 e^{+ijk}
function accKernel(k, p) {
  let out = C(0, 0);
  const mf = Math.max(1, p.mFront || 1), mb = Math.max(1, p.mBack || 1);
  for (let j = 1; j <= mf; j++) {
    const w = (p.fa || 0) * Math.pow(p.alpha || 0, j - 1);
    out = cAdd(out, C(w * Math.cos(j * k), -w * Math.sin(j * k)));
  }
  for (let j = 1; j <= mb; j++) {
    const w = (p.fb || 0) * Math.pow(p.beta || 0, j - 1);
    out = cAdd(out, C(w * Math.cos(j * k), w * Math.sin(j * k)));
  }
  return out;
}

function qrAtK(k, p) {
  const { fs, fv, fvl } = p.part;
  const emik = C(Math.cos(k), -Math.sin(k));
  const epik = C(Math.cos(k), Math.sin(k));
  const hs = (p.rearP || 0) * fs;
  const hv = (p.rearP || 0) * fvl;
  return {
    Q: cAdd(C(fv - hv, 0), cAdd(cScale(emik, fvl), cScale(epik, hv))),
    R: cAdd(cScale(cSub(emik, C(1, 0)), fs), cScale(cSub(epik, C(1, 0)), hs)),
  };
}

// 环道模态特征方程 F(lambda) = 0
function charF(lam, k, p) {
  const { tau0, tauA } = p;
  const eA = cExp(cScale(lam, -tauA));
  const e0 = cExp(cScale(lam, -tau0));
  const lam2 = cMul(lam, lam);
  const Pk = cSub(C(1, 0), cMul(eA, accKernel(k, p)));
  const { Q: Qk, R: Rk } = qrAtK(k, p);
  // F = lam^2 P - e^{-lam tau0}(Q lam + R)
  return cSub(cMul(lam2, Pk), cMul(e0, cAdd(cMul(Qk, lam), Rk)));
}

function rootsNoDelay(k, p) {
  const Pk = cSub(C(1, 0), accKernel(k, p));
  const { Q: Qk, R: Rk } = qrAtK(k, p);
  const disc = cSqrt(cAdd(cMul(Qk, Qk), cScale(cMul(Pk, Rk), 4)));
  const twoP = cScale(Pk, 2);
  return [cDiv(cAdd(Qk, disc), twoP), cDiv(cSub(Qk, disc), twoP)];
}

function newtonRoot(lam0, k, p, iters = 40) {
  let lam = lam0;
  for (let i = 0; i < iters; i++) {
    const F = charF(lam, k, p);
    const h = 1e-6 * (1 + cAbs(lam));
    const dF = cScale(cSub(charF(C(lam.re + h, lam.im), k, p), F), 1 / h);
    if (cAbs(dF) < 1e-14) break;
    const step = cDiv(F, dF);
    lam = cSub(lam, step);
    if (cAbs(step) < 1e-13) break;
  }
  return lam;
}

// 用 tau 连续化求最右根（已与暴力多起点搜索核对一致）
function rightmostAtK(k, p) {
  const base = { ...p, tau0: 0, tauA: 0 };
  let best = null;
  for (const r0 of rootsNoDelay(k, base)) {
    let lam = r0;
    const NS = p.tau0 > 0 || p.tauA > 0 ? 6 : 0;
    for (let s = 1; s <= NS; s++) {
      const frac = s / NS;
      lam = newtonRoot(lam, k, { ...p, tau0: p.tau0 * frac, tauA: p.tauA * frac }, 18);
    }
    if (!best || lam.re > best.re) best = lam;
  }
  return best;
}

function ringSpectrum(p, N) {
  const out = [];
  for (let m = 1; m < N; m++) {
    const k = (2 * Math.PI * m) / N;
    if (k > Math.PI + 1e-9) continue;
    out.push({ m, k, lam: rightmostAtK(k, p) });
  }
  return out;
}

function kernelTotal(p) {
  let s = 0;
  for (let j = 1; j <= Math.max(1, p.mFront || 1); j++) s += (p.fa || 0) * Math.pow(p.alpha || 0, j - 1);
  for (let j = 1; j <= Math.max(1, p.mBack || 1); j++) s += (p.fb || 0) * Math.pow(p.beta || 0, j - 1);
  return s;
}

function longWaveMargin(part, p) {
  const { fs, fv, fvl } = part;
  if ((p.rearP || 0) > 0 && Math.abs(kernelTotal(p)) < 1e-12) {
    const hs = p.rearP * fs, hv = p.rearP * fvl, mu = fv + fvl;
    return -((fs - hs) * (fs - hs) + mu * (fvl - hv) * (fs - hs) - 0.5 * mu * mu * (fs + hs));
  }
  return 0.5 * (fv * fv - fvl * fvl) - (1 - kernelTotal(p)) * fs;
}

/* ---------- IDM 加速度 ---------- */
function idmAcc(v, s, dv, P) {
  const sStar = P.s0 + Math.max(0, v * P.T + (v * dv) / (2 * Math.sqrt(P.a * P.b)));
  const sSafe = Math.max(s, 0.05);
  return P.a * (1 - Math.pow(Math.max(v, 0) / P.v0, P.delta) - Math.pow(sStar / sSafe, 2));
}

function idmGapAtSpeed(v, P) {
  const z2 = Math.max(1e-10, 1 - Math.pow(v / P.v0, P.delta));
  return (P.s0 + v * P.T) / Math.sqrt(z2);
}

function commonEquilibriumSpeed(meanGap, models) {
  const avgGap = (v) => models.reduce((q, P) => q + idmGapAtSpeed(v, P), 0) / models.length;
  if (avgGap(0) > meanGap) return 0;
  let lo = 0, hi = Math.min(...models.map((P) => P.v0)) * 0.999;
  for (let i = 0; i < 70; i++) {
    const mid = 0.5 * (lo + hi);
    if (avgGap(mid) < meanGap) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

function makeArrangement(N, count, mode, groups = 4, seed = 17) {
  const out = Array(N).fill(false);
  count = Math.max(0, Math.min(N, Math.round(count)));
  if (!count) return out;
  if (mode === "成组") {
    const ng = Math.max(1, Math.min(count, Math.round(groups)));
    let placed = 0;
    for (let g = 0; g < ng; g++) {
      const size = Math.floor(count / ng) + (g < count % ng ? 1 : 0);
      const start = Math.floor((g * N) / ng);
      for (let j = 0; j < size; j++) out[(start + j) % N] = true;
      placed += size;
    }
    return out;
  }
  if (mode === "随机") {
    const ids = Array.from({ length: N }, (_, i) => i);
    let x = seed >>> 0;
    for (let i = N - 1; i > 0; i--) {
      x = (1664525 * x + 1013904223) >>> 0;
      const j = x % (i + 1);
      [ids[i], ids[j]] = [ids[j], ids[i]];
    }
    for (let i = 0; i < count; i++) out[ids[i]] = true;
    return out;
  }
  if (mode === "连续") {
    for (let i = 0; i < count; i++) out[i] = true;
    return out;
  }
  for (let j = 0; j < count; j++) out[Math.floor(((j + 0.5) * N) / count) % N] = true;
  return out;
}

/* ---------- 仿真核心 ---------- */
function createSim(cfg) {
  const {
    N, L, l, P, tau0, tauA, dt, pert,
    fa = 0, fb = 0, alpha = 0, beta = 0, mFront = 1, mBack = 1,
    rearP = 0, hetero = false, mixed = false, heavyShare = 0.25,
    avShare = 0.3, arrangement = "均匀", avGroups = 4, kappa0 = 0.15,
  } = cfg;
  const se = L / N - l;
  const heavy = hetero ? makeArrangement(N, N * heavyShare, arrangement, avGroups, 29) : Array(N).fill(false);
  const av = mixed ? makeArrangement(N, N * avShare, arrangement, avGroups, 71) : Array(N).fill(false);
  const models = heavy.map((on) => on ? ({ ...P, v0: P.v0 * 0.78, T: P.T * 1.32, a: P.a * 0.62, b: P.b * 0.82 }) : P);
  const ve = hetero ? commonEquilibriumSpeed(se, models) : idmEquilibriumSpeed(se, P);
  const eqGaps = models.map((Pi) => idmGapAtSpeed(ve, Pi));
  const parts = models.map((Pi) => idmPartials(ve, Pi));
  const maxTau = Math.max(tau0, tauA);
  const H = Math.max(2, Math.ceil(maxTau / dt) + 3);

  const x = new Float64Array(N);
  const v = new Float64Array(N);
  const acc = new Float64Array(N);
  const hx = new Float64Array(H * N);
  const hv = new Float64Array(H * N);
  const ha = new Float64Array(H * N);

  for (let n = 0; n < N; n++) {
    x[n] = n === 0 ? 0 : x[n - 1] + l + eqGaps[n - 1];
    v[n] = ve;
    acc[n] = 0;
  }
  // 扰动：单车速度脉冲
  v[0] = ve + pert;

  for (let h = 0; h < H; h++)
    for (let n = 0; n < N; n++) {
      hx[h * N + n] = x[n];
      hv[h * N + n] = v[n];
      ha[h * N + n] = 0;
    }

  return {
    N, L, l, P, models, parts, heavy, av, fa, fb, alpha, beta, mFront, mBack,
    rearP, mixed, kappa0, tau0, tauA, dt, se, ve,
    x, v, acc, hx, hv, ha, H, head: 0, t: 0,
    lead(n) { return (n + 1) % N; },
    // 从历史取值（线性插值）
    delayed(arr, tau, n) {
      if (tau <= 0) return null;
      const back = tau / this.dt;
      const i0 = Math.floor(back), fr = back - i0;
      const a = (this.head - i0 + this.H * 4) % this.H;
      const b = (this.head - i0 - 1 + this.H * 4) % this.H;
      return arr[a * N + n] * (1 - fr) + arr[b * N + n] * fr;
    },
    step() {
      const { N, L, l, tau0, tauA, dt } = this;
      const useD = tau0 > 0;
      const raw = new Float64Array(N);
      for (let n = 0; n < N; n++) {
        const ld = (n + 1) % N;
        let xn, vn, xl, vl;
        if (useD) {
          xn = this.delayed(this.hx, tau0, n);
          vn = this.delayed(this.hv, tau0, n);
          xl = this.delayed(this.hx, tau0, ld);
          vl = this.delayed(this.hv, tau0, ld);
        } else {
          xn = this.x[n]; vn = this.v[n]; xl = this.x[ld]; vl = this.v[ld];
        }
        if (ld === 0) xl += L;
        const s = xl - xn - l;
        let A = idmAcc(vn, s, vn - vl, this.models[n]);
        if (this.rearP > 0) {
          const rd = (n - 1 + N) % N;
          let xr, vr, xs = xn;
          if (useD) {
            xr = this.delayed(this.hx, tau0, rd);
            vr = this.delayed(this.hv, tau0, rd);
          } else {
            xr = this.x[rd]; vr = this.v[rd];
          }
          if (n === 0) xs += L;
          const sr = xs - xr - l;
          A += this.rearP * this.parts[n].fs * (sr - eqGaps[n]) + this.rearP * this.parts[n].fvl * (vr - vn);
        }
        raw[n] = Math.max(-12, Math.min(6, A));
      }
      // a = raw + K a：固定点迭代到机器精度；有通信延迟的项直接由历史缓冲给出。
      let guess = Float64Array.from(this.acc), next = new Float64Array(N);
      for (let it = 0; it < 18; it++) {
        let err = 0;
        for (let n = 0; n < N; n++) {
          let A = raw[n];
          for (let j = 1; j <= Math.max(1, this.mFront); j++) {
            const q = (n + j) % N;
            let gain = this.fa * Math.pow(this.alpha, j - 1);
            if (this.mixed) gain = this.av[n] ? (this.av[q] ? gain : Math.min(gain, this.kappa0 * Math.pow(this.alpha, j - 1))) : 0;
            if (gain) A += gain * (tauA > 0 ? this.delayed(this.ha, tauA, q) : guess[q]);
          }
          for (let j = 1; j <= Math.max(1, this.mBack); j++) {
            const q = (n - j + N) % N;
            const gain = this.fb * Math.pow(this.beta, j - 1);
            if (gain) A += gain * (tauA > 0 ? this.delayed(this.ha, tauA, q) : guess[q]);
          }
          next[n] = Math.max(-20, Math.min(10, A));
          err = Math.max(err, Math.abs(next[n] - guess[n]));
        }
        [guess, next] = [next, guess];
        if (err < 1e-9) break;
      }
      this.acc.set(guess);
      for (let n = 0; n < N; n++) {
        const A = this.acc[n];
        let nv = this.v[n] + A * dt;
        if (nv < 0) nv = 0;
        this.x[n] += this.v[n] * dt + 0.5 * A * dt * dt;
        this.v[n] = nv;
      }
      this.head = (this.head + 1) % this.H;
      const h = this.head;
      for (let n = 0; n < N; n++) {
        this.hx[h * N + n] = this.x[n];
        this.hv[h * N + n] = this.v[n];
        this.ha[h * N + n] = this.acc[n];
      }
      this.t += dt;
    },
    applyRuntime(next) {
      for (const key of ["fa", "fb", "alpha", "beta", "mFront", "mBack", "rearP", "kappa0"]) {
        if (next[key] !== undefined && Number.isFinite(next[key])) this[key] = next[key];
      }
    },
    gaps() {
      const g = new Float64Array(this.N);
      for (let n = 0; n < this.N; n++) {
        const ld = (n + 1) % this.N;
        let xl = this.x[ld];
        if (ld === 0) xl += this.L;
        g[n] = xl - this.x[n] - this.l;
      }
      return g;
    },
    modeAmp(m) {
      const { N, v } = this;
      let mean = 0;
      for (let n = 0; n < N; n++) mean += v[n];
      mean /= N;
      let re = 0, im = 0;
      for (let n = 0; n < N; n++) {
        const th = (-2 * Math.PI * m * n) / N;
        re += (v[n] - mean) * Math.cos(th);
        im += (v[n] - mean) * Math.sin(th);
      }
      return { re: re / N, im: im / N, mean };
    },
    velStd() {
      const { N, v } = this;
      let m = 0;
      for (let n = 0; n < N; n++) m += v[n];
      m /= N;
      let s = 0;
      for (let n = 0; n < N; n++) s += (v[n] - m) * (v[n] - m);
      return Math.sqrt(s / N);
    },
    gapStd() {
      const g = this.gaps();
      let m = 0; for (const q of g) m += q; m /= g.length;
      let z = 0; for (const q of g) z += (q - m) * (q - m);
      return Math.sqrt(z / g.length);
    },
  };
}

/* ---------- 最小二乘 ---------- */
function lsSlope(pts) {
  const n = pts.length;
  if (n < 8) return null;
  let sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (const [X, Y] of pts) { sx += X; sy += Y; sxx += X * X; sxy += X * Y; }
  const den = n * sxx - sx * sx;
  if (Math.abs(den) < 1e-12) return null;
  return (n * sxy - sx * sy) / den;
}

/* ---------- 配色 ---------- */
const CLR = {
  paper: "#E7EBE6",
  panel: "#F2F5F0",
  ink: "#151E27",
  soft: "#5C6B78",
  rule: "#C3CBC1",
  blue: "#3E5F8A",
  teal: "#1E6E5B",
  brick: "#A83A2B",
  amber: "#B07D24",
};
// 围绕平衡速度的发散色标
function speedColor(dv, span) {
  const t = Math.max(-1, Math.min(1, dv / span));
  if (t < 0) {
    const u = -t;
    return `rgb(${Math.round(231 + (168 - 231) * u)},${Math.round(235 + (58 - 235) * u)},${Math.round(230 + (43 - 230) * u)})`;
  }
  return `rgb(${Math.round(231 + (62 - 231) * t)},${Math.round(235 + (95 - 235) * t)},${Math.round(230 + (138 - 230) * t)})`;
}

/* ---------- 小型绘图 ---------- */
function plot(ctx, W, Ht, opts) {
  const { xs, series, xlab, ylab, xr, yr, marks = [], hline = null } = opts;
  ctx.clearRect(0, 0, W, Ht);
  ctx.fillStyle = CLR.panel;
  ctx.fillRect(0, 0, W, Ht);
  const P = { l: 44, r: 10, t: 12, b: 28 };
  const iw = W - P.l - P.r, ih = Ht - P.t - P.b;
  const X = (x) => P.l + ((x - xr[0]) / (xr[1] - xr[0])) * iw;
  const Y = (y) => P.t + ih - ((y - yr[0]) / (yr[1] - yr[0])) * ih;
  // 网格
  ctx.strokeStyle = CLR.rule;
  ctx.lineWidth = 0.5;
  ctx.font = "9px ui-monospace, monospace";
  ctx.fillStyle = CLR.soft;
  for (let i = 0; i <= 4; i++) {
    const yy = yr[0] + ((yr[1] - yr[0]) * i) / 4;
    ctx.beginPath(); ctx.moveTo(P.l, Y(yy)); ctx.lineTo(W - P.r, Y(yy)); ctx.stroke();
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    ctx.fillText(Math.abs(yy) < 1e4 ? yy.toFixed(Math.abs(yy) < 1 ? 2 : 1) : yy.toExponential(0), P.l - 4, Y(yy));
  }
  for (let i = 0; i <= 4; i++) {
    const xx = xr[0] + ((xr[1] - xr[0]) * i) / 4;
    ctx.beginPath(); ctx.moveTo(X(xx), P.t); ctx.lineTo(X(xx), Ht - P.b); ctx.stroke();
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(xx.toFixed(xr[1] - xr[0] < 5 ? 2 : 0), X(xx), Ht - P.b + 4);
  }
  if (hline !== null && hline >= yr[0] && hline <= yr[1]) {
    ctx.strokeStyle = CLR.ink; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(P.l, Y(hline)); ctx.lineTo(W - P.r, Y(hline)); ctx.stroke();
  }
  for (const s of series) {
    ctx.strokeStyle = s.color; ctx.lineWidth = s.w || 1.6;
    ctx.beginPath();
    let started = false;
    for (let i = 0; i < xs.length; i++) {
      const yv = s.ys[i];
      if (!isFinite(yv)) { started = false; continue; }
      const px = X(xs[i]), py = Y(Math.max(yr[0], Math.min(yr[1], yv)));
      if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
    }
    ctx.stroke();
    if (s.fillNeg) {
      ctx.fillStyle = "rgba(168,58,43,0.16)";
      ctx.beginPath();
      let open = false;
      for (let i = 0; i < xs.length; i++) {
        if (s.ys[i] < 0) {
          if (!open) { ctx.moveTo(X(xs[i]), Y(0)); open = true; }
          ctx.lineTo(X(xs[i]), Y(Math.max(yr[0], s.ys[i])));
        } else if (open) { ctx.lineTo(X(xs[i]), Y(0)); ctx.closePath(); ctx.fill(); open = false; }
      }
      if (open) { ctx.lineTo(X(xs[xs.length - 1]), Y(0)); ctx.closePath(); ctx.fill(); }
    }
  }
  for (const mk of marks) {
    ctx.strokeStyle = mk.color; ctx.lineWidth = 1.4;
    ctx.setLineDash([3, 3]);
    ctx.beginPath(); ctx.moveTo(X(mk.x), P.t); ctx.lineTo(X(mk.x), Ht - P.b); ctx.stroke();
    ctx.setLineDash([]);
    if (mk.label) {
      ctx.fillStyle = mk.color; ctx.font = "9px ui-monospace, monospace";
      ctx.textAlign = "left"; ctx.textBaseline = "top";
      ctx.fillText(mk.label, X(mk.x) + 3, P.t + 2);
    }
  }
  ctx.fillStyle = CLR.soft; ctx.font = "9px ui-monospace, monospace";
  ctx.textAlign = "right"; ctx.textBaseline = "bottom";
  ctx.fillText(xlab, W - P.r, Ht - 2);
  ctx.save(); ctx.translate(9, P.t + 2); ctx.rotate(-Math.PI / 2);
  ctx.textAlign = "right"; ctx.textBaseline = "top";
  ctx.fillText(ylab, 0, 0); ctx.restore();
}

/* ---------- 复用组件（定义在主组件外，避免拖动滑块时丢焦点） ---------- */
function Slider({ label, value, set, min, max, step, unit, disabled }) {
  return (
    <div style={{ opacity: disabled ? 0.35 : 1, marginBottom: 9 }}>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "baseline", marginBottom: 2 }}>
        <span style={{ fontSize: 11, color: CLR.soft, letterSpacing: "0.02em" }}>{label}</span>
        <span style={{ fontSize: 11.5, fontFamily: "ui-monospace, monospace", color: CLR.ink }}>
          {value.toFixed(step < 0.01 ? 3 : step < 1 ? 2 : 0)}
          <span style={{ color: CLR.soft, marginLeft: 3 }}>{unit}</span>
        </span>
      </div>
      <input
        type="range" min={min} max={max} step={step} value={value} disabled={disabled}
        onChange={(e) => set(parseFloat(e.target.value))}
        style={{ width: "100%", accentColor: CLR.blue, height: 3 }}
      />
    </div>
  );
}

function Group({ title, children }) {
  return (
    <div style={{ marginBottom: 16 }}>
      <div style={{
        fontSize: 10, letterSpacing: "0.16em", textTransform: "uppercase",
        color: CLR.soft, borderBottom: `1px solid ${CLR.rule}`, paddingBottom: 4, marginBottom: 9,
      }}>{title}</div>
      {children}
    </div>
  );
}

function Card({ title, note, children }) {
  return (
    <div style={{
      background: CLR.panel, border: `1px solid ${CLR.rule}`, padding: "10px 12px 8px",
      display: "flex", flexDirection: "column", minWidth: 0,
    }}>
      <div style={{ display: "flex", justifyContent: "space-between", alignItems: "baseline", marginBottom: 6, gap: 8 }}>
        <span style={{ fontSize: 11, fontWeight: 600, color: CLR.ink, letterSpacing: "0.02em" }}>{title}</span>
        {note && <span style={{ fontSize: 10, color: CLR.soft, fontFamily: "ui-monospace, monospace", textAlign: "right" }}>{note}</span>}
      </div>
      {children}
    </div>
  );
}

function MiniModal({ title, subtitle, onClose, children }) {
  return (
    <div className="tsl-modal-backdrop" role="presentation" onMouseDown={onClose}>
      <section className="tsl-mini-modal" role="dialog" aria-modal="true" aria-label={title}
        onMouseDown={(event) => event.stopPropagation()}>
        <div className="tsl-modal-head">
          <div><b>{title}</b>{subtitle && <span>{subtitle}</span>}</div>
          <button className="tsl-icon-close" type="button" onClick={onClose} aria-label="关闭参数窗口">×</button>
        </div>
        {children}
      </section>
    </div>
  );
}

const CSS = `
.tsl-shell{max-width:1540px;margin:0 auto;padding:16px 18px 36px}
.tsl-toolbar{display:flex;gap:8px;align-items:center;flex-wrap:wrap}
.tsl-case-strip{max-width:1504px;margin:13px auto 0;padding:9px 12px;background:#eaf3ef;border:1px solid ${CLR.teal};display:flex;align-items:center;justify-content:space-between;gap:12px;font-size:11px;line-height:1.45}
.tsl-case-strip b{display:block;color:${CLR.teal};font-size:12px}.tsl-case-strip span{color:${CLR.ink}}
.tsl-tabs{display:flex;gap:4px;overflow-x:auto;border-bottom:1px solid ${CLR.rule};padding:0 18px;background:${CLR.panel}}
.tsl-tab{border:0;border-bottom:3px solid transparent;background:transparent;color:${CLR.soft};padding:11px 15px 9px;cursor:pointer;white-space:nowrap;font-size:12px}
.tsl-tab.on{color:${CLR.ink};border-bottom-color:${CLR.blue};font-weight:600}
.tsl-panel{padding-top:14px}
.tsl-live-grid{display:grid;grid-template-columns:minmax(520px,1fr) minmax(300px,360px);gap:14px;align-items:start}
.tsl-ring-wrap{position:relative;min-width:0}
.tsl-ring-canvas{width:100%;max-height:720px;aspect-ratio:1/1;display:block;background:${CLR.panel}}
.tsl-focus .tsl-live-grid{grid-template-columns:1fr}
.tsl-focus .tsl-live-side{display:none}
.tsl-focus .tsl-ring-wrap{max-width:min(900px,92vw);margin:0 auto}
.tsl-live-side{display:grid;gap:12px}
.tsl-scenario-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:7px}
.tsl-scenario-tile{min-height:70px;border:1px solid ${CLR.rule};background:${CLR.paper};padding:7px 5px;color:${CLR.soft};cursor:pointer;text-align:center;display:grid;place-items:center;gap:3px}
.tsl-scenario-tile .symbol{font:700 18px/1 ui-monospace,monospace;color:${CLR.blue}}
.tsl-scenario-tile b{font-size:10px;color:${CLR.ink};line-height:1.15}.tsl-scenario-tile span{font-size:8.5px;line-height:1.15}
.tsl-scenario-tile.on{border-color:${CLR.teal};background:#eaf3ef;box-shadow:inset 0 0 0 1px ${CLR.teal}}
.tsl-scenario-tile.on .symbol{color:${CLR.teal}}
.tsl-core-controls{display:grid;grid-template-columns:1fr 1fr;gap:8px}
.tsl-core-button{border:1px solid ${CLR.rule};background:${CLR.paper};color:${CLR.ink};padding:10px;text-align:left;cursor:pointer;min-width:0}
.tsl-core-button:hover{border-color:${CLR.blue};background:#edf1f4}
.tsl-core-button .symbol{display:block;color:${CLR.blue};font:700 15px/1 ui-monospace,monospace;margin-bottom:6px}
.tsl-core-button b{display:block;font-size:11.5px;margin-bottom:4px}.tsl-core-button small{display:block;color:${CLR.soft};font-size:9.5px;line-height:1.35}
.tsl-modal-backdrop{position:fixed;inset:0;z-index:1200;background:rgba(22,35,29,.28);display:grid;place-items:center;padding:18px}
.tsl-mini-modal{width:min(420px,100%);max-height:calc(100vh - 36px);overflow:auto;background:${CLR.panel};border:1px solid ${CLR.rule};box-shadow:0 18px 55px rgba(20,37,29,.28);padding:14px}
.tsl-modal-head{display:flex;justify-content:space-between;gap:12px;align-items:flex-start;border-bottom:1px solid ${CLR.rule};padding-bottom:9px;margin-bottom:12px}.tsl-modal-head b{display:block;font-size:14px}.tsl-modal-head span{display:block;font-size:10px;color:${CLR.soft};margin-top:3px;line-height:1.4}
.tsl-icon-close{border:0;background:transparent;color:${CLR.soft};font-size:24px;line-height:1;cursor:pointer;padding:0 2px}
.tsl-modal-toggle{display:flex;gap:8px;align-items:center;font-size:11.5px;margin-bottom:12px}
.tsl-statusline{display:grid;grid-template-columns:repeat(3,1fr);gap:8px;margin-top:10px}
.tsl-status{background:${CLR.paper};padding:8px 9px;min-width:0}
.tsl-status b{display:block;font-family:ui-monospace,monospace;font-size:14px;margin-top:2px}
.tsl-status span{font-size:10px;color:${CLR.soft}}
.tsl-live-note{font-size:10.5px;line-height:1.55;color:${CLR.soft};margin-top:8px}
.tsl-event{display:grid;grid-template-columns:72px 1fr;gap:7px;padding:6px 0;border-top:1px solid ${CLR.rule};font-size:10.5px}
.tsl-event:first-child{border-top:0}
.tsl-event time{font-family:ui-monospace,monospace;color:${CLR.soft}}
.tsl-propagation{display:grid;grid-template-columns:minmax(0,1.15fr) minmax(0,.85fr);gap:12px;align-items:start}
.tsl-prop-stack{display:grid;gap:12px}
.tsl-canvas-caption{display:flex;justify-content:space-between;gap:10px;font-size:10px;color:${CLR.soft};font-family:ui-monospace,monospace;margin-top:5px}
.tsl-settings{display:grid;grid-template-columns:repeat(3,minmax(250px,1fr));gap:12px;align-items:start}
.tsl-main{display:none}
.tsl-main.theory{display:block}
.tsl-main.theory>.tsl-side,.tsl-main.theory .tsl-hero,.tsl-main.theory .tsl-diagnostics{display:none}
.tsl-main.settings{display:block;padding:14px 18px 36px;max-width:1540px;margin:0 auto}
.tsl-main.settings>.tsl-side{display:grid;grid-template-columns:repeat(3,minmax(250px,1fr));gap:12px 22px;padding:0;background:transparent!important;border:0!important}
.tsl-main.settings>div:last-child{display:none}
.tsl-hero{display:none}
.tsl-side{min-width:0}
.tsl-charts{display:grid;grid-template-columns:repeat(2,1fr);gap:12px;margin-top:12px}
.tsl-diagnostics{display:grid;grid-template-columns:1fr 1fr;gap:12px;margin-top:12px}
.tsl-expgrid{display:grid;grid-template-columns:minmax(220px,.9fr) 1.6fr auto;gap:10px;align-items:center}
.tsl-nonlinear{display:grid;gap:12px}
.tsl-nl-picker{display:grid;grid-template-columns:repeat(5,minmax(120px,1fr));gap:7px}
.tsl-nl-choice{border:1px solid ${CLR.rule};background:${CLR.panel};color:${CLR.soft};padding:9px 10px;text-align:left;cursor:pointer;min-height:60px}
.tsl-nl-choice b{display:block;color:${CLR.ink};font-size:11.5px;margin-bottom:3px}
.tsl-nl-choice span{display:block;font-size:9.5px;line-height:1.35}
.tsl-nl-choice.on{border-color:${CLR.blue};box-shadow:inset 0 0 0 1px ${CLR.blue};background:#edf1f4}
.tsl-nl-grid{display:grid;grid-template-columns:minmax(0,1.45fr) minmax(280px,.55fr);gap:12px;align-items:start}
.tsl-nl-metric{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin:9px 0}
.tsl-nl-metric>div{background:${CLR.paper};padding:8px 9px;min-width:0}
.tsl-nl-metric span{display:block;font-size:9.5px;color:${CLR.soft}}
.tsl-nl-metric b{display:block;margin-top:2px;font-family:ui-monospace,monospace;font-size:13px}
.tsl-nl-badge{display:inline-block;padding:3px 7px;background:${CLR.teal};color:#fff;font-size:9.5px;letter-spacing:.04em}
.tsl-nl-badge.warn{background:${CLR.amber}}
.tsl-applications{display:grid;gap:12px}
.tsl-app-picker{display:grid;grid-template-columns:repeat(2,minmax(240px,1fr));gap:8px}
.tsl-app-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:12px}
.tsl-app-bottom{display:grid;grid-template-columns:minmax(0,1.35fr) minmax(280px,.65fr);gap:12px;align-items:start}
.tsl-app-benefits{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin-top:8px}
.tsl-app-benefits>div{background:${CLR.paper};padding:8px 9px;min-width:0}
.tsl-app-benefits span{display:block;font-size:9.5px;color:${CLR.soft}}
.tsl-app-benefits b{display:block;margin-top:2px;font-family:ui-monospace,monospace;font-size:13px}
.tsl-select{width:100%;padding:7px 8px;border:1px solid ${CLR.rule};background:${CLR.panel};color:${CLR.ink};font-size:11.5px}
.tsl-chip{border:1px solid ${CLR.rule};background:transparent;color:${CLR.soft};padding:4px 8px;font-size:10.5px;cursor:pointer}
.tsl-chip.on{background:${CLR.ink};border-color:${CLR.ink};color:white}
@media(max-width:1120px){
  .tsl-live-grid{grid-template-columns:minmax(440px,1fr) 310px}
  .tsl-settings{grid-template-columns:repeat(2,minmax(240px,1fr))}
  .tsl-propagation{grid-template-columns:1fr}
  .tsl-nl-picker{grid-template-columns:repeat(3,minmax(120px,1fr))}
  .tsl-nl-grid{grid-template-columns:1fr}
  .tsl-app-bottom{grid-template-columns:1fr}
  .tsl-main.settings>.tsl-side{grid-template-columns:repeat(2,minmax(240px,1fr))}
}
@media(max-width:760px){
  .tsl-shell{padding:10px 10px 28px}
  .tsl-live-grid{grid-template-columns:1fr}
  .tsl-settings{grid-template-columns:1fr}
  .tsl-charts{grid-template-columns:1fr}
  .tsl-diagnostics{grid-template-columns:1fr}
  .tsl-expgrid{grid-template-columns:1fr}
  .tsl-nl-picker{display:flex;overflow-x:auto;padding-bottom:4px}
  .tsl-nl-choice{min-width:155px}
  .tsl-app-picker,.tsl-app-grid{grid-template-columns:1fr}
  .tsl-statusline{grid-template-columns:repeat(2,1fr)}
  .tsl-ring-canvas{max-height:none}
  .tsl-scenario-grid{grid-template-columns:repeat(2,minmax(0,1fr))}
  .tsl-case-strip{margin:10px 12px 0;align-items:flex-start;flex-direction:column}
  .tsl-main.settings{padding:10px}
  .tsl-main.settings>.tsl-side{grid-template-columns:1fr}
}
input[type=range]{cursor:pointer}
.tsl-btn:focus-visible,.tsl-tab:focus-visible,input:focus-visible,select:focus-visible{outline:2px solid ${CLR.blue};outline-offset:2px}
`;

const EXPERIMENTS = {
  E1: { title: "中性稳定线", section: "§4–5 基准与 IDM", does: "扫描平衡速度，寻找扰动增长率变号的临界点", metric: "不稳定速度带" },
  E2: { title: "线性增长率", section: "§4–11 通用", does: "在线性阶段拟合 ln|A(t)| 的斜率", metric: "Re λ：理论 vs 实测" },
  E3: { title: "最不稳定波数", section: "§8–9 环形 Fourier", does: "对车辆速度做空间 FFT，识别主导模态", metric: "m*：理论 vs 实测" },
  E4: { title: "逐车传递比", section: "§7 单向前馈", does: "扫描正弦频率，计算单车最大传递增益", metric: "maxω |G(iω)|" },
  E5: { title: "有限 N 效应", section: "附录：有限环道", does: "比较当前 N 与 N→∞ 的稳定裕度", metric: "最低可用波数与判定" },
  E6: { title: "异质排序效应", section: "§12 异质车队", does: "保持配方不变，比较均匀、随机和连续排列", metric: "沿车队中途峰值" },
  E7: { title: "AV 临界渗透率", section: "§12.10 混合交通", does: "由异质长波预算求最小 AV 比例", metric: "p_min" },
  E8: { title: "延迟势阱", section: "§11 时间延迟", does: "比较 Ψ(0) 与全频段 min Ψ(ω)", metric: "有限频率稳定裕度" },
  E9: { title: "后视闭合点", section: "§10 后向间距", does: "扫描后视权重，寻找不稳定带闭合点", metric: "p_c" },
  E10: { title: "走停波波速", section: "§10 波速", does: "由时空模态相位斜率测波速", metric: "c₁：理论 vs 实测" },
};

/* ---------- E11–E12：非线性现象的可复现实验数据 ----------
   全 IDM、N=120、RK4、dt=0.05 s；曲线为离线长时积分的紧凑采样。
   孤波/kink 面板是约化方程模板，明确不冒充基准 IDM 的直接数值解。 */
const NONLINEAR_CASES = {
  steepening: {
    label: "波形陡化", kicker: "高次谐波增长", evidence: "同一 m=4 波形的波肩变尖，高次谐波能量相对基波由 0 增至 0.190。",
    verdict: "完整 IDM 证据", caveat: "外形变陡说明非线性已介入，但不能单凭外形判成孤波。",
  },
  saturation: {
    label: "非线性饱和", kicker: "增长率随振幅下降", evidence: "vₑ=8 m/s 时 σᵥ 从 0.0141 增至 2.962 m/s，末 300 s 均值 2.554 m/s；对数曲线后段明显弯折。",
    verdict: "完整 IDM 证据", caveat: "后段不能再用一个固定 Re λ 拟合；1800 s 末仍缓慢增长，因此称“进入饱和区”，不声称已达严格定常波。",
  },
  soliton: {
    label: "孤立波与 kink", kicker: "幅宽与波速缩放", evidence: "KdV 模板满足 L∝|A|⁻¹ᐟ²；mKdV kink 模板满足 L∝|A|⁻¹。拖动振幅可直接观察幅宽收缩。",
    verdict: "约化方程模板", caveat: "基准 IDM 的中性点不严格满足 V″=0；kink 只用于说明判据，不能当作 IDM 已出现 mKdV kink 的证明。",
  },
  trigger: {
    label: "有限幅触发与迟滞", kicker: "同一稳定工况、不同初幅", evidence: "vₑ=20 m/s 的小脉冲迅速恢复，大脉冲在 1800 s 仍保留显著波动，显示有限幅暂态的强非线性。",
    verdict: "触发证据；迟滞未证实", caveat: "这是单向幅值扫描。真正迟滞还必须从均匀态和拥堵态做双向参数扫描，并得到不同转变点。",
  },
  recovery: {
    label: "振幅依赖恢复", kicker: "恢复时间 tᵣ(A₀)", evidence: "恢复定义为 σᵥ<0.02 m/s 连续保持 300 s：A₀=0.5、2、4、6 m/s 分别需 4、164、1374、4700 s。",
    verdict: "完整 IDM 证据", caveat: "A₀≥8 m/s 在 5400 s 内尚未满足恢复阈值，图中按右删失显示，而不是宣称永不恢复。",
  },
};

const NONLINEAR_DATA = {"time":[0,30,60,90,120,150,180,210,240,270,300,330,360,390,420,450,480,510,540,570,600,630,660,690,720,750,780,810,840,870,900,930,960,990,1020,1050,1080,1110,1140,1170,1200,1230,1260,1290,1320,1350,1380,1410,1440,1470,1500,1530,1560,1590,1620,1650,1680,1710,1740,1770,1800],"unstableStd":[0.01414214,0.00725894,0.00819213,0.00926034,0.01046784,0.0118328,0.01337573,0.01511985,0.01709138,0.01931995,0.02183907,0.0246866,0.02790533,0.03154364,0.03565615,0.04030461,0.04555877,0.05149741,0.05820949,0.06579544,0.07436861,0.08405685,0.0950043,0.10737334,0.12134672,0.13712984,0.15495324,0.1750751,0.19778384,0.22340054,0.25228111,0.28481781,0.32143977,0.36261192,0.40883144,0.46062084,0.5185162,0.58304909,0.6547206,0.73396585,0.82110863,0.91630665,1.01949126,1.13030842,1.24807173,1.37174105,1.49993933,1.63101469,1.76314438,1.89446486,2.02320433,2.14779358,2.26693963,2.37965845,2.48527314,2.5833885,2.67385225,2.75671089,2.83216541,2.90053016,2.96219714],"stableStd":[0.0141421356,0.0076495468,0.0062374603,0.0050542748,0.0040985837,0.0033234319,0.0026948903,0.0021852209,0.0017719425,0.0014368251,0.0011650866,0.0009447404,0.000766067,0.0006211852,0.000503704,0.0004084413,0.0003311951,0.0002685581,0.0002177672,0.0001765822,0.0001431862,0.0001161062,0.0000941477,0.0000763421,0.0000619039,0.0000501964,0.000040703,0.0000330051,0.000026763,0.0000217015,0.0000175972,0.0000142692,0.0000115705,0.0000093823,0.0000076078,0.000006169,0.0000050023,0.0000040562,0.0000032891,0.0000026671,0.0000021627,0.0000017536,0.000001422,0.0000011531,0.000000935,0.0000007582,0.0000006148,0.0000004985,0.0000004042,0.0000003278,0.0000002658,0.0000002155,0.0000001748,0.0000001417,0.0000001149,0.0000000932,0.0000000756,0.0000000613,0.0000000497,0.0000000403,0.0000000327],"profileTimes":[0,450,900,1350,1800],"profiles":[[8,8.00416,8.00813,8.01176,8.01486,8.01732,8.01902,8.01989,8.01989,8.01902,8.01732,8.01486,8.01176,8.00813,8.00416,8,7.99584,7.99187,7.98824,7.98514,7.98268,7.98098,7.98011,7.98011,7.98098,7.98268,7.98514,7.98824,7.99187,7.99584],[7.95062,7.94582,7.94336,7.94335,7.94581,7.95064,7.95763,7.96648,7.97683,7.98821,8.00012,8.01204,8.02344,8.0338,8.04266,8.04964,8.05442,8.05679,8.05667,8.05405,8.04907,8.04195,8.03302,8.02266,8.01134,7.99954,7.98777,7.97654,7.96632,7.95756],[8.02632,7.95439,7.88464,7.81975,7.76218,7.71419,7.6777,7.65428,7.64508,7.65074,7.67137,7.70645,7.7548,7.8146,7.88339,7.95813,8.03536,8.11136,8.18237,8.24486,8.29573,8.33255,8.35369,8.3584,8.34682,8.31988,8.27918,8.22686,8.16542,8.09762],[9.59856,9.32548,9.01648,8.68609,8.3457,8.00453,7.67037,7.35025,7.05091,6.7791,6.54181,6.34637,6.20044,6.11194,6.08864,6.13773,6.26504,6.47397,6.76428,7.13062,7.56105,8.03575,8.52651,8.99785,9.4108,9.72972,9.93016,10.00425,9.96031,9.81743],[11.63957,11.07524,10.45917,9.81939,9.17131,8.52445,7.88577,7.26142,6.6576,6.08117,5.54007,5.04392,4.60455,4.23663,3.95824,3.7911,3.76006,3.89124,4.20854,4.72899,5.45806,6.3856,7.4818,8.69073,9.92056,11.03648,11.87841,12.32549,12.3682,12.10007]],"harmonicRatio":[0,0.00341,0.02134,0.10912,0.19027],"pulseAmplitude":[0.5,1,2,4,6,8,10,12,14,16,18],"pulsePeak":[0.045453,0.090906,0.181812,0.363624,0.545436,0.776178,1.12649,1.648254,2.238683,2.956423,3.921149],"pulseLate":[0.00088,0.001966,0.004996,0.017766,0.055737,0.177045,0.461179,0.666703,0.640622,0.774663,1.35145],"recoveryTime":[4,26,164,1374,4700,5400,5400,5400,5400,5400,5400],"recoveryCensored":[false,false,false,false,false,true,true,true,true,true,true]};

const APPLICATION_DATA = typeof window !== "undefined" ? window.APPLICATION_CASES_DATA : null;

function decodeQ8(encoded) {
  const raw = atob(encoded), out = new Uint8Array(raw.length);
  for (let i = 0; i < raw.length; i++) out[i] = raw.charCodeAt(i);
  return out;
}

function decodeQ16(encoded) {
  const raw = atob(encoded), out = new Uint16Array(raw.length / 2);
  for (let i = 0; i < out.length; i++) out[i] = raw.charCodeAt(2 * i) | (raw.charCodeAt(2 * i + 1) << 8);
  return out;
}

function transferMag(part, fa, w) {
  const num = C(part.fs - fa * w * w, part.fvl * w);
  const den = C(part.fs - w * w, -part.fv * w);
  return cAbs(cDiv(num, den));
}

function unstableIntervals(samples) {
  const bands = [];
  let start = null;
  for (let i = 0; i < samples.length; i++) {
    if (samples[i].bad && start === null) start = samples[i].v;
    if ((!samples[i].bad || i === samples.length - 1) && start !== null) {
      bands.push([start, samples[Math.max(0, i - (samples[i].bad ? 0 : 1))].v]); start = null;
    }
  }
  return bands;
}

/* ---------- 混合交通的实际增益与稳定性 ---------- */
function mixedGainProfile(N, fa, kappa0, avShare, arrangement, avGroups) {
  const av = makeArrangement(N, N * avShare, arrangement, avGroups, 71);
  const gains = new Float64Array(N);
  let avCount = 0, aaCount = 0, ahCount = 0, sum = 0;
  for (let n = 0; n < N; n++) {
    if (!av[n]) continue;
    avCount++;
    const q = (n + 1) % N;
    const gain = av[q] ? fa : Math.min(fa, kappa0);
    if (av[q]) aaCount++; else ahCount++;
    gains[n] = gain;
    sum += gain;
  }
  return { av, gains, avCount, aaCount, ahCount, meanGain: sum / N };
}

function mixedTransferMargin(part, gains, w) {
  if (w === 0) return 0;
  let logProduct = 0;
  for (const gain of gains) {
    logProduct += Math.log(Math.max(1e-14, transferMag(part, gain, w)));
  }
  // 正值表示一圈传播后衰减，负值表示一圈传播后放大。
  return -logProduct / gains.length;
}

function buildAvStabilityPlane({ P, N, fa, kappa0, arrangement, avGroups, currentShare, currentSpeed }) {
  const nx = 61, ny = 61, vMax = P.v0 * 0.999;
  const profiles = Array.from({ length: nx }, (_, ix) =>
    mixedGainProfile(N, fa, kappa0, ix / (nx - 1), arrangement, avGroups)
  );
  const cells = new Array(nx * ny);
  for (let iy = 0; iy < ny; iy++) {
    const v = (iy / (ny - 1)) * vMax;
    const part = idmPartials(v, P);
    const frequencyLogs = [];
    for (let iw = 1; iw <= 90; iw++) {
      const w = iw * (6 / 90);
      frequencyLogs.push({
        g0: Math.log(Math.max(1e-14, transferMag(part, 0, w))),
        gaa: Math.log(Math.max(1e-14, transferMag(part, fa, w))),
        gah: Math.log(Math.max(1e-14, transferMag(part, Math.min(fa, kappa0), w))),
      });
    }
    for (let ix = 0; ix < nx; ix++) {
      const share = ix / (nx - 1), profile = profiles[ix];
      const pMean = {
        part, fa: profile.meanGain, fb: 0, alpha: 0, beta: 0,
        mFront: 1, mBack: 1, rearP: 0, tau0: 0, tauA: 0,
      };
      const longMargin = 2 * longWaveMargin(part, pMean);
      let transferMargin = Infinity;
      const n0 = N - profile.aaCount - profile.ahCount;
      for (const logs of frequencyLogs) {
        const margin = -(n0 * logs.g0 + profile.aaCount * logs.gaa + profile.ahCount * logs.gah) / N;
        transferMargin = Math.min(transferMargin, margin);
      }
      const longNorm = longMargin / Math.max(2 * part.fs, 1e-8);
      const score = Math.min(longNorm, transferMargin * 10);
      cells[iy * nx + ix] = {
        share, v, score, longMargin, transferMargin,
        stable: longMargin >= -1e-9 && transferMargin >= -1e-9,
        meanGain: profile.meanGain,
        aaCount: profile.aaCount,
        ahCount: profile.ahCount,
      };
    }
  }
  const ix = Math.max(0, Math.min(nx - 1, Math.round(currentShare * (nx - 1))));
  const iy = Math.max(0, Math.min(ny - 1, Math.round((currentSpeed / vMax) * (ny - 1))));
  return { nx, ny, vMax, cells, current: cells[iy * nx + ix] };
}

function analyzeTheorySnapshot({ se, P, eff, N, useHetero, useMixed, avShare, arrangement, avGroups, kappa0 }) {
  if (se <= P.s0 + 0.01) return null;
  const ve = idmEquilibriumSpeed(se, P);
  if (ve <= 0.02) return null;
  const part = idmPartials(ve, P);
  const p = { part, ...eff };
  const mixedSimple = useMixed && !useHetero && eff.fb === 0 && eff.rearP === 0 &&
    eff.mFront === 1 && eff.mBack === 1 && eff.tau0 === 0 && eff.tauA === 0;
  const mixed = useMixed
    ? mixedGainProfile(N, eff.fa, kappa0, avShare, arrangement, avGroups)
    : null;

  // 混合交通不能把全车都当作 κ=fa。有限环谱只作为等效均匀参考；
  // 稳定标签使用实际 {0, κ0, κ} 增益序列的环周传递乘积。
  const pSpectrum = mixedSimple
    ? { ...p, fa: mixed.meanGain, alpha: 0, mFront: 1 }
    : p;
  const exactPsi = eff.fb === 0 && eff.rearP === 0 && eff.mFront === 1 && !useHetero && !useMixed;
  const ws = [], psis = [];
  let psi0;
  if (mixedSimple) {
    psi0 = 2 * longWaveMargin(part, pSpectrum);
    for (let i = 0; i <= 400; i++) {
      const w = (i / 400) * 6;
      ws.push(w);
      psis.push(mixedTransferMargin(part, mixed.gains, w));
    }
  } else {
    psi0 = 2 * longWaveMargin(part, p);
    for (let i = 0; i <= 400; i++) {
      const w = (i / 400) * 6;
      ws.push(w);
      psis.push(exactPsi ? psiOmega(w, p) : psi0 + (1 - Math.min(0.98, Math.abs(kernelTotal(p)))) * w * w);
    }
  }
  const psiMin = Math.min(...psis);
  const spec = ringSpectrum(pSpectrum, N);
  let best = { m: 0, k: 0, lam: C(-1, 0) };
  for (const s of spec) if (s.lam.re > best.lam.re) best = s;

  const band = [];
  for (let i = 1; i <= 240; i++) {
    const vv = (i / 240) * P.v0 * 0.999;
    const pt = idmPartials(vv, P);
    const pp = { part: pt, ...eff };
    let mn;
    if (mixedSimple) {
      mn = 2 * longWaveMargin(pt, { ...pp, fa: mixed.meanGain, alpha: 0, mFront: 1 });
    } else {
      mn = 2 * longWaveMargin(pt, pp);
      if (exactPsi) {
        mn = Infinity;
        for (let j = 0; j <= 120; j++) mn = Math.min(mn, psiOmega((j / 120) * 6, pp));
      }
    }
    band.push({ v: vv, val: mn });
  }

  const inconclusive = useHetero || (useMixed && !mixedSimple);
  const isUnstable = mixedSimple
    ? psi0 < -1e-9 || psiMin < -1e-9
    : best.lam.re > 1e-6 || psiMin < -1e-9;
  return {
    ve, part, psi0, ws, psis, psiMin, spec, best, band, p, exactPsi,
    mixed, mixedSimple, isUnstable, inconclusive,
    statusBasis: mixedSimple ? "实际 AV/HDV 环周传递" : inconclusive ? "需异质矩阵复核" : exactPsi ? "精确解析判据" : "均匀参考谱",
  };
}

/* ============================================================
   书稿直达链接：PDF 中的 ?case=… 会定位到对应实验分区，并给出复现口径。
   ============================================================ */
const BOOK_CASES = {
  L01: { tab: "live", title: "局部稳定、临界阻尼与欠阻尼恢复", hint: "保持 S0；用完整参数把平衡工作点与 IDM 偏导带入局部二阶根，再在探针页观察恢复。" },
  W01: { tab: "live", title: "三辆 IDM 的两个 RK4 更新步", hint: "设置 N=3 后打开完整参数；书中 W01 给出与环道仿真相同的同步 RK4 子步。" },
  M01: { tab: "live", title: "局部欠阻尼恢复与 RK4 超调", hint: "本案例的解析根与 RK4 更新见书中 M01；实验台用于观察同类扰动的探针响应。" },
  M02: { tab: "live", title: "有限环道与离散波长", hint: "依次比较 N=12、20、30，并在理论与验证区检查最低允许波数。" },
  E1: { tab: "live", title: "中性稳定线扫描", hint: "保持 S0，扫描密度或环道长度；目标是寻找当前理论由稳定变为失稳的边界。" },
  E2: { tab: "live", title: "线性增长率与角频率拟合", hint: "保持 S0、小幅扰动；运行后在理论与验证区比对解析谱与实测增长率。" },
  E3: { tab: "live", title: "最不稳定离散波数", hint: "保持 S0，运行后在传播与探针区查看时空图，并用理论谱核对主导模态。" },
  E4: { tab: "live", title: "前车加速度前馈的逐车传递比", hint: "启用前车加速度图标，取 κ=0.15；在理论与验证区比较频率响应。" },
  E5: { tab: "live", title: "有限 N 效应", hint: "保持 S0，逐次调整车辆数 N；危险长波只会在足够大的环道中出现。" },
  E6: { tab: "live", title: "异质车流的排列效应", hint: "启用 HV 图标，固定配方后依次选择均匀、随机、连续和成组排列。" },
  E7: { tab: "live", title: "AV 渗透率、速度与编组二维稳定域", hint: "启用 AV 图标；改变渗透率、编组数与 κ₀，并在理论区读取二维稳定平面。" },
  E8: { tab: "live", title: "反应与通信延迟的频域扫描", hint: "启用时间延迟图标，扫描 τ₀ 与 τₐ；低频门槛不变，有限频率响应会改变。" },
  E9: { tab: "live", title: "后向间距权重与稳定带闭合", hint: "启用后向间距图标，扫描 p 并比较理论稳定带。" },
  E10: { tab: "live", title: "后向反馈下的长波传播速度", hint: "启用后向间距图标；在传播页量取色带斜率并与理论波速比较。" },
  S0: { tab: "live", title: "原始 IDM 的 1800 s 稳定／失稳轨迹", hint: "保持 S0；分别设置书中稳定与失稳工作点，运行后转到传播与探针查看位置--时间轨迹、速度和间距。" },
  S1: { tab: "live", title: "前车加速度前馈的轨迹对照", hint: "启用前车加速度图标并设置 κ=0.15；与 S0 的同一扰动作单因素比较。" },
  S2: { tab: "live", title: "双向加速度耦合的轨迹对照", hint: "依次启用前车和后车加速度图标，取 (κ,fᵦ)=(0.30,0.10)。" },
  S3: { tab: "live", title: "多前／后车反馈的空间分配", hint: "启用多车信息图标，并在弹窗中设置 α、β、m₊ 与 m₋。" },
  S4: { tab: "live", title: "后向间距反馈的轨迹对照", hint: "启用后向间距图标，比较 p=0 与书中稳定化权重下的时空图。" },
  S5: { tab: "live", title: "时间延迟诱发走停波", hint: "启用时间延迟图标，保持其他参数不变，只扫描 τ₀ 或 τₐ。" },
  S6: { tab: "live", title: "异质车流的三种排序轨迹", hint: "启用 HV 图标；固定配方，比较不同排序下的中途放大和探针曲线。" },
  S7: { tab: "live", title: "AV/HDV 混合编组的稳定平面", hint: "启用 AV 图标，扫描渗透率与平衡速度，并改变编组数。" },
  M03: { tab: "live", title: "原始 IDM 单模态包络", hint: "使用 S0 的小幅单模态扰动；稳定与失稳工作点只改变平衡速度。" },
  M04: { tab: "live", title: "运行中开启前车加速度反馈", hint: "先让 S0 形成波动，再通过 κ 图标在线设置 κ=0.25；不要重置车辆状态。" },
  M05: { tab: "live", title: "时间步长收敛检查", hint: "对同一工况分别取 Δt 和 Δt/2，比较增长率、波速和临界点。" },
  E3B: { tab: "live", title: "双向加速度隐式耦合的环道谱", hint: "启用前/后车加速度，理论判据需要检查全部离散波数而非只看长波。" },
  E3M: { tab: "live", title: "多前／后车信息的空间分配", hint: "启用多车信息；固定加速度核总量后，比较不同空间衰减的有限波长结果。" },
  E11: { tab: "nonlinear", title: "有限幅扰动的波形陡化与非线性饱和", hint: "已定位到非线性现象页；用选择器切换剖面陡化、饱和、孤波/kink 和触发恢复证据。" },
  E12: { tab: "nonlinear", title: "稳定工况下的有限幅制动恢复", hint: "已定位到非线性现象页；选择有限幅恢复并查看恢复时间与残余波动。" },
  M06: { tab: "live", title: "IDM 基本图下的 LWR 激波", hint: "书中 M06 与 W02 给出有限体积更新；在线实验台用于把微观位置--时间图与宏观守恒波联系起来。" },
  W02: { tab: "live", title: "四格 LWR 有限体积更新", hint: "书中 W02 展开四格的通量和两个时间步；可对照平台中的宏观稳定性实验。" },
  W03: { tab: "live", title: "三格 ARZ 守恒量、通量与松弛更新", hint: "书中 W03 展开 ARZ 两个更新步；可对照 E13--E14 的宏观稳定性结果。" },
  E13: { tab: "live", title: "IDM 一致 ARZ 的稳定／失稳谱实验", hint: "书中 E13 给出完整 ARZ 谱推进与时空图；本平台用于核验共享 IDM 工作点的微观证据。" },
  E14: { tab: "live", title: "IDM 与 ARZ 的速度时空图对照", hint: "书中 E14 的宏微观图使用相同 IDM 平衡态与长波；平台可复核微观端的位置--时间轨迹。" },
  D01: { tab: "data", title: "NGSIM US--101 自然扰动传播证据", hint: "已定位到真实数据页；阅读片段筛选、增益、时滞及其不确定性。" },
  D02: { tab: "data", title: "稳定性标定、留出验证与谱约束优化", hint: "真实数据页给出 D01 证据；书中 D02 再将其用于标定、递推验证与反事实优化。" },
  A01: { tab: "applications", application: "open", title: "开放车队中的协同 AV 稳定性收益", hint: "已定位到案例 A；比较 HDV 与协同 AV 的两张位置--时间图和指定车辆速度曲线。" },
  A02: { tab: "applications", application: "bottleneck", title: "低速瓶颈下的事件触发型固定 VSL", hint: "已定位到案例 B；VSL 只在慢车事件发生后开启，并以稳定裕度选择 114 km/h。" },
};

const BOOK_CASES_EN = {
  L01: { title: "Local Stability, Critical Damping, and Underdamped Recovery", hint: "Keep S0 active. Use the IDM parameters to evaluate the local second-order roots at the equilibrium point, then inspect recovery on the probe page." },
  W01: { title: "Two RK4 Steps for Three IDM Vehicles", hint: "Set N=3 and open the scenario parameters. W01 uses the same synchronous RK4 substeps as the ring-road simulation." },
  M01: { title: "Local Underdamped Recovery and RK4 Overshoot", hint: "The analytical roots and RK4 update appear in M01; use the laboratory to inspect the probe response to a comparable disturbance." },
  M02: { title: "Finite Ring Road and Discrete Wavelengths", hint: "Compare N=12, 20, and 30 in sequence, then inspect the lowest admissible wavenumber under Theory and Validation." },
  E1: { title: "Neutral-Stability Line Scan", hint: "Keep S0 active and scan density or ring length to locate where the current theoretical verdict changes from stable to unstable." },
  E2: { title: "Fitting Linear Growth Rate and Angular Frequency", hint: "Keep S0 active with a small disturbance, then compare the analytical spectrum with the measured growth rate under Theory and Validation." },
  E3: { title: "Most Unstable Discrete Wavenumber", hint: "Keep S0 active, run the simulation, inspect the spatiotemporal diagram under Propagation and Probe, and compare the dominant mode with the theoretical spectrum." },
  E4: { title: "Vehicle-to-Vehicle Transfer Ratio with Leader-Acceleration Feedforward", hint: "Enable the leader-acceleration tile and set κ=0.15, then compare the frequency responses under Theory and Validation." },
  E5: { title: "Finite-N Effect", hint: "Keep S0 active and vary the vehicle count N. Dangerous long waves appear only when the ring is sufficiently large." },
  E6: { title: "Arrangement Effects in Heterogeneous Traffic", hint: "Enable the HV tile, keep the mixture fixed, and compare uniform, random, contiguous, and grouped arrangements." },
  E7: { title: "AV Penetration, Speed, and Platooning Stability Plane", hint: "Enable the AV tile, vary penetration, platoon count, and κ₀, then inspect the two-dimensional stability plane under Theory and Validation." },
  E8: { title: "Frequency-Domain Scan of Reaction and Communication Delays", hint: "Enable the delay tile and scan τ₀ and τₐ. The low-frequency threshold remains unchanged while the finite-frequency response varies." },
  E9: { title: "Rear-Gap Weight and Stability-Band Closure", hint: "Enable the rear-gap tile, scan p, and compare the resulting theoretical stability bands." },
  E10: { title: "Long-Wave Propagation Speed with Rear Feedback", hint: "Enable the rear-gap tile, measure the slope of the color bands on the propagation page, and compare it with the theoretical wave speed." },
  S0: { title: "Original IDM: 1800 s Stable and Unstable Trajectories", hint: "Keep S0 active, enter the stable and unstable operating points from the book, run the simulation, and inspect position-time trajectories, speed, and spacing." },
  S1: { title: "Trajectory Comparison with Leader-Acceleration Feedforward", hint: "Enable the leader-acceleration tile and set κ=0.15, then compare it with S0 under the same disturbance." },
  S2: { title: "Trajectory Comparison with Bidirectional Acceleration Coupling", hint: "Enable the leader- and follower-acceleration tiles and set (κ, fᵦ)=(0.30, 0.10)." },
  S3: { title: "Spatial Allocation of Multi-Leader and Multi-Follower Feedback", hint: "Enable the multi-vehicle-information tile and set α, β, m₊, and m₋ in its dialog." },
  S4: { title: "Trajectory Comparison with Rear-Gap Feedback", hint: "Enable the rear-gap tile and compare p=0 with the stabilizing weight used in the book." },
  S5: { title: "Stop-and-Go Waves Induced by Time Delay", hint: "Enable the delay tile, keep the other parameters fixed, and scan only τ₀ or τₐ." },
  S6: { title: "Three Arrangement Trajectories for Heterogeneous Traffic", hint: "Enable the HV tile, keep the mixture fixed, and compare transient amplification and probe traces across arrangements." },
  S7: { title: "Stability Plane for Mixed AV/HDV Platoons", hint: "Enable the AV tile, scan penetration and equilibrium speed, and vary the number of platoons." },
  M03: { title: "Single-Mode Envelope of the Original IDM", hint: "Use the small single-mode disturbance from S0. The stable and unstable operating points differ only in equilibrium speed." },
  M04: { title: "Enabling Leader-Acceleration Feedback During a Run", hint: "Let waves develop under S0, then set κ=0.25 online through the κ tile without resetting the vehicle state." },
  M05: { title: "Time-Step Convergence Check", hint: "Run the same case with Δt and Δt/2, then compare growth rate, wave speed, and the critical point." },
  E3B: { title: "Ring Spectrum with Implicit Bidirectional Acceleration Coupling", hint: "Enable leader- and follower-acceleration feedback. The theoretical criterion must inspect every discrete wavenumber, not only the long-wave limit." },
  E3M: { title: "Spatial Allocation of Multi-Leader and Multi-Follower Information", hint: "Enable multi-vehicle information, keep the total acceleration-kernel weight fixed, and compare finite-wavelength results across spatial decays." },
  E11: { title: "Waveform Steepening and Nonlinear Saturation under Finite Disturbances", hint: "The Nonlinear Phenomena page is open. Use its selector to compare steepening, saturation, solitary-wave or kink evidence, and triggered recovery." },
  E12: { title: "Finite-Braking Recovery in a Stable Regime", hint: "The Nonlinear Phenomena page is open. Select finite-amplitude recovery and inspect recovery time and residual oscillation." },
  M06: { title: "LWR Shock under the IDM Fundamental Diagram", hint: "M06 and W02 provide the finite-volume update. Use the laboratory to connect microscopic position-time diagrams with macroscopic conservation waves." },
  W02: { title: "Two Finite-Volume Updates on Four LWR Cells", hint: "W02 expands the fluxes and two time steps for four cells; compare them with the laboratory's macroscopic stability experiments." },
  W03: { title: "Conserved Variables, Fluxes, and Relaxation on Three ARZ Cells", hint: "W03 expands two ARZ updates for three cells; compare them with the macroscopic stability results in E13 and E14." },
  E13: { title: "Stable and Unstable Spectra of an IDM-Consistent ARZ Model", hint: "E13 presents the full ARZ spectral evolution and spatiotemporal diagram; use this laboratory to verify microscopic evidence at the shared IDM operating point." },
  E14: { title: "Comparison of IDM and ARZ Speed Fields", hint: "E14 uses the same IDM equilibrium and long-wave disturbance in both models; use the laboratory to reproduce the microscopic position-time trajectories." },
  D01: { title: "Natural Disturbance Propagation in NGSIM US-101", hint: "The Empirical Data page is open. Review episode screening, gain, delay, and their uncertainty." },
  D02: { title: "Stability Calibration, Holdout Validation, and Spectrum-Constrained Optimization", hint: "The Empirical Data page presents the D01 evidence; D02 then uses it for calibration, recursive validation, and counterfactual optimization." },
  A01: { title: "Stability Benefits of Cooperative AVs in an Open Platoon", hint: "Application A is open. Compare the HDV and cooperative-AV position-time diagrams and the selected-vehicle speed traces." },
  A02: { title: "Event-Triggered Fixed VSL at a Slow-Vehicle Bottleneck", hint: "Application B is open. The VSL activates only after the slow-vehicle event and selects 114 km/h using the stability margin." },
};

/* ============================================================
   主组件
   ============================================================ */
function App() {
  const linkedCaseId = typeof window === "undefined" ? "" : new URLSearchParams(window.location.search).get("case") || "";
  const linkedCase = BOOK_CASES[linkedCaseId] || null;
  const linkedCaseCopy = linkedCase && window.TSL_I18N?.language === "en" ? BOOK_CASES_EN[linkedCaseId] : linkedCase;
  const [P, setP] = useState({ v0: 33.3, T: 1.5, s0: 2.0, a: 1.0, b: 1.5, delta: 4 });
  const [N, setN] = useState(60);
  const [rho, setRho] = useState(28);      // veh/km
  const [fixedLength, setFixedLength] = useState(false);
  const [ringLength, setRingLength] = useState(2.14);
  const [l, setL] = useState(5);
  const [kappa, setKappa] = useState(0);
  const [fb, setFb] = useState(0.08);
  const [alpha, setAlpha] = useState(0.6);
  const [beta, setBeta] = useState(0.35);
  const [mFront, setMFront] = useState(3);
  const [mBack, setMBack] = useState(2);
  const [rearP, setRearP] = useState(0.15);
  const [tau0, setTau0] = useState(0);
  const [tauA, setTauA] = useState(0);
  const [pert, setPert] = useState(0.05);
  const [dt, setDt] = useState(0.02);
  const [speedMul, setSpeedMul] = useState(4);
  const [running, setRunning] = useState(false);
  const [useFF, setUseFF] = useState(false);
  const [useBackAcc, setUseBackAcc] = useState(false);
  const [useMulti, setUseMulti] = useState(false);
  const [useRearGap, setUseRearGap] = useState(false);
  const [useDelay, setUseDelay] = useState(false);
  const [useHetero, setUseHetero] = useState(false);
  const [useMixed, setUseMixed] = useState(false);
  const [heavyShare, setHeavyShare] = useState(0.25);
  const [avShare, setAvShare] = useState(0.3);
  const [arrangement, setArrangement] = useState("均匀");
  const [avGroups, setAvGroups] = useState(4);
  const [kappa0, setKappa0] = useState(0.15);
  const [selectedVehicle, setSelectedVehicle] = useState(0);
  const [fieldMode, setFieldMode] = useState("速度");
  const [activeTab, setActiveTab] = useState(() => linkedCase ? linkedCase.tab : "live");
  const [parameterModal, setParameterModal] = useState(null);
  const [feedbackModal, setFeedbackModal] = useState(null);
  const [nonlinearCase, setNonlinearCase] = useState(() => linkedCaseId === "E12" ? "recovery" : "steepening");
  const [nonlinearTimeIndex, setNonlinearTimeIndex] = useState(4);
  const [nonlinearAmplitude, setNonlinearAmplitude] = useState(1.0);
  const [applicationCase, setApplicationCase] = useState(() => linkedCase && linkedCase.application ? linkedCase.application : "open");
  const [applicationProbeIndex, setApplicationProbeIndex] = useState(3);
  const [ringFocus, setRingFocus] = useState(false);
  const [rampDuration, setRampDuration] = useState(120);
  const [appliedKappa, setAppliedKappa] = useState(0);
  const [appliedRuntime, setAppliedRuntime] = useState({
    fb: 0, alpha: 0, beta: 0, mFront: 1, mBack: 1,
    rearP: 0, tau0: 0, tauA: 0,
  });
  const [autoControl, setAutoControl] = useState(false);
  const [interventions, setInterventions] = useState([]);
  const [experiment, setExperiment] = useState("E1");
  const [experimentResult, setExperimentResult] = useState(null);
  const [avPlaneHover, setAvPlaneHover] = useState(null);
  const [nonce, setNonce] = useState(0);
  const [live, setLive] = useState({ t: 0, std: 0, gapStd: 0, vbar: 0, gmin: 0, gmax: 0 });
  const [meas, setMeas] = useState(null);

  const Lring = fixedLength ? ringLength * 1000 : (N / rho) * 1000;
  const rhoNow = (N / Lring) * 1000;
  const se = Lring / N - l;

  const eff = useMemo(
    () => ({
      fa: useFF ? kappa : 0,
      fb: useBackAcc ? fb : 0,
      alpha: useMulti ? alpha : 0,
      beta: useMulti ? beta : 0,
      mFront: useMulti ? mFront : 1,
      mBack: useMulti ? mBack : 1,
      rearP: useRearGap ? rearP : 0,
      tau0: useDelay ? tau0 : 0,
      tauA: useDelay ? tauA : 0,
    }),
    [useFF, kappa, useBackAcc, fb, useMulti, alpha, beta, mFront, mBack, useRearGap, rearP, useDelay, tau0, tauA]
  );

  /* ---- 理论量 ---- */
  const theory = useMemo(() => analyzeTheorySnapshot({
    se, P, eff, N, useHetero, useMixed, avShare, arrangement, avGroups, kappa0,
  }), [se, P, eff, N, useHetero, useMixed, avShare, arrangement, avGroups, kappa0]);

  const currentEff = useMemo(() => ({ ...appliedRuntime, fa: appliedKappa }), [appliedRuntime, appliedKappa]);
  const currentTheory = useMemo(() => analyzeTheorySnapshot({
    se, P, eff: currentEff, N, useHetero, useMixed, avShare, arrangement, avGroups, kappa0,
  }), [se, P, currentEff, N, useHetero, useMixed, avShare, arrangement, avGroups, kappa0]);

  const avPlane = useMemo(() => buildAvStabilityPlane({
    P, N, fa: useFF ? kappa : 0, kappa0, arrangement, avGroups,
    currentShare: avShare, currentSpeed: theory ? theory.ve : 0,
  }), [P, N, useFF, kappa, kappa0, arrangement, avGroups, avShare, theory]);

  /* ---- 仿真实例 ---- */
  const simRef = useRef(null);
  const speedSTCanvas = useRef(null);
  const gapSTCanvas = useRef(null);
  const trajectoryCanvas = useRef(null);
  const ringCanvas = useRef(null);
  const legacySTCanvas = useRef(null), legacyRingCanvas = useRef(null);
  const legacyVehicleC = useRef(null), legacyEnvelopeC = useRef(null);
  const nonlinearCanvas = useRef(null);
  const applicationBaseCanvas = useRef(null), applicationControlCanvas = useRef(null), applicationProbeCanvas = useRef(null), applicationVslCanvas = useRef(null);
  const logRef = useRef({ lnA: [], phase: [], lnStd: [], vehicle: [], envelope: [], lastPhase: 0, wraps: 0, trajLast: null });
  const colRef = useRef(0);
  const trajRowRef = useRef(0);
  const frameRef = useRef(0);
  const rampRef = useRef({ active: false, from: 0, to: 0, start: 0, end: 0 });
  const autoControlRef = useRef(false);
  const autoTriggeredRef = useRef(false);
  const skipAutoResetRef = useRef(false);

  const resetSim = useCallback(() => {
    if (se <= P.s0 + 0.01) { simRef.current = null; return; }
    simRef.current = createSim({
      N, L: Lring, l, P, ...eff, dt, pert,
      hetero: useHetero, mixed: useMixed, heavyShare, avShare, arrangement, avGroups, kappa0,
    });
    logRef.current = { lnA: [], phase: [], lnStd: [], vehicle: [], envelope: Array(N).fill(0), lastPhase: null, wraps: 0, trajLast: null };
    colRef.current = 0;
    trajRowRef.current = 0;
    frameRef.current = 0;
    rampRef.current = { active: false, from: eff.fa, to: eff.fa, start: 0, end: 0 };
    autoTriggeredRef.current = false;
    setAppliedKappa(eff.fa);
    setAppliedRuntime({
      fb: eff.fb, alpha: eff.alpha, beta: eff.beta,
      mFront: eff.mFront, mBack: eff.mBack, rearP: eff.rearP,
      tau0: eff.tau0, tauA: eff.tauA,
    });
    setInterventions([]);
    setMeas(null);
    setLive({ t: 0, std: 0, gapStd: 0, vbar: simRef.current.ve, gmin: simRef.current.se, gmax: simRef.current.se });
    for (const ref of [speedSTCanvas, gapSTCanvas, trajectoryCanvas]) {
      const c = ref.current;
      if (!c) continue;
      const ctx = c.getContext("2d");
      ctx.fillStyle = CLR.panel;
      ctx.fillRect(0, 0, c.width, c.height);
    }
    setNonce((x) => x + 1);
  }, [N, Lring, l, P, eff, dt, pert, se, useHetero, useMixed, heavyShare, avShare, arrangement, avGroups, kappa0]);

  useEffect(() => {
    if (skipAutoResetRef.current) { skipAutoResetRef.current = false; return; }
    resetSim();
    /* eslint-disable-next-line */
  }, [N, Lring, l, P, dt, pert, useHetero, useMixed, heavyShare, avShare, arrangement, avGroups, kappa0, tau0, tauA, useDelay]);
  useEffect(() => { if (selectedVehicle >= N) setSelectedVehicle(N - 1); }, [N, selectedVehicle]);
  useEffect(() => { autoControlRef.current = autoControl; }, [autoControl]);

  const markIntervention = useCallback(() => {
    for (const ref of [speedSTCanvas, gapSTCanvas]) {
      const c = ref.current;
      if (!c) continue;
      const ctx = c.getContext("2d"), x = Math.min(colRef.current, c.width - 1);
      ctx.fillStyle = CLR.amber;
      ctx.fillRect(x, 0, 2, c.height);
    }
    const c = trajectoryCanvas.current;
    if (c) {
      const ctx = c.getContext("2d"), y = Math.min(trajRowRef.current, c.height - 1);
      ctx.fillStyle = CLR.amber;
      ctx.fillRect(0, y, c.width, 2);
    }
  }, []);

  const beginRamp = useCallback((sim, target, duration, source = "手动") => {
    if (!sim) return;
    const to = Math.max(0, Math.min(0.95, target));
    const dur = Math.max(0, duration);
    const from = sim.fa;
    sim.applyRuntime({
      fb: useBackAcc ? fb : 0,
      alpha: useMulti ? alpha : 0,
      beta: useMulti ? beta : 0,
      mFront: useMulti ? mFront : 1,
      mBack: useMulti ? mBack : 1,
      rearP: useRearGap ? rearP : 0,
      kappa0,
    });
    setAppliedRuntime({
      fb: useBackAcc ? fb : 0,
      alpha: useMulti ? alpha : 0,
      beta: useMulti ? beta : 0,
      mFront: useMulti ? mFront : 1,
      mBack: useMulti ? mBack : 1,
      rearP: useRearGap ? rearP : 0,
      tau0: sim.tau0,
      tauA: sim.tauA,
    });
    rampRef.current = { active: dur > 0.001, from, to, start: sim.t, end: sim.t + dur };
    if (dur <= 0.001) sim.fa = to;
    setAppliedKappa(sim.fa);
    setUseFF(to > 0);
    setInterventions((items) => [...items.slice(-5), { t: sim.t, from, to, duration: dur, source }]);
    if (source === "自动") autoTriggeredRef.current = true;
    markIntervention();
  }, [useBackAcc, fb, useMulti, alpha, beta, mFront, mBack, useRearGap, rearP, kappa0, markIntervention]);

  const applyOnline = useCallback(() => {
    beginRamp(simRef.current, useFF ? kappa : 0, rampDuration, "手动");
  }, [beginRamp, useFF, kappa, rampDuration]);

  const prepareJamControlDemo = useCallback(() => {
    setRunning(false);
    const demoP = { v0: 33.3, T: 1.5, s0: 2.0, a: 1.0, b: 1.5, delta: 4 };
    const demoN = 60, demoRho = 50, demoL = (demoN / demoRho) * 1000, demoPert = -2.5, demoDt = 0.02;
    skipAutoResetRef.current = true;
    setP(demoP); setN(demoN); setRho(demoRho); setFixedLength(false); setL(5); setPert(0.4); setDt(demoDt);
    setUseBackAcc(false); setUseMulti(false); setUseRearGap(false); setUseDelay(false); setUseHetero(false); setUseMixed(false);
    setUseFF(true);
    setKappa(0.45);
    setRampDuration(120);
    setAutoControl(false);
    const sim = createSim({
      N: demoN, L: demoL, l: 5, P: demoP,
      fa: 0, fb: 0, alpha: 0, beta: 0, mFront: 1, mBack: 1, rearP: 0,
      tau0: 0, tauA: 0, dt: demoDt, pert: demoPert,
      hetero: false, mixed: false, heavyShare: 0, avShare: 0, arrangement: "均匀", avGroups: 1, kappa0: 0,
    });
    simRef.current = sim;
    logRef.current = { lnA: [], phase: [], lnStd: [], vehicle: [], envelope: Array(demoN).fill(0), lastPhase: null, wraps: 0, trajLast: null };
    colRef.current = 0;
    trajRowRef.current = 0;
    frameRef.current = 0;
    rampRef.current = { active: false, from: 0, to: 0, start: 0, end: 0 };
    autoTriggeredRef.current = false;
    setAppliedKappa(0);
    setAppliedRuntime({ fb: 0, alpha: 0, beta: 0, mFront: 1, mBack: 1, rearP: 0, tau0: 0, tauA: 0 });
    setInterventions([]);
    setMeas(null);
    setLive({ t: 0, std: sim.velStd(), gapStd: sim.gapStd(), vbar: sim.ve, gmin: sim.se, gmax: sim.se });
    for (const ref of [speedSTCanvas, gapSTCanvas, trajectoryCanvas]) {
      const c = ref.current;
      if (!c) continue;
      const ctx = c.getContext("2d");
      ctx.fillStyle = CLR.panel;
      ctx.fillRect(0, 0, c.width, c.height);
    }
    setNonce((x) => x + 1);
    setActiveTab("live");
  }, []);

  /* ---- 主循环 ---- */
  useEffect(() => {
    if (!running) return;
    let raf;
    const loop = () => {
      const sim = simRef.current;
      if (sim && theory) {
        const stepsPerFrame = Math.max(1, Math.round((speedMul * 0.05) / dt));
        const sampleEvery = Math.max(1, Math.round((1800 / 960) / dt));
        for (let i = 0; i < stepsPerFrame; i++) {
          if (autoControlRef.current && !autoTriggeredRef.current && sim.t >= 600) {
            beginRamp(sim, Math.max(0.45, kappa), rampDuration, "自动");
          }
          const ramp = rampRef.current;
          if (ramp.active) {
            if (sim.t >= ramp.end) {
              sim.fa = ramp.to;
              ramp.active = false;
            } else {
              const q = Math.max(0, Math.min(1, (sim.t - ramp.start) / Math.max(1e-9, ramp.end - ramp.start)));
              sim.fa = ramp.from + (ramp.to - ramp.from) * q;
            }
          }
          sim.step();
          if (Math.round(sim.t / dt) % sampleEvery === 0) {
            const m = theory.best.m || 1;
            const A = sim.modeAmp(m);
            const mag = Math.hypot(A.re, A.im);
            const lg = logRef.current;
            if (mag > 1e-14) {
              lg.lnA.push([sim.t, Math.log(mag)]);
              let ph = Math.atan2(A.im, A.re);
              if (lg.lastPhase !== null) {
                let d = ph - lg.lastPhase;
                while (d > Math.PI) { d -= 2 * Math.PI; lg.wraps -= 1; }
                while (d < -Math.PI) { d += 2 * Math.PI; lg.wraps += 1; }
              }
              lg.lastPhase = ph;
              lg.phase.push([sim.t, ph + 2 * Math.PI * lg.wraps]);
            }
            logRef.current.lnStd.push([sim.t, Math.log(Math.max(sim.velStd(), 1e-14))]);
            const gaps = sim.gaps();
            const sv = Math.min(sim.N - 1, selectedVehicle);
            lg.vehicle.push([sim.t, sim.v[sv], gaps[sv]]);
            if (lg.vehicle.length > 1200) lg.vehicle.shift();
            for (let n = 0; n < sim.N; n++) lg.envelope[n] = Math.max(lg.envelope[n] || 0, Math.abs(sim.v[n] - sim.ve));
            drawSTColumns(sim);
            drawTrajectorySample(sim);
          }
        }
        drawRing(sim);
        frameRef.current = (frameRef.current + 1) % 5;
        if (frameRef.current === 0) {
          const g = sim.gaps();
          let gmin = 1e9, gmax = -1e9, vb = 0;
          for (let n = 0; n < sim.N; n++) { gmin = Math.min(gmin, g[n]); gmax = Math.max(gmax, g[n]); vb += sim.v[n]; }
          setAppliedKappa(sim.fa);
          setLive({ t: sim.t, std: sim.velStd(), gapStd: sim.gapStd(), vbar: vb / sim.N, gmin, gmax });
          computeMeasured();
        }
      }
      raf = requestAnimationFrame(loop);
    };
    raf = requestAnimationFrame(loop);
    return () => cancelAnimationFrame(raf);
    /* eslint-disable-next-line */
  }, [running, theory, speedMul, dt, selectedVehicle, beginRamp, kappa, rampDuration]);

  /* ---- 测量：从对数幅值与相位拟合 ---- */
  const computeMeasured = useCallback(() => {
    const lg = logRef.current;
    const sim = simRef.current;
    if (!sim || !theory || lg.lnA.length < 20) return;
    // 拟合窗口：跳过初始瞬态（单车脉冲先在各模态间重分配），停在饱和之前
    const satLevel = Math.log(Math.max(0.02 * theory.ve, 0.05));
    let iSat = lg.lnA.length - 1;
    for (let i = 0; i < lg.lnA.length; i++) {
      if (lg.lnA[i][1] > satLevel) { iSat = i; break; }
    }
    const tEnd = lg.lnA[iSat][0];
    const tStart = Math.max(0.28 * tEnd, 15);
    let i0 = 0;
    while (i0 < iSat && lg.lnA[i0][0] < tStart) i0++;
    if (iSat - i0 < 12) return;
    const win = lg.lnA.slice(i0, iSat + 1);
    const winP = lg.phase.slice(i0, iSat + 1);
    const gr = lsSlope(win);
    const om = lsSlope(winP);
    if (gr === null) return;
    // 空间主导模态
    let bm = 1, bmag = -1;
    for (let m = 1; m <= Math.floor(sim.N / 2); m++) {
      const A = sim.modeAmp(m);
      const mg = Math.hypot(A.re, A.im);
      if (mg > bmag) { bmag = mg; bm = m; }
    }
    const k = (2 * Math.PI * (theory.best.m || 1)) / sim.N;
    setMeas({ growth: gr, omega: om, c1: om !== null ? -om / k : null, domM: bm, nWin: win.length });
  }, [theory, pert]);

  /* ---- 1800 s 速度/间距时空图 ---- */
  const drawSTColumns = (sim) => {
    const gaps = sim.gaps();
    const col = colRef.current;
    const drawOne = (ref, field) => {
      const c = ref.current;
      if (!c) return;
      const ctx = c.getContext("2d"), W = c.width, H = c.height;
      if (col >= W) {
        ctx.drawImage(c, -1, 0);
        ctx.fillStyle = CLR.panel;
        ctx.fillRect(W - 1, 0, 1, H);
      }
      const cx = Math.min(col, W - 1), rowH = H / sim.N;
      const span = field === "gap" ? Math.max(1, sim.se * 0.35) : Math.max(0.6, sim.ve * 0.25);
      for (let n = 0; n < sim.N; n++) {
        const eqGap = sim.models ? idmGapAtSpeed(sim.ve, sim.models[n]) : sim.se;
        const value = field === "gap" ? gaps[n] - eqGap : sim.v[n] - sim.ve;
        ctx.fillStyle = speedColor(value, span);
        ctx.fillRect(cx, n * rowH, 1, Math.ceil(rowH));
      }
    };
    drawOne(speedSTCanvas, "speed");
    drawOne(gapSTCanvas, "gap");
    colRef.current = col + 1;
  };

  /* ---- 位置–时间轨迹；轨迹颜色映射车辆速度 ---- */
  const drawTrajectorySample = (sim) => {
    if (colRef.current % 2 !== 0) return;
    const c = trajectoryCanvas.current;
    if (!c) return;
    const ctx = c.getContext("2d"), W = c.width, H = c.height;
    let row = trajRowRef.current;
    if (row >= H) {
      ctx.drawImage(c, 0, -1);
      ctx.fillStyle = CLR.panel;
      ctx.fillRect(0, H - 1, W, 1);
      row = H - 1;
    }
    const xs = Array.from(sim.x, (x) => ((x % sim.L) + sim.L) % sim.L / sim.L * (W - 1));
    const last = logRef.current.trajLast;
    const span = Math.max(0.6, sim.ve * 0.25);
    for (let n = 0; n < sim.N; n++) {
      ctx.strokeStyle = speedColor(sim.v[n] - sim.ve, span);
      ctx.fillStyle = ctx.strokeStyle;
      ctx.lineWidth = 1;
      if (last && Math.abs(xs[n] - last[n]) < W * 0.45) {
        ctx.beginPath();
        ctx.moveTo(last[n], Math.max(0, row - 1));
        ctx.lineTo(xs[n], row);
        ctx.stroke();
      } else {
        ctx.fillRect(xs[n], row, 1.5, 1.5);
      }
    }
    logRef.current.trajLast = xs;
    trajRowRef.current += 1;
  };

  /* ---- 环道动画 ---- */
  const drawRing = (sim) => {
    const c = ringCanvas.current;
    if (!c) return;
    const ctx = c.getContext("2d");
    const W = c.width, H = c.height;
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = CLR.panel;
    ctx.fillRect(0, 0, W, H);
    const cx = W / 2, cy = H / 2, R = Math.min(W, H) / 2 - 46;
    ctx.strokeStyle = CLR.rule; ctx.lineWidth = Math.max(18, Math.min(W, H) * 0.055);
    ctx.beginPath(); ctx.arc(cx, cy, R, 0, 2 * Math.PI); ctx.stroke();
    ctx.strokeStyle = "rgba(21,30,39,0.22)"; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.arc(cx, cy, R - Math.min(W, H) * 0.028, 0, 2 * Math.PI); ctx.stroke();
    ctx.beginPath(); ctx.arc(cx, cy, R + Math.min(W, H) * 0.028, 0, 2 * Math.PI); ctx.stroke();
    const span = Math.max(0.6, sim.ve * 0.25);
    const carW = Math.max(5, Math.min(9, W / 85));
    const carL = Math.max(10, Math.min(18, H / 38));
    for (let n = 0; n < sim.N; n++) {
      const th = ((sim.x[n] % sim.L) / sim.L) * 2 * Math.PI - Math.PI / 2;
      ctx.save();
      ctx.translate(cx + R * Math.cos(th), cy + R * Math.sin(th));
      ctx.rotate(th + Math.PI / 2);
      ctx.fillStyle = speedColor(sim.v[n] - sim.ve, span);
      ctx.fillRect(-carW / 2, -carL / 2, carW, carL);
      ctx.strokeStyle = sim.av && sim.av[n] ? CLR.teal : sim.heavy && sim.heavy[n] ? CLR.amber : "rgba(21,30,39,0.35)";
      ctx.lineWidth = sim.av && sim.av[n] ? 1.6 : sim.heavy && sim.heavy[n] ? 1.3 : 0.6;
      ctx.strokeRect(-carW / 2, -carL / 2, carW, carL);
      ctx.restore();
    }
    ctx.fillStyle = CLR.soft;
    ctx.font = `${Math.max(11, Math.round(W / 52))}px ui-monospace, monospace`;
    ctx.textAlign = "center";
    ctx.fillText(`${(sim.L / 1000).toFixed(2)} km · ${sim.N} 辆`, cx, cy - 16);
    ctx.fillText(`t = ${sim.t.toFixed(1)} s`, cx, cy + 5);
    ctx.fillText(`κ = ${sim.fa.toFixed(3)} · σᵥ = ${sim.velStd().toFixed(3)}`, cx, cy + 26);
  };

  useEffect(() => {
    if (simRef.current) drawRing(simRef.current);
    /* eslint-disable-next-line */
  }, [activeTab, nonce, ringFocus]);

  /* ---- 静态图重绘 ---- */
  const psiC = useRef(null), specC = useRef(null), bandC = useRef(null), growC = useRef(null), avPlaneC = useRef(null);
  const vehicleC = useRef(null), envelopeC = useRef(null);
  useEffect(() => {
    if (!theory) return;
    // Psi(omega)
    {
      const c = psiC.current; if (c) {
        const ctx = c.getContext("2d");
        const psiMax = Math.max(...theory.psis);
        const yl = Math.min(theory.mixedSimple ? -0.005 : -0.02, theory.psiMin * 1.3);
        const yh = Math.max(theory.mixedSimple ? 0.005 : 0.05, psiMax * 1.2, theory.mixedSimple ? 0 : theory.psi0 * 2.2);
        plot(ctx, c.width, c.height, {
          xs: theory.ws, xr: [0, 6], yr: [yl, yh],
          series: [{ ys: theory.psis, color: CLR.blue, fillNeg: true }],
          xlab: "ω [rad/s]", ylab: theory.mixedSimple ? "−mean ln|Gₙ|" : "Ψ(ω)", hline: 0,
        });
      }
    }
    // Re lambda(k)
    {
      const c = specC.current; if (c) {
        const ctx = c.getContext("2d");
        const ks = theory.spec.map((s) => s.k);
        const rs = theory.spec.map((s) => s.lam.re);
        const mx = Math.max(1e-4, ...rs), mn = Math.min(-1e-3, ...rs.map((r) => r));
        plot(ctx, c.width, c.height, {
          xs: ks, xr: [0, Math.PI], yr: [Math.max(mn, -0.06), Math.max(mx * 1.4, 0.005)],
          series: [{ ys: rs, color: CLR.brick }],
          xlab: "k [rad/veh]", ylab: "Re λ [1/s]", hline: 0,
          marks: theory.best.lam.re > 0 ? [{ x: theory.best.k, color: CLR.teal, label: `m*=${theory.best.m}` }] : [],
        });
      }
    }
    // 不稳定带
    {
      const c = bandC.current; if (c) {
        const ctx = c.getContext("2d");
        const xs = theory.band.map((d) => d.v);
        const ys = theory.band.map((d) => d.val);
        const mn = Math.min(...ys), mx = Math.max(...ys);
        plot(ctx, c.width, c.height, {
          xs, xr: [0, P.v0], yr: [Math.max(mn * 1.2, -0.25), Math.min(mx * 1.1, 0.45)],
          series: [{ ys, color: CLR.teal, fillNeg: true }],
          xlab: "平衡速度 vₑ [m/s]", ylab: theory.mixedSimple ? "长波裕度" : "min Ψ", hline: 0,
          marks: [{ x: theory.ve, color: CLR.ink, label: `vₑ=${theory.ve.toFixed(2)}` }],
        });
      }
    }
    // 增长曲线
    drawGrowth();
    /* eslint-disable-next-line */
  }, [theory, nonce, activeTab]);

  useEffect(() => {
    const c = avPlaneC.current;
    if (!c || !avPlane) return;
    const ctx = c.getContext("2d"), W = c.width, H = c.height;
    const pad = { l: 62, r: 20, t: 20, b: 42 };
    const iw = W - pad.l - pad.r, ih = H - pad.t - pad.b;
    const cw = iw / avPlane.nx, ch = ih / avPlane.ny;
    const blend = (a, b, q) => Math.round(a + (b - a) * q);
    ctx.clearRect(0, 0, W, H);
    ctx.fillStyle = CLR.panel;
    ctx.fillRect(0, 0, W, H);
    for (let iy = 0; iy < avPlane.ny; iy++) {
      for (let ix = 0; ix < avPlane.nx; ix++) {
        const cell = avPlane.cells[iy * avPlane.nx + ix];
        const q = 0.28 + 0.72 * Math.min(1, Math.sqrt(Math.abs(cell.score) / 0.08));
        const target = cell.stable ? [30, 110, 91] : [168, 58, 43];
        ctx.fillStyle = `rgb(${blend(242, target[0], q)},${blend(245, target[1], q)},${blend(240, target[2], q)})`;
        ctx.fillRect(pad.l + ix * cw, pad.t + ih - (iy + 1) * ch, Math.ceil(cw + 0.5), Math.ceil(ch + 0.5));
      }
    }
    // 稳定边界：相邻网格判定改变的位置。
    ctx.strokeStyle = "rgba(21,30,39,0.86)";
    ctx.lineWidth = 1.15;
    ctx.beginPath();
    for (let iy = 0; iy < avPlane.ny; iy++) {
      for (let ix = 0; ix < avPlane.nx; ix++) {
        const here = avPlane.cells[iy * avPlane.nx + ix].stable;
        const x = pad.l + ix * cw, y = pad.t + ih - (iy + 1) * ch;
        if (ix + 1 < avPlane.nx && here !== avPlane.cells[iy * avPlane.nx + ix + 1].stable) {
          ctx.moveTo(x + cw, y); ctx.lineTo(x + cw, y + ch);
        }
        if (iy + 1 < avPlane.ny && here !== avPlane.cells[(iy + 1) * avPlane.nx + ix].stable) {
          ctx.moveTo(x, y); ctx.lineTo(x + cw, y);
        }
      }
    }
    ctx.stroke();
    ctx.strokeStyle = CLR.ink;
    ctx.lineWidth = 1;
    ctx.strokeRect(pad.l, pad.t, iw, ih);
    ctx.font = "11px ui-monospace, monospace";
    ctx.fillStyle = CLR.soft;
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (let i = 0; i <= 5; i++) {
      const x = pad.l + iw * i / 5;
      ctx.fillText(`${i * 20}%`, x, pad.t + ih + 7);
    }
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (let i = 0; i <= 5; i++) {
      const value = avPlane.vMax * i / 5, y = pad.t + ih - ih * i / 5;
      ctx.fillText(value.toFixed(1), pad.l - 7, y);
    }
    ctx.textAlign = "right"; ctx.textBaseline = "bottom";
    ctx.fillText("AV 渗透率", W - pad.r, H - 3);
    ctx.save();
    ctx.translate(11, pad.t);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "right"; ctx.textBaseline = "top";
    ctx.fillText("平衡速度 vₑ [m/s]", 0, 0);
    ctx.restore();
    // 当前工作点。
    if (theory) {
      const x = pad.l + avShare * iw;
      const y = pad.t + ih - (theory.ve / avPlane.vMax) * ih;
      ctx.beginPath(); ctx.arc(x, y, 5.5, 0, 2 * Math.PI);
      ctx.fillStyle = "#FFFFFF"; ctx.fill();
      ctx.strokeStyle = CLR.ink; ctx.lineWidth = 2; ctx.stroke();
    }
    // 紧凑图例。
    ctx.textAlign = "left"; ctx.textBaseline = "middle"; ctx.font = "10.5px ui-monospace, monospace";
    ctx.fillStyle = CLR.teal; ctx.fillRect(pad.l + 8, pad.t + 8, 11, 11);
    ctx.fillStyle = CLR.ink; ctx.fillText("稳定", pad.l + 24, pad.t + 13.5);
    ctx.fillStyle = CLR.brick; ctx.fillRect(pad.l + 70, pad.t + 8, 11, 11);
    ctx.fillStyle = CLR.ink; ctx.fillText("失稳", pad.l + 86, pad.t + 13.5);
  }, [avPlane, avShare, theory, activeTab, experiment, useMixed]);

  /* ---- E11–E12 五类非线性现象 ---- */
  useEffect(() => {
    const c = nonlinearCanvas.current;
    if (!c || activeTab !== "nonlinear") return;
    const ctx = c.getContext("2d");
    if (nonlinearCase === "steepening") {
      const p0 = NONLINEAR_DATA.profiles[0];
      const p1 = NONLINEAR_DATA.profiles[nonlinearTimeIndex];
      const xs = Array.from({ length: 120 }, (_, i) => i);
      const expand = (p) => xs.map((i) => p[i % 30] - 8);
      const ys0 = expand(p0), ys1 = expand(p1);
      const span = Math.max(0.06, ...ys1.map(Math.abs)) * 1.12;
      plot(ctx, c.width, c.height, {
        xs, xr: [0, 119], yr: [-span, span],
        series: [{ ys: ys0, color: CLR.soft, w: 1.1 }, { ys: ys1, color: CLR.brick, w: 2.0 }],
        xlab: "车辆编号 n", ylab: "vₙ − v̄ [m/s]", hline: 0,
      });
    } else if (nonlinearCase === "saturation") {
      plot(ctx, c.width, c.height, {
        xs: NONLINEAR_DATA.time, xr: [0, 1800], yr: [-8, 1],
        series: [
          { ys: NONLINEAR_DATA.unstableStd.map((x) => Math.log10(Math.max(x, 1e-8))), color: CLR.brick, w: 2 },
          { ys: NONLINEAR_DATA.stableStd.map((x) => Math.log10(Math.max(x, 1e-8))), color: CLR.blue, w: 1.6 },
        ],
        xlab: "时间 t [s]", ylab: "log₁₀ σᵥ", marks: [{ x: 1500, color: CLR.amber, label: "末 300 s" }],
      });
    } else if (nonlinearCase === "soliton") {
      const xs = Array.from({ length: 321 }, (_, i) => -8 + i * 0.05);
      const A = nonlinearAmplitude, width = 2.2 / Math.sqrt(A);
      const kdv = xs.map((x) => A / Math.pow(Math.cosh(x / width), 2));
      const kink = xs.map((x) => Math.tanh(x / (1.4 / A)));
      plot(ctx, c.width, c.height, {
        xs, xr: [-8, 8], yr: [-1.15, Math.max(1.15, A * 1.12)],
        series: [{ ys: kdv, color: CLR.blue, w: 2 }, { ys: kink, color: CLR.amber, w: 1.7 }],
        xlab: "慢坐标 X", ylab: "归一化 R", hline: 0,
      });
    } else if (nonlinearCase === "trigger") {
      plot(ctx, c.width, c.height, {
        xs: NONLINEAR_DATA.pulseAmplitude, xr: [0, 18], yr: [0, 4.2],
        series: [
          { ys: NONLINEAR_DATA.pulsePeak, color: CLR.soft, w: 1.4 },
          { ys: NONLINEAR_DATA.pulseLate, color: CLR.brick, w: 2 },
        ],
        xlab: "单车制动幅值 A₀ [m/s]", ylab: "σᵥ [m/s]", hline: 0.02,
        marks: [{ x: 8, color: CLR.amber, label: "长暂态区" }],
      });
    } else if (nonlinearCase === "recovery") {
      plot(ctx, c.width, c.height, {
        xs: NONLINEAR_DATA.pulseAmplitude, xr: [0, 18], yr: [0, 5700],
        series: [{ ys: NONLINEAR_DATA.recoveryTime, color: CLR.teal, w: 2 }],
        xlab: "单车制动幅值 A₀ [m/s]", ylab: "恢复时间 tᵣ [s]",
        marks: [{ x: 8, color: CLR.brick, label: "≥8 m/s：右删失" }],
      });
    }
  }, [activeTab, nonlinearCase, nonlinearTimeIndex, nonlinearAmplitude]);

  /* ---- 交通管理应用：开放车队与移动瓶颈 ---- */
  useEffect(() => {
    if (activeTab !== "applications" || !APPLICATION_DATA) return;
    const pair = APPLICATION_DATA[applicationCase];
    if (!pair) return;
    const positionMin = Math.min(pair.baseline.positionMin, pair.control.positionMin);
    const positionMax = Math.max(pair.baseline.positionMax, pair.control.positionMax);
    const drawTrajectory = (ref, data, controlled) => {
      const c = ref.current; if (!c) return;
      const ctx = c.getContext("2d"), W = c.width, H = c.height;
      const P0 = { l: 55, r: 12, t: 12, b: 29 };
      const iw = W - P0.l - P0.r, ih = H - P0.t - P0.b;
      const q = decodeQ8(data.speedQ8), aq = decodeQ8(data.activeQ8), xq = decodeQ16(data.positionQ16), nt = data.heatTime.length, n = data.n;
      ctx.clearRect(0, 0, W, H); ctx.fillStyle = CLR.panel; ctx.fillRect(0, 0, W, H);
      const decodePosition = (it, j) => data.positionMin + (data.positionMax - data.positionMin) * xq[it * n + j] / 65535;
      const tmax = data.heatTime[data.heatTime.length - 1];
      if (controlled && isFinite(data.events?.zone_start) && isFinite(data.events?.zone_end)) {
        const t0 = data.events.control_start, t1 = data.events.control_end;
        const x0 = P0.l + t0 / tmax * iw, x1 = P0.l + t1 / tmax * iw;
        const y0 = P0.t + ih - (data.events.zone_end - positionMin) / (positionMax - positionMin) * ih;
        const y1 = P0.t + ih - (data.events.zone_start - positionMin) / (positionMax - positionMin) * ih;
        ctx.fillStyle = CLR.teal; ctx.globalAlpha = 0.12;
        ctx.fillRect(x0, y0, x1 - x0, y1 - y0);
        ctx.globalAlpha = 1;
      }
      ctx.lineWidth = 0.78;
      for (let j = 0; j < n; j++) {
        for (let it = 0; it < nt - 1; it++) {
          if (aq[it * n + j] < 128 || aq[(it + 1) * n + j] < 128) continue;
          const v = data.vmin + (data.vmax - data.vmin) * (q[it * n + j] + q[(it + 1) * n + j]) / 510;
          const x0 = P0.l + data.heatTime[it] / tmax * iw, x1 = P0.l + data.heatTime[it + 1] / tmax * iw;
          const p0 = decodePosition(it, j), p1 = decodePosition(it + 1, j);
          const y0 = P0.t + ih - (p0 - positionMin) / (positionMax - positionMin) * ih;
          const y1 = P0.t + ih - (p1 - positionMin) / (positionMax - positionMin) * ih;
          ctx.strokeStyle = speedColor(v - 0.5 * (data.vmin + data.vmax), 0.5 * (data.vmax - data.vmin));
          ctx.beginPath(); ctx.moveTo(x0, y0); ctx.lineTo(x1, y1); ctx.stroke();
        }
      }
      const eventEntries = Object.entries(data.events || {});
      for (const [key, value] of eventEntries) {
        if (!isFinite(value) || key.startsWith("zone_") || (!key.includes("start") && !key.includes("end")) || value > tmax) continue;
        const x = P0.l + value / tmax * iw;
        ctx.strokeStyle = key.includes("control") ? CLR.teal : CLR.amber;
        ctx.setLineDash(key.includes("end") ? [4, 3] : []); ctx.lineWidth = 1.4;
        ctx.beginPath(); ctx.moveTo(x, P0.t); ctx.lineTo(x, P0.t + ih); ctx.stroke(); ctx.setLineDash([]);
      }
      if (Array.isArray(data.queueTail)) {
        ctx.strokeStyle = CLR.ink; ctx.lineWidth = 1.35; ctx.setLineDash([6, 4]); ctx.beginPath();
        let started = false;
        for (let i = 0; i < data.queueTail.length; i++) {
          const tail = data.queueTail[i];
          if (!isFinite(tail)) { started = false; continue; }
          const x = P0.l + data.time[i] / tmax * iw;
          const y = P0.t + ih - (tail - positionMin) / (positionMax - positionMin) * ih;
          if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
        }
        ctx.stroke(); ctx.setLineDash([]);
        ctx.fillStyle = CLR.ink; ctx.font = "10px ui-monospace, monospace"; ctx.textAlign = "right"; ctx.textBaseline = "top";
        ctx.fillText(`排队尾部 · 形成段 cᵩ=${data.metrics.queue_tail_wave_speed_mps.toFixed(2)} m/s`, P0.l + iw - 6, P0.t + 5);
      }
      ctx.strokeStyle = CLR.ink; ctx.lineWidth = 1; ctx.strokeRect(P0.l, P0.t, iw, ih);
      ctx.fillStyle = CLR.soft; ctx.font = "10px ui-monospace, monospace";
      ctx.textAlign = "center"; ctx.textBaseline = "top";
      for (let k = 0; k <= 4; k++) { const x = P0.l + k * iw / 4; ctx.fillText(String(Math.round(k * tmax / 4)), x, P0.t + ih + 5); }
      ctx.textAlign = "right"; ctx.textBaseline = "middle";
      for (let k = 0; k <= 4; k++) { const y = P0.t + ih - k * ih / 4; const km = (positionMin + k * (positionMax - positionMin) / 4) / 1000; ctx.fillText(km.toFixed(km >= 10 ? 0 : 1), P0.l - 5, y); }
      ctx.textAlign = "right"; ctx.textBaseline = "bottom"; ctx.fillText("时间 t [s]", W - P0.r, H - 2);
      ctx.save(); ctx.translate(9, P0.t); ctx.rotate(-Math.PI / 2); ctx.textAlign = "right"; ctx.textBaseline = "top"; ctx.fillText("道路位置 x [km]", 0, 0); ctx.restore();
      ctx.fillStyle = controlled ? CLR.teal : CLR.brick; ctx.textAlign = "left"; ctx.textBaseline = "top";
      ctx.fillText(controlled ? "控制" : "无控制", P0.l + 6, P0.t + 5);
    };
    drawTrajectory(applicationBaseCanvas, pair.baseline, false);
    drawTrajectory(applicationControlCanvas, pair.control, true);

    const c = applicationProbeCanvas.current;
    if (c) {
      const b = pair.baseline, u = pair.control;
      const idx = Math.min(applicationProbeIndex, b.probeIds.length - 1);
      const ymax = Math.max(...b.probes[idx], ...u.probes[idx]);
      const ymin = Math.min(...b.probes[idx], ...u.probes[idx]);
      const pad = Math.max(0.3, (ymax - ymin) * 0.12);
      const marks = [];
      for (const [key, value] of Object.entries(b.events || {})) {
        if (isFinite(value) && !key.includes("control") && !key.startsWith("zone_") && (key.includes("start") || key.includes("end"))) marks.push({ x: value, color: CLR.amber, label: key.includes("start") ? "扰动开始" : "扰动结束" });
      }
      if (applicationCase === "bottleneck") {
        marks.push({ x: u.events.control_start, color: CLR.teal, label: "VSL 开始" });
      }
      plot(c.getContext("2d"), c.width, c.height, {
        xs: b.time, xr: [0, b.time[b.time.length - 1]], yr: [ymin - pad, ymax + pad],
        series: [{ ys: b.probes[idx], color: CLR.brick, w: 1.5 }, { ys: u.probes[idx], color: CLR.teal, w: 2 }],
        xlab: "时间 t [s]", ylab: "速度 v [m/s]", marks,
      });
    }

    if (applicationCase === "bottleneck" && applicationVslCanvas.current) {
      const design = APPLICATION_DATA.metrics.case2.vsl_design;
      const rows = design.scan;
      plot(applicationVslCanvas.current.getContext("2d"), applicationVslCanvas.current.width, applicationVslCanvas.current.height, {
        xs: rows.map(r => r.limit_kph), xr: design.admissible_range, yr: [-0.0017, 0.0027],
        series: [
          { ys: rows.map(r => r.phi), color: CLR.blue, w: 2 },
          { ys: rows.map(() => 0), color: CLR.brick, w: 1.2 },
        ],
        xlab: "候选 VSL [km/h]", ylab: "长波稳定裕度 Φ [s⁻²]",
        marks: [{ x: design.selected.limit_kph, color: CLR.teal, label: `推荐 ${design.selected.limit_kph.toFixed(0)} km/h` }],
      });
    }
  }, [activeTab, applicationCase, applicationProbeIndex]);

  const inspectAvPlane = useCallback((e) => {
    const c = e.currentTarget, rect = c.getBoundingClientRect();
    const x = (e.clientX - rect.left) * c.width / rect.width;
    const y = (e.clientY - rect.top) * c.height / rect.height;
    const pad = { l: 62, r: 20, t: 20, b: 42 };
    const iw = c.width - pad.l - pad.r, ih = c.height - pad.t - pad.b;
    if (x < pad.l || x > pad.l + iw || y < pad.t || y > pad.t + ih) { setAvPlaneHover(null); return; }
    const ix = Math.max(0, Math.min(avPlane.nx - 1, Math.floor(((x - pad.l) / iw) * avPlane.nx)));
    const iy = Math.max(0, Math.min(avPlane.ny - 1, Math.floor(((pad.t + ih - y) / ih) * avPlane.ny)));
    setAvPlaneHover(avPlane.cells[iy * avPlane.nx + ix]);
  }, [avPlane]);

  const avPlaneDetail = avPlaneHover || avPlane.current;

  const drawGrowth = useCallback(() => {
    const c = growC.current; if (!c || !theory) return;
    const lg = logRef.current;
    const ctx = c.getContext("2d");
    if (lg.lnA.length < 3) {
      ctx.fillStyle = CLR.panel; ctx.fillRect(0, 0, c.width, c.height);
      ctx.fillStyle = CLR.soft; ctx.font = "11px ui-monospace, monospace";
      ctx.textAlign = "center";
      ctx.fillText("按「运行」开始采样", c.width / 2, c.height / 2);
      return;
    }
    const xs = lg.lnA.map((d) => d[0]);
    const ys = lg.lnA.map((d) => d[1] / Math.LN10);
    const tmax = Math.max(10, xs[xs.length - 1]);
    plot(ctx, c.width, c.height, {
      xs, xr: [0, tmax], yr: [-8, 1],
      series: [
        { ys, color: CLR.blue, w: 1.8 },
        {
          ys: xs.map((t) => (Math.log(pert) + theory.best.lam.re * t) / Math.LN10),
          color: CLR.brick, w: 1.2,
        },
      ],
      xlab: "t [s]", ylab: "log₁₀|A_m*|",
    });
  }, [theory, pert]);

  const drawDiagnostics = useCallback(() => {
    const sim = simRef.current, lg = logRef.current;
    if (!sim) return;
    if (vehicleC.current) {
      const c = vehicleC.current, ctx = c.getContext("2d");
      if (lg.vehicle.length < 2) {
        ctx.fillStyle = CLR.panel; ctx.fillRect(0, 0, c.width, c.height);
        ctx.fillStyle = CLR.soft; ctx.font = "11px ui-monospace, monospace"; ctx.textAlign = "center";
        ctx.fillText("运行后显示所选车辆的速度与间距", c.width / 2, c.height / 2);
      } else {
        const xs = lg.vehicle.map((d) => d[0]);
        const eqg = idmGapAtSpeed(sim.ve, sim.models[Math.min(selectedVehicle, sim.N - 1)]);
        const vv = lg.vehicle.map((d) => 100 * (d[1] / Math.max(sim.ve, 0.1) - 1));
        const gg = lg.vehicle.map((d) => 100 * (d[2] / Math.max(eqg, 0.1) - 1));
        const span = Math.max(1, ...vv.map(Math.abs), ...gg.map(Math.abs));
        plot(ctx, c.width, c.height, {
          xs, xr: [xs[0], Math.max(xs[0] + 1, xs[xs.length - 1])], yr: [-span * 1.15, span * 1.15],
          series: [{ ys: vv, color: CLR.blue }, { ys: gg, color: CLR.amber }],
          xlab: "t [s]", ylab: "相对平衡 [%]", hline: 0,
        });
      }
    }
    if (envelopeC.current) {
      const c = envelopeC.current, ctx = c.getContext("2d");
      const xs = Array.from({ length: sim.N }, (_, i) => i);
      const ys = (lg.envelope || []).slice(0, sim.N);
      plot(ctx, c.width, c.height, {
        xs, xr: [0, Math.max(1, sim.N - 1)], yr: [0, Math.max(0.02, ...ys) * 1.15],
        series: [{ ys, color: CLR.teal }], xlab: "车号 n", ylab: "max |δv| [m/s]",
      });
    }
  }, [selectedVehicle]);

  useEffect(() => { drawGrowth(); drawDiagnostics(); }, [live.t, drawGrowth, drawDiagnostics]);

  const runExperiment = useCallback(() => {
    if (!theory) return;
    const rows = [];
    let headline = "", interpretation = "";
    if (experiment === "E1") {
      const bands = unstableIntervals(theory.band.map((d) => ({ v: d.v, bad: d.val < 0 })));
      headline = bands.length ? bands.map((b) => `${b[0].toFixed(2)}–${b[1].toFixed(2)} m/s`).join("；") : "全速域稳定";
      rows.push({ k: "扫描结果", v: headline }, { k: "文档基准", v: "0.91–18.61 m/s", note: "默认 IDM、N→∞" });
      interpretation = "负稳定裕度区间即微小扰动会放大的速度带；改变 IDM 或反馈参数后边界会同步更新。";
    } else if (experiment === "E2") {
      headline = meas ? `误差 ${Math.abs((meas.growth - theory.best.lam.re) / Math.max(1e-9, Math.abs(theory.best.lam.re))) * 100 < 999 ? (Math.abs((meas.growth - theory.best.lam.re) / Math.max(1e-9, Math.abs(theory.best.lam.re))) * 100).toFixed(1) : ">999"}%` : "请先运行约 30–60 s";
      rows.push({ k: "解析 Re λ", v: `${theory.best.lam.re.toFixed(6)} 1/s` }, { k: "仿真实测", v: meas ? `${meas.growth.toFixed(6)} 1/s` : "等待线性拟合窗口" });
      interpretation = "小扰动阶段两者应接近；进入非线性饱和后不再用指数增长拟合。";
    } else if (experiment === "E3") {
      headline = `预测 m*=${theory.best.m}${meas ? `，实测 m=${meas.domM}` : ""}`;
      rows.push({ k: "预测波数", v: `k*=${theory.best.k.toFixed(4)} rad/veh` }, { k: "主导模态", v: meas ? `${meas.domM}` : "运行后测量" });
      interpretation = "有限环只允许离散波数 2πm/N，因此 N 会决定能否采样到最危险长波。";
    } else if (experiment === "E4") {
      let mx = 0, wm = 0;
      for (let i = 0; i <= 500; i++) { const w = i * 0.012; const g = transferMag(theory.part, eff.fa, w); if (g > mx) { mx = g; wm = w; } }
      headline = mx <= 1 + 1e-5 ? "逐车不放大" : "逐车放大";
      rows.push({ k: "max |G|", v: mx.toFixed(5) }, { k: "最危险频率", v: `${wm.toFixed(3)} rad/s` }, { k: "判据", v: "max |G| ≤ 1" });
      interpretation = "这是开链单向级联的精确检查；启用后车或多向耦合时应以环道谱为准。";
    } else if (experiment === "E5") {
      const infMargin = longWaveMargin(theory.part, theory.p);
      headline = theory.best.lam.re <= 0 ? `N=${N} 环稳定` : `N=${N} 环失稳`;
      rows.push({ k: "最小波数", v: `${(2 * Math.PI / N).toFixed(4)} rad/veh` }, { k: "有限环 max Re λ", v: theory.best.lam.re.toFixed(6) }, { k: "N→∞ 长波裕度", v: infMargin.toFixed(6) });
      interpretation = "短环道可能跳过最危险长波，从而系统性高估稳定性；文档建议 N 至少上百。";
    } else if (experiment === "E6") {
      const modes = ["均匀", "随机", "连续"];
      for (const mode of modes) {
        const avs = makeArrangement(N, N * avShare, mode, avGroups, 71);
        let peak = 1;
        for (let iw = 1; iw <= 150; iw++) {
          const w = iw * 0.02; let logp = 0, pm = 0;
          for (let n = 0; n < N; n++) {
            const q = (n + 1) % N;
            const gain = avs[n] ? (avs[q] ? kappa : Math.min(kappa, kappa0)) : 0;
            logp += Math.log(Math.max(1e-12, transferMag(theory.part, gain, w)));
            pm = Math.max(pm, logp);
          }
          peak = Math.max(peak, Math.exp(Math.min(pm, 20)));
        }
        rows.push({ k: `${mode}排列`, v: `${peak.toFixed(3)}×` });
      }
      headline = rows.reduce((a, b) => parseFloat(a.v) <= parseFloat(b.v) ? a : b).k + "中途峰值最低";
      interpretation = "配方相同不代表传播过程相同；前缀连乘决定中途峰值，均匀分散通常能截短连续放大段。";
    } else if (experiment === "E7") {
      const gain = useFF ? kappa : 0;
      const criticalMean = Math.max(0, 1 - 0.5 * (theory.part.fv ** 2 - theory.part.fvl ** 2) / theory.part.fs);
      let pmin = null;
      if (gain > 0) {
        for (let i = 0; i <= 1000; i++) {
          const share = i / 1000;
          const profile = mixedGainProfile(N, gain, kappa0, share, arrangement, avGroups);
          if (profile.meanGain >= criticalMean - 1e-12) { pmin = share; break; }
        }
      }
      const profile = theory.mixed || mixedGainProfile(N, gain, kappa0, avShare, arrangement, avGroups);
      headline = pmin === null ? "当前增益下不存在可行渗透率" : `当前排列 p_min≈${(100 * pmin).toFixed(1)}%`;
      rows.push(
        { k: "实际 AV 数", v: `${profile.avCount}/${N}` },
        { k: "实际平均增益", v: profile.meanGain.toFixed(4) },
        { k: "临界平均增益", v: criticalMean.toFixed(4) },
        { k: "AV–AV / AV–HDV", v: `${profile.aaCount} / ${profile.ahCount}` }
      );
      interpretation = "按当前 N、排列、κ 与 κ₀ 逐个生成 AV/HDV 邻接关系，再扫描渗透率；不再把所有车辆误当成具有同一个 κ。";
    } else if (experiment === "E8") {
      headline = theory.psiMin >= 0 ? "全频段稳定" : "有限频率出现失稳势阱";
      rows.push({ k: "Ψ(0)", v: theory.psi0.toFixed(6), note: "长波门槛，与延迟无关" }, { k: "min Ψ(ω)", v: theory.psiMin.toFixed(6), note: `τ₀=${eff.tau0.toFixed(2)} s, τₐ=${eff.tauA.toFixed(2)} s` });
      interpretation = "若 Ψ(0)≥0 但 min Ψ<0，失稳完全由延迟在有限频率处挖出的势阱造成。";
    } else if (experiment === "E9") {
      let worst = 0;
      for (let i = 1; i < 300; i++) {
        const pt = idmPartials((i / 300) * P.v0 * 0.999, P), mu = pt.fv + pt.fvl;
        worst = Math.max(worst, 2 * (pt.fs + mu * pt.fvl) / (mu * mu));
      }
      let lo = 0, hi = 0.95;
      for (let i = 0; i < 60; i++) { const p = 0.5 * (lo + hi); if ((1 + p) / ((1 - p) ** 2) < worst) lo = p; else hi = p; }
      headline = `全速域闭合点 p_c=${(0.5 * (lo + hi)).toFixed(3)}`;
      rows.push({ k: "当前 p", v: eff.rearP.toFixed(3) }, { k: "最坏放宽需求", v: worst.toFixed(3) }, { k: "文档基准", v: "p_c=0.223" });
      interpretation = "后向间距反馈直接降低波速，稳定带从两端向中间收缩；它与后车加速度反馈不是同一机制。";
    } else if (experiment === "E10") {
      const vp = -theory.part.fs / (theory.part.fv + theory.part.fvl);
      const c1 = (1 - eff.rearP) * vp;
      headline = meas && meas.c1 !== null ? `理论 ${c1.toFixed(4)} / 实测 ${Math.abs(meas.c1).toFixed(4)}` : `理论 c₁=${c1.toFixed(4)}`;
      rows.push({ k: "V′(sₑ)", v: vp.toFixed(5) }, { k: "后视修正", v: `(1−p)=${(1 - eff.rearP).toFixed(3)}` }, { k: "实测", v: meas && meas.c1 !== null ? Math.abs(meas.c1).toFixed(5) : "运行后拟合相位" });
      interpretation = "加速度反馈不改变长波波速；后向间距反馈将其乘以 (1−p)，对称极限 p→1 时传播趋于停止。";
    }
    setExperimentResult({ key: experiment, headline, rows, interpretation });
  }, [experiment, theory, meas, N, avShare, avGroups, arrangement, useFF, kappa, kappa0, eff, P]);

  /* ---- 判定 ---- */
  const targetVerdict = theory ? (theory.inconclusive ? "待复核" : theory.isUnstable ? "失稳" : "稳定") : "—";
  const currentVerdict = currentTheory ? (currentTheory.inconclusive ? "待复核" : currentTheory.isUnstable ? "失稳" : "稳定") : "—";
  const waveState = live.std >= 0.5 ? "stop-and-go 已形成" : live.std >= 0.18 ? "扰动正在放大" : live.std >= 0.03 ? "弱扰动" : "近似均匀流";
  const runtimeKeys = ["fa", "fb", "alpha", "beta", "mFront", "mBack", "rearP", "tau0", "tauA"];
  const parametersPending = runtimeKeys.some((key) => Math.abs((currentEff[key] || 0) - (eff[key] || 0)) > 0.005);
  const verdictColor = (value) => value === "失稳" ? CLR.brick : value === "稳定" ? CLR.teal : CLR.soft;

  const num = (x, d = 4) => (x === null || x === undefined || !isFinite(x) ? "—" : x.toFixed(d));

  const cvs = { width: "100%", display: "block", background: CLR.panel };
  const nonlinearMeta = NONLINEAR_CASES[nonlinearCase];
  const nonlinearMetrics = nonlinearCase === "steepening" ? [
    ["剖面时刻", `${NONLINEAR_DATA.profileTimes[nonlinearTimeIndex]} s`],
    ["高次谐波比", NONLINEAR_DATA.harmonicRatio[nonlinearTimeIndex].toFixed(3)],
    ["注入模态", "m = 4"], ["平衡速度", "8.0 m/s"],
  ] : nonlinearCase === "saturation" ? [
    ["初始 σᵥ", "0.0141 m/s"], ["1800 s σᵥ", "2.962 m/s"],
    ["末 300 s 均值", "2.554 m/s"], ["稳定对照末值", "3.27×10⁻⁸ m/s"],
  ] : nonlinearCase === "soliton" ? [
    ["模板振幅 A", nonlinearAmplitude.toFixed(2)], ["KdV 相对幅宽", (2.2 / Math.sqrt(nonlinearAmplitude)).toFixed(2)],
    ["KdV 缩放", "L ∝ |A|⁻¹ᐟ²"], ["mKdV 缩放", "L ∝ |A|⁻¹"],
  ] : nonlinearCase === "trigger" ? [
    ["A₀=0.5 末窗", "0.00088 m/s"], ["A₀=18 末窗", "1.351 m/s"],
    ["A₀=18 峰值", "3.921 m/s"], ["迟滞结论", "尚未证实"],
  ] : [
    ["A₀=0.5", "tᵣ = 4 s"], ["A₀=4", "tᵣ = 1374 s"],
    ["A₀=6", "tᵣ = 4700 s"], ["A₀≥8", "tᵣ > 5400 s"],
  ];
  const applicationPair = APPLICATION_DATA ? APPLICATION_DATA[applicationCase] : null;
  const applicationMeta = applicationCase === "open" ? {
    title: "案例 A · 开放车队的正弦扰动与协同 AV",
    baseline: "HDV / IDM", control: "多前车前馈 + 后车速度反馈",
    setup: "60 车开放车队，vₑ=8 m/s；头车施加 0.6 m/s、周期 48 s、持续 5 周期的正弦速度扰动。",
    conclusion: "HDV 末车振幅大于头车输入，说明扰动沿上游放大；协同 AV 把传播增益压到 1 以下并更快恢复。",
  } : {
    title: "案例 B · 连续开放车流、低速瓶颈与固定空间 VSL",
    baseline: "法定上限 120 km/h · 小 T 失稳", control: "最大稳定裕度 VSL 114 km/h",
    setup: "12 km 单车道高速公路，法定最高限速 120 km/h；取短车头时距 T=0.7 s、净间距 18 m，事件前均匀速度约 75 km/h。0–1200 s 持续到达（1087 辆新车进入）；慢车在 x≈6 km 处于 300–420 s 降至约 65 km/h。VSL 固定在 x=5–6 km，并严格在慢车事件于 t=300 s 发生后才开启。",
    conclusion: "120 km/h 上限下，短时距工作点 Φ=−0.001415 s⁻²，处于长波失稳侧。t<300 s 时控制关闭；事件触发后才在 114–120 km/h 内执行由稳定裕度确定的 114 km/h 命令。该命令以 30 s 斜坡施加，再用连续流仿真事后评价排队、安全、延误和能耗。",
  };
  const applicationBenefits = !applicationPair ? [] : applicationCase === "open" ? [
    ["末车扰动振幅", "0.927 → 0.099 m/s", "下降 89.3%"],
    ["峰值速度离差", "0.549 → 0.193 m/s", "下降 64.9%"],
    ["恢复时刻", "454 → 392 s", "提前 62 s"],
    ["最小 TTC", "65.3 → 116.2 s", "提高 78.0%"],
    ["Wu 模型扰动附加电耗", "0.718 → 0.076 Wh/车", "下降 89.4%"],
  ] : [
    ["形成阶段排队尾波速度", "−2.69 → −0.28 m/s", "逆向传播减缓 89.4%"],
    ["低速车辆峰值", "121 → 102 辆", "下降 15.7%"],
    ["峰值速度离差", "8.912 → 4.107 m/s", "下降 53.9%"],
    ["最小 TTC", "2.34 → 6.70 s", "提高 186%"],
    ["累计延误", "6.866 → 4.799 veh·h", "下降 30.1%"],
    ["Wu 模型单位里程电耗", "17.521 → 17.328 kWh/100km", "下降 1.10%"],
  ];

  const feedbackTiles = [
    { id: "front", symbol: "κ", title: "前车加速度", note: "前馈", on: useFF, enable: () => setUseFF(true) },
    { id: "back", symbol: "fᵦ", title: "后车加速度", note: "双向", on: useBackAcc, enable: () => setUseBackAcc(true) },
    { id: "multi", symbol: "m", title: "多车信息", note: "核衰减", on: useMulti, enable: () => setUseMulti(true) },
    { id: "rear", symbol: "p", title: "后向间距", note: "后视", on: useRearGap, enable: () => setUseRearGap(true) },
    { id: "delay", symbol: "τ", title: "时间延迟", note: "感知 / 通信", on: useDelay, enable: () => setUseDelay(true) },
    { id: "hetero", symbol: "HV", title: "异质车流", note: "重型 / 保守", on: useHetero, enable: () => setUseHetero(true) },
    { id: "mixed", symbol: "AV", title: "混合交通", note: "AV 编组", on: useMixed, enable: () => setUseMixed(true) },
  ];

  const openFeedback = (tile) => {
    tile.enable();
    setFeedbackModal(tile.id);
  };
  const activeFeedback = feedbackTiles.find((tile) => tile.id === feedbackModal);

  return (
    <div style={{
      background: CLR.paper, color: CLR.ink, minHeight: "100vh", padding: "0",
      fontFamily: "'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', system-ui, sans-serif",
    }}>
      <style>{CSS}</style>
      {/* 顶栏 */}
      <div style={{
        display: "flex", alignItems: "center", justifyContent: "space-between",
        padding: "12px 18px", borderBottom: `1px solid ${CLR.rule}`, background: CLR.panel, flexWrap: "wrap", gap: 10,
      }}>
        <div>
          <div style={{ fontSize: 15, fontWeight: 700, letterSpacing: "-0.01em" }}>跟驰稳定性实验台</div>
          <div style={{ fontSize: 10.5, color: CLR.soft, marginTop: 1 }}>
            环道微观仿真 × 解析判据 · 理论与实测同屏对照
          </div>
        </div>
        <div style={{ display: "flex", gap: 18, alignItems: "center", fontFamily: "ui-monospace, monospace", fontSize: 11 }}>
          <div><span style={{ color: CLR.soft }}>ρ </span>{rhoNow.toFixed(1)} veh/km</div>
          <div><span style={{ color: CLR.soft }}>sₑ </span>{se.toFixed(2)} m</div>
          <div><span style={{ color: CLR.soft }}>vₑ </span>{theory ? theory.ve.toFixed(2) : "—"} m/s</div>
          <div style={{ display: "flex", gap: 5, alignItems: "center" }}>
            <div style={{
              padding: "3px 9px", fontWeight: 700, fontSize: 10.5, letterSpacing: "0.05em",
              color: "#fff", background: verdictColor(currentVerdict),
            }} title={`当前已应用参数；${currentTheory ? currentTheory.statusBasis : "无可用判据"}`}>当前 {currentVerdict}</div>
            <div style={{
              padding: "3px 9px", fontWeight: 700, fontSize: 10.5, letterSpacing: "0.05em",
              color: "#fff", background: verdictColor(targetVerdict), opacity: parametersPending ? 1 : 0.58,
            }} title={`控制面板目标参数；${theory ? theory.statusBasis : "无可用判据"}`}>目标 {targetVerdict}</div>
          </div>
        </div>
      </div>

      <nav className="tsl-tabs" aria-label="实验台分区">
        {[
          ["live", "实时环道"],
          ["propagation", "传播与探针"],
          ["theory", "理论与验证"],
          ["nonlinear", "非线性现象"],
          ["applications", "应用案例"],
          ["data", "真实数据"],
        ].map(([key, label]) => (
          <button key={key} className={`tsl-tab ${activeTab === key ? "on" : ""}`}
            aria-selected={activeTab === key} onClick={() => setActiveTab(key)}>{label}</button>
        ))}
      </nav>

      {linkedCase && (
        <div className="tsl-case-strip" role="status">
          <div><b>{window.TSL_I18N?.language === "en" ? "Book case" : "书中"} {linkedCaseId} · {linkedCaseCopy.title}</b><span>{linkedCaseCopy.hint}</span></div>
          <button className="tsl-chip" type="button" onClick={() => {
            setActiveTab(linkedCase.tab);
            if (linkedCase.tab === "live") setParameterModal("scenario");
          }}>{linkedCase.tab === "live" ? "调整场景参数" : "查看案例"}</button>
        </div>
      )}

      <div className={`tsl-shell ${ringFocus ? "tsl-focus" : ""}`}
        style={{ display: activeTab === "live" || activeTab === "propagation" ? "block" : "none" }}>
        <section className="tsl-panel" style={{ display: activeTab === "live" ? "block" : "none" }} aria-label="实时环道与在线控制">
          <div className="tsl-live-grid">
            <div className="tsl-ring-wrap">
              <Card title="环形道路实时仿真" note={`速度颜色：红=慢 · 蓝=快 · 平衡=${theory ? theory.ve.toFixed(2) : "—"} m/s`}>
                <canvas ref={ringCanvas} className="tsl-ring-canvas" width={720} height={720}
                  role="img" aria-label="环形道路车辆运行状态，车辆颜色映射速度" />
                <div className="tsl-statusline" aria-live="polite">
                  <div className="tsl-status"><span>仿真时间</span><b>{live.t.toFixed(1)} s</b></div>
                  <div className="tsl-status"><span>速度标准差 σᵥ</span><b>{live.std.toFixed(3)} m/s</b></div>
                  <div className="tsl-status"><span>当前前馈 κ</span><b>{appliedKappa.toFixed(3)}</b></div>
                  <div className="tsl-status"><span>平均速度</span><b>{live.vbar.toFixed(2)} m/s</b></div>
                  <div className="tsl-status"><span>净间距范围</span><b>{live.gmin.toFixed(1)}–{live.gmax.toFixed(1)} m</b></div>
                  <div className="tsl-status"><span>交通状态</span><b style={{ color: live.std >= 0.18 ? CLR.brick : CLR.teal }}>{waveState}</b></div>
                </div>
              </Card>
            </div>

            <aside className="tsl-live-side">
              <Card title="场景与反馈" note="点击图标开启并调整">
                <div className="tsl-scenario-grid">
                  {feedbackTiles.map((tile) => (
                    <button key={tile.id} type="button" className={`tsl-scenario-tile ${tile.on ? "on" : ""}`}
                      onClick={() => openFeedback(tile)} aria-pressed={tile.on} title={`调整${tile.title}参数`}>
                      <span className="symbol">{tile.symbol}</span><b>{tile.title}</b><span>{tile.on ? "已启用 · 点击调整" : tile.note}</span>
                    </button>
                  ))}
                </div>
              </Card>

              <Card title="基础参数" note="无需离开仿真页面">
                <div className="tsl-core-controls">
                  <button className="tsl-core-button" type="button" onClick={() => setParameterModal("idm")}>
                    <span className="symbol">IDM</span>
                    <b>IDM 基本参数</b>
                    <small>T={P.T.toFixed(2)} s · a={P.a.toFixed(2)} · b={P.b.toFixed(2)}</small>
                  </button>
                  <button className="tsl-core-button" type="button" onClick={() => setParameterModal("scenario")}>
                    <span className="symbol">RING</span>
                    <b>场景参数</b>
                    <small>N={N} · ρ={rhoNow.toFixed(1)} veh/km · L={(Lring / 1000).toFixed(2)} km</small>
                  </button>
                </div>
                <div className="tsl-live-note">IDM 或环道几何改变后按新的均匀平衡态重建仿真；反馈参数仍可通过上方图标不中断地在线施加。</div>
              </Card>

              <Card title="运行控制" note="不中断仿真">
                <div className="tsl-toolbar">
                  <button className="tsl-btn" onClick={() => setRunning((r) => !r)} style={{
                    flex: 1, minWidth: 100, padding: "9px 12px", border: "none", cursor: "pointer",
                    background: running ? CLR.amber : CLR.blue, color: "#fff", fontWeight: 600,
                  }}>{running ? "暂停仿真" : live.t > 0 ? "继续仿真" : "开始仿真"}</button>
                  <button className="tsl-btn" onClick={() => { setRunning(false); resetSim(); }} style={{
                    padding: "9px 12px", border: `1px solid ${CLR.rule}`, background: "transparent", color: CLR.ink, cursor: "pointer",
                  }}>重置</button>
                  <button className="tsl-btn" onClick={() => setRingFocus((x) => !x)} style={{
                    padding: "9px 12px", border: `1px solid ${CLR.rule}`, background: "transparent", color: CLR.ink, cursor: "pointer",
                  }}>{ringFocus ? "退出专注" : "放大环道"}</button>
                </div>
                <Slider label="播放倍速" value={speedMul} set={setSpeedMul} min={1} max={80} step={1} unit="×" />
                <div className="tsl-live-note">播放倍速只改变观看速度，不改变积分步长。环道放大模式会隐藏右侧面板，仿真仍继续。</div>
              </Card>

              <Card title="在线干预" note="参数变化不会重置车辆状态">
                <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", gap: 8, marginBottom: 8, fontSize: 11.5 }}>
                  <span>前车加速度反馈：<b style={{ color: useFF ? CLR.teal : CLR.soft }}>{useFF ? `κ = ${kappa.toFixed(3)}` : "未启用"}</b></span>
                  <button className="tsl-chip" type="button" onClick={() => { setUseFF(true); setFeedbackModal("front"); }}>调整 κ</button>
                </div>
                <Slider label="平滑过渡时间" value={rampDuration} set={setRampDuration} min={0} max={300} step={5} unit="s" />
                <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr", gap: 7 }}>
                  <button className="tsl-btn" onClick={applyOnline} style={{
                    padding: "9px 10px", border: "none", background: CLR.teal, color: "#fff", cursor: "pointer", fontWeight: 600,
                  }}>应用在线干预</button>
                  <button className="tsl-btn" onClick={() => beginRamp(simRef.current, 0, rampDuration, "手动关闭")} style={{
                    padding: "9px 10px", border: `1px solid ${CLR.rule}`, background: "transparent", color: CLR.ink, cursor: "pointer",
                  }}>逐步关闭反馈</button>
                </div>
                <div className="tsl-live-note">
                  当前 κ={appliedKappa.toFixed(3)}，目标 κ={useFF ? kappa.toFixed(3) : "0.000"}。
                  当前理论={currentVerdict}，目标理论={targetVerdict}。
                  {parametersPending ? " 目标参数尚未完全应用。" : " 当前参数已与目标一致。"}
                  {useMixed && currentTheory && currentTheory.mixed &&
                    ` 当前混合车流实际平均增益=${currentTheory.mixed.meanGain.toFixed(4)}（AV–AV ${currentTheory.mixed.aaCount} 对，AV–HDV ${currentTheory.mixed.ahCount} 对）。`}
                  在 stop-and-go 已形成后再应用，可直接观察波动是否衰减。
                </div>
              </Card>

              <Card title="拥堵形成 → 控制消散演示" note="推荐入口">
                <button className="tsl-btn" onClick={prepareJamControlDemo} style={{
                  width: "100%", padding: "9px 10px", border: `1px solid ${CLR.blue}`, background: "transparent", color: CLR.blue, cursor: "pointer", fontWeight: 600,
                }}>准备“不稳定 IDM”演示</button>
                <label style={{ display: "flex", gap: 8, alignItems: "center", fontSize: 11.5, marginTop: 10 }}>
                  <input type="checkbox" checked={autoControl} onChange={(e) => { setAutoControl(e.target.checked); autoTriggeredRef.current = false; }} style={{ accentColor: CLR.blue }} />
                  600 s 时自动将 κ 平滑提高到 {Math.max(0.45, kappa).toFixed(2)}
                </label>
                <div className="tsl-live-note">预设采用 50 veh/km，并给 0 号车一个 −2.5 m/s 的制动脉冲。手动流程：开始仿真 → 等待 σᵥ 明显增大 → 点击“应用在线干预”。自动流程会在 600 s 触发。</div>
              </Card>

              <Card title="干预记录" note={interventions.length ? `${interventions.length} 次` : "尚未干预"}>
                {interventions.length === 0 ? <div className="tsl-live-note">每次在线改变都会保留仿真时间、增益变化和过渡时长。</div> :
                  interventions.slice().reverse().map((e, i) => (
                    <div className="tsl-event" key={`${e.t}-${i}`}>
                      <time>t={e.t.toFixed(1)} s</time>
                      <span>{e.source}：κ {e.from.toFixed(3)} → {e.to.toFixed(3)}（{e.duration.toFixed(0)} s）</span>
                    </div>
                  ))}
              </Card>
            </aside>
          </div>
        </section>

        <section className="tsl-panel" style={{ display: activeTab === "propagation" ? "block" : "none" }} aria-label="扰动传播与车辆探针">
          <div className="tsl-propagation">
            <div className="tsl-prop-stack">
              <Card title="位置–时间轨迹图" note="轨迹颜色映射速度 · 时间窗 1800 s">
                <canvas ref={trajectoryCanvas} width={960} height={480} style={{ ...cvs, imageRendering: "auto" }}
                  role="img" aria-label="车辆位置随时间变化的轨迹，轨迹颜色表示车辆速度" />
                <div className="tsl-canvas-caption"><span>环道位置 x/L →</span><span>时间向下 · 黄线=在线干预</span></div>
              </Card>
              <Card title="速度时空图" note="车号 × 时间 · 1800 s">
                <canvas ref={speedSTCanvas} width={960} height={270} style={{ ...cvs, imageRendering: "pixelated" }}
                  role="img" aria-label="各车辆速度随时间传播的时空图" />
                <div className="tsl-canvas-caption"><span>车号 n ↓ · 红=低于平衡速度</span><span>时间 → · 蓝=高于平衡速度</span></div>
              </Card>
              <Card title="间距时空图" note="车号 × 时间 · 1800 s">
                <canvas ref={gapSTCanvas} width={960} height={270} style={{ ...cvs, imageRendering: "pixelated" }}
                  role="img" aria-label="各车辆净间距随时间传播的时空图" />
                <div className="tsl-canvas-caption"><span>车号 n ↓ · 红=间距减小</span><span>时间 → · 蓝=间距增大</span></div>
              </Card>
            </div>
            <div className="tsl-prop-stack">
              <Card title={`车辆 ${Math.min(selectedVehicle, N - 1)} 探针`} note="速度与间距相对平衡值">
                <Slider label="观测车辆" value={selectedVehicle} set={(x) => setSelectedVehicle(Math.round(x))} min={0} max={Math.max(0, N - 1)} step={1} unit="号" />
                <canvas ref={vehicleC} width={720} height={300} style={cvs} role="img" aria-label="指定车辆速度和净间距随时间变化曲线" />
                <div className="tsl-canvas-caption"><span>蓝=速度偏差</span><span>黄=间距偏差</span></div>
              </Card>
              <Card title="扰动包络" note="每辆车经历的最大速度偏差">
                <canvas ref={envelopeC} width={720} height={260} style={cvs} role="img" aria-label="沿车辆编号的最大速度扰动包络" />
              </Card>
              <Card title="当前读数" note={`t=${live.t.toFixed(1)} s`}>
                <div className="tsl-statusline">
                  <div className="tsl-status"><span>σᵥ</span><b>{live.std.toFixed(3)}</b></div>
                  <div className="tsl-status"><span>σₛ</span><b>{live.gapStd.toFixed(3)}</b></div>
                  <div className="tsl-status"><span>κ</span><b>{appliedKappa.toFixed(3)}</b></div>
                </div>
                <div className="tsl-live-note">黄线标出干预时刻。曲线和时空图不会在在线调参时清空，因此可以直接比较干预前后的放大与衰减。</div>
              </Card>
            </div>
          </div>
        </section>
      </div>

      <div className="tsl-shell" style={{ display: activeTab === "nonlinear" ? "block" : "none" }}>
        <section className="tsl-panel tsl-nonlinear" aria-label="五类典型非线性稳定性现象">
          <div className="tsl-nl-picker" role="tablist" aria-label="选择非线性现象">
            {Object.entries(NONLINEAR_CASES).map(([key, item]) => (
              <button key={key} className={`tsl-nl-choice ${nonlinearCase === key ? "on" : ""}`}
                onClick={() => setNonlinearCase(key)} aria-selected={nonlinearCase === key}>
                <b>{item.label}</b><span>{item.kicker}</span>
              </button>
            ))}
          </div>

          <div className="tsl-nl-grid">
            <Card title={`E11–E12 · ${nonlinearMeta.label}`} note={nonlinearMeta.kicker}>
              <canvas ref={nonlinearCanvas} width={980} height={430} style={cvs}
                role="img" aria-label={`${nonlinearMeta.label}的可观测证据图`} />
              <div className="tsl-canvas-caption">
                <span>{nonlinearCase === "steepening" ? "灰=初始剖面 · 红=所选时刻" : nonlinearCase === "saturation" ? "红=失稳工况 · 蓝=稳定对照" : nonlinearCase === "soliton" ? "蓝=KdV 孤波 · 黄=mKdV kink" : nonlinearCase === "trigger" ? "灰=全过程峰值 · 红=末窗残余" : "A₀≥8 m/s 按 5400 s 上界显示右删失"}</span>
                <span>基准参数：N=120 · RK4 · Δt=0.05 s</span>
              </div>
            </Card>

            <div className="tsl-prop-stack">
              <Card title="交互控制" note="不改变实时环道状态">
                {nonlinearCase === "steepening" && <Slider label="速度剖面时刻" value={nonlinearTimeIndex} set={(x) => setNonlinearTimeIndex(Math.round(x))} min={0} max={4} step={1} unit={`→ ${NONLINEAR_DATA.profileTimes[nonlinearTimeIndex]} s`} />}
                {nonlinearCase === "soliton" && <Slider label="模板振幅 A" value={nonlinearAmplitude} set={setNonlinearAmplitude} min={0.4} max={2.2} step={0.05} unit="" />}
                {(nonlinearCase === "trigger" || nonlinearCase === "recovery") && <div className="tsl-live-note">横轴扫描单车制动幅值 A₀。点由完整 IDM 长时积分得到；恢复曲线将 5400 s 内未达到阈值的样本作为右删失数据。</div>}
                {nonlinearCase === "saturation" && <div className="tsl-live-note">纵轴使用 log₁₀σᵥ：直线段对应近似固定增长率，弯折说明瞬时增长率随振幅改变。</div>}
                <div className="tsl-nl-metric">
                  {nonlinearMetrics.map(([k, v]) => <div key={k}><span>{k}</span><b>{v}</b></div>)}
                </div>
              </Card>

              <Card title="如何判读" note="现象 → 证据 → 限制">
                <div><span className={`tsl-nl-badge ${nonlinearMeta.verdict.includes("未证实") || nonlinearMeta.verdict.includes("模板") ? "warn" : ""}`}>{nonlinearMeta.verdict}</span></div>
                <div style={{ fontSize: 11.5, lineHeight: 1.65, marginTop: 9 }}>{nonlinearMeta.evidence}</div>
                <div className="tsl-live-note" style={{ borderLeft: `2px solid ${nonlinearMeta.verdict.includes("完整") ? CLR.teal : CLR.amber}`, paddingLeft: 8 }}>{nonlinearMeta.caveat}</div>
              </Card>

              <Card title="统一实验口径" note="便于论文复现">
                <div className="tsl-live-note" style={{ marginTop: 0 }}>
                  E11：vₑ=8 m/s、m=4、A₀=0.02 m/s、1800 s；稳定对照 vₑ=20 m/s。<br />
                  E12：vₑ=20 m/s，0 号车单次制动 A₀=0.5–18 m/s；恢复定义为 σᵥ&lt;0.02 m/s 连续 300 s。
                </div>
              </Card>
            </div>
          </div>
        </section>
      </div>

      <div className="tsl-shell" style={{ display: activeTab === "applications" ? "block" : "none" }}>
        <section className="tsl-panel tsl-applications" aria-label="稳定性分析的交通管理应用案例">
          <div className="tsl-app-picker" role="tablist" aria-label="选择应用案例">
            <button className={`tsl-nl-choice ${applicationCase === "open" ? "on" : ""}`}
              onClick={() => { setApplicationCase("open"); setApplicationProbeIndex(3); }} aria-selected={applicationCase === "open"}>
              <b>案例 A · 开放车队协同稳定</b><span>头车有限周期正弦扰动 · HDV 与多前车/后车协同 AV 对照</span>
            </button>
            <button className={`tsl-nl-choice ${applicationCase === "bottleneck" ? "on" : ""}`}
              onClick={() => { setApplicationCase("bottleneck"); setApplicationProbeIndex(3); }} aria-selected={applicationCase === "bottleneck"}>
              <b>案例 B · 低速车辆与固定 VSL</b><span>连续到达车流 · 事件点上游固定 1 km 限速区</span>
            </button>
          </div>

          {!APPLICATION_DATA ? <Card title="应用数据未加载"><div role="alert" className="tsl-live-note">请使用完整版或双击版实验台打开此模块。</div></Card> : <>
            <Card title={applicationMeta.title} note="完整 IDM · 开放边界 · 1200 s">
              <div style={{ fontSize: 11.5, lineHeight: 1.6 }}>{applicationMeta.setup}</div>
              <div className="tsl-live-note" style={{ borderLeft: `2px solid ${CLR.blue}`, paddingLeft: 8 }}>{applicationMeta.conclusion}</div>
            </Card>

            <div className="tsl-app-grid">
              <Card title={applicationMeta.baseline} note="位置–时间轨迹 · 颜色=速度">
                <canvas ref={applicationBaseCanvas} width={720} height={350} style={{ ...cvs, imageRendering: "auto" }}
                  role="img" aria-label={`${applicationMeta.baseline}条件下按速度着色的道路位置时间轨迹图`} />
                <div className="tsl-canvas-caption"><span>每条线=一辆车的道路轨迹</span><span>{applicationCase === "bottleneck" ? "黑色虚线=向上游移动的排队尾部" : "轨迹颜色映射瞬时速度"}</span></div>
              </Card>
              <Card title={applicationMeta.control} note="位置–时间轨迹 · 颜色=速度">
                <canvas ref={applicationControlCanvas} width={720} height={350} style={{ ...cvs, imageRendering: "auto" }}
                  role="img" aria-label={`${applicationMeta.control}条件下按速度着色的道路位置时间轨迹图`} />
                <div className="tsl-canvas-caption"><span>与左图共用位置范围和速度色标</span><span>{applicationCase === "bottleneck" ? "绿色带=固定 VSL 区；黑色虚线=排队尾波" : "事件时刻以竖线标记"}</span></div>
              </Card>
            </div>

            {applicationCase === "bottleneck" && <Card title="120 km/h 上限下为何取 114 km/h" note="允许范围内直接最大化稳定裕度 Φ">
              <canvas ref={applicationVslCanvas} width={1160} height={300} style={cvs}
                role="img" aria-label="候选可变限速的长波稳定裕度与纯稳定性目标扫描" />
              <div className="tsl-canvas-caption"><span>红线 Φ=0 为中性边界；118–120 km/h 仍在失稳侧</span><span>在一公里控制区的运营约束 114–120 km/h 内最大化 Φ，得到 114 km/h；安全与运行指标不进入在线目标</span></div>
            </Card>}

            <div className="tsl-app-bottom">
              <Card title={`车辆 ${applicationPair.baseline.probeIds[Math.min(applicationProbeIndex, applicationPair.baseline.probeIds.length - 1)]} 速度曲线`} note="红=无控制 · 绿=稳定性控制">
                <Slider label="指定车辆" value={applicationProbeIndex} set={(x) => setApplicationProbeIndex(Math.round(x))}
                  min={0} max={applicationPair.baseline.probeIds.length - 1} step={1}
                  unit={`→ ${applicationPair.baseline.probeIds[Math.min(applicationProbeIndex, applicationPair.baseline.probeIds.length - 1)]} 号`} />
                <canvas ref={applicationProbeCanvas} width={920} height={310} style={cvs}
                  role="img" aria-label="指定车辆在无控制和稳定性控制下的速度时间曲线对比" />
              </Card>

              <div className="tsl-prop-stack">
                <Card title="安全—效率—能耗读数" note="事后仿真评价 · 不参与在线选速">
                  <div className="tsl-app-benefits">
                    {applicationBenefits.map(([k, v, d]) => <div key={k}><span>{k}</span><b>{v}</b><span style={{ color: d.includes("增加") ? CLR.amber : CLR.teal, marginTop: 3 }}>{d}</span></div>)}
                  </div>
                </Card>
                <Card title="管理含义" note={applicationCase === "open" ? "稳定性收益" : "控制收益与代价"}>
                  <div className="tsl-live-note" style={{ marginTop: 0 }}>
                    {applicationCase === "open"
                      ? "协同 AV 不改变头车扰动，也几乎不改变平均速度；收益来自削弱逐车放大、提高 TTC，并把 Wu 模型估计的扰动附加电耗降低 89.4%。总电耗只下降 0.15%，因为大部分能量用于维持基准巡航。"
                      : "法定上限仍是 120 km/h；真正导致失稳的是短时距 T=0.7 s 与高密度工作点。慢车在 t=300 s 触发扰动，VSL 此前保持关闭并从该时刻起响应。无控制排队尾部以 −2.69 m/s 向上游传播；事件触发型 114 km/h VSL 使尾波速度降到 −0.28 m/s、速度离差下降 53.9%、延误下降 30.1%。除 Φ 外均为事后读数。"}
                  </div>
                </Card>
              </div>
            </div>
          </>}
        </section>
      </div>

      <div className="tsl-shell" style={{ display: activeTab === "data" ? "block" : "none" }}>
        <section className="tsl-panel" aria-label="NGSIM 真实轨迹稳定性证据">
          <div style={{ maxWidth: 1540, margin: "0 auto", padding: "14px 18px 36px", display: "grid", gap: 12 }}>
            <Card title="D01 · NGSIM US-101 真实轨迹" note="观测证据 · 不是对整体串行失稳的单独证明">
              <div style={{ display: "grid", gridTemplateColumns: "minmax(0,1.75fr) minmax(300px,.75fr)", gap: 14, alignItems: "start" }}>
                <img src="ngsim_us101_d01.png" alt="NGSIM US-101 车辆轨迹、车队速度和逐车扰动传递增益"
                  style={{ width: "100%", display: "block", border: `1px solid ${CLR.rule}`, background: "#fff" }} />
                <div style={{ display: "grid", gap: 9 }}>
                  {[
                    ["数据规模", "39,187 条 · 119 辆车", "US-101 · 2 号车道 · 179.5 s"],
                    ["有效跟驰片段", "78 组", "每组连续时间不少于 20 s"],
                    ["时间车头时距", "1.77 s", "中位数；四分位区间 1.44–2.18 s"],
                    ["速度扰动增益", "1.033", "中位数；bootstrap 95% CI 0.954–1.077"],
                    ["增益大于 1", "55.1%", "bootstrap 95% CI 44.9%–65.4%"],
                    ["响应时滞 / 相关", "1.85 s / 0.794", "配对片段的中位数"],
                  ].map(([k, v, d]) => (
                    <div key={k} style={{ padding: "9px 11px", background: CLR.mist, borderLeft: `3px solid ${CLR.blue}` }}>
                      <span style={{ display: "block", color: CLR.soft, fontSize: 10 }}>{k}</span>
                      <b style={{ display: "block", fontFamily: "ui-monospace, monospace", fontSize: 17, margin: "2px 0" }}>{v}</b>
                      <span style={{ display: "block", color: CLR.soft, fontSize: 10, lineHeight: 1.45 }}>{d}</span>
                    </div>
                  ))}
                </div>
              </div>
            </Card>

            <div style={{ display: "grid", gridTemplateColumns: "repeat(2,minmax(0,1fr))", gap: 12 }}>
              <Card title="估计方法" note="可复现处理链">
                <div className="tsl-live-note" style={{ marginTop: 0, lineHeight: 1.75 }}>
                  对速度序列采用 21 帧三阶 Savitzky–Golay 平滑，并去除片段线性趋势；以跟驰车与前车速度标准差之比定义经验传递增益，时滞由互相关峰值估计。筛选仅保留连续不少于 20 s、前车速度标准差不少于 0.25 m/s 的同车道片段。
                </div>
              </Card>
              <Card title="怎样解释结果" note="统计不确定性必须保留">
                <div className="tsl-live-note" style={{ marginTop: 0, lineHeight: 1.75 }}>
                  样本中略多于一半的局部片段出现增益大于 1，但中位增益置信区间跨过 1。因此它支持“扰动可能逐车放大”的现象性证据，却不足以单独判定整段车流串行失稳；还需控制车道变换、共同外部激励与测量噪声。
                </div>
                <a href="https://doi.org/10.21949/1504477" target="_blank" rel="noreferrer"
                  style={{ display: "inline-block", marginTop: 10, color: CLR.blue, fontSize: 11 }}>官方数据 DOI ↗</a>
              </Card>
            </div>
          </div>
        </section>
      </div>

      <div className={`tsl-main ${activeTab === "theory" ? "theory" : ""}`}>
        {/* 左侧控制栏 */}
        <div className="tsl-side" style={{ borderRight: `1px solid ${CLR.rule}`, background: CLR.paper }}>
          <Group title="场景">
            {[
              ["S0 基准 IDM", true, null, null],
              ["S1 前车加速度前馈", useFF, setUseFF, "κ"],
              ["S2 后车加速度", useBackAcc, setUseBackAcc, "fᵦ"],
              ["S3 多前车 / 多后车", useMulti, setUseMulti, "αβm"],
              ["S4 后向间距反馈", useRearGap, setUseRearGap, "p"],
              ["S5 时间延迟", useDelay, setUseDelay, "τ"],
              ["S6 异质车流", useHetero, setUseHetero, "HV"],
              ["S7 混合交通 / AV 编组", useMixed, setUseMixed, "AV"],
            ].map(([lab, on, set, badge]) => (
              <label key={lab} style={{
                display: "flex", alignItems: "center", gap: 8, marginBottom: 6,
                cursor: set ? "pointer" : "default", fontSize: 12,
                color: on ? CLR.ink : CLR.soft,
              }}>
                <input type="checkbox" checked={on} disabled={!set}
                  onChange={(e) => set && set(e.target.checked)}
                  style={{ accentColor: CLR.blue }} />
                <span style={{ flex: 1 }}>{lab}</span>
                {badge && <span style={{ fontFamily: "ui-monospace, monospace", fontSize: 10, color: CLR.soft }}>{badge}</span>}
              </label>
            ))}
          </Group>

          <Group title="环道">
            <div style={{ display: "flex", gap: 5, marginBottom: 10 }}>
              <button className={`tsl-chip ${!fixedLength ? "on" : ""}`} onClick={() => setFixedLength(false)}>密度控制</button>
              <button className={`tsl-chip ${fixedLength ? "on" : ""}`} onClick={() => setFixedLength(true)}>周长控制</button>
            </div>
            <Slider label="车辆数 N" value={N} set={(x) => setN(Math.round(x))} min={12} max={160} step={1} unit="辆" />
            <Slider label="密度 ρ" value={fixedLength ? rhoNow : rho} set={(x) => { setRho(x); setFixedLength(false); }} min={6} max={110} step={0.5} unit="veh/km" disabled={fixedLength} />
            <Slider label="环道周长 L" value={Lring / 1000} set={(x) => { setRingLength(x); setFixedLength(true); }} min={0.2} max={8} step={0.01} unit="km" disabled={!fixedLength} />
            <Slider label="车长 l" value={l} set={setL} min={3} max={16} step={0.5} unit="m" />
            <div style={{ fontSize: 10, color: CLR.soft, fontFamily: "ui-monospace, monospace", marginTop: 2 }}>
              ρ = {rhoNow.toFixed(2)} veh/km · 平均净间距 {se.toFixed(2)} m
            </div>
          </Group>

          <Group title="IDM 参数">
            <Slider label="期望速度 v₀" value={P.v0} set={(x) => setP({ ...P, v0: x })} min={10} max={40} step={0.1} unit="m/s" />
            <Slider label="安全时距 T" value={P.T} set={(x) => setP({ ...P, T: x })} min={0.6} max={3} step={0.05} unit="s" />
            <Slider label="最小间距 s₀" value={P.s0} set={(x) => setP({ ...P, s0: x })} min={0.5} max={8} step={0.1} unit="m" />
            <Slider label="最大加速度 a" value={P.a} set={(x) => setP({ ...P, a: x })} min={0.2} max={3} step={0.05} unit="m/s²" />
            <Slider label="舒适减速度 b" value={P.b} set={(x) => setP({ ...P, b: x })} min={0.5} max={4} step={0.05} unit="m/s²" />
            <Slider label="加速度指数 δ" value={P.delta} set={(x) => setP({ ...P, delta: x })} min={1} max={8} step={1} unit="" />
          </Group>

          <Group title="反馈与延迟">
            <Slider label="前馈增益 κ" value={kappa} set={setKappa} min={0} max={0.95} step={0.005} unit="" disabled={!useFF} />
            <Slider label="后车加速度 fᵦ" value={fb} set={setFb} min={-0.5} max={0.55} step={0.005} unit="" disabled={!useBackAcc} />
            <Slider label="前向衰减 α" value={alpha} set={setAlpha} min={0} max={0.9} step={0.01} unit="" disabled={!useMulti || !useFF} />
            <Slider label="后向衰减 β" value={beta} set={setBeta} min={0} max={0.9} step={0.01} unit="" disabled={!useMulti || !useBackAcc} />
            <Slider label="前视车辆数 m₊" value={mFront} set={(x) => setMFront(Math.round(x))} min={1} max={8} step={1} unit="辆" disabled={!useMulti || !useFF} />
            <Slider label="后视车辆数 m₋" value={mBack} set={(x) => setMBack(Math.round(x))} min={1} max={8} step={1} unit="辆" disabled={!useMulti || !useBackAcc} />
            <Slider label="后向间距权重 p" value={rearP} set={setRearP} min={0} max={0.8} step={0.005} unit="" disabled={!useRearGap} />
            <Slider label="反应延迟 τ₀" value={tau0} set={setTau0} min={0} max={1.6} step={0.02} unit="s" disabled={!useDelay} />
            <Slider label="通信延迟 τₐ" value={tauA} set={setTauA} min={0} max={1.6} step={0.02} unit="s" disabled={!useDelay || (!useFF && !useBackAcc)} />
            {Math.abs(kernelTotal(eff)) >= 0.98 && <div style={{ fontSize: 10, lineHeight: 1.45, color: CLR.brick, borderLeft: `2px solid ${CLR.brick}`, paddingLeft: 7 }}>
              当前加速度核总量接近或超过 1，隐式方程可能奇异；这是模型失稳，不是绘图错误。
            </div>}
          </Group>

          <Group title="异质与混合交通">
            <Slider label="重型 / 保守车辆占比" value={heavyShare} set={setHeavyShare} min={0} max={0.8} step={0.01} unit="" disabled={!useHetero} />
            <Slider label="AV 渗透率" value={avShare} set={setAvShare} min={0} max={1} step={0.01} unit="" disabled={!useMixed} />
            <Slider label="AV–HV 降级增益 κ₀" value={kappa0} set={setKappa0} min={0} max={0.5} step={0.01} unit="" disabled={!useMixed} />
            <Slider label="AV 编组数" value={avGroups} set={(x) => setAvGroups(Math.round(x))} min={1} max={12} step={1} unit="组" disabled={!useMixed || arrangement !== "成组"} />
            <select className="tsl-select" value={arrangement} disabled={!useHetero && !useMixed} onChange={(e) => setArrangement(e.target.value)} aria-label="车辆排列方式">
              <option value="均匀">均匀</option><option value="随机">随机</option><option value="连续">连续</option><option value="成组">成组</option>
            </select>
          </Group>

          <Group title="扰动与数值">
            <Slider label="初始扰动幅值" value={pert} set={setPert} min={0.005} max={1} step={0.005} unit="m/s" />
            <Slider label="步长 dt" value={dt} set={setDt} min={0.005} max={0.05} step={0.005} unit="s" />
            <Slider label="播放倍速" value={speedMul} set={setSpeedMul} min={1} max={40} step={1} unit="×" />
            <Slider label="观测车辆" value={selectedVehicle} set={(x) => setSelectedVehicle(Math.round(x))} min={0} max={Math.max(0, N - 1)} step={1} unit="号" />
          </Group>

          <div style={{ display: "flex", gap: 6 }}>
            <button className="tsl-btn" onClick={() => setRunning((r) => !r)} style={{
              flex: 1, padding: "8px 0", border: "none", cursor: "pointer",
              background: running ? CLR.amber : CLR.blue, color: "#fff",
              fontSize: 12, fontWeight: 600, letterSpacing: "0.04em",
            }}>{running ? "暂停" : "运行"}</button>
            <button className="tsl-btn" onClick={() => { setRunning(false); resetSim(); }} style={{
              flex: 1, padding: "8px 0", border: `1px solid ${CLR.rule}`, cursor: "pointer",
              background: "transparent", color: CLR.ink, fontSize: 12, fontWeight: 600,
            }}>重置</button>
          </div>
        </div>

        {/* 右侧主区 */}
        <div style={{ padding: 14 }}>
          {/* 主图：时空图 + 环道 */}
          <div className="tsl-hero">
            <Card title="时空图  车号 × 时间"
              note={fieldMode === "速度" ? `速度偏差（红=慢 蓝=快）` : "净间距偏差（红=小 蓝=大）"}>
              <div style={{ display: "flex", gap: 5, marginBottom: 7 }}>
                {["速度", "间距"].map((x) => <button key={x} className={`tsl-chip ${fieldMode === x ? "on" : ""}`} onClick={() => { setFieldMode(x); resetSim(); }}>{x}</button>)}
              </div>
              <canvas ref={legacySTCanvas} width={880} height={300} style={{ ...cvs, imageRendering: "pixelated" }} />
              <div style={{ display: "flex", justifyContent: "space-between", fontSize: 9.5, color: CLR.soft, marginTop: 4, fontFamily: "ui-monospace, monospace" }}>
                <span>↑ 车号 0（波向上游即向下传播）</span>
                <span>时间 →</span>
              </div>
            </Card>
            <Card title="环道" note={`${live.t.toFixed(0)} s`}>
              <canvas ref={legacyRingCanvas} width={230} height={230} style={cvs} />
              <div style={{ fontSize: 10, fontFamily: "ui-monospace, monospace", color: CLR.soft, marginTop: 6, lineHeight: 1.6 }}>
                <div>σᵥ = {live.std.toFixed(4)} m/s</div>
                <div>σₛ = {live.gapStd.toFixed(4)} m</div>
                <div>v̄ = {live.vbar.toFixed(3)} m/s</div>
                <div>间距 {live.gmin.toFixed(2)} – {live.gmax.toFixed(2)} m</div>
                {(useHetero || useMixed) && <div style={{ marginTop: 3 }}>描边：绿=AV · 黄=异质车型</div>}
              </div>
            </Card>
          </div>

          {/* 按章节组织的可复现实验 */}
          <div style={{ marginBottom: 12 }}>
            <Card title="章节验证实验 E1–E10" note={EXPERIMENTS[experiment].section}>
              <div className="tsl-expgrid">
                <select className="tsl-select" value={experiment} onChange={(e) => { setExperiment(e.target.value); setExperimentResult(null); }} aria-label="选择验证实验">
                  {Object.entries(EXPERIMENTS).map(([k, e]) => <option key={k} value={k}>{k} · {e.title}</option>)}
                </select>
                <div style={{ fontSize: 11, color: CLR.soft, lineHeight: 1.45 }}>
                  {EXPERIMENTS[experiment].does}<br /><span style={{ color: CLR.ink }}>测量：{EXPERIMENTS[experiment].metric}</span>
                </div>
                <button className="tsl-btn" onClick={runExperiment} style={{ border: "none", background: CLR.teal, color: "#fff", padding: "8px 14px", cursor: "pointer", fontWeight: 600 }}>运行验证</button>
              </div>
              {experimentResult && experimentResult.key === experiment && (
                <div style={{ marginTop: 10, borderTop: `1px solid ${CLR.rule}`, paddingTop: 9 }}>
                  <div style={{ fontSize: 13, fontWeight: 700, color: experimentResult.headline.includes("失稳") || experimentResult.headline.includes("放大") ? CLR.brick : CLR.teal }}>{experimentResult.headline}</div>
                  <div style={{ display: "grid", gridTemplateColumns: "repeat(auto-fit,minmax(160px,1fr))", gap: 8, marginTop: 8 }}>
                    {experimentResult.rows.map((r, i) => <div key={i} style={{ background: CLR.paper, border: `1px solid ${CLR.rule}`, padding: "7px 8px", fontSize: 10.5 }}>
                      <div style={{ color: CLR.soft }}>{r.k}</div><div style={{ fontFamily: "ui-monospace, monospace", marginTop: 2 }}>{r.v}</div>{r.note && <div style={{ color: CLR.soft, marginTop: 2 }}>{r.note}</div>}
                    </div>)}
                  </div>
                  <div style={{ fontSize: 10.5, color: CLR.soft, lineHeight: 1.55, marginTop: 8 }}>{experimentResult.interpretation}</div>
                </div>
              )}
            </Card>
          </div>

          {(experiment === "E7" || useMixed) && avPlaneDetail && (
            <div style={{ marginBottom: 12 }}>
              <Card
                title="E7 AV 渗透率–速度稳定性二维平面"
                note={`κ=${(useFF ? kappa : 0).toFixed(3)} · κ₀=${kappa0.toFixed(3)} · ${arrangement}排列 · N=${N}`}>
                <canvas
                  ref={avPlaneC}
                  width={920}
                  height={430}
                  style={{ ...cvs, cursor: "crosshair" }}
                  onMouseMove={inspectAvPlane}
                  onMouseLeave={() => setAvPlaneHover(null)}
                  role="img"
                  aria-label="AV 渗透率和平衡速度组成的二维稳定性平面，绿色稳定，红色失稳"
                />
                <div className="tsl-canvas-caption">
                  <span>绿=稳定 · 红=失稳 · 黑线=稳定边界</span>
                  <span>白点=当前工作点 · 移动鼠标读取参数</span>
                </div>
                <div className="tsl-live-note" style={{ color: avPlaneDetail.stable ? CLR.teal : CLR.brick }}>
                  {avPlaneHover ? "指向位置" : "当前工作点"}：AV={Math.round(avPlaneDetail.share * 100)}%，
                  vₑ={avPlaneDetail.v.toFixed(2)} m/s，{avPlaneDetail.stable ? "稳定" : "失稳"}；
                  平均增益={avPlaneDetail.meanGain.toFixed(4)}，长波裕度={avPlaneDetail.longMargin.toFixed(5)}，
                  最小环周传递裕度={avPlaneDetail.transferMargin.toFixed(5)}。
                </div>
                {!useFF && <div className="tsl-live-note" style={{ color: CLR.brick }}>当前未启用 S1，因此平面按 κ=0 计算；启用前车加速度反馈后才会出现 AV 稳定区域。</div>}
              </Card>
            </div>
          )}

          {/* 验证读数 */}
          <Card title="验证读数  理论 vs 实测"
            note={meas ? `拟合窗口 ${meas.nWin} 点` : "运行后自动拟合线性增长段"}>
            <div style={{ overflowX: "auto" }}>
              <table style={{ width: "100%", borderCollapse: "collapse", fontSize: 11.5, fontFamily: "ui-monospace, monospace" }}>
                <thead>
                  <tr style={{ color: CLR.soft, fontSize: 10, letterSpacing: "0.06em" }}>
                    {["量", "解析预测", "仿真实测", "相对误差", "对应实验"].map((h) => (
                      <th key={h} style={{ textAlign: "left", padding: "4px 8px 6px 0", borderBottom: `1px solid ${CLR.rule}`, fontWeight: 500 }}>{h}</th>
                    ))}
                  </tr>
                </thead>
                <tbody>
                  {[
                    {
                      k: "增长率 Re λ [1/s]",
                      th: theory ? theory.best.lam.re : null,
                      ms: meas ? meas.growth : null,
                      e: "E2",
                    },
                    {
                      k: "振荡频率 Im λ [rad/s]",
                      th: theory ? Math.abs(theory.best.lam.im) : null,
                      ms: meas ? Math.abs(meas.omega) : null,
                      e: "E2",
                    },
                    {
                      k: "波速 c₁ [veh⁻¹]",
                      th: theory && theory.best.k > 0 ? Math.abs(theory.best.lam.im) / theory.best.k : null,
                      ms: meas && meas.c1 !== null ? Math.abs(meas.c1) : null,
                      e: "E10",
                    },
                    {
                      k: "最不稳定模态 m*",
                      th: theory ? theory.best.m : null,
                      ms: meas ? meas.domM : null,
                      e: "E3", int: true,
                    },
                  ].map((r) => {
                    const err = r.th && r.ms && Math.abs(r.th) > 1e-9
                      ? Math.abs((r.ms - r.th) / r.th) : null;
                    return (
                      <tr key={r.k}>
                        <td style={{ padding: "5px 8px 5px 0", color: CLR.soft }}>{r.k}</td>
                        <td style={{ padding: "5px 8px 5px 0" }}>{r.int ? (r.th ?? "—") : num(r.th, 5)}</td>
                        <td style={{ padding: "5px 8px 5px 0" }}>{r.int ? (r.ms ?? "—") : num(r.ms, 5)}</td>
                        <td style={{ padding: "5px 8px 5px 0", color: err === null ? CLR.soft : err < 0.05 ? CLR.teal : err < 0.2 ? CLR.amber : CLR.brick }}>
                          {err === null ? "—" : (err * 100).toFixed(1) + "%"}
                        </td>
                        <td style={{ padding: "5px 0", color: CLR.soft }}>{r.e}</td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
            <div style={{ fontSize: 10.5, color: CLR.soft, marginTop: 8, lineHeight: 1.6 }}>
              增长率由模态 m* 的复振幅 A(t) 拟合：Re λ = d ln|A|/dt，Im λ = d arg A/dt，波速 c₁ = −Im λ / k。
              仅在扰动仍处线性阶段（|A| 未饱和）时有效。
            </div>
          </Card>

          <div className="tsl-diagnostics">
            <Card title={`车辆 ${Math.min(selectedVehicle, N - 1)} 时序`} note="蓝=速度偏差 · 黄=间距偏差">
              <canvas ref={legacyVehicleC} width={520} height={190} style={cvs} />
            </Card>
            <Card title="扰动包络" note="各车经历的最大速度偏差">
              <canvas ref={legacyEnvelopeC} width={520} height={190} style={cvs} />
            </Card>
          </div>

          {/* 分析图 */}
          <div className="tsl-charts">
            <Card title="E1 全速域不稳定带" note={theory && theory.mixedSimple ? "按实际平均增益计算长波边界" : theory && theory.exactPsi ? "min Ψ < 0 即失稳" : "长波裕度 < 0 即失稳"}>
              <canvas ref={bandC} width={300} height={180} style={cvs} />
            </Card>
            <Card title="E2 模态振幅增长" note={theory && theory.mixedSimple ? "蓝=实测 红=平均增益参考" : "蓝=实测 红=理论斜率"}>
              <canvas ref={growC} width={300} height={180} style={cvs} />
            </Card>
            <Card title="E3 环道增长率谱" note={theory ? (theory.mixedSimple ? `平均增益等效参考 · N=${N}` : `N=${N} 个模态`) : ""}>
              <canvas ref={specC} width={300} height={180} style={cvs} />
            </Card>
            <Card
              title={theory && theory.mixedSimple ? "混合环周传递裕度" : theory && theory.exactPsi ? "频域判据 Ψ(ω)" : "等效长波裕度"}
              note={theory && theory.mixedSimple ? "−mean ln|Gₙ| ≥ 0 即不放大" : theory && theory.exactPsi ? "全频段非负即线稳定" : "复杂耦合以环道谱复核"}>
              <canvas ref={psiC} width={300} height={180} style={cvs} />
            </Card>
          </div>

          {/* 目标工作点的解析量 */}
          {theory && (
            <div style={{ marginTop: 12 }}>
              <Card title="目标工作点的解析量" note={`${theory.statusBasis} · 目标${targetVerdict}`}>
                <div style={{
                  display: "grid", gridTemplateColumns: "repeat(auto-fit, minmax(132px, 1fr))",
                  gap: "8px 16px", fontFamily: "ui-monospace, monospace", fontSize: 11.5,
                }}>
                  {[
                    ["f_s", theory.part.fs, 4],
                    ["f_v", theory.part.fv, 4],
                    ["f_vl", theory.part.fvl, 4],
                    ["V′(sₑ)", -theory.part.fs / (theory.part.fv + theory.part.fvl), 4],
                    [theory.mixedSimple ? "混合长波裕度" : theory.exactPsi ? "Ψ(0)" : "长波裕度", theory.psi0, 5],
                    [theory.mixedSimple ? "min 传递裕度" : theory.exactPsi ? "min Ψ(ω)" : "等效最小裕度", theory.psiMin, 5],
                    [theory.mixedSimple ? "等效 max Re λ" : "max Re λ", theory.best.lam.re, 6],
                    ...(theory.mixed ? [["实际平均增益", theory.mixed.meanGain, 4]] : []),
                  ].map(([k, v, d]) => (
                    <div key={k}>
                      <div style={{ fontSize: 10, color: CLR.soft }}>{k}</div>
                      <div style={{ color: k.startsWith("min") || k.startsWith("Ψ") ? (v < 0 ? CLR.brick : CLR.teal) : CLR.ink }}>
                        {v.toFixed(d)}
                      </div>
                    </div>
                  ))}
                </div>
                <div style={{ fontSize: 10.5, color: CLR.soft, marginTop: 10, lineHeight: 1.7, borderTop: `1px solid ${CLR.rule}`, paddingTop: 8 }}>
                  {theory.mixedSimple
                    ? `S7 已按实际排列使用逐车增益：HDV=0，AV–HDV=min(κ,κ₀)，AV–AV=κ；当前目标中 AV–AV ${theory.mixed.aaCount} 对、AV–HDV ${theory.mixed.ahCount} 对。稳定标签来自完整环周传递乘积；增长率谱仅作平均增益参考。`
                    : eff.rearP > 0
                    ? "后向间距反馈进入 Q、R 并改变波速；加速度核进入 P，环道谱同时检查有限波长。"
                    : `长波门槛只取决于加速度核总量 m₀=${kernelTotal(theory.p).toFixed(3)}，与延迟和分配距离无关。`}
                  {useDelay && tau0 > 0 && " 反应延迟只在有限频率处把 Ψ 压下去——看 Ψ(ω) 图上的势阱。"}
                  {useHetero && " 异质车型场景的 Fourier 谱仍是均匀参考；最终判断需结合逐车传递与前缀放大实验 E6。"}
                </div>
              </Card>
            </div>
          )}
        </div>
      </div>

      {parameterModal === "idm" && (
        <MiniModal title="IDM 基本参数" subtitle="控制原始 IDM 的自由流、跟驰和制动响应。修改后将从新的均匀平衡态重新开始。" onClose={() => setParameterModal(null)}>
          <Slider label="期望速度 v₀" value={P.v0} set={(x) => setP({ ...P, v0: x })} min={10} max={40} step={0.1} unit="m/s" />
          <Slider label="安全时距 T" value={P.T} set={(x) => setP({ ...P, T: x })} min={0.6} max={3} step={0.05} unit="s" />
          <Slider label="最小间距 s₀" value={P.s0} set={(x) => setP({ ...P, s0: x })} min={0.5} max={8} step={0.1} unit="m" />
          <Slider label="最大加速度 aₘₐₓ" value={P.a} set={(x) => setP({ ...P, a: x })} min={0.2} max={3} step={0.05} unit="m/s²" />
          <Slider label="舒适减速度 b" value={P.b} set={(x) => setP({ ...P, b: x })} min={0.5} max={4} step={0.05} unit="m/s²" />
          <Slider label="加速度指数 δ" value={P.delta} set={(x) => setP({ ...P, delta: x })} min={1} max={8} step={1} unit="" />
          <div className="tsl-live-note">当前工作点：sₑ={se.toFixed(2)} m，vₑ={theory ? theory.ve.toFixed(2) : "—"} m/s。减小 T 或改变 a、b 后，理论稳定域和仿真初始平衡态会同步更新。</div>
          <div style={{ display: "flex", justifyContent: "space-between", gap: 8, marginTop: 12 }}>
            <button className="tsl-chip" type="button" onClick={() => setP({ v0: 33.3, T: 1.5, s0: 2.0, a: 1.0, b: 1.5, delta: 4 })}>恢复默认 IDM</button>
            <button className="tsl-btn" type="button" onClick={() => setParameterModal(null)} style={{ border: "none", background: CLR.teal, color: "#fff", padding: "7px 12px", cursor: "pointer" }}>完成</button>
          </div>
        </MiniModal>
      )}

      {parameterModal === "scenario" && (
        <MiniModal title="场景参数" subtitle="设置环道规模、交通密度、扰动和数值积分参数。修改后将从新的均匀平衡态重新开始。" onClose={() => setParameterModal(null)}>
          <div style={{ display: "flex", gap: 5, marginBottom: 10 }}>
            <button className={`tsl-chip ${!fixedLength ? "on" : ""}`} type="button" onClick={() => setFixedLength(false)}>密度控制</button>
            <button className={`tsl-chip ${fixedLength ? "on" : ""}`} type="button" onClick={() => setFixedLength(true)}>周长控制</button>
          </div>
          <Slider label="车辆数 N" value={N} set={(x) => setN(Math.round(x))} min={12} max={160} step={1} unit="辆" />
          <Slider label="密度 ρ" value={fixedLength ? rhoNow : rho} set={(x) => { setRho(x); setFixedLength(false); }} min={6} max={110} step={0.5} unit="veh/km" disabled={fixedLength} />
          <Slider label="环道周长 L" value={Lring / 1000} set={(x) => { setRingLength(x); setFixedLength(true); }} min={0.2} max={8} step={0.01} unit="km" disabled={!fixedLength} />
          <Slider label="车长 l" value={l} set={setL} min={3} max={16} step={0.5} unit="m" />
          <Slider label="初始扰动幅值" value={pert} set={setPert} min={0.005} max={1} step={0.005} unit="m/s" />
          <Slider label="积分步长 dt" value={dt} set={setDt} min={0.005} max={0.05} step={0.005} unit="s" />
          <Slider label="观测车辆" value={selectedVehicle} set={(x) => setSelectedVehicle(Math.round(x))} min={0} max={Math.max(0, N - 1)} step={1} unit="号" />
          <div className="tsl-live-note">ρ={rhoNow.toFixed(2)} veh/km · L={(Lring / 1000).toFixed(2)} km · 平均净间距={se.toFixed(2)} m。密度控制时 L 随 N、ρ 联动；周长控制时 ρ 随 N、L 联动。</div>
          <div style={{ display: "flex", justifyContent: "flex-end", marginTop: 12 }}>
            <button className="tsl-btn" type="button" onClick={() => setParameterModal(null)} style={{ border: "none", background: CLR.teal, color: "#fff", padding: "7px 12px", cursor: "pointer" }}>完成</button>
          </div>
        </MiniModal>
      )}

      {feedbackModal && activeFeedback && (
        <MiniModal title={activeFeedback.title} subtitle="打开后即可在线修改；运行中的车辆状态不会被重置。" onClose={() => setFeedbackModal(null)}>
          <label className="tsl-modal-toggle">
            <input type="checkbox" checked={activeFeedback.on}
              onChange={(event) => {
                const enabled = event.target.checked;
                if (feedbackModal === "front") setUseFF(enabled);
                if (feedbackModal === "back") setUseBackAcc(enabled);
                if (feedbackModal === "multi") setUseMulti(enabled);
                if (feedbackModal === "rear") setUseRearGap(enabled);
                if (feedbackModal === "delay") setUseDelay(enabled);
                if (feedbackModal === "hetero") setUseHetero(enabled);
                if (feedbackModal === "mixed") setUseMixed(enabled);
              }} style={{ accentColor: CLR.teal }} />
            启用此机制
          </label>

          {feedbackModal === "front" && <>
            <Slider label="前馈增益 κ" value={kappa} set={setKappa} min={0} max={0.95} step={0.005} unit="" disabled={!useFF} />
            <div className="tsl-live-note">前车加速度通过通信或感知前馈进入控制律，主要削弱车队放大，而不改变单车局部闭环根。</div>
          </>}
          {feedbackModal === "back" && <>
            <Slider label="后车加速度增益 fᵦ" value={fb} set={setFb} min={-0.5} max={0.55} step={0.005} unit="" disabled={!useBackAcc} />
            <div className="tsl-live-note">该项形成双向加速度耦合；环道谱会同时检查有限波长模态。</div>
          </>}
          {feedbackModal === "multi" && <>
            <Slider label="前向衰减 α" value={alpha} set={setAlpha} min={0} max={0.9} step={0.01} unit="" disabled={!useMulti} />
            <Slider label="后向衰减 β" value={beta} set={setBeta} min={0} max={0.9} step={0.01} unit="" disabled={!useMulti} />
            <Slider label="前视车辆数 m₊" value={mFront} set={(x) => setMFront(Math.round(x))} min={1} max={8} step={1} unit="辆" disabled={!useMulti} />
            <Slider label="后视车辆数 m₋" value={mBack} set={(x) => setMBack(Math.round(x))} min={1} max={8} step={1} unit="辆" disabled={!useMulti} />
          </>}
          {feedbackModal === "rear" && <>
            <Slider label="后向间距权重 p" value={rearP} set={setRearP} min={0} max={0.8} step={0.005} unit="" disabled={!useRearGap} />
            <div className="tsl-live-note">后向间距反馈既改变稳定带，也会把长波传播速度乘以 (1−p)。</div>
          </>}
          {feedbackModal === "delay" && <>
            <Slider label="反应延迟 τ₀" value={tau0} set={setTau0} min={0} max={1.6} step={0.02} unit="s" disabled={!useDelay} />
            <Slider label="通信延迟 τₐ" value={tauA} set={setTauA} min={0} max={1.6} step={0.02} unit="s" disabled={!useDelay || (!useFF && !useBackAcc)} />
          </>}
          {feedbackModal === "hetero" && <>
            <Slider label="重型 / 保守车辆占比" value={heavyShare} set={setHeavyShare} min={0} max={0.8} step={0.01} unit="" disabled={!useHetero} />
            <select className="tsl-select" value={arrangement} disabled={!useHetero && !useMixed} onChange={(event) => setArrangement(event.target.value)}>
              <option value="均匀">均匀</option><option value="随机">随机</option><option value="连续">连续</option><option value="成组">成组</option>
            </select>
          </>}
          {feedbackModal === "mixed" && <>
            <Slider label="AV 渗透率" value={avShare} set={setAvShare} min={0} max={1} step={0.01} unit="" disabled={!useMixed} />
            <Slider label="AV–HV 降级增益 κ₀" value={kappa0} set={setKappa0} min={0} max={0.5} step={0.01} unit="" disabled={!useMixed} />
            <Slider label="AV 编组数" value={avGroups} set={(x) => setAvGroups(Math.round(x))} min={1} max={12} step={1} unit="组" disabled={!useMixed || arrangement !== "成组"} />
          </>}
          <div style={{ display: "flex", justifyContent: "flex-end", gap: 8, marginTop: 12 }}>
            <button className="tsl-chip" type="button" onClick={() => { setParameterModal("idm"); setFeedbackModal(null); }}>IDM 基本参数</button>
            <button className="tsl-chip" type="button" onClick={() => { setParameterModal("scenario"); setFeedbackModal(null); }}>场景参数</button>
            <button className="tsl-btn" type="button" onClick={() => setFeedbackModal(null)} style={{ border: "none", background: CLR.teal, color: "#fff", padding: "7px 12px", cursor: "pointer" }}>完成</button>
          </div>
        </MiniModal>
      )}
    </div>
  );
}

ReactDOM.createRoot(document.getElementById("root")).render(<App />);
