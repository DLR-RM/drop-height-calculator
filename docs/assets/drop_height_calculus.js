
  // Copyright (c) 2026 German Aerospace Center (DLR) 
  // Author: Anton Shu 
  
  const el = (id) => document.getElementById(id);

  const R = {
    K_s: { min: 10,  max: 1e5 },
    I_m: { min: 1e-6, max: 5e-3 },
    I_l: { min: 1e-4, max: 2e0 },
    I_eff: { min: 1e-4, max: 2e0 },
    h_norm: { min: 0.01, max: 10 }
  };

  const defaults = {
    K_s: 12000, I_m: 5e-4, N: 9, Delta_t_ms: 0.5,
    I_l: 0.01, Icoll: 1.14, l_eff: 0.1, tau_max: 150, dq_max: 15, p_norm_manual: 0.5, h_norm_manual: 1.0
  };

  const actuatorPresets = {
    "ral26":                   { K_s: 570, I_m: 1.54e-4, N: 25, tau_max: 120 },
    "humanoids15":             { K_s: 2300, I_m: 7.5e-5, N: 80, tau_max: 137 },    
  };

  const useCasePresets = {
    "foot":       { I_l: 0.0082, Icoll: 1.14, l_eff: 0.1, dq_max: 15, p_norm_manual: 0.5, h_norm_manual: 1.0, Delta_t_ms: 2 },
    "knee":       { I_l: 0.13, Icoll: 4, l_eff: 0.1, dq_max: 15, p_norm_manual: 0.5, h_norm_manual: 1.0, Delta_t_ms: 2 },
    "light":      { I_l: 0.0022, Icoll: 0.7, l_eff: 0.1, dq_max: 15, h_norm_manual: 1.0, Delta_t_ms: 1 },
    "agile_case": { I_l: 0.005, Icoll: 0.5, l_eff: 0.05, dq_max: 20, h_norm_manual: 1.0, Delta_t_ms: 1 }
  };

  const presetNames = {
    "jd12": "JD12",
    "ral26": "RAL26",
    "sensodrive_ambitious": "Sensodrive Ambitious",
    "sensodrive_conservative": "Sensodrive Conservative",
    "agile": "Agile"
  };

  const fmt = (x) => {
    if (!isFinite(x)) return String(x);
    const ax = Math.abs(x);
    if (ax !== 0 && (ax < 1e-3 || ax >= 1e4)) return x.toExponential(3);
    return x.toLocaleString(undefined, { maximumSignificantDigits: 6 });
  };

  function sliderToLog(v, min, max) {
    const t = v / 1000, logMin = Math.log10(min), logMax = Math.log10(max);
    return Math.pow(10, logMin + t * (logMax - logMin));
  }
  function logToSlider(x, min, max) {
    const logMin = Math.log10(min), logMax = Math.log10(max);
    const t = (Math.log10(x) - logMin) / (logMax - logMin);
    return Math.max(0, Math.min(1000, Math.round(t * 1000)));
  }

  function linspace(a, b, n) {
    const arr = new Array(n), step = (b - a) / (n - 1);
    for (let i = 0; i < n; i++) arr[i] = a + step * i;
    return arr;
  }
  function logspace(a, b, n) {
    const la = Math.log10(a), lb = Math.log10(b), arr = new Array(n), step = (lb - la) / (n - 1);
    for (let i = 0; i < n; i++) arr[i] = Math.pow(10, la + step * i);
    return arr;
  }

  function feasibilityMetric(p) {
    const Delta_t = p.Delta_t_ms / 1000;
    const I_eff = p.I_eff;
    if (!isFinite(I_eff) || I_eff <= 0) return Infinity;
    const omega_r = computeomega_r(p);
    return (Delta_t * omega_r / Math.PI);
  }
  
  function computebeta(p) {
    const pnorm = computePnorm(p, p.normalizePMode);
    const h =Math.sqrt(2)*Math.sqrt(p.Icoll*p.I_l*p.K_s*9.81*p.l_eff)*p.tau_max/(4*p.Icoll*p.K_s*9.81*p.l_eff*pnorm)
    return (!isFinite(h) || h < 0) ? NaN : h;
  }

    function computeomega_r(p) {
    const I_eff = p.I_eff;
    const omega_r = Math.sqrt(p.K_s * (I_eff + p.I_l) / ((p.I_l) *(I_eff)));
    const h = omega_r;
    return (!isFinite(h) || h < 0) ? NaN : h;
  }

  function computecsysfreq(p) {
    const h = computeomega_r(p) / (Math.PI * 2);
    return (!isFinite(h) || h < 0) ? NaN : h;
  }
  
  function computecDt(p) {
    const h = Math.PI/(computeomega_r(p));
    return (!isFinite(h) || h < 0) ? NaN : h;
  }
  
  function computeEnergyRatio(p) {
    const pnorm = computePnorm(p, p.normalizePMode);
    if (pnorm <= 0 || p.K_s <= 0) return NaN;
    return (p.tau_max * p.tau_max * p.I_l) / (p.K_s * pnorm * pnorm);
  }

  const setKsFromEnergyRatio = (params, r_E) => {
    const pnorm = computePnorm(params, params.normalizePMode);
    if (pnorm > 0 && r_E > 0) {
      params.K_s = (params.tau_max * params.tau_max * params.I_l) / (r_E * pnorm * pnorm);
    }
  };
  
  function computePext(p) {
    const Delta_t = p.Delta_t_ms / 1000;
    const I_eff = p.I_eff;
    const omega_r = computeomega_r(p);
    const afact = Delta_t * omega_r / Math.PI;

    if (afact < 1) {
      const secArg = Math.PI * afact / 2;
      const c = Math.cos(secArg);
      
      if (!isFinite(c) || Math.abs(c) < 1e-12) return Infinity;
      const sec = 1 / c;
      return sec*((p.I_l*omega_r*omega_r - p.K_s*afact*afact)*p.I_eff - p.I_l*p.K_s*afact*afact)*p.tau_max/(omega_r*p.K_s*p.I_eff);    
    } 
    else {
      const get_pext_t = (t) => {
        const tau_max = p.tau_max, J_l = p.I_l, K_s = p.K_s, Pi = Math.PI, I_eff_local = p.I_eff;
        const Jm_N2_Jl = I_eff_local + J_l;
        const term1 = I_eff_local * (-Delta_t * Delta_t * K_s + J_l * Pi * Pi) - J_l * K_s * Delta_t * Delta_t;
        const sqrt_term = Math.sqrt(Jm_N2_Jl);
        
        const num = 2 * tau_max * term1 * sqrt_term;
        const arg1 = (Math.sqrt(K_s) * sqrt_term * t) / (Math.sqrt(J_l) * Math.sqrt(I_eff_local));
        const arg2 = (Pi * t) / Delta_t;

        const den1 = Math.sqrt(K_s) * Pi * I_eff_local;
        const den2 = Math.sin(arg1) * Pi * Math.sqrt(I_eff_local) * Math.sqrt(J_l) - Delta_t * Math.sqrt(K_s) * sqrt_term * Math.sin(arg2);

        return num / (den1 * den2);
      };

      const get_dp_ext_dt = (t) => {
        const tau_max = p.tau_max, J_l = p.I_l, K_s = p.K_s, Pi = Math.PI, I_eff_local = p.I_eff;
        const Jm_N2_Jl = I_eff_local + J_l;
        const term1 = I_eff_local * (-Delta_t * Delta_t * K_s + J_l * Pi * Pi) - J_l * K_s * Delta_t * Delta_t;

        const arg1 = (Math.sqrt(K_s) * Math.sqrt(Jm_N2_Jl) * t) / (Math.sqrt(J_l) * Math.sqrt(I_eff_local));
        const arg2 = (Pi * t) / Delta_t;

        const num = 2 * tau_max * term1 * Jm_N2_Jl * (-Math.cos(arg1) + Math.cos(arg2));
        const den_part = Math.sin(arg1) * Pi * Math.sqrt(I_eff_local) * Math.sqrt(J_l) - Delta_t * Math.sqrt(K_s) * Math.sqrt(Jm_N2_Jl) * Math.sin(arg2);
        const den = I_eff_local * den_part * den_part;

        return num / den;
      };

      let min_pext = Infinity;
      const steps = 200;
      const dt_step = Delta_t / steps;

      for (let i = 1; i < steps; i++) {
        let t0 = i * dt_step;
        let t1 = (i + 1) * dt_step;
        let d0 = get_dp_ext_dt(t0);
        let d1 = get_dp_ext_dt(t1);

        if (d0 < 0 && d1 >= 0) {
          let low = t0, high = t1;
          for (let iter = 0; iter < 20; iter++) {
            let mid = (low + high) / 2;
            if (get_dp_ext_dt(mid) < 0) {
              low = mid;
            } else {
              high = mid;
            }
          }
          
          let pext_val = Math.abs(get_pext_t((low + high) / 2));
          if (pext_val < min_pext) {
            min_pext = pext_val;
          }
        }
      }

      if (min_pext === Infinity) {
        for (let i = 1; i < steps; i++) {
          let val = Math.abs(get_pext_t(i * dt_step));
          if (val < min_pext) min_pext = val;
        }
      }

      return min_pext;
    }
  }

  function computePextAbs(p) {
    const p_ext = computePext(p);
    return isFinite(p_ext) ? Math.abs(p_ext) : p_ext;
  }
  
  function computeH(p) {
    const p_ext = computePext(p);
    if (!isFinite(p_ext)) return Infinity;
    const h = p_ext * p_ext / (8 * p.Icoll * p.Icoll * 9.81 * p.l_eff * p.l_eff);
    return (!isFinite(h) || h < 0) ? NaN : h;
  }

  function computePnorm(p, mode) {
    if (mode === undefined) mode = p.normalizePMode || 'a';
    if (mode === 'a') {
      return p.p_norm_manual;
    } else if (mode === 'b') {
      return p.I_l * p.dq_max;
    } else if (mode === 'c') {
      return p.l_eff * p.Icoll * Math.sqrt(2 * 9.81 * p.h_norm_manual);
    }
    return p.p_norm_manual;
  }
  
  function computeHnorm(p, mode) {
    if (mode === undefined) mode = p.normalizeHMode || 'a';
    if (mode === 'c') {
      return p.h_norm_manual;
    }
    const p_n = computePnorm(p, mode);
    if (!isFinite(p_n)) return Infinity;
    const h = p_n * p_n / (8 * p.Icoll * p.Icoll * 9.81 * p.l_eff * p.l_eff);
    return (!isFinite(h) || h < 0) ? NaN : h;
  }

  function bindRange(id, paramKey, isLog) {
    const minEl = el(id + "_min"), maxEl = el(id + "_max"), slider = el(id), numEl = el(id + "_num");
    
    const updateBounds = () => {
      let nMin = parseFloat(minEl.value), nMax = parseFloat(maxEl.value);
      if (isNaN(nMin) || isNaN(nMax) || nMin >= nMax || (isLog && nMin <= 0)) return;
      
      let currentReal = isLog ? sliderToLog(+slider.value, R[paramKey].min, R[paramKey].max) : parseFloat(slider.value);
      currentReal = Math.max(nMin, Math.min(nMax, currentReal));
      
      if (isLog) {
        R[paramKey].min = nMin; R[paramKey].max = nMax;
        slider.value = logToSlider(currentReal, nMin, nMax);
      } else {
        slider.min = nMin; slider.max = nMax; slider.value = currentReal;
      }
      triggerUpdate();
    };

    if (minEl) minEl.addEventListener("change", updateBounds);
    if (maxEl) maxEl.addEventListener("change", updateBounds);

    if (numEl) {
      numEl.addEventListener("change", () => {
        let val = parseFloat(numEl.value);
        if (isNaN(val)) return;
        const nMin = isLog ? R[paramKey].min : parseFloat(minEl.value);
        const nMax = isLog ? R[paramKey].max : parseFloat(maxEl.value);
        val = Math.max(nMin, Math.min(nMax, val));
        numEl.value = parseFloat(Number(val).toPrecision(6));
        
        if (isLog) slider.value = logToSlider(val, nMin, nMax);
        else slider.value = val;
        
        triggerUpdate();
      });
    }
  }

  function triggerUpdate() {
    window.cancelAnimationFrame(window.__raf);
    window.__raf = window.requestAnimationFrame(updateAll);
  }
  
  let __syncingInertia = false;
  
  function syncJeffFromJmAndN() {
    if (__syncingInertia) return;
    __syncingInertia = true;
    try {
      const Jm = sliderToLog(+el("Jm").value, R.I_m.min, R.I_m.max);
      const N = Math.max(1e-12, +el("N").value);
      const Jeff = Jm * N * N;
      if (isFinite(Jeff) && Jeff > 0) {
        const v = logToSlider(Jeff, R.I_eff.min, R.I_eff.max);
        if (+el("Jeff").value !== v) el("Jeff").value = v;
        const JeffNum = el("Jeff_num");
        if (JeffNum && document.activeElement !== JeffNum) JeffNum.value = parseFloat(Number(Jeff).toPrecision(6));
        const JeffVal = el("Jeff_val");
        if (JeffVal) JeffVal.textContent = fmt(Jeff);
      }
    } finally {
      __syncingInertia = false;
    }
  }
  
  function syncJmFromJeffAndN() {
    if (__syncingInertia) return;
    __syncingInertia = true;
    try {
      const Jeff = sliderToLog(+el("Jeff").value, R.I_eff.min, R.I_eff.max);
      const N = Math.max(1e-12, +el("N").value);
      const Jm = Jeff / (N * N);
      if (isFinite(Jm) && Jm > 0) {
        const v = logToSlider(Jm, R.I_m.min, R.I_m.max);
        if (+el("Jm").value !== v) el("Jm").value = v;
        const JmNum = el("Jm_num");
        if (JmNum && document.activeElement !== JmNum) JmNum.value = parseFloat(Number(Jm).toPrecision(6));
        const JmVal = el("Jm_val");
        if (JmVal) JmVal.textContent = fmt(Jm);
      }
    } finally {
      __syncingInertia = false;
    }
  }

  el("Ks").value = logToSlider(defaults.K_s, R.K_s.min, R.K_s.max);
  el("Jm").value = logToSlider(defaults.I_m, R.I_m.min, R.I_m.max);
  el("Jl").value = logToSlider(defaults.I_l, R.I_l.min, R.I_l.max);
  el("Jeff").value = logToSlider(defaults.I_m * defaults.N * defaults.N, R.I_eff.min, R.I_eff.max);
  el("N").value = defaults.N; 
  el("Dt").value = defaults.Delta_t_ms;
  el("Icoll").value = defaults.Icoll; el("leff").value = defaults.l_eff; 
  el("taumax").value = defaults.tau_max; el("dqmax").value = defaults.dq_max;
  el("pnorm_manual").value = defaults.p_norm_manual;
  el("hnorm_manual").value = defaults.h_norm_manual;
  syncJeffFromJmAndN();

  bindRange("Ks", "K_s", true); bindRange("Jeff", "I_eff", true);
  bindRange("Jm", "I_m", true); bindRange("Jl", "I_l", true);
  bindRange("N", null, false); bindRange("Dt", null, false);
  bindRange("Icoll", null, false);
  bindRange("leff", null, false); bindRange("taumax", null, false);
  bindRange("dqmax", null, false); bindRange("pnorm_manual", null, false);
  bindRange("hnorm_manual", "h_norm", false);

  function plotInteractiveCurve(divId, xLabel, yLabel, xs, xIsLog, baseParams, paramSetter, computeY, xOpt = null, xSelProp = null, xOpt2 = null, extraOpts = {}) {
    const logY = el("logY").checked;
    const includeInvalidY = el("includeInvalidY") ? el("includeInvalidY").checked : false;
    
    const showPlotBeta = baseParams.showPlotBeta;
    const showPlotGamma = baseParams.showPlotGamma;
    
    const testedRatios = Array.isArray(extraOpts.testedRatios) ? extraOpts.testedRatios : [];
    const showTestedRatios = extraOpts.showTestedRatios === true && testedRatios.length > 0;

    const candidateName = el("customPresetName") ? el("customPresetName").value.trim() || "candidate1" : "candidate1";

    if (!showPlotBeta) xOpt = null;
    if (!showPlotGamma) xOpt2 = null;
    
    const getXSel = (p) => {
      if (typeof xSelProp === "function") return xSelProp(p);
      if (typeof xSelProp === "string") return p[xSelProp];
      return xSelProp;
    };

    const traces = [];
    const yScaleBounds = [];

    const addActuatorTraces = (pParams, name, lineColor, markerColor, isMain = false) => {
      const ysValid = [], ysInvalid = [], xV = [], yV = [];
      for (let i = 0; i < xs.length; i++) {
        const p = { ...pParams };
        if (typeof paramSetter === "function") paramSetter(p, xs[i]); 
        else p[paramSetter] = xs[i];
        
        const y = computeY(p);
        const feas = feasibilityMetric(p) < 1;
        
        if (isFinite(y) && !isNaN(y)) {
          ysValid.push(feas ? y : null);
          ysInvalid.push(!feas ? y : null);
          
          if (!logY || y > 0) {
            if (feas) {
              xV.push(xs[i]);
              yV.push(y);
            }
            if (feas || includeInvalidY) {
              yScaleBounds.push(y);
            }
          }
        } else {
          ysValid.push(null); ysInvalid.push(null);
        }
      }

      traces.push({
        x: xs, y: ysValid, type: "scatter", mode: "lines",
        line: { width: 3, color: lineColor },
        name: isMain ? candidateName : name,
        showlegend: true
      });
      traces.push({
        x: xs, y: ysInvalid, type: "scatter", mode: "lines",
        line: { width: 3, color: lineColor, dash: "dot" },
        name: isMain ? `${candidateName} (Long)` : `${name} (Long)`,
        showlegend: false
      });

      const actXSel = getXSel(pParams);
      if (actXSel !== null && isFinite(actXSel)) {
        const pSel = { ...pParams };
        if (typeof paramSetter === "function") paramSetter(pSel, actXSel);
        else pSel[paramSetter] = actXSel;
        let ySel = computeY(pSel);
        if (isFinite(ySel) && (!logY || ySel > 0)) {
          traces.push({
            x: [actXSel], y: [ySel], type: "scatter", mode: "markers",
            marker: { size: 12, color: markerColor, symbol: "x", line: { width: 2, color: markerColor } },
            name: isMain ? `${candidateName} Selected` : `${name} Selected`,
            showlegend: false
          });
        }
      }

      if (baseParams.showPlotGamma && typeof extraOpts.getGammaX === "function") {
        const xGam = extraOpts.getGammaX(pParams);
        if (xGam !== null && isFinite(xGam)) {
          const pGam = { ...pParams };
          if (typeof paramSetter === "function") paramSetter(pGam, xGam);
          else pGam[paramSetter] = xGam;
          const yGam = computeY(pGam);
          if (isFinite(yGam) && (!logY || yGam > 0)) {
            traces.push({
              x: [xGam], y: [yGam], type: "scatter", mode: "markers",
              marker: { size: 10, color: markerColor, symbol: "diamond" },
              name: isMain ? `${candidateName} γ` : `${name} γ`,
              showlegend: false
            });

          }
        }
      }
    };

    addActuatorTraces(baseParams, "Current", "#00658b", "#00658b", true);

    if (Array.isArray(extraOpts.comparisons)) {
      extraOpts.comparisons.forEach(c => {
        addActuatorTraces(c.params, c.name, c.color, c.color, false);
      });
    }

    traces.push({
      x: [null], y: [null], type: "scatter", mode: "lines",
      line: { width: 2.5, color: "#888888" },
      name: "short contact",
      showlegend: true
    });
    traces.push({
      x: [null], y: [null], type: "scatter", mode: "lines",
      line: { width: 2.5, color: "#888888", dash: "dot" },
      name: "long contact",
      showlegend: true
    });
	

	traces.push({
          x: [null], y: [null], type: "scatter", mode: "markers", 
          marker: { size: 12, color: "#888888", symbol: "x" }, name: "selected"
        });
    const hasGamma = showPlotGamma && (typeof extraOpts.getGammaX === "function" || (xOpt2 !== null && isFinite(xOpt2)));
    if (hasGamma) {
      traces.push({
        x: [null], y: [null], type: "scatter", mode: "markers", 
        marker: { size: 12, color: "#888888", symbol: "diamond" }, name: "$\\Large\\gamma_{k\\mathrm{max}}$"
      });
    }
    const layout = {
      margin: { l: 60, r: 20, t: 80, b: 52 },
      showlegend: true, legend: { orientation: "h", x: 0, y: 1.25},
      xaxis: { title: xLabel, type: xIsLog ? "log" : "linear" },
      yaxis: { title: yLabel, type: logY ? "log" : "linear" },
      paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",
      shapes: []
    };

    if (extraOpts.forceEqualAspect) {
      layout.yaxis.scaleanchor = "x";
      layout.yaxis.scaleratio = 1;
    }

    if (xOpt !== null && isFinite(xOpt)) {
      const pOpt = { ...baseParams };
      if (typeof paramSetter === "function") paramSetter(pOpt, xOpt); else pOpt[paramSetter] = xOpt;
      const yOpt = computeY(pOpt);
      
      if (isFinite(yOpt)) {
        traces.push({
          x: [xOpt], y: [yOpt], type: "scatter", mode: "markers", 
          marker: { size: 12, color: "#a855f7", symbol: "star" }, name: "β"
        });
        
        layout.shapes.push({
          type: 'line',
          x0: xOpt, x1: xOpt,
          y0: 0, y1: 1, yref: 'paper',
          line: { color: '#a855f7', width: 2, dash: 'dash' }
        });
      }
    }

    if (xOpt2 !== null && isFinite(xOpt2)) {
      const pOpt2 = { ...baseParams };
      if (typeof paramSetter === "function") paramSetter(pOpt2, xOpt2); else pOpt2[paramSetter] = xOpt2;
      const yOpt2 = computeY(pOpt2);
      
      if (isFinite(yOpt2)) {
        traces.push({
          x: [xOpt2], y: [yOpt2], type: "scatter", mode: "markers", 
          marker: { size: 12, color: "#ea580c", symbol: "diamond" }, name: "γ",
        showlegend: false
        });
		
        layout.shapes.push({
          type: 'line',
          x0: xOpt2, x1: xOpt2,
          y0: 0, y1: 1, yref: 'paper',
          line: { color: '#ea580c', width: 2, dash: 'dash' }
        });
      }
    }

    const PAD = 0.05;
    if (xs.length > 0) {
      const xmin = xs[0], xmax = xs[xs.length - 1];
      if (xIsLog) {
        const lo = Math.log10(xmin), hi = Math.log10(xmax), span = (hi - lo) || 1;
        layout.xaxis.range = [lo - PAD*span, hi + PAD*span];
      } else {
        const span = (xmax - xmin) || 1; layout.xaxis.range = [xmin - PAD*span, xmax + PAD*span];
      }
    }
    
    if (yScaleBounds.length > 0) {
      const ymin = Math.min(...yScaleBounds), ymax = Math.max(...yScaleBounds);
      if (logY) {
        const lo = Math.log10(ymin), hi = Math.log10(ymax), span = (hi - lo) || 1;
        layout.yaxis.range = [lo - PAD*span, hi + PAD*span];
      } else {
        const span = (ymax - ymin) || 1; layout.yaxis.range = [ymin - PAD*span, ymax + PAD*span];
      }
    }
    if (showTestedRatios) {
      const testedColors = ["#0f766e", "#7c3aed", "#dc2626"];
      testedRatios.forEach((ratio, idx) => {
        if (!(isFinite(ratio) && ratio > 0)) return;
        const pTest = { ...baseParams };
        if (typeof paramSetter === "function") paramSetter(pTest, ratio); else pTest[paramSetter] = ratio;
        const yTest = computeY(pTest);
        const color = testedColors[idx % testedColors.length];
        layout.shapes.push({
          type: "line",
          x0: ratio, x1: ratio,
          y0: 0, y1: 1, yref: "paper",
          line: { color: color, width: 2, dash: "dot" }
        });
        if (isFinite(yTest) && (!logY || yTest > 0)) {
          traces.push({
            x: [ratio], y: [yTest], type: "scatter", mode: "markers",
            marker: { size: 10, color: color, symbol: "circle-open", line: { width: 2, color: color } },
            name: `tested ${ratio.toFixed(2)}`
          });
        }
      });
    }

    Plotly.react(divId, traces, layout, { responsive: true, displayModeBar: false });
  }

  function setReferenceNeededPlot(divId, isNeeded) {
    const plotDiv = el(divId);
    if (!plotDiv) return;
    const wrapper = el("wrap_" + divId);
    if (wrapper) wrapper.classList.toggle("reference-needed", isNeeded);
    if (isNeeded) {
      const layout = {
        margin: { l: 60, r: 20, t: 80, b: 52 },
        xaxis: { visible: false },
        yaxis: { visible: false },
        annotations: [],
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)"
      };
      Plotly.react(divId, [], layout, { responsive: true, displayModeBar: false });
    }
  }

  function renderReferenceText(id, active, val, unit) {
    const out = el(id);
    if (!out) return;
    if (!active) {
      out.textContent = "—";
      out.style.color = "#6b7280";
    } else if (!isFinite(val)) {
      out.textContent = "∞ (singularity)";
      out.style.color = "#b45309";
    } else if (isNaN(val)) {
      out.textContent = "— (invalid)";
      out.style.color = "#b91c1c";
    } else {
      out.textContent = `${fmt(val)} ${unit}`;
      out.style.color = "#111827";
    }
  }

  function readParams() {
    const N = +el("N").value;
    const Jm = sliderToLog(+el("Jm").value, R.I_m.min, R.I_m.max);
    const Jeff = sliderToLog(+el("Jeff").value, R.I_eff.min, R.I_eff.max);
    const Icoll = +el("Icoll").value;
    const l_eff = +el("leff").value;
    
    let I_l;
    if (el("computeJl") && el("computeJl").checked) {
      I_l = Icoll * (l_eff * l_eff);
    } else {
      I_l = sliderToLog(+el("Jl").value, R.I_l.min, R.I_l.max);
    }

    return {
      normalizeHActive: el("normalizeH").checked,
      normalizePActive: el("normalizeP").checked,
      normalizeHMode: document.querySelector('input[name="normH_option"]:checked')?.value || 'a',
      normalizePMode: document.querySelector('input[name="normP_option"]:checked')?.value || 'a',
      showPlotBeta: el("showPlotBeta").checked, 
      showPlotGamma: el("showPlotGamma").checked, 
      showTestedRatios: el("showTestedRatios").checked,
      p_norm_manual: +el("pnorm_manual").value,
      h_norm_manual: +el("hnorm_manual").value,
      K_s: sliderToLog(+el("Ks").value, R.K_s.min, R.K_s.max),
      I_m: Jm, I_eff: Jeff, N, Delta_t_ms: +el("Dt").value,
      I_l: I_l,
      Icoll: Icoll, l_eff: l_eff, 
      tau_max: +el("taumax").value, dq_max: +el("dqmax").value
    };
  }

  function updateAll() {
    const p = readParams();

    const selVal = el("actuatorPresetSelector").value;
    if (selVal !== "custom" && actuatorPresets[selVal]) {
        actuatorPresets[selVal].K_s = p.K_s;
        actuatorPresets[selVal].I_m = p.I_m;
        actuatorPresets[selVal].N = p.N;
        actuatorPresets[selVal].tau_max = p.tau_max;
        actuatorPresets[selVal].I_eff = p.I_eff;
    }

    const ucSelVal = el("useCasePresetSelector").value;
    if (ucSelVal !== "custom" && useCasePresets[ucSelVal]) {
        useCasePresets[ucSelVal].I_l = p.I_l;
        useCasePresets[ucSelVal].Icoll = p.Icoll;
        useCasePresets[ucSelVal].l_eff = p.l_eff;
        useCasePresets[ucSelVal].dq_max = p.dq_max;
        useCasePresets[ucSelVal].p_norm_manual = p.p_norm_manual;
        useCasePresets[ucSelVal].h_norm_manual = p.h_norm_manual;
        useCasePresets[ucSelVal].Delta_t_ms = p.Delta_t_ms;
    }

    const normalizeH = p.normalizeHActive;
    const normalizeP = p.normalizePActive;
    const computeJl = el("computeJl") ? el("computeJl").checked : false;

    // Disabling logic based on active normalization configuration
    const isPnormNeeded = (normalizeH && p.normalizeHMode === 'a') || (normalizeP && p.normalizePMode === 'a');
    const isDqmaxNeeded = (normalizeH && p.normalizeHMode === 'b') || (normalizeP && p.normalizePMode === 'b');
    const isHnormNeeded = (normalizeH && p.normalizeHMode === 'c') || (normalizeP && p.normalizePMode === 'c');

    ["Jl", "Jl_num", "Jl_min", "Jl_max"].forEach(id => {
      if (el(id)) el(id).disabled = computeJl;
    });
    const jlControl = el("Jl").closest('.control');
    if (jlControl) {
        jlControl.style.opacity = computeJl ? "0.6" : "1";
    }

    ["dqmax", "dqmax_num", "dqmax_min", "dqmax_max"].forEach(id => {
      if (el(id)) el(id).disabled = !isDqmaxNeeded;
    });
    const dqmaxControl = el("dqmax").closest('.control');
    if (dqmaxControl) dqmaxControl.style.opacity = !isDqmaxNeeded ? "0.6" : "1";

    ["pnorm_manual", "pnorm_manual_num", "pnorm_manual_min", "pnorm_manual_max"].forEach(id => {
      if (el(id)) el(id).disabled = !isPnormNeeded;
    });
    const pnormControl = el("pnorm_manual").closest('.control');
    if (pnormControl) pnormControl.style.opacity = !isPnormNeeded ? "0.6" : "1";

    ["hnorm_manual", "hnorm_manual_num", "hnorm_manual_min", "hnorm_manual_max"].forEach(id => {
      if (el(id)) el(id).disabled = !isHnormNeeded;
    });
    const hnormControl = el("hnorm_manual").closest('.control');
    if (hnormControl) hnormControl.style.opacity = !isHnormNeeded ? "0.6" : "1";

    const displayToggle = (id, condition) => { const e = el(id); if (e) e.style.display = condition ? "block" : "none"; };

    el("Ks_val").textContent = fmt(p.K_s);
    el("Jm_val").textContent = fmt(p.I_m);
    el("Jeff_val").textContent = fmt(p.I_eff);
    el("N_val").textContent = fmt(p.N);
    el("Dt_val").textContent = fmt(p.Delta_t_ms);
    el("Jl_val").textContent = fmt(p.I_l);
    el("Icoll_val").textContent = fmt(p.Icoll);
    el("leff_val").textContent = fmt(p.l_eff);
    el("taumax_val").textContent = fmt(p.tau_max);
    el("dqmax_val").textContent = fmt(p.dq_max);
    el("pnorm_manual_val").textContent = fmt(p.p_norm_manual);
    el("hnorm_manual_val").textContent = fmt(p.h_norm_manual);

    const setNumVal = (id, val) => { const elNum = el(id + "_num"); if (elNum && document.activeElement !== elNum) elNum.value = parseFloat(Number(val).toPrecision(6)); };
    setNumVal("Ks", p.K_s); setNumVal("Jm", p.I_m); setNumVal("Jeff", p.I_eff);
    setNumVal("N", p.N); setNumVal("Dt", p.Delta_t_ms); setNumVal("Jl", p.I_l); setNumVal("Icoll", p.Icoll); setNumVal("leff", p.l_eff); 
    setNumVal("taumax", p.tau_max); setNumVal("dqmax", p.dq_max); setNumVal("pnorm_manual", p.p_norm_manual);
    setNumVal("hnorm_manual", p.h_norm_manual);

    if (isFinite(p.I_m) && p.I_m > 0) {
        const v = logToSlider(p.I_m, R.I_m.min, R.I_m.max);
        if (+el("Jm").value !== v) el("Jm").value = v;
    }
    if (isFinite(p.I_eff) && p.I_eff > 0) {
        const v = logToSlider(p.I_eff, R.I_eff.min, R.I_eff.max);
        if (+el("Jeff").value !== v) el("Jeff").value = v;
    }
    if (computeJl && isFinite(p.I_l) && p.I_l > 0) {
        const v = logToSlider(p.I_l, R.I_l.min, R.I_l.max);
        if (+el("Jl").value !== v) el("Jl").value = v;
    }

    const beta = computebeta(p); 
    const gamma = 1.176927892;
    const Dtcrit = computecDt(p);

    const renderText = (id, val, unit) => {
      const out = el(id);
      if (!isFinite(val)) { out.textContent = "∞ (singularity)"; out.style.color = "#b45309"; }
      else if (isNaN(val)) { out.textContent = "— (invalid)"; out.style.color = "#b91c1c"; }
      else { out.textContent = `${fmt(val)} ${unit}`; out.style.color = "#111827"; }
    };
    
    const showPnormReference = normalizeP || (normalizeH && (p.normalizeHMode === 'a' || p.normalizeHMode === 'b'));
    const showHnormReference = normalizeH || (normalizeP && p.normalizePMode === 'c');
    renderText("pext_big", computePext(p), "Nms");
    renderText("h_big", computeH(p), "m");
    renderReferenceText("pnorm_big", showPnormReference, computePnorm(p, p.normalizePMode), "Nms");
    renderReferenceText("hnorm_big", showHnormReference, computeHnorm(p, p.normalizeHMode), "m");
    renderText("h_Dt", Dtcrit * 1000, "ms");
    renderText("h_sysfreq", computecsysfreq(p), "Hz");
    renderText("beta_big", beta, "");

    // Only show reference values after the corresponding normalization is active.
    renderReferenceText("normP_value", normalizeP, computePnorm(p, p.normalizePMode), "Nms");
    renderReferenceText("normH_value", normalizeH, computeHnorm(p, p.normalizeHMode), "m");

    const feasVal = feasibilityMetric(p);
	el("period_explenation_no_svg").style.display = (feasVal >= 1) ? "block" : "none";
	el("period_explenation_yes_svg").style.display = (feasVal < 1) ? "block" : "none";
    if (feasVal >= 1) el("feasValue").textContent = feasVal.toFixed(3);

    // Update dynamic System Characteristics Card SVG and text
    const sysSvgContainer = el("sys_svg_container");
    const sysText = el("sys_text");
    if (sysSvgContainer && sysText) {
      if (feasVal < 1) {
        sysSvgContainer.innerHTML = '<svg width="45" height="45" viewBox="0 0 100 100"><circle cx="50" cy="50" r="45" fill="#ecfdf5" stroke="#10b981" stroke-width="4"/><path d="M30 50 l15 15 l25 -25" stroke="#10b981" stroke-width="8" stroke-linecap="round" stroke-linejoin="round" fill="none"/></svg>';
        sysText.textContent = "$\\Delta_t$ is smaller than $T_{1/2}$ and the impact is a short-term impact";
        sysText.style.color = "#10b981";
      } else {
        sysSvgContainer.innerHTML = '<svg width="45" height="45" viewBox="0 0 100 100"><rect x="5" y="5" width="90" height="90" rx="20" fill="#fef2f2" stroke="#ef4444" stroke-width="4"/><path d="M35 35 l30 30 m0 -30 l-30 30" stroke="#ef4444" stroke-width="8" stroke-linecap="round" stroke-linejoin="round"/></svg>';
        sysText.textContent = "$\\Delta_t$ is bigger than $T_{1/2}$ and the impact is a long-term impact";
        sysText.style.color = "#ef4444";
      }
    }

    const deltaT_sec = p.Delta_t_ms / 1000;
    const gammaRatio = deltaT_sec / (Math.sqrt(p.I_l / p.K_s) * Math.PI);
    if (gammaRatio > 0.25) el("gammaRatioValue").textContent = gammaRatio.toFixed(3);

    if (window.MathJax) MathJax.typesetPromise();

    const xLog = el("logX").checked;

    const Jeff_min_val = parseFloat(el("Jeff_min").value) || R.I_eff.min;
    const Jeff_max_val = parseFloat(el("Jeff_max").value) || R.I_eff.max;
    const ratio_min = Jeff_min_val / p.I_l;
    const ratio_max = Jeff_max_val / p.I_l;

    const xsRatio = xLog ? logspace(ratio_min, ratio_max, 1000) : linspace(ratio_min, ratio_max, 1000);
    
    const N_min = parseFloat(el("N_min").value);
    const N_max = parseFloat(el("N_max").value);
    const xsN    = xLog ? logspace(N_min, N_max, 1000) : linspace(N_min, N_max, 1000);
    
    const Dt_min = parseFloat(el("Dt_min").value) || 0;
    const Dt_max = parseFloat(el("Dt_max").value) || 50;
    const xsDt   = linspace(Dt_min, Dt_max, 1000);
    
    const xsJeff = xLog ? logspace(R.I_eff.min, R.I_eff.max, 1000) : linspace(R.I_eff.min, R.I_eff.max, 1000);
    const xsKs   = xLog ? logspace(R.K_s.min, R.K_s.max, 1000) : linspace(R.K_s.min, R.K_s.max, 1000);
    const xsJm   = xLog ? logspace(R.I_m.min, R.I_m.max, 1000) : linspace(R.I_m.min, R.I_m.max, 1000);
    const xsJl   = xLog ? logspace(R.I_l.min, R.I_l.max, 1000) : linspace(R.I_l.min, R.I_l.max, 1000);

    const Icoll_min = parseFloat(el("Icoll_min").value);
    const Icoll_max = Math.max(Icoll_min + 0.1, parseFloat(el("Icoll_max").value));
    const xsIcoll = xLog ? logspace(Icoll_min, Icoll_max, 1000) : linspace(Icoll_min, Icoll_max, 1000);
    
    const leff_min = parseFloat(el("leff_min").value);
    const leff_max = Math.max(leff_min + 0.01, parseFloat(el("leff_max").value));
    const xsleff = xLog ? logspace(leff_min, leff_max, 1000) : linspace(leff_min, leff_max, 1000);

	const taumax_min = parseFloat(el("taumax_min").value) || 1;
	const taumax_max = parseFloat(el("taumax_max").value) || 400;
	const xsTaumax = xLog ? logspace(taumax_min, taumax_max, 1000) : linspace(taumax_min, taumax_max, 1000);

    const pnorm_current = computePnorm(p, p.normalizePMode);
    let r_E_min = 1e-4;
    let r_E_max = 1e4;
    if (pnorm_current > 0 && p.tau_max > 0 && p.I_l > 0) {
      r_E_min = (p.tau_max * p.tau_max * p.I_l) / (R.K_s.max * pnorm_current * pnorm_current);
      r_E_max = (p.tau_max * p.tau_max * p.I_l) / (R.K_s.min * pnorm_current * pnorm_current);
    }
    const xsEnergyRatio = xLog ? logspace(r_E_min, r_E_max, 1000) : linspace(r_E_min, r_E_max, 1000);

    const comp1Val = el("comparePreset1").value;
    const comp2Val = el("comparePreset2").value;
    const compParamsList = [];

    if (comp1Val !== "none" && actuatorPresets[comp1Val]) {
      const pr = actuatorPresets[comp1Val];
      compParamsList.push({
        params: {
          ...p,
          K_s: pr.K_s, I_m: pr.I_m, N: pr.N, tau_max: pr.tau_max,
          I_eff: pr.I_eff || (pr.I_m * pr.N * pr.N)
        },
        color: "#82a043", name: presetNames[comp1Val] || comp1Val
      });
    }

    if (comp2Val !== "none" && actuatorPresets[comp2Val]) {
      const pr = actuatorPresets[comp2Val];
      compParamsList.push({
        params: {
          ...p,
          K_s: pr.K_s, I_m: pr.I_m, N: pr.N, tau_max: pr.tau_max,
          I_eff: pr.I_eff || (pr.I_m * pr.N * pr.N)
        },
        color: "#f8de53", name: presetNames[comp2Val] || comp2Val
      });
    }

    const basePlotOpts = { comparisons: compParamsList };
    const basePlotOptsHN = { comparisons: compParamsList, getGammaX: (pParams) => Math.sqrt((gamma * pParams.I_l) / pParams.I_m) };
    const basePlotOptsHJeff = { comparisons: compParamsList, getGammaX: (pParams) => gamma * pParams.I_l };
    const setJeffFromRatio = (params, ratio) => { params.I_eff = params.I_l * ratio; };

    const isHOpen = el("details_h").open;
    const isPextOpen = el("details_pext").open;
    const isFreqOpen = el("details_freq").open;

    // --- Drop Height ($h$) Plots ---
    if (isHOpen) {
      const cur_paramH = normalizeH ? (c) => computeH(c)/computeHnorm(c, p.normalizeHMode) : computeH;
	  const cur_pLabelH = normalizeH ? "$h_\\rm{max} / h_\\rm{ref}$" : "$h_\\rm{max}$ (m)";
	  
      plotInteractiveCurve("plotHN_h", "$N$", cur_pLabelH, xsN, xLog, p, (p, v) => { p.N = v; p.I_eff = p.I_m * v * v; }, cur_paramH, null, "N", null, basePlotOptsHN);
      plotInteractiveCurve("plotHJeff_h", "$I_\\rm{a}$", cur_pLabelH, xsJeff, xLog, p, "I_eff", cur_paramH, null, "I_eff", null, basePlotOptsHJeff);
      plotInteractiveCurve("plotJm_h", "Motor Inertia $I_{m}$", cur_pLabelH, xsJm, xLog, p, (p, v) => { p.I_m = v; p.I_eff = v * p.N * p.N; }, cur_paramH, null, "I_m", null, basePlotOpts);
      plotInteractiveCurve("plotKs_h", "$K$", cur_pLabelH, xsKs, xLog, p, "K_s", cur_paramH, null, "K_s", null, basePlotOpts);
      plotInteractiveCurve("plotJl_h", "$I_l$", cur_pLabelH, xsJl, xLog, p, "I_l", cur_paramH, null, "$I_l$", null, basePlotOpts);
      plotInteractiveCurve("plotDt_h", "$\\Delta t$", cur_pLabelH, xsDt, false, p, "Delta_t_ms", cur_paramH, null, "Delta_t_ms", null, basePlotOpts);
      plotInteractiveCurve("plotIcoll_h", "$m_f$", cur_pLabelH, xsIcoll, xLog, p, "Icoll", cur_paramH, null, "Icoll", null, basePlotOpts);
      plotInteractiveCurve("plotleff_h", "$l$", cur_pLabelH, xsleff, xLog, p, "l_eff", cur_paramH, null, "l_eff", null, basePlotOpts);
		plotInteractiveCurve("plotTaumax_h", "$\\tau_\rm{max}$", cur_pLabelH, xsTaumax, xLog, p, "tau_max", cur_paramH, null, "tau_max", null, basePlotOpts);
      if (normalizeH) {
        setReferenceNeededPlot("plotRatio_h", false);
        setReferenceNeededPlot("plotEnergyRatio_h", false);
        plotInteractiveCurve("plotRatio_h", "$\\gamma = I_{a}/{I}_{l}$", cur_pLabelH, xsRatio, xLog, p, setJeffFromRatio, cur_paramH, beta, (params) => params.I_eff / params.I_l, gamma, Object.assign({ forceEqualAspect: normalizeH, testedRatios: [0.68, 0.45, 0.33], showTestedRatios: p.showTestedRatios }, basePlotOpts));
        plotInteractiveCurve("plotEnergyRatio_h", "Energy Ratio $E_{pot}/{E}_{kin}$", cur_pLabelH, xsEnergyRatio, xLog, p, setKsFromEnergyRatio, cur_paramH, null, computeEnergyRatio, null, Object.assign({ forceEqualAspect: normalizeH }, basePlotOpts));
      } else {
        setReferenceNeededPlot("plotRatio_h", false);
        setReferenceNeededPlot("plotEnergyRatio_h", true);
		plotInteractiveCurve("plotRatio_h", "$\\gamma = I_{a}/{I}_{l}$", cur_pLabelH, xsRatio, xLog, p, setJeffFromRatio, cur_paramH, beta, (params) => params.I_eff / params.I_l, gamma, Object.assign({ forceEqualAspect: normalizeH, testedRatios: [0.68, 0.45, 0.33], showTestedRatios: p.showTestedRatios }, basePlotOpts));
      }
    }

    // --- Admissible Impact ($p_{ext}$) Plots ---
    if (isPextOpen) {
      const cur_paramPext = normalizeP ? (c) => computePextAbs(c)/computePnorm(c, p.normalizePMode) : computePextAbs;
	const cur_pLabelPext = normalizeP ? "$p_{\\rm{m}\\rm{a}\\rm{x}} / \\it{p}_{\\rm{r}\\rm{e}\\rm{f}}$" : "$p_{\\rm{m}\\rm{a}\\rm{x}} $";

      plotInteractiveCurve("plotHN_pext", "$N$", cur_pLabelPext, xsN, xLog, p, (p, v) => { p.N = v; p.I_eff = p.I_m * v * v; }, cur_paramPext, null, "N", null, basePlotOptsHN);
      plotInteractiveCurve("plotHJeff_pext", "$I_{a}$", cur_pLabelPext, xsJeff, xLog, p, "I_eff", cur_paramPext, null, "I_eff", null, basePlotOptsHJeff);
      plotInteractiveCurve("plotJm_pext", "Motor Inertia $I_{m}$", cur_pLabelPext, xsJm, xLog, p, (p, v) => { p.I_m = v; p.I_eff = v * p.N * p.N; }, cur_paramPext, null, "I_m", null, basePlotOpts);
      plotInteractiveCurve("plotKs_pext", "Stiffness $K$", cur_pLabelPext, xsKs, xLog, p, "K_s", cur_paramPext, null, "K_s", null, basePlotOpts);
      plotInteractiveCurve("plotJl_pext", "Link Inertia $I_l$", cur_pLabelPext, xsJl, xLog, p, "I_l", cur_paramPext, null, "I_l", null, basePlotOpts);
      plotInteractiveCurve("plotDt_pext", "Δt", cur_pLabelPext, xsDt, false, p, "Delta_t_ms", cur_paramPext, null, "Delta_t_ms", null, basePlotOpts);
      plotInteractiveCurve("plotTaumax_pext", "Max collision torque $\\tau_\\rm{max}$", cur_pLabelPext, xsTaumax, xLog, p, "tau_max", cur_paramPext, null, "tau_max", null, basePlotOpts);
      displayToggle("cardPlotHIcoll_pext", computeJl);
      displayToggle("cardPlotHleff_pext", computeJl);
	  
      if (computeJl) {
        plotInteractiveCurve("plotIcoll_pext", "Foot mass $m_f$", cur_pLabelPext, xsIcoll, xLog, p, "Icoll", cur_paramPext, null, "Icoll", null, basePlotOpts);
        plotInteractiveCurve("plotleff_pext", "Lever length $l$", cur_pLabelPext, xsleff, xLog, p, "l_eff", cur_paramPext, null, "l_eff", null, basePlotOpts);
      }
      if (normalizeP) {
        setReferenceNeededPlot("plotRatio_pext", false);
        setReferenceNeededPlot("plotEnergyRatio_pext", false);
        plotInteractiveCurve("plotRatio_pext", "$\\gamma = I_{a} / I_{l}$", cur_pLabelPext, xsRatio, xLog, p, setJeffFromRatio, cur_paramPext, beta, (params) => params.I_eff / params.I_l, gamma, Object.assign({ forceEqualAspect: normalizeP, testedRatios: [0.68, 0.45, 0.33], showTestedRatios: p.showTestedRatios }, basePlotOpts));
        plotInteractiveCurve("plotEnergyRatio_pext", "Energy Ratio $E_{pot}/E_{kin}$", cur_pLabelPext, xsEnergyRatio, xLog, p, setKsFromEnergyRatio, cur_paramPext, null, computeEnergyRatio, null, Object.assign({ forceEqualAspect: normalizeP }, basePlotOpts));
      } else {
        setReferenceNeededPlot("plotRatio_pext", false);
        setReferenceNeededPlot("plotEnergyRatio_pext", true);
		plotInteractiveCurve("plotRatio_pext", "$\\gamma = I_{a} / I_{l}$", cur_pLabelPext, xsRatio, xLog, p, setJeffFromRatio, cur_paramPext, beta, (params) => params.I_eff / params.I_l, gamma, Object.assign({ forceEqualAspect: normalizeP, testedRatios: [0.68, 0.45, 0.33], showTestedRatios: p.showTestedRatios }, basePlotOpts));
      }
    }

    // --- Eigenfrequency ($f$) Plots ---
    if (isFreqOpen) {
      plotInteractiveCurve("plotFKs", "$K$", "Frequency (Hz)", xsKs, xLog, p, "K_s", computecsysfreq, null, "K_s", null, basePlotOpts);
      plotInteractiveCurve("plotFJl", "$I_l$", "Frequency (Hz)", xsJl, xLog, p, "I_l", computecsysfreq, null, "I_l", null, basePlotOpts);
      plotInteractiveCurve("plotFJm", "$I_{m}$", "Frequency (Hz)", xsJm, xLog, p, (p, v) => { p.I_m = v; p.I_eff = v * p.N * p.N; }, computecsysfreq, null, "I_m", null, basePlotOpts);
      plotInteractiveCurve("plotFJeff", "$I_{a}$", "Frequency (Hz)", xsJeff, xLog, p, "I_eff", computecsysfreq, null, "I_eff", null, basePlotOpts);
    }
  }

  // --- Synchronization & Toggle logic for open expandable panels ---
  let lastModified = 'H';

  function syncIfBothOpen(changedBy) {
    const hOpen = el("details_h").open;
    const pOpen = el("details_pext").open;
    if (hOpen && pOpen) {
      if (changedBy === 'H') {
        el("normalizeP").checked = el("normalizeH").checked;
        const selectedRadio = document.querySelector('input[name="normH_option"]:checked');
        if (selectedRadio) {
          const targetRadio = document.querySelector(`input[name="normP_option"][value="${selectedRadio.value}"]`);
          if (targetRadio) targetRadio.checked = true;
        }
      } else if (changedBy === 'P') {
        el("normalizeH").checked = el("normalizeP").checked;
        const selectedRadio = document.querySelector('input[name="normP_option"]:checked');
        if (selectedRadio) {
          const targetRadio = document.querySelector(`input[name="normH_option"][value="${selectedRadio.value}"]`);
          if (targetRadio) targetRadio.checked = true;
        }
      }
    }
    el("normH_options_container").style.display = el("normalizeH").checked ? "inline-flex" : "none";
    el("normP_options_container").style.display = el("normalizeP").checked ? "inline-flex" : "none";
  }

  // Bind Normalization option event listeners
  el("normalizeH").addEventListener("change", () => {
    lastModified = 'H';
    syncIfBothOpen('H');
    triggerUpdate();
  });
  document.querySelectorAll('input[name="normH_option"]').forEach(radio => {
    radio.addEventListener("change", () => {
      lastModified = 'H';
      syncIfBothOpen('H');
      triggerUpdate();
    });
  });

  el("normalizeP").addEventListener("change", () => {
    lastModified = 'P';
    syncIfBothOpen('P');
    triggerUpdate();
  });
  document.querySelectorAll('input[name="normP_option"]').forEach(radio => {
    radio.addEventListener("change", () => {
      lastModified = 'P';
      syncIfBothOpen('P');
      triggerUpdate();
    });
  });

  el("actuatorPresetSelector").addEventListener("change", (e) => {
    const v = e.target.value;
    if (actuatorPresets[v]) {
      const pr = actuatorPresets[v];
      el("Ks").value = logToSlider(pr.K_s, R.K_s.min, R.K_s.max);
      el("Jm").value = logToSlider(pr.I_m, R.I_m.min, R.I_m.max);
      el("N").value = pr.N;
      el("taumax").value = pr.tau_max;
      syncJeffFromJmAndN();
      triggerUpdate();
    }
  });

  el("useCasePresetSelector").addEventListener("change", (e) => {
    const v = e.target.value;
    if (useCasePresets[v]) {
      const pr = useCasePresets[v];
      el("Jl").value = logToSlider(pr.I_l, R.I_l.min, R.I_l.max);
      el("Icoll").value = pr.Icoll;
      el("leff").value = pr.l_eff;
      el("dqmax").value = pr.dq_max;
      if (pr.p_norm_manual !== undefined) el("pnorm_manual").value = pr.p_norm_manual;
      if (pr.h_norm_manual !== undefined) el("hnorm_manual").value = pr.h_norm_manual;
      if (pr.Delta_t_ms !== undefined) el("Dt").value = pr.Delta_t_ms;
      triggerUpdate();
    }
  });

  el("comparePreset1").addEventListener("change", triggerUpdate);
  el("comparePreset2").addEventListener("change", triggerUpdate);

  el("details_h").addEventListener("toggle", () => {
    syncIfBothOpen(lastModified);
    triggerUpdate();
  });
  el("details_pext").addEventListener("toggle", () => {
    syncIfBothOpen(lastModified);
    triggerUpdate();
  });
  el("details_freq").addEventListener("toggle", triggerUpdate);

  document.querySelectorAll("input[type=range], input[type=checkbox]").forEach(input => {
    input.addEventListener("input", () => {
      if (["Jl", "Icoll", "leff", "dqmax", "pnorm_manual", "hnorm_manual", "computeJl", "Dt"].includes(input.id)) {
          el("useCasePresetSelector").value = "custom";
      }

      if (input.id === "Jm" || input.id === "N") syncJeffFromJmAndN();
      if (input.id === "Jeff") syncJmFromJeffAndN();
      triggerUpdate();
    });
  });

  el("customPresetName").addEventListener("input", triggerUpdate);

  el("savePresetBtn").addEventListener("click", () => {
    const name = el("customPresetName").value.trim();
    if (!name) return;
    
    const id = "custom_" + name.toLowerCase().replace(/[^a-z0-9_]/g, '_');
    
    const p = readParams();
    actuatorPresets[id] = {
      K_s: p.K_s,
      I_m: p.I_m,
      N: p.N,
      tau_max: p.tau_max,
      I_eff: p.I_eff
    };
    presetNames[id] = name;
    
    ["actuatorPresetSelector", "comparePreset1", "comparePreset2"].forEach(selId => {
      const sel = el(selId);
      let exists = Array.from(sel.options).some(opt => opt.value === id);
      if (!exists) {
        const opt = document.createElement("option");
        opt.value = id;
        opt.textContent = name;
        sel.appendChild(opt);
      }
    });
    
    el("actuatorPresetSelector").value = id;
    triggerUpdate();
  });

  el("saveUseCasePresetBtn").addEventListener("click", () => {
    const name = el("customUseCasePresetName").value.trim();
    if (!name) return;
    
    const id = "custom_uc_" + name.toLowerCase().replace(/[^a-z0-9_]/g, '_');
    
    const p = readParams();
    useCasePresets[id] = {
      I_l: p.I_l,
      Icoll: p.Icoll,
      l_eff: p.l_eff,
      dq_max: p.dq_max,
      p_norm_manual: p.p_norm_manual,
      h_norm_manual: p.h_norm_manual,
      Delta_t_ms: p.Delta_t_ms
    };
    
    const sel = el("useCasePresetSelector");
    let exists = Array.from(sel.options).some(opt => opt.value === id);
    if (!exists) {
      const opt = document.createElement("option");
      opt.value = id;
      opt.textContent = name;
      sel.appendChild(opt);
    }
    
    sel.value = id;
    triggerUpdate();
  });

  window.addEventListener("resize", triggerUpdate);
  syncIfBothOpen(lastModified);
  triggerUpdate();