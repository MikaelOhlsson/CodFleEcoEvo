// Client-side animation for the Simulation tab.
//
// The server sends the simulated trajectory once per simulation (message
// "anim-data"). Playing, pausing and scrubbing then happen entirely in the
// browser, so they respond immediately no matter how busy the server is. The
// time-series plot is a static image rendered by the server; the server also
// sends where its time axis sits on the image (message "ts-map"), so that a
// time cursor can be drawn over it here.

$(function () {
  var data = null;       // trajectory and display settings from the server
  var tsMap = null;      // position of the time-series panels on the image
  var tNow = 0;          // time currently displayed
  var playing = false;
  var lastStamp = null;  // timestamp of the previous animation frame
  var DURATION = 20;     // seconds to play a full run at speed 1x
  var SLIDER_MAX = 1000;

  function el(id) { return document.getElementById(id); }

  // Linear interpolation of the trajectory at time t
  function stateAt(t) {
    var T = data.time, n = T.length, lo = 0, hi = n - 1;
    if (t <= T[0]) hi = 0;
    if (t >= T[n - 1]) lo = n - 1;
    while (hi - lo > 1) {
      var mid = (lo + hi) >> 1;
      if (T[mid] <= t) lo = mid; else hi = mid;
    }
    var w = hi > lo ? (t - T[lo]) / (T[hi] - T[lo]) : 0;
    function f(v) { return v[lo] + w * (v[hi] - v[lo]); }
    return { NC: f(data.NC), muC: f(data.muC), NF: f(data.NF) };
  }

  function niceStep(span, n) {
    var step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var r = span / n / step;
    if (r >= 7.5) step *= 10; else if (r >= 3.5) step *= 5;
    else if (r >= 1.5) step *= 2;
    return step;
  }
  function ticks(lo, hi, n) {
    var step = niceStep(hi - lo, n), out = [];
    var digits = Math.max(0, -Math.floor(Math.log10(step) + 1e-9));
    for (var v = Math.ceil(lo / step - 1e-9) * step; v <= hi + 1e-9 * step;
         v += step) {
      out.push({ v: v, label: (Math.abs(v) < 1e-9 * step ? 0 : v).toFixed(digits) });
    }
    return out;
  }
  function fmt(v) {
    return (v !== 0 && Math.abs(v) < 0.01) ? v.toExponential(1) : v.toFixed(2);
  }
  function dnorm(z, mu, sd) {
    var u = (z - mu) / sd;
    return Math.exp(-0.5 * u * u) / (sd * Math.sqrt(2 * Math.PI));
  }

  function drawGauss(s) {
    var svg = el("anim-gauss");
    if (!svg || svg.clientWidth === 0) return; // tab not visible
    var W = svg.clientWidth, H = svg.clientHeight;
    var m = { l: 72, r: 16, t: 12, b: 52 };
    var pw = W - m.l - m.r, ph = H - m.t - m.b;
    var z0 = data.zlim[0], z1 = data.zlim[1], ymax = data.ymax * 1.05;
    function X(z) { return m.l + (z - z0) / (z1 - z0) * pw; }
    function Y(y) { return m.t + ph - y / ymax * ph; }
    var out = [];
    var font = 'font-size="14" fill="#4d4d4d"';

    // Grid, ticks and tick labels
    ticks(z0, z1, 8).forEach(function (t) {
      var x = X(t.v).toFixed(1);
      out.push('<line x1="' + x + '" x2="' + x + '" y1="' + m.t + '" y2="' +
               (m.t + ph) + '" stroke="#ebebeb"/>');
      out.push('<text x="' + x + '" y="' + (m.t + ph + 18) +
               '" text-anchor="middle" ' + font + '>' + t.label + '</text>');
    });
    ticks(0, ymax, 5).forEach(function (t) {
      var y = Y(t.v).toFixed(1);
      out.push('<line x1="' + m.l + '" x2="' + (m.l + pw) + '" y1="' + y +
               '" y2="' + y + '" stroke="#ebebeb"/>');
      out.push('<text x="' + (m.l - 6) + '" y="' + y +
               '" dy="0.35em" text-anchor="end" ' + font + '>' + t.label +
               '</text>');
    });

    // Size threshold z*
    var xs = X(data.zstar).toFixed(1);
    out.push('<line x1="' + xs + '" x2="' + xs + '" y1="' + m.t + '" y2="' +
             (m.t + ph) + '" stroke="grey" stroke-dasharray="6,4"/>');

    // Population curves: area under each curve equals population density
    var species = [
      { N: s.NC, mu: s.muC, sd: data.sC, col: data.col_cod },
      { N: s.NF, mu: 0, sd: data.sF, col: data.col_flounder }
    ];
    species.forEach(function (sp) {
      var pts = [];
      for (var i = 0; i <= 400; i++) {
        var z = z0 + (z1 - z0) * i / 400;
        var y = Math.min(sp.N * dnorm(z, sp.mu, sp.sd), ymax);
        pts.push(X(z).toFixed(1) + "," + Y(y).toFixed(1));
      }
      var line = "M" + pts.join("L");
      out.push('<path d="' + line + "L" + X(z1).toFixed(1) + "," + Y(0) +
               "L" + X(z0).toFixed(1) + "," + Y(0) + 'Z" fill="' + sp.col +
               '" fill-opacity="0.25" stroke="none"/>');
      out.push('<path d="' + line + '" fill="none" stroke="' + sp.col +
               '" stroke-width="2"/>');
    });

    // Panel border, axis titles, legend
    out.push('<rect x="' + m.l + '" y="' + m.t + '" width="' + pw +
             '" height="' + ph + '" fill="none" stroke="#333"/>');
    out.push('<text x="' + (m.l + pw / 2) + '" y="' + (H - 10) +
             '" text-anchor="middle" font-size="16">Body size, ' +
             '<tspan font-style="italic">z</tspan></text>');
    out.push('<text transform="translate(18,' + (m.t + ph / 2) +
             ') rotate(-90)" text-anchor="middle" font-size="16">' +
             'Density &#215; trait distribution</text>');
    [["cod", data.col_cod], ["flounder", data.col_flounder]]
      .forEach(function (lg, i) {
        var y = m.t + 16 + 20 * i, x = m.l + pw - 90;
        out.push('<line x1="' + x + '" x2="' + (x + 22) + '" y1="' + y +
                 '" y2="' + y + '" stroke="' + lg[1] + '" stroke-width="3"/>');
        out.push('<text x="' + (x + 28) + '" y="' + y + '" dy="0.35em" ' +
                 'font-size="14">' + lg[0] + '</text>');
      });

    svg.innerHTML = out.join("");
  }

  function drawCursor() {
    var c1 = el("ts-cursor-1"), c2 = el("ts-cursor-2");
    if (!c1 || !c2) return;
    if (!tsMap || !data) { c1.style.display = c2.style.display = "none"; return; }
    var x = tsMap.x[0] + (tNow - tsMap.t[0]) / (tsMap.t[1] - tsMap.t[0]) *
            (tsMap.x[1] - tsMap.x[0]);
    [[c1, tsMap.y1], [c2, tsMap.y2]].forEach(function (cy) {
      var st = cy[0].style;
      st.display = "block";
      st.left = (100 * x) + "%";
      st.top = (100 * (1 - cy[1][1])) + "%";
      st.height = (100 * (cy[1][1] - cy[1][0])) + "%";
    });
  }

  function render() {
    if (!data) return;
    var s = stateAt(tNow);
    var fishing = tNow >= data.fishing[0] && tNow < data.fishing[1];
    var title = el("anim-title");
    title.innerHTML = "Time " + tNow.toFixed(1) + " &#183; " +
      (fishing ? "fishing (&#951; = " + data.eta.toFixed(2) + ")"
               : "no fishing");
    title.style.color = fishing ? data.col_fishing : "";
    el("anim-subtitle").innerHTML =
      "<i>N</i><sub>C</sub> = " + fmt(s.NC) +
      ", &#956;<sub>C</sub> = " + s.muC.toFixed(2) +
      "; <i>N</i><sub>F</sub> = " + fmt(s.NF);
    el("anim-time").value = Math.round(tNow / data.tmax * SLIDER_MAX);
    drawGauss(s);
    drawCursor();
  }

  function setPlaying(on) {
    playing = on;
    lastStamp = null;
    el("anim-play").innerHTML = on
      ? '<i class="fas fa-pause" role="presentation"></i> Pause'
      : '<i class="fas fa-play" role="presentation"></i> Play';
    if (on) requestAnimationFrame(step);
  }

  function step(stamp) {
    if (!playing || !data) return;
    if (lastStamp !== null) {
      var speed = parseFloat(el("anim-speed").value);
      tNow += (stamp - lastStamp) / 1000 * speed * data.tmax / DURATION;
      if (tNow >= data.tmax) { tNow = data.tmax; setPlaying(false); }
    }
    lastStamp = stamp;
    render();
    if (playing) requestAnimationFrame(step);
  }

  $(document).on("click", "#anim-play", function () {
    if (!data) return;
    if (!playing && tNow >= data.tmax) tNow = 0; // replay from the start
    setPlaying(!playing);
  });

  $(document).on("input", "#anim-time", function () {
    if (!data) return;
    tNow = this.value / SLIDER_MAX * data.tmax;
    render();
  });

  Shiny.addCustomMessageHandler("anim-data", function (msg) {
    var failed = !!msg.failed;
    el("anim-failed").style.display = failed ? "block" : "none";
    el("anim-display").style.display = failed ? "none" : "block";
    if (failed) { data = null; setPlaying(false); return; }
    data = msg;
    tNow = Math.min(tNow, data.tmax);
    render();
  });

  Shiny.addCustomMessageHandler("ts-map", function (msg) {
    tsMap = msg;
    drawCursor();
  });

  // Redraw when the tab becomes visible or the window changes size
  $(document).on("shown.bs.tab", render);
  $(window).on("resize", render);
});
