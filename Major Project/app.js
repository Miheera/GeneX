// ============================================================
// Cancer Gene Network Explorer — app.js
// Connects to FastAPI at http://localhost:8000
// ============================================================

const API_BASE = 'http://localhost:8000';

// ── Colors ──────────────────────────────────────────────────
const COLOR_NORMAL      = '#4A90D9';
const COLOR_PERTURBED   = '#E85D30';
const COLOR_HIGHLIGHT   = '#FF6B35';
const COLOR_LINK        = '#94a3b8';

// Tolerance chart palette (one per strategy)
const TOL_COLORS = {
  degree:      '#2563eb',
  betweenness: '#7c3aed',
  random:      '#059669',
};

// ── State ────────────────────────────────────────────────────
let graphData        = { nodes: [], edges: [] };
let selectedNode     = null;
let perturbationType = 'inhibition';
let d3Simulation     = null;

// Tolerance state
let tolStrategy  = 'degree';
let tolMetric    = 'efficiency';
let tolDataCache = {};          // keyed by strategy

// Druggability state
let topN = 20;

// ── DOM refs ─────────────────────────────────────────────────
const svgEl              = document.getElementById('network-svg');
const geneNameEl         = document.getElementById('gene-name');
const geneDegreeEl       = document.getElementById('gene-degree');
const geneScoreEl        = document.getElementById('gene-score');
const slider             = document.getElementById('strength-slider');
const sliderValueEl      = document.getElementById('slider-value');
const btnSimulate        = document.getElementById('btn-simulate');
const btnReset           = document.getElementById('btn-reset');
const metricBaselineEl   = document.getElementById('metric-baseline');
const metricCurrentEl    = document.getElementById('metric-current');
const metricChangeEl     = document.getElementById('metric-change');
const metricImpactEl     = document.getElementById('metric-impact');
const impactChartEl      = document.getElementById('impact-chart-container');
const cancerTypeSelect   = document.getElementById('cancer-type');

// ── Tooltip ──────────────────────────────────────────────────
const tooltip = document.createElement('div');
tooltip.id = 'node-tooltip';
Object.assign(tooltip.style, {
  position:        'fixed',
  background:      'rgba(15,23,42,0.92)',
  color:           '#f1f5f9',
  padding:         '8px 12px',
  borderRadius:    '6px',
  fontSize:        '12px',
  lineHeight:      '1.6',
  pointerEvents:   'none',
  zIndex:          '9999',
  display:         'none',
  backdropFilter:  'blur(4px)',
  border:          '1px solid rgba(255,255,255,0.1)',
  boxShadow:       '0 4px 12px rgba(0,0,0,0.3)',
  maxWidth:        '160px',
});
document.body.appendChild(tooltip);

// ── Error message element ─────────────────────────────────────
const errorMsg = document.createElement('p');
errorMsg.id = 'simulate-error';
Object.assign(errorMsg.style, {
  color:      '#ef4444',
  fontSize:   '0.8rem',
  marginTop:  '6px',
  display:    'none',
});

// ── Slider live update ────────────────────────────────────────
slider.addEventListener('input', () => {
  sliderValueEl.textContent = `${slider.value}%`;
});

// ── Cancer type change → reload graph ─────────────────────────
cancerTypeSelect.addEventListener('change', () => {
  const type = cancerTypeSelect.value;
  loadGraph(`${type}_cancer`);
});

// ============================================================
// TAB NAVIGATION
// ============================================================
function setupTabs() {
  const tabs = document.querySelectorAll('.tab-btn');
  const panes = document.querySelectorAll('.view-pane');

  tabs.forEach(tab => {
    tab.addEventListener('click', () => {
      const target = tab.getAttribute('aria-controls');

      tabs.forEach(t => { t.classList.remove('active'); t.setAttribute('aria-selected', 'false'); });
      panes.forEach(p => { p.classList.remove('active'); p.setAttribute('aria-hidden', 'true'); });

      tab.classList.add('active');
      tab.setAttribute('aria-selected', 'true');
      const pane = document.getElementById(target);
      pane.classList.add('active');
      pane.removeAttribute('aria-hidden');

      // Lazy-load on first visit
      if (target === 'view-tolerance' && Object.keys(tolDataCache).length === 0) {
        loadTolerance(tolStrategy);
      }
      if (target === 'view-expression' && !expressionCorrelationData) {
        loadExpressionCorrelation();
      }
    });
  });
}

// ============================================================
// PERTURBATION TYPE TOGGLE (existing)
// ============================================================
function setupPerturbationToggle() {
  const btns = document.querySelectorAll('.perturb-btn');
  btns.forEach(btn => {
    btn.addEventListener('click', () => {
      // Only handle perturbation buttons (not strategy/metric toggles)
      if (!btn.closest('#slider-group') && btn.closest('.control-group')) {
        btns.forEach(b => { if (b.closest('.control-group') && !b.closest('#strategy-toggle') && !b.closest('#metric-toggle')) b.classList.remove('active'); });
        btn.classList.add('active');
        perturbationType = btn.dataset.type;
        const sliderGroup = document.getElementById('slider-group');
        sliderGroup.style.display = perturbationType === 'knockout' ? 'none' : 'block';
      }
    });
  });

  // Handle separately via delegation for each toggle group
  const perturbSection = document.querySelector('.control-group .perturb-toggle');
  if (perturbSection) {
    perturbSection.querySelectorAll('.perturb-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        perturbSection.querySelectorAll('.perturb-btn').forEach(b => b.classList.remove('active'));
        btn.classList.add('active');
        perturbationType = btn.dataset.type;
        const sliderGroup = document.getElementById('slider-group');
        sliderGroup.style.display = perturbationType === 'knockout' ? 'none' : 'block';
      });
    });
  }
}

// ============================================================
// GRAPH LOADING (existing, unchanged)
// ============================================================
async function loadGraph(graphId = 'breast_cancer') {
  clearGraph();
  resetMetrics();

  try {
    const res  = await fetch(`${API_BASE}/api/graph/${graphId}`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();

    graphData = {
      nodes: data.nodes.map(n => ({ ...n })),
      links: (data.edges || data.links || []).map(e => ({ ...e })),
    };

    renderGraph(graphData);
  } catch (err) {
    console.warn('API unavailable – using mock data.', err);
    graphData = generateMockData();
    renderGraph(graphData);
  }
}

// ============================================================
// D3 GRAPH RENDER (existing, unchanged)
// ============================================================
function renderGraph(data) {
  const svgSel = d3.select('#network-svg');
  svgSel.selectAll('*').remove();

  const rect   = svgEl.getBoundingClientRect();
  const W      = rect.width  || 800;
  const H      = rect.height || 600;

  const degrees  = data.nodes.map(n => n.degree || 1);
  const minDeg   = d3.min(degrees);
  const maxDeg   = d3.max(degrees);
  const rScale   = d3.scaleLinear().domain([minDeg, maxDeg]).range([8, 24]).clamp(true);

  const weights  = data.links.map(l => l.weight || 0.5);
  const minW     = d3.min(weights) ?? 0;
  const maxW     = d3.max(weights) ?? 1;
  const opScale  = d3.scaleLinear().domain([minW, maxW]).range([0.2, 0.8]).clamp(true);

  const root = svgSel.append('g').attr('class', 'zoom-root');

  const zoom = d3.zoom()
    .scaleExtent([0.3, 5])
    .on('zoom', (event) => root.attr('transform', event.transform));
  svgSel.call(zoom);
  svgSel.call(zoom.transform, d3.zoomIdentity.translate(W / 2, H / 2).scale(1));

  if (d3Simulation) d3Simulation.stop();

  d3Simulation = d3.forceSimulation(data.nodes)
    .force('link',      d3.forceLink(data.links).id(d => d.id).distance(120))
    .force('charge',    d3.forceManyBody().strength(-500))
    .force('center',    d3.forceCenter(0, 0))
    .force('collision', d3.forceCollide(28));

  const linkSel = root.append('g').attr('class', 'links')
    .selectAll('line')
    .data(data.links)
    .join('line')
    .attr('stroke', COLOR_LINK)
    .attr('stroke-width', 1.5)
    .attr('stroke-opacity', d => opScale(d.weight ?? 0.5));

  const nodeSel = root.append('g').attr('class', 'nodes')
    .selectAll('g')
    .data(data.nodes)
    .join('g')
    .attr('class', 'node-group')
    .style('cursor', 'pointer')
    .call(
      d3.drag()
        .on('start', dragStart)
        .on('drag',  dragged)
        .on('end',   dragEnd)
    );

  nodeSel.append('circle')
    .attr('class', 'node-ring')
    .attr('r', d => rScale(d.degree || 1) + 6)
    .attr('fill', 'none')
    .attr('stroke', COLOR_HIGHLIGHT)
    .attr('stroke-width', 2)
    .attr('opacity', 0);

  nodeSel.append('circle')
    .attr('class', 'node-circle')
    .attr('r', d => rScale(d.degree || 1))
    .attr('fill', COLOR_NORMAL)
    .attr('stroke', 'rgba(255,255,255,0.6)')
    .attr('stroke-width', 1.5);

  nodeSel.append('text')
    .attr('class', 'node-label')
    .text(d => d.label || d.id)
    .attr('dx', d => rScale(d.degree || 1) + 4)
    .attr('dy', '0.35em')
    .attr('font-size', '10px')
    .attr('fill', '#64748b')
    .attr('pointer-events', 'none');

  nodeSel
    .on('mouseover', (event, d) => {
      tooltip.innerHTML = `<strong>${d.label || d.id}</strong><br/>Degree: ${d.degree ?? '—'}<br/>Score: ${d.disease_score?.toFixed ? d.disease_score.toFixed(3) : (d.score ?? '—')}`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', (event) => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top  = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => { tooltip.style.display = 'none'; });

  nodeSel.on('click', (event, d) => {
    event.stopPropagation();
    selectNode(d, nodeSel);
  });

  svgSel.on('click', () => {
    selectedNode = null;
    nodeSel.selectAll('.node-ring').attr('opacity', 0);
    nodeSel.selectAll('.node-circle').attr('fill', COLOR_NORMAL);
    geneNameEl.textContent   = 'None';
    geneDegreeEl.textContent = '—';
    geneScoreEl.textContent  = '—';
  });

  d3Simulation.on('tick', () => {
    linkSel
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y);

    nodeSel.attr('transform', d => `translate(${d.x},${d.y})`);
  });

  window._nodeSel = nodeSel;
  window._linkSel = linkSel;
  window._rScale  = rScale;
}

// ── Node selection ────────────────────────────────────────────
function selectNode(d, nodeSel) {
  selectedNode = d;
  hideError();

  geneNameEl.textContent   = d.label || d.id;
  geneDegreeEl.textContent = d.degree ?? '—';
  geneScoreEl.textContent  = d.disease_score != null
    ? (+d.disease_score).toFixed(3)
    : (d.score ?? '—');

  nodeSel.selectAll('.node-ring').attr('opacity', 0);
  nodeSel.selectAll('.node-circle').attr('fill', COLOR_NORMAL);

  nodeSel.filter(n => n.id === d.id).selectAll('.node-ring').attr('opacity', 1);
  nodeSel.filter(n => n.id === d.id).selectAll('.node-circle').attr('fill', COLOR_HIGHLIGHT);
}

// ── Drag helpers ──────────────────────────────────────────────
function dragStart(event, d) {
  if (!event.active) d3Simulation.alphaTarget(0.3).restart();
  d.fx = d.x; d.fy = d.y;
}
function dragged(event, d) { d.fx = event.x; d.fy = event.y; }
function dragEnd(event, d) {
  if (!event.active) d3Simulation.alphaTarget(0);
  d.fx = null; d.fy = null;
}

// ============================================================
// SIMULATE BUTTON (existing, unchanged)
// ============================================================
btnSimulate.addEventListener('click', async () => {
  if (!selectedNode) { showError('Please select a gene by clicking a node'); return; }
  hideError();
  const strength = parseFloat(slider.value) / 100;
  btnSimulate.disabled    = true;
  btnSimulate.textContent = 'Computing…';

  const body = { gene_id: selectedNode.id, perturbation_type: perturbationType, strength };

  try {
    const res = await fetch(`${API_BASE}/api/perturb`, {
      method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body),
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const result = await res.json();
    updateMetrics(result);
    applyVisualPerturbation(selectedNode, perturbationType, strength);
    const propagated = Boolean(result.propagation_model);
    renderImpactChart(propagated ? result.top_impacted_genes : result.top_impacted || [], propagated);
  } catch (err) {
    console.error('Perturbation API request failed.', err);
    clearPerturbationResults();
    showError('Perturbation could not be computed. Check the API connection and try again.');
  } finally {
    btnSimulate.disabled    = false;
    btnSimulate.textContent = 'Simulate Perturbation';
  }
});

function clearPerturbationResults() {
  document.getElementById('metric-baseline-label').textContent = 'Baseline Efficiency';
  document.getElementById('metric-current-label').textContent = 'Current Efficiency';
  document.getElementById('metric-change-label').textContent = 'Percentage Change';
  document.getElementById('metric-impact-label').textContent = 'Structural Impact';
  document.getElementById('perturbation-results-note').textContent = 'No result available for this request.';
  document.getElementById('impact-chart-title').textContent = 'Perturbation Results';
  metricBaselineEl.textContent = '—';
  metricCurrentEl.textContent = '—';
  metricChangeEl.textContent = '—';
  metricImpactEl.textContent = '—';
  impactChartEl.innerHTML = '<p class="chart-placeholder">No result available.</p>';
}

function updateMetrics(result) {
  const labels = {
    baseline: document.getElementById('metric-baseline-label'),
    current: document.getElementById('metric-current-label'),
    change: document.getElementById('metric-change-label'),
    impact: document.getElementById('metric-impact-label'),
  };
  const note = document.getElementById('perturbation-results-note');
  const chartTitle = document.getElementById('impact-chart-title');
  if (result.propagation_model) {
    labels.baseline.textContent = 'Target';
    labels.current.textContent = 'Perturbation';
    labels.change.textContent = 'Magnitude';
    labels.impact.textContent = 'Propagation Model';
    metricBaselineEl.textContent = result.target_gene;
    metricCurrentEl.textContent = result.perturbation_type;
    metricChangeEl.textContent = `${(result.perturbation_strength * 100).toFixed(0)}%`;
    metricImpactEl.textContent = 'Confidence-weighted RWR';
    metricChangeEl.style.color = result.perturbation_type === 'inhibition' ? '#ef4444' : '#10b981';
    note.textContent = `${result.perturbation_type === 'inhibition' ? 'Negative' : 'Positive'} target-centered signal on a fixed graph. STRING confidence is unchanged. Scores describe topological proximity, not calibrated biological response.`;
    chartTitle.textContent = 'Top Propagated Nodes';
    return;
  }

  labels.baseline.textContent = 'Efficiency Before';
  labels.current.textContent = 'Efficiency After';
  labels.change.textContent = 'ΔE%';
  labels.impact.textContent = 'Structural Impact';
  metricBaselineEl.textContent = (result.efficiency_baseline ?? 0).toFixed(3);
  metricCurrentEl.textContent = (result.efficiency_perturbed ?? 0).toFixed(3);
  const delta = result.efficiency_delta_pct ?? 0;
  metricChangeEl.textContent = `${delta >= 0 ? '+' : ''}${delta.toFixed(2)}%`;
  metricImpactEl.textContent = Number(result.impact_score ?? 0).toFixed(3);
  metricChangeEl.style.color = delta < 0 ? '#ef4444' : '#10b981';
  note.textContent = `Virtual Topological Knockout: ${result.connected_components ?? '—'} components; LCC ${result.lcc_size ?? '—'} nodes (${((result.lcc_fraction ?? 0) * 100).toFixed(1)}%). Structural ΔE% is not a biological dose-response prediction.`;
  chartTitle.textContent = 'Largest Betweenness Changes';
}

function applyVisualPerturbation(node, type, strength) {
  const rScale = window._rScale;
  if (!rScale) return;
  const baseR = rScale(node.degree || 1);

  if (type === 'knockout') {
    const nodeSel = window._nodeSel;
    if (!nodeSel) return;
    nodeSel.filter(n => n.id === node.id).transition().duration(600).style('opacity', 0).remove();
    window._linkSel
      .filter(l => l.source.id === node.id || l.target.id === node.id)
      .transition().duration(600).style('opacity', 0).remove();
    graphData.nodes = graphData.nodes.filter(n => n.id !== node.id);
    graphData.links = graphData.links.filter(l => (l.source.id ?? l.source) !== node.id && (l.target.id ?? l.target) !== node.id);
    d3Simulation.nodes(graphData.nodes);
    d3Simulation.force('link').links(graphData.links);
    d3Simulation.alpha(0.3).restart();
    selectedNode = null;
    geneNameEl.textContent = 'None'; geneDegreeEl.textContent = '—'; geneScoreEl.textContent = '—';

  } else if (type === 'inhibition') {
    const newR = baseR * (1 + strength * 0.5);
    const liveGroup = d3.select('#network-svg').selectAll('.node-group').filter(n => n.id === node.id);
    liveGroup.style('opacity', 1);
    liveGroup.selectAll('.node-circle').transition().duration(500).attr('r', newR).attr('fill', COLOR_PERTURBED);
    liveGroup.selectAll('.node-ring').transition().duration(500).attr('r', newR + 6);

  } else if (type === 'overexpression') {
    const newR = baseR * (1 + strength * 0.5);
    const liveGroup = d3.select('#network-svg').selectAll('.node-group').filter(n => n.id === node.id);
    liveGroup.style('opacity', 1);
    liveGroup.selectAll('.node-circle').transition().duration(500).attr('r', newR).attr('fill', COLOR_PERTURBED);
    liveGroup.selectAll('.node-ring').transition().duration(500).attr('r', newR + 6);
  }
}

// ============================================================
// IMPACT BAR CHART (existing, unchanged)
// ============================================================
function renderImpactChart(topImpacted, isPropagation = false) {
  impactChartEl.innerHTML = '';
  if (!topImpacted || topImpacted.length === 0) {
    impactChartEl.innerHTML = '<p style="color:#94a3b8;font-size:0.8rem;padding:8px 0;">No impact data available.</p>';
    return;
  }
  const containerWidth = impactChartEl.getBoundingClientRect().width || 260;
  const labelW = 60, barW = containerWidth - labelW - 24, barH = 18, gap = 10;
  const svgH = topImpacted.length * (barH + gap);
  const scoreOf = d => isPropagation ? d.propagation_score || 0 : d.delta_bc || 0;
  const maxDelta = d3.max(topImpacted, d => Math.abs(scoreOf(d))) || 1;
  const xScale = d3.scaleLinear().domain([0, maxDelta]).range([0, barW]).clamp(true);

  const chartSvg = d3.select(impactChartEl).append('svg')
    .attr('width', '100%').attr('height', svgH).attr('viewBox', `0 0 ${containerWidth} ${svgH}`);

  const row = chartSvg.selectAll('g.bar-entry').data(topImpacted).join('g')
    .attr('class', 'bar-entry').attr('transform', (d, i) => `translate(0, ${i * (barH + gap)})`);

  row.append('text').text(d => d.gene || d.gene_id || '')
    .attr('x', labelW - 6).attr('y', barH / 2 + 1)
    .attr('text-anchor', 'end').attr('dominant-baseline', 'middle')
    .attr('font-size', '11px').attr('fill', '#475569').attr('font-family', 'Inter, sans-serif');

  row.append('rect').attr('x', labelW).attr('y', 0).attr('width', barW).attr('height', barH)
    .attr('rx', 4).attr('fill', '#f1f5f9');

  row.append('rect').attr('class', 'impact-bar').attr('x', labelW).attr('y', 0)
    .attr('width', 0).attr('height', barH).attr('rx', 4).attr('fill', COLOR_PERTURBED)
    .transition().duration(700).delay((d, i) => i * 60)
    .attr('width', d => xScale(Math.abs(d.delta_bc || 0)));

  row.append('text').text(d => { const v = scoreOf(d); return v.toFixed(3); })
    .attr('x', labelW + barW + 4).attr('y', barH / 2 + 1)
    .attr('dominant-baseline', 'middle').attr('font-size', '10px')
    .attr('fill', '#94a3b8').attr('font-family', 'Inter, sans-serif');
}

// ============================================================
// RESET BUTTON (existing, unchanged)
// ============================================================
btnReset.addEventListener('click', () => {
  selectedNode = null; slider.value = 50; sliderValueEl.textContent = '50%';
  hideError(); resetMetrics();
  document.getElementById('slider-group').style.display = 'block';
  const cancerType = cancerTypeSelect.value;
  loadGraph(`${cancerType}_cancer`);
});

function resetMetrics() {
  document.getElementById('metric-baseline-label').textContent = 'Baseline Efficiency';
  document.getElementById('metric-current-label').textContent = 'Current Efficiency';
  document.getElementById('metric-change-label').textContent = 'Percentage Change';
  document.getElementById('metric-impact-label').textContent = 'Structural Impact';
  document.getElementById('perturbation-results-note').textContent = 'Run a perturbation to view results.';
  document.getElementById('impact-chart-title').textContent = 'Top Impacted Genes';
  metricBaselineEl.textContent = '—';
  metricCurrentEl.textContent  = '—';
  metricChangeEl.textContent   = '—';
  metricChangeEl.style.color   = '';
  metricImpactEl.textContent   = '—';
  geneNameEl.textContent       = 'None';
  geneDegreeEl.textContent     = '—';
  geneScoreEl.textContent      = '—';
  impactChartEl.innerHTML = '<p class="chart-placeholder">No perturbation results yet.</p>';
}

function clearGraph() {
  if (d3Simulation) d3Simulation.stop();
  d3.select('#network-svg').selectAll('*').remove();
  window._nodeSel = null; window._linkSel = null;
}

function showError(msg) { errorMsg.textContent = msg; errorMsg.style.display = 'block'; }
function hideError() { errorMsg.style.display = 'none'; }

// ============================================================
// ★ ATTACK TOLERANCE (new)
// ============================================================
function setupToleranceControls() {
  // Strategy toggle
  document.getElementById('strategy-toggle').querySelectorAll('.perturb-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('strategy-toggle').querySelectorAll('.perturb-btn')
        .forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      tolStrategy = btn.dataset.strategy;
      loadTolerance(tolStrategy);
    });
  });

  // Metric toggle (efficiency / lcc_fraction)
  document.getElementById('metric-toggle').querySelectorAll('.perturb-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('metric-toggle').querySelectorAll('.perturb-btn')
        .forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      tolMetric = btn.dataset.metric;
      // Re-render from cache if available
      if (tolDataCache[tolStrategy]) {
        renderToleranceChart(tolDataCache[tolStrategy], tolStrategy, tolMetric);
      }
    });
  });
}

async function loadTolerance(strategy) {
  const loadingEl = document.getElementById('tol-loading');
  const errorEl   = document.getElementById('tol-error');
  const chartEl   = document.getElementById('tol-chart');

  loadingEl.style.display = 'flex';
  errorEl.style.display   = 'none';
  chartEl.innerHTML       = '';

  // Update title
  const stratLabel = { degree: 'Degree', betweenness: 'Betweenness', random: 'Random' }[strategy];
  document.getElementById('tol-chart-title').textContent = `Attack Tolerance — ${stratLabel} Removal`;
  document.getElementById('tol-chart-subtitle').textContent =
    tolMetric === 'efficiency'
      ? 'Global efficiency vs. fraction of nodes removed'
      : 'Largest connected component fraction vs. fraction of nodes removed';

  // Use cache if available
  if (tolDataCache[strategy]) {
    loadingEl.style.display = 'none';
    renderToleranceChart(tolDataCache[strategy], strategy, tolMetric);
    return;
  }

  try {
    const res = await fetch(`${API_BASE}/api/analysis/attack-tolerance?strategy=${strategy}`);
    if (!res.ok) throw new Error(`Server returned HTTP ${res.status}`);
    const data = await res.json();
    tolDataCache[strategy] = data;
    loadingEl.style.display = 'none';
    renderToleranceChart(data, strategy, tolMetric);
  } catch (err) {
    loadingEl.style.display = 'none';
    errorEl.style.display   = 'block';
    errorEl.textContent     = `Failed to load attack tolerance data: ${err.message}`;
    console.error('Attack tolerance fetch error:', err);
  }
}

function renderToleranceChart(data, strategy, metric) {
  const container = document.getElementById('tol-chart');
  container.innerHTML = '';

  if (!data || data.length === 0) {
    container.innerHTML = '<p style="color:#94a3b8;padding:32px 0;text-align:center;">No data returned from server.</p>';
    return;
  }

  // Update subtitle dynamically
  document.getElementById('tol-chart-subtitle').textContent =
    metric === 'efficiency'
      ? 'Global efficiency vs. fraction of nodes removed'
      : 'Largest connected component fraction vs. fraction of nodes removed';

  const margin = { top: 20, right: 40, bottom: 56, left: 60 };
  const rect   = container.getBoundingClientRect();
  const W      = Math.max(rect.width  || 700, 400);
  const H      = Math.max(rect.height || 360, 300);
  const innerW = W - margin.left - margin.right;
  const innerH = H - margin.top  - margin.bottom;

  const color  = TOL_COLORS[strategy] || '#2563eb';
  const isRand = strategy === 'random';
  const yKey   = metric;
  const stdKey = metric === 'efficiency' ? 'efficiency_std' : 'lcc_fraction_std';

  const xScale = d3.scaleLinear()
    .domain([0, d3.max(data, d => +d.fraction_removed) || 1])
    .range([0, innerW]);

  const yMax = isRand && data[0][stdKey] != null
    ? d3.max(data, d => +d[yKey] + (+d[stdKey] || 0))
    : d3.max(data, d => +d[yKey]);

  const yScale = d3.scaleLinear()
    .domain([0, Math.min((yMax || 1) * 1.05, 1.0)])
    .range([innerH, 0]);

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`)
    .style('overflow', 'visible');

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

  // Grid lines
  g.append('g').attr('class', 'grid')
    .call(d3.axisLeft(yScale).ticks(5).tickSize(-innerW).tickFormat(''))
    .selectAll('line').attr('stroke', '#e2e8f0').attr('stroke-dasharray', '3,3');
  g.select('.grid .domain').remove();

  // Shaded band (random ± 1 std)
  if (isRand && data[0][stdKey] != null) {
    const area = d3.area()
      .x(d => xScale(+d.fraction_removed))
      .y0(d => yScale(Math.max(0, +d[yKey] - (+d[stdKey] || 0))))
      .y1(d => yScale(Math.min(1, +d[yKey] + (+d[stdKey] || 0))))
      .curve(d3.curveMonotoneX);

    g.append('path').datum(data).attr('class', 'std-band')
      .attr('fill', color).attr('fill-opacity', 0.12)
      .attr('d', area);
  }

  // Main line
  const line = d3.line()
    .x(d => xScale(+d.fraction_removed))
    .y(d => yScale(+d[yKey]))
    .curve(d3.curveMonotoneX);

  const path = g.append('path').datum(data)
    .attr('fill', 'none')
    .attr('stroke', color)
    .attr('stroke-width', 2.5)
    .attr('d', line);

  // Animate line draw
  const totalLen = path.node().getTotalLength();
  path.attr('stroke-dasharray', totalLen)
      .attr('stroke-dashoffset', totalLen)
      .transition().duration(900).ease(d3.easeCubicOut)
      .attr('stroke-dashoffset', 0);

  // Dots
  g.selectAll('.dot').data(data).join('circle')
    .attr('class', 'dot')
    .attr('cx', d => xScale(+d.fraction_removed))
    .attr('cy', d => yScale(+d[yKey]))
    .attr('r', 3.5)
    .attr('fill', color)
    .attr('stroke', 'white')
    .attr('stroke-width', 1.5)
    .style('opacity', 0)
    .transition().delay(800).duration(200)
    .style('opacity', 1);

  // Tooltip on dots
  g.selectAll('.dot-hit').data(data).join('circle')
    .attr('class', 'dot-hit')
    .attr('cx', d => xScale(+d.fraction_removed))
    .attr('cy', d => yScale(+d[yKey]))
    .attr('r', 10)
    .attr('fill', 'transparent')
    .on('mouseover', (event, d) => {
      const val = (+d[yKey]).toFixed(4);
      const std = d[stdKey] != null ? ` ± ${(+d[stdKey]).toFixed(4)}` : '';
      tooltip.innerHTML = `<strong>Removed: ${(+d.fraction_removed * 100).toFixed(0)}%</strong><br/>${metric === 'efficiency' ? 'Efficiency' : 'LCC'}: ${val}${std}`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', event => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top  = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => { tooltip.style.display = 'none'; });

  // X axis
  g.append('g').attr('transform', `translate(0,${innerH})`)
    .call(d3.axisBottom(xScale).ticks(6).tickFormat(d => `${Math.round(d * 100)}%`))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');

  // Y axis
  g.append('g').call(d3.axisLeft(yScale).ticks(5).tickFormat(d3.format('.2f')))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');

  g.select('.domain').attr('stroke', '#e2e8f0');

  // Axis labels
  g.append('text').attr('text-anchor', 'middle')
    .attr('x', innerW / 2).attr('y', innerH + 46)
    .attr('fill', '#64748b').attr('font-size', '12px')
    .text('Fraction of nodes removed');

  g.append('text').attr('text-anchor', 'middle')
    .attr('transform', `translate(-44, ${innerH / 2}) rotate(-90)`)
    .attr('fill', '#64748b').attr('font-size', '12px')
    .text(metric === 'efficiency' ? 'Global Efficiency' : 'LCC Fraction');

  // Legend
  const legendEl = document.getElementById('tol-legend');
  const stratLabel = { degree: 'Degree', betweenness: 'Betweenness', random: 'Random' }[strategy];
  legendEl.innerHTML = `
    <div class="legend-line">
      <div class="legend-line-swatch" style="background:${color};"></div>
      <span>${stratLabel}</span>
    </div>
    ${isRand && data[0][stdKey] != null ? `
    <div class="legend-line">
      <div class="legend-line-swatch" style="background:${color};opacity:0.25;height:12px;border-radius:3px;"></div>
      <span>± 1 std</span>
    </div>` : ''}
  `;
}

// ============================================================
// ★ DRUGGABILITY / OPEN TARGETS (new)
// ============================================================
function setupDruggabilityControls() {
  const topnSlider = document.getElementById('topn-slider');
  const topnLabel  = document.getElementById('topn-value');

  topnSlider.addEventListener('input', () => {
    topN = parseInt(topnSlider.value);
    topnLabel.textContent = topN;
  });

  document.getElementById('btn-load-druggable').addEventListener('click', () => {
    loadDruggability(topN);
  });

  // Sub-tabs
  document.querySelectorAll('.drug-subtab').forEach(tab => {
    tab.addEventListener('click', () => {
      document.querySelectorAll('.drug-subtab').forEach(t => t.classList.remove('active'));
      document.querySelectorAll('.drug-table-panel').forEach(p => p.style.display = 'none');
      tab.classList.add('active');
      document.getElementById(tab.dataset.target).style.display = 'block';
    });
  });
}

async function loadDruggability(topN) {
  const loadingEl  = document.getElementById('drug-loading');
  const errorEl    = document.getElementById('drug-error');
  const contentEl  = document.getElementById('drug-content');
  const noteEl     = document.getElementById('drug-methodology-note');

  loadingEl.style.display = 'flex';
  errorEl.style.display   = 'none';
  contentEl.style.display = 'none';
  noteEl.style.display    = 'none';

  try {
    const res = await fetch(`${API_BASE}/api/analysis/open-targets?top_n=${topN}`);
    if (!res.ok) throw new Error(`Server returned HTTP ${res.status}`);
    const data = await res.json();

    loadingEl.style.display = 'none';
    contentEl.style.display = 'block';

    // Methodology note banner
    const noteText = data.methodology_note || (data.druggable?.[0]?.methodology_note) || '';
    if (noteText) {
      noteEl.innerHTML = `<strong>⚠ Methodology Note:</strong> ${noteText}`;
      noteEl.style.display = 'block';
    }

    // Update sub-tab counts
    document.getElementById('count-druggable').textContent = (data.druggable || []).length;
    document.getElementById('count-open').textContent      = (data.open_targets || []).length;

    // Max structural score for bar normalisation
    const allScores = [
      ...(data.druggable    || []).map(d => +d.structural_score || 0),
      ...(data.open_targets || []).map(d => +d.structural_score || 0),
    ];
    const maxScore = d3.max(allScores) || 1;

    renderDruggableTable(data.druggable    || [], 'druggable-table-body', true,  maxScore);
    renderDruggableTable(data.open_targets || [], 'open-table-body',      false, maxScore);

  } catch (err) {
    loadingEl.style.display = 'none';
    errorEl.style.display   = 'block';
    errorEl.textContent     = `Failed to load druggability data: ${err.message}`;
    console.error('Druggability fetch error:', err);
  }
}

function renderDruggableTable(rows, containerId, showDrugCols, maxScore) {
  const container = document.getElementById(containerId);
  container.innerHTML = '';

  if (!rows.length) {
    container.innerHTML = `
      <div class="empty-state">
        <svg width="32" height="32" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.5"><circle cx="12" cy="12" r="10"/><path d="M8 12h8M12 8v8"/></svg>
        <p>No genes to display.</p>
      </div>`;
    return;
  }

  const table = document.createElement('table');
  table.className = 'gene-table';

  // Header
  const headerCols = showDrugCols
    ? ['Rank', 'Gene', 'Structural Score', 'Drugs', 'Approved Drug', 'Sources']
    : ['Rank', 'Gene', 'Structural Score', 'Sources'];

  table.innerHTML = `<thead><tr>${headerCols.map(c => `<th>${c}</th>`).join('')}</tr></thead>`;

  const tbody = document.createElement('tbody');

  rows.forEach(row => {
    const tr = document.createElement('tr');

    const drugs = row.drug_names
      ? row.drug_names.split('|').filter(Boolean)
      : [];

    const sources = row.sources
      ? row.sources.split('|').filter(Boolean)
      : [];

    // Interpret "True"/"False" string from the CSV-backed API
    const hasApproved = String(row.has_approved_drug).toLowerCase() === 'true';

    const scoreWidth = maxScore > 0 ? ((+row.structural_score || 0) / maxScore * 100).toFixed(1) : 0;

    const rankCell   = `<td><span class="rank-badge">#${row.rank}</span></td>`;
    const geneCell   = `<td><span class="gene-symbol-cell">${row.gene_symbol}</span></td>`;
    const scoreCell  = `<td class="score-bar-cell">
      <div class="score-bar-wrap">
        <div class="score-bar-track"><div class="score-bar-fill" style="width:${scoreWidth}%"></div></div>
        <span class="score-bar-label">${(+row.structural_score || 0).toFixed(2)}</span>
      </div>
    </td>`;

    const sourcesCell = `<td><div class="source-chips">${
      sources.length
        ? sources.map(s => `<span class="source-chip">${s}</span>`).join('')
        : '<span class="source-chip">—</span>'
    }</div></td>`;

    if (showDrugCols) {
      const drugsCell = `<td><div class="drug-chips">${
        drugs.length
          ? drugs.map(d => `<span class="drug-chip${hasApproved ? ' approved' : ''}">${d}</span>`).join('')
          : '<span style="color:#94a3b8;font-size:0.8rem;">—</span>'
      }</div></td>`;
      const approvedCell = `<td>${
        hasApproved
          ? `<span class="badge-approved"><svg width="11" height="11" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2.5"><path d="M20 6L9 17l-5-5"/></svg> Yes</span>`
          : `<span class="badge-none">No</span>`
      }</td>`;

      tr.innerHTML = rankCell + geneCell + scoreCell + drugsCell + approvedCell + sourcesCell;
    } else {
      tr.innerHTML = rankCell + geneCell + scoreCell + sourcesCell;
    }

    tbody.appendChild(tr);
  });

  table.appendChild(tbody);
  container.appendChild(table);
}

// ============================================================
// MOCK DATA (existing, unchanged)
// ============================================================
function generateMockData() {
  const genes = [
    { id: 'BRCA1', label: 'BRCA1', degree: 9,  disease_score: 0.921 },
    { id: 'TP53',  label: 'TP53',  degree: 11, disease_score: 0.887 },
    { id: 'EGFR',  label: 'EGFR',  degree: 7,  disease_score: 0.754 },
    { id: 'MYC',   label: 'MYC',   degree: 6,  disease_score: 0.712 },
    { id: 'PTEN',  label: 'PTEN',  degree: 5,  disease_score: 0.680 },
    { id: 'MTOR',  label: 'MTOR',  degree: 8,  disease_score: 0.645 },
    { id: 'ATM',   label: 'ATM',   degree: 4,  disease_score: 0.609 },
    { id: 'MAPK1', label: 'MAPK1', degree: 3,  disease_score: 0.573 },
    { id: 'AKT1',  label: 'AKT1',  degree: 6,  disease_score: 0.540 },
    { id: 'PIK3CA',label: 'PIK3CA',degree: 5,  disease_score: 0.508 },
    { id: 'VEGFA', label: 'VEGFA', degree: 4,  disease_score: 0.472 },
    { id: 'JAK2',  label: 'JAK2',  degree: 3,  disease_score: 0.435 },
  ];

  const edgePairs = [
    ['BRCA1','TP53'],['BRCA1','PTEN'],['BRCA1','ATM'],
    ['TP53','EGFR'],['TP53','MYC'],  ['TP53','MTOR'],
    ['EGFR','AKT1'],['EGFR','MAPK1'],['EGFR','PIK3CA'],
    ['MYC','VEGFA'],['MYC','JAK2'], ['PTEN','AKT1'],
    ['MTOR','AKT1'],['PIK3CA','AKT1'],['VEGFA','JAK2'],
    ['ATM','BRCA1'],['MAPK1','AKT1'],
  ];

  const links = edgePairs.map(([s, t]) => ({
    source: s, target: t,
    weight: parseFloat((Math.random() * 0.8 + 0.2).toFixed(2)),
  }));

  return { nodes: genes, links };
}

// ============================================================
// ★ THRESHOLD SENSITIVITY (new)
// ============================================================
let threshData = null;
let threshMetric = 'edge_count';

const THRESH_METRIC_LABELS = {
  edge_count:  'Edge Count',
  node_count:  'Node Count',
  density:     'Network Density',
  mean_degree: 'Mean Degree',
};

const THRESH_METRIC_FMT = {
  edge_count:  d3.format('d'),
  node_count:  d3.format('d'),
  density:     d3.format('.4f'),
  mean_degree: d3.format('.2f'),
};

function setupThresholdControls() {
  document.getElementById('thresh-metric-toggle').querySelectorAll('.perturb-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('thresh-metric-toggle').querySelectorAll('.perturb-btn')
        .forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      threshMetric = btn.dataset.tmetric;
      if (threshData) renderThresholdCharts(threshData);
    });
  });
}

async function loadThreshold() {
  const loadingEl = document.getElementById('thresh-loading');
  const errorEl   = document.getElementById('thresh-error');
  const contentEl = document.getElementById('thresh-content');

  loadingEl.style.display = 'flex';
  errorEl.style.display   = 'none';
  contentEl.style.display = 'none';

  try {
    const res = await fetch(`${API_BASE}/api/analysis/threshold-sensitivity`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const raw = await res.json();

    // Parse all numeric-looking string fields
    threshData = raw.map(row => ({
      threshold:                  +row.threshold,
      node_count:                 +row.node_count,
      edge_count:                 +row.edge_count,
      mean_degree:                parseFloat(row.mean_degree) || 0,
      density:                    parseFloat(row.density) || 0,
      n_connected_components:     +row.n_connected_components,
      is_fully_connected:         String(row.is_fully_connected).toLowerCase() === 'true',
      top10_betweenness_genes:    row.top10_betweenness_genes || '',
      jaccard_vs_shipping_threshold: row.jaccard_vs_shipping_threshold !== '' && row.jaccard_vs_shipping_threshold != null
                                     ? parseFloat(row.jaccard_vs_shipping_threshold)
                                     : null,
      status:        row.status || 'unknown',
      error_message: row.error_message || '',
    }));

    loadingEl.style.display = 'none';
    contentEl.style.display = 'block';
    renderThresholdCharts(threshData);
  } catch (err) {
    loadingEl.style.display = 'none';
    errorEl.style.display   = 'block';
    errorEl.textContent     = `Failed to load threshold data: ${err.message}`;
    console.error('Threshold fetch error:', err);
  }
}

function renderThresholdCharts(data) {
  renderThreshMainChart(data, threshMetric);
  renderThreshJaccardChart(data);
  renderThreshStatusTable(data);
}

function renderThreshMainChart(data, metric) {
  const container = document.getElementById('thresh-main-chart');
  container.innerHTML = '';

  const successData = data.filter(d => d.status === 'success');
  const allData     = data;

  const margin = { top: 16, right: 32, bottom: 50, left: 64 };
  const rect   = container.getBoundingClientRect();
  const W      = Math.max(rect.width  || 680, 400);
  const H      = 220;
  const innerW = W - margin.left - margin.right;
  const innerH = H - margin.top  - margin.bottom;

  const color = '#2563eb';

  const xScale = d3.scaleLinear()
    .domain(d3.extent(allData, d => d.threshold))
    .range([0, innerW]).nice();

  const yMax = d3.max(successData, d => d[metric]) || 1;
  const yScale = d3.scaleLinear().domain([0, yMax * 1.08]).range([innerH, 0]);

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`);

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

  // Grid
  g.append('g').attr('class', 'grid')
    .call(d3.axisLeft(yScale).ticks(4).tickSize(-innerW).tickFormat(''))
    .selectAll('line').attr('stroke', '#e2e8f0').attr('stroke-dasharray', '3,3');
  g.select('.grid .domain').remove();

  // Degraded / non-success vertical bands
  allData.filter(d => d.status !== 'success').forEach(d => {
    g.append('rect')
      .attr('x', xScale(d.threshold) - 10)
      .attr('y', 0)
      .attr('width', 20)
      .attr('height', innerH)
      .attr('fill', d.status === 'failed' ? 'rgba(239,68,68,0.08)' : 'rgba(245,158,11,0.08)');
  });

  // Line (success points only)
  if (successData.length > 1) {
    const line = d3.line()
      .x(d => xScale(d.threshold))
      .y(d => yScale(d[metric]))
      .curve(d3.curveMonotoneX);

    const path = g.append('path').datum(successData)
      .attr('fill', 'none').attr('stroke', color).attr('stroke-width', 2.5)
      .attr('d', line);

    const len = path.node().getTotalLength();
    path.attr('stroke-dasharray', len).attr('stroke-dashoffset', len)
      .transition().duration(800).ease(d3.easeCubicOut)
      .attr('stroke-dashoffset', 0);
  }

  // Dots — colour by status
  const dotColor = d => d.status === 'success' ? color : d.status === 'degraded' ? '#f59e0b' : '#ef4444';

  g.selectAll('.thresh-dot').data(allData).join('circle')
    .attr('class', 'thresh-dot')
    .attr('cx', d => xScale(d.threshold))
    .attr('cy', d => d.status === 'success' ? yScale(d[metric]) : innerH * 0.5)
    .attr('r', 5).attr('fill', dotColor).attr('stroke', 'white').attr('stroke-width', 1.5)
    .on('mouseover', (event, d) => {
      tooltip.innerHTML = `<strong>t = ${d.threshold}</strong><br/>${THRESH_METRIC_LABELS[metric]}: ${d.status === 'success' ? THRESH_METRIC_FMT[metric](d[metric]) : '—'}<br/>Status: ${d.status}`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', event => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top  = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => { tooltip.style.display = 'none'; });

  // Axes
  g.append('g').attr('transform', `translate(0,${innerH})`)
    .call(d3.axisBottom(xScale).ticks(allData.length).tickFormat(d3.format('d')))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');

  g.append('g').call(d3.axisLeft(yScale).ticks(4).tickFormat(THRESH_METRIC_FMT[metric]))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');

  g.append('text').attr('text-anchor', 'middle')
    .attr('x', innerW / 2).attr('y', innerH + 42)
    .attr('fill', '#64748b').attr('font-size', '12px')
    .text('STRING Confidence Threshold');

  g.append('text').attr('text-anchor', 'middle')
    .attr('transform', `translate(-50,${innerH / 2}) rotate(-90)`)
    .attr('fill', '#64748b').attr('font-size', '12px')
    .text(THRESH_METRIC_LABELS[metric]);
}

function renderThreshJaccardChart(data) {
  const container = document.getElementById('thresh-jaccard-chart');
  container.innerHTML = '';

  const jData = data.filter(d => d.jaccard_vs_shipping_threshold != null);
  if (!jData.length) {
    container.innerHTML = '<p style="color:#94a3b8;font-size:0.82rem;padding:16px 0;">Jaccard data not available for any threshold.</p>';
    return;
  }

  const margin = { top: 12, right: 28, bottom: 48, left: 56 };
  const rect   = container.getBoundingClientRect();
  const W      = Math.max(rect.width || 300, 200);
  const H      = 180;
  const innerW = W - margin.left - margin.right;
  const innerH = H - margin.top  - margin.bottom;
  const color  = '#7c3aed';

  const xScale = d3.scaleLinear()
    .domain(d3.extent(data, d => d.threshold)).range([0, innerW]).nice();

  const yScale = d3.scaleLinear().domain([0, 1]).range([innerH, 0]);

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H);

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

  g.append('g').attr('class', 'grid')
    .call(d3.axisLeft(yScale).ticks(4).tickSize(-innerW).tickFormat(''))
    .selectAll('line').attr('stroke', '#e2e8f0').attr('stroke-dasharray', '3,3');
  g.select('.grid .domain').remove();

  if (jData.length > 1) {
    const line = d3.line()
      .x(d => xScale(d.threshold)).y(d => yScale(d.jaccard_vs_shipping_threshold))
      .curve(d3.curveMonotoneX);
    const path = g.append('path').datum(jData)
      .attr('fill', 'none').attr('stroke', color).attr('stroke-width', 2)
      .attr('d', line);
    const len = path.node().getTotalLength();
    path.attr('stroke-dasharray', len).attr('stroke-dashoffset', len)
      .transition().duration(700).attr('stroke-dashoffset', 0);
  }

  g.selectAll('.jacc-dot').data(jData).join('circle')
    .attr('class', 'jacc-dot')
    .attr('cx', d => xScale(d.threshold)).attr('cy', d => yScale(d.jaccard_vs_shipping_threshold))
    .attr('r', 4).attr('fill', color).attr('stroke', 'white').attr('stroke-width', 1.5)
    .on('mouseover', (event, d) => {
      tooltip.innerHTML = `<strong>t = ${d.threshold}</strong><br/>Jaccard: ${d.jaccard_vs_shipping_threshold.toFixed(3)}`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', event => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top  = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => { tooltip.style.display = 'none'; });

  g.append('g').attr('transform', `translate(0,${innerH})`)
    .call(d3.axisBottom(xScale).ticks(jData.length).tickFormat(d3.format('d')))
    .selectAll('text').attr('font-size', '10px').attr('fill', '#64748b');

  g.append('g').call(d3.axisLeft(yScale).ticks(4).tickFormat(d3.format('.2f')))
    .selectAll('text').attr('font-size', '10px').attr('fill', '#64748b');

  g.append('text').attr('text-anchor', 'middle')
    .attr('x', innerW / 2).attr('y', innerH + 40)
    .attr('fill', '#64748b').attr('font-size', '11px').text('Threshold');

  g.append('text').attr('text-anchor', 'middle')
    .attr('transform', `translate(-42,${innerH / 2}) rotate(-90)`)
    .attr('fill', '#64748b').attr('font-size', '11px').text('Jaccard');
}

function renderThreshStatusTable(data) {
  const container = document.getElementById('thresh-status-table');
  container.innerHTML = '';

  const table = document.createElement('table');
  table.className = 'thresh-status-tbl';

  table.innerHTML = `<thead><tr><th>Threshold</th><th>Status</th><th>Nodes</th><th>Edges</th></tr></thead>`;
  const tbody = document.createElement('tbody');

  data.forEach(row => {
    const tr = document.createElement('tr');
    const statusClass = { success: 'success', degraded: 'degraded', empty: 'empty', failed: 'failed' }[row.status] || 'empty';
    const statusIcon  = { success: '✓', degraded: '⚠', empty: '○', failed: '✕' }[row.status] || '?';
    tr.innerHTML = `
      <td style="font-weight:600;font-variant-numeric:tabular-nums;">${row.threshold}</td>
      <td>
        <span class="status-badge ${statusClass}">${statusIcon} ${row.status}</span>
        ${row.error_message ? `<div class="thresh-error-msg">${row.error_message}</div>` : ''}
      </td>
      <td>${row.node_count || '—'}</td>
      <td>${row.edge_count || '—'}</td>
    `;
    tbody.appendChild(tr);
  });

  table.appendChild(tbody);
  container.appendChild(table);
}

// ============================================================
// ★ COSMIC BENCHMARK (new)
// ============================================================
let cosmicData   = null;
let cosmicMetric = 'f1';

const COSMIC_COLORS = {
  betweenness:  '#2563eb',
  degree:       '#7c3aed',
  impact_score: '#059669',
};

const COSMIC_LABELS = {
  betweenness:  'Betweenness',
  degree:       'Degree',
  impact_score: 'Impact Score',
};

function setupCosmicControls() {
  document.getElementById('cosmic-metric-toggle').querySelectorAll('.perturb-btn').forEach(btn => {
    btn.addEventListener('click', () => {
      document.getElementById('cosmic-metric-toggle').querySelectorAll('.perturb-btn')
        .forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      cosmicMetric = btn.dataset.cmetric;
      document.getElementById('cosmic-chart-title').textContent =
        `COSMIC Benchmark — ${cosmicMetric === 'f1' ? 'F1 Score' : cosmicMetric.charAt(0).toUpperCase() + cosmicMetric.slice(1)}`;
      if (cosmicData) renderCosmicChart(cosmicData, cosmicMetric);
    });
  });
}

async function loadCosmic() {
  const loadingEl = document.getElementById('cosmic-loading');
  const errorEl   = document.getElementById('cosmic-error');
  const contentEl = document.getElementById('cosmic-content');

  loadingEl.style.display = 'flex';
  errorEl.style.display   = 'none';
  contentEl.style.display = 'none';

  try {
    const res = await fetch(`${API_BASE}/api/analysis/cosmic-benchmark`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    cosmicData = await res.json();

    loadingEl.style.display = 'none';
    contentEl.style.display = 'block';

    // Methodology note from first row
    const note = cosmicData[0]?.methodology_note || '';
    if (note) {
      const noteEl = document.getElementById('cosmic-methodology-note');
      noteEl.innerHTML = `<strong>ℹ Methodology:</strong> ${note}`;
      noteEl.style.display = 'block';
    }

    // Positives bar
    const ref = cosmicData[0] || {};
    const posBar = document.getElementById('cosmic-positives-bar');
    posBar.innerHTML = `
      <div class="cosmic-stat">
        <span class="cosmic-stat-value">${ref.n_positives_total ?? '—'}</span>
        <span class="cosmic-stat-label">COSMIC breast<br>drivers (total)</span>
      </div>
      <div class="cosmic-divider"></div>
      <div class="cosmic-stat">
        <span class="cosmic-stat-value in-graph">${ref.n_positives_in_graph ?? '—'}</span>
        <span class="cosmic-stat-label">Present in<br>this network</span>
      </div>
      <div class="cosmic-divider"></div>
      <div class="cosmic-stat">
        <span class="cosmic-stat-value excluded">${ref.n_positives_excluded_from_graph ?? '—'}</span>
        <span class="cosmic-stat-label">Excluded<br>(below threshold)</span>
      </div>
      <span style="font-size:0.8rem;color:var(--text-muted);margin-left:8px;">of ${ref.n_positives_total ?? '?'} known COSMIC breast cancer genes, ${ref.n_positives_in_graph ?? '?'} are present in this network (${ref.n_positives_excluded_from_graph ?? '?'} excluded below STRING threshold)</span>
    `;

    renderCosmicChart(cosmicData, cosmicMetric);
  } catch (err) {
    loadingEl.style.display = 'none';
    errorEl.style.display   = 'block';
    errorEl.textContent     = `Failed to load COSMIC benchmark data: ${err.message}`;
    console.error('COSMIC fetch error:', err);
  }
}

function renderCosmicChart(data, metric) {
  const container = document.getElementById('cosmic-chart');
  container.innerHTML = '';

  // Organise data
  const methods  = ['betweenness', 'degree', 'impact_score'];
  const topNVals = [10, 20, 50];

  // Build lookup: data[method][topN] = value
  const lookup = {};
  data.forEach(row => {
    if (!lookup[row.rank_by]) lookup[row.rank_by] = {};
    lookup[row.rank_by][+row.top_n] = +row[metric];
  });

  const margin = { top: 20, right: 32, bottom: 60, left: 64 };
  const rect   = container.getBoundingClientRect();
  const W      = Math.max(rect.width || 680, 400);
  const H      = 300;
  const innerW = W - margin.left - margin.right;
  const innerH = H - margin.top  - margin.bottom;

  // X: grouped by top_n
  const x0 = d3.scaleBand().domain(topNVals.map(String)).range([0, innerW]).padding(0.25);
  const x1 = d3.scaleBand().domain(methods).range([0, x0.bandwidth()]).padding(0.08);

  const allVals = data.map(d => +d[metric]);
  const yScale  = d3.scaleLinear().domain([0, Math.max(d3.max(allVals) || 1, 0.1) * 1.15]).range([innerH, 0]);

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H)
    .attr('viewBox', `0 0 ${W} ${H}`);

  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

  // Grid
  g.append('g').attr('class', 'grid')
    .call(d3.axisLeft(yScale).ticks(5).tickSize(-innerW).tickFormat(''))
    .selectAll('line').attr('stroke', '#e2e8f0').attr('stroke-dasharray', '3,3');
  g.select('.grid .domain').remove();

  // Grouped bars
  const groups = g.selectAll('.bar-group').data(topNVals.map(String)).join('g')
    .attr('class', 'bar-group')
    .attr('transform', n => `translate(${x0(n)},0)`);

  groups.selectAll('rect').data(n => methods.map(m => ({ method: m, value: (lookup[m] || {})[+n] ?? 0, n })))
    .join('rect')
    .attr('x', d => x1(d.method))
    .attr('y', innerH)               // start at baseline
    .attr('width', x1.bandwidth())
    .attr('height', 0)              // animate from 0
    .attr('rx', 3)
    .attr('fill', d => COSMIC_COLORS[d.method])
    .attr('opacity', 0.85)
    .transition().duration(700).delay((d, i) => i * 80)
    .attr('y', d => yScale(d.value))
    .attr('height', d => innerH - yScale(d.value));

  // Value labels on bars (after transition)
  groups.selectAll('.bar-val').data(n => methods.map(m => ({ method: m, value: (lookup[m] || {})[+n] ?? 0, n })))
    .join('text')
    .attr('class', 'bar-val')
    .attr('x', d => x1(d.method) + x1.bandwidth() / 2)
    .attr('y', innerH)
    .attr('text-anchor', 'middle')
    .attr('font-size', '9px')
    .attr('fill', '#64748b')
    .attr('font-family', 'Inter,sans-serif')
    .transition().duration(700).delay((d, i) => i * 80)
    .attr('y', d => yScale(d.value) - 4)
    .text(d => d.value > 0 ? d3.format('.2f')(d.value) : '');

  // Tooltip on bars
  groups.selectAll('.bar-hit').data(n => methods.map(m => ({ method: m, value: (lookup[m] || {})[+n] ?? 0, n })))
    .join('rect')
    .attr('class', 'bar-hit')
    .attr('x', d => x1(d.method))
    .attr('y', 0)
    .attr('width', x1.bandwidth())
    .attr('height', innerH)
    .attr('fill', 'transparent')
    .on('mouseover', (event, d) => {
      tooltip.innerHTML = `<strong>${COSMIC_LABELS[d.method]}</strong><br/>Top-${d.n}: ${d3.format('.4f')(d.value)}`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', event => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top  = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => { tooltip.style.display = 'none'; });

  // X axis
  g.append('g').attr('transform', `translate(0,${innerH})`)
    .call(d3.axisBottom(x0).tickFormat(n => `Top ${n}`))
    .selectAll('text').attr('font-size', '12px').attr('fill', '#64748b');

  // Y axis
  g.append('g').call(d3.axisLeft(yScale).ticks(5).tickFormat(d3.format('.2f')))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');

  g.append('text').attr('text-anchor', 'middle')
    .attr('x', innerW / 2).attr('y', innerH + 52)
    .attr('fill', '#64748b').attr('font-size', '12px').text('Top-N Cutoff');

  const metricLabel = metric === 'f1' ? 'F1 Score' : metric.charAt(0).toUpperCase() + metric.slice(1);
  g.append('text').attr('text-anchor', 'middle')
    .attr('transform', `translate(-50,${innerH / 2}) rotate(-90)`)
    .attr('fill', '#64748b').attr('font-size', '12px').text(metricLabel);

  // Legend
  const legendEl = document.getElementById('cosmic-legend');
  legendEl.innerHTML = methods.map(m =>
    `<div class="legend-line">
       <div class="legend-line-swatch" style="background:${COSMIC_COLORS[m]};height:12px;border-radius:2px;"></div>
       <span>${COSMIC_LABELS[m]}</span>
     </div>`
  ).join('');
}

// ============================================================
// ★ STRUCTURAL IMPORTANCE VS. TCGA EXPRESSION
// ============================================================
let expressionCorrelationData = null;

function correlationDescription(rho) {
  const strength = Math.abs(rho) < 0.2 ? 'negligible' : Math.abs(rho) < 0.4 ? 'weak'
    : Math.abs(rho) < 0.6 ? 'moderate' : Math.abs(rho) < 0.8 ? 'strong' : 'very strong';
  const direction = rho > 0 ? 'positive' : rho < 0 ? 'negative' : 'no';
  return `${strength.charAt(0).toUpperCase() + strength.slice(1)} ${direction} monotonic association.`;
}

async function loadExpressionCorrelation() {
  const loadingEl = document.getElementById('expression-loading');
  const errorEl = document.getElementById('expression-error');
  const contentEl = document.getElementById('expression-content');

  loadingEl.style.display = 'flex';
  errorEl.style.display = 'none';
  contentEl.style.display = 'none';
  try {
    const res = await fetch(`${API_BASE}/api/analysis/structural-expression`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    expressionCorrelationData = await res.json();

    document.getElementById('expression-rho').textContent =
      d3.format('+.3f')(+expressionCorrelationData.spearman_rho);
    document.getElementById('expression-pvalue').textContent =
      d3.format('.4f')(+expressionCorrelationData.spearman_pvalue);
    document.getElementById('expression-gene-count').textContent = expressionCorrelationData.n_genes ?? '—';
    document.getElementById('expression-interpretation').textContent =
      correlationDescription(+expressionCorrelationData.spearman_rho);
    const raw = expressionCorrelationData.raw_expression_correlation || {};
    document.getElementById('expression-raw-correlation').textContent =
      `Raw-expression comparison: ρ ${d3.format('+.3f')(+raw.spearman_rho)}, p ${d3.format('.4f')(+raw.spearman_pvalue)}.`;
    const significanceCaveat = +expressionCorrelationData.spearman_pvalue > 0.05
      ? ' This z-score association is not statistically significant (p > 0.05).'
      : '';
    document.getElementById('expression-methodology-note').innerHTML =
      `<strong>ℹ Methodology:</strong> ${expressionCorrelationData.methodology_note || ''}${significanceCaveat}`;
    renderExpressionScatter(expressionCorrelationData.points || []);

    loadingEl.style.display = 'none';
    contentEl.style.display = 'block';
  } catch (err) {
    loadingEl.style.display = 'none';
    errorEl.style.display = 'block';
    errorEl.textContent = `Failed to load structural-expression analysis: ${err.message}`;
    console.error('Structural-expression fetch error:', err);
  }
}

function renderExpressionScatter(points) {
  const container = document.getElementById('expression-scatter');
  container.innerHTML = '';
  const data = points.map(point => ({
    ...point,
    median_zscore: +point.median_zscore,
    structural_score: +point.structural_score,
  })).filter(point => Number.isFinite(point.median_zscore) && Number.isFinite(point.structural_score));
  if (!data.length) {
    container.innerHTML = '<p class="empty-state">No valid expression-score points were returned.</p>';
    return;
  }

  const margin = { top: 20, right: 28, bottom: 62, left: 70 };
  const rect = container.getBoundingClientRect();
  const W = Math.max(rect.width || 680, 400);
  const H = 380;
  const innerW = W - margin.left - margin.right;
  const innerH = H - margin.top - margin.bottom;
  const zscoreExtent = d3.extent(data, d => d.median_zscore);
  const scoreMax = d3.max(data, d => d.structural_score) || 1;
  const zscorePadding = Math.max((zscoreExtent[1] - zscoreExtent[0]) * 0.1, 0.1);
  const xScale = d3.scaleLinear()
    .domain([zscoreExtent[0] - zscorePadding, zscoreExtent[1] + zscorePadding]).nice().range([0, innerW]);
  const yScale = d3.scaleLinear().domain([0, scoreMax * 1.1]).nice().range([innerH, 0]);

  const svg = d3.select(container).append('svg')
    .attr('width', W).attr('height', H).attr('viewBox', `0 0 ${W} ${H}`);
  const g = svg.append('g').attr('transform', `translate(${margin.left},${margin.top})`);

  g.append('g').attr('class', 'grid')
    .call(d3.axisLeft(yScale).ticks(5).tickSize(-innerW).tickFormat(''))
    .selectAll('line').attr('stroke', '#e2e8f0').attr('stroke-dasharray', '3,3');
  g.select('.grid .domain').remove();

  g.selectAll('circle').data(data).join('circle')
    .attr('cx', d => xScale(d.median_zscore))
    .attr('cy', d => yScale(d.structural_score))
    .attr('r', 0).attr('fill', '#2563eb').attr('fill-opacity', 0.78)
    .attr('stroke', '#1d4ed8').attr('stroke-width', 1)
    .on('mouseover', (event, d) => {
      tooltip.innerHTML = `<strong>${d.gene_symbol}</strong><br/>Median z-score: ${d3.format('.3f')(d.median_zscore)}<br/>Raw median expression: ${d3.format(',.2f')(+d.median_expression)}<br/>Knockout impact: ${d3.format('.2f')(d.structural_score)}<br/>Overexpressed: ${d3.format('.2f')(+d.pct_overexpressed)}%`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', event => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => { tooltip.style.display = 'none'; })
    .transition().duration(500).delay((d, i) => i * 20).attr('r', 5);

  g.append('g').attr('transform', `translate(0,${innerH})`)
    .call(d3.axisBottom(xScale).ticks(6))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');
  g.append('g').call(d3.axisLeft(yScale).ticks(5))
    .selectAll('text').attr('font-size', '11px').attr('fill', '#64748b');
  g.append('text').attr('x', innerW / 2).attr('y', innerH + 52).attr('text-anchor', 'middle')
    .attr('fill', '#64748b').attr('font-size', '12px').text('Median TCGA BRCA expression z-score');
  g.append('text').attr('transform', `translate(-52,${innerH / 2}) rotate(-90)`).attr('text-anchor', 'middle')
    .attr('fill', '#64748b').attr('font-size', '12px').text('Knockout impact score (% efficiency change)');
}


document.addEventListener('DOMContentLoaded', () => {
  // Wire error message element after button group
  const buttonGroup = document.querySelector('.button-group');
  if (buttonGroup) buttonGroup.after(errorMsg);

  // Set defaults
  slider.value = 50;
  sliderValueEl.textContent = '50%';

  // Setup all controls
  setupTabs();
  setupToleranceControls();
  setupDruggabilityControls();
  setupThresholdControls();
  setupCosmicControls();

  // Setup perturbation toggle (network view)
  const perturbButtons = document.querySelectorAll('#view-network .perturb-btn');
  perturbButtons.forEach(btn => {
    btn.addEventListener('click', () => {
      perturbButtons.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      perturbationType = btn.dataset.type;
      const sliderGroup = document.getElementById('slider-group');
      sliderGroup.style.display = perturbationType === 'knockout' ? 'none' : 'block';
    });
  });

  // Load default graph
  loadGraph('breast_cancer');

  // Auto-load on page load
  loadDruggability(20);
  loadThreshold();
  loadCosmic();
});

// Resize handler
window.addEventListener('resize', () => {
  if (!d3Simulation) return;
  const rect = svgEl.getBoundingClientRect();
  d3Simulation
    .force('center', d3.forceCenter(rect.width / 2, rect.height / 2))
    .alpha(0.2)
    .restart();
});
