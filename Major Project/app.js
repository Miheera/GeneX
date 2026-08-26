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

// ── State ────────────────────────────────────────────────────
let graphData        = { nodes: [], edges: [] };
let selectedNode     = null;
let perturbationType = 'inhibition';   // default
let d3Simulation     = null;

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
// Inserted after button group (done after DOM ready via insertAfter helper)

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
// PERTURBATION TYPE TOGGLE (injected by JS, see index.html too)
// ============================================================
function setupPerturbationToggle() {
  const btns = document.querySelectorAll('.perturb-btn');
  btns.forEach(btn => {
    btn.addEventListener('click', () => {
      btns.forEach(b => b.classList.remove('active'));
      btn.classList.add('active');
      perturbationType = btn.dataset.type;

      // Show/hide slider for knockout
      const sliderGroup = document.getElementById('slider-group');
      if (perturbationType === 'knockout') {
        sliderGroup.style.display = 'none';
      } else {
        sliderGroup.style.display = 'block';
      }
    });
  });
}

// ============================================================
// GRAPH LOADING
// ============================================================
async function loadGraph(graphId = 'breast_cancer') {
  clearGraph();
  resetMetrics();

  try {
    const res  = await fetch(`${API_BASE}/api/graph/${graphId}`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();

    // Normalise: API returns "edges", D3 needs "links"
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
// D3 GRAPH RENDER
// ============================================================
function renderGraph(data) {
  const svgSel = d3.select('#network-svg');
  svgSel.selectAll('*').remove();

  const rect   = svgEl.getBoundingClientRect();
  const W      = rect.width  || 800;
  const H      = rect.height || 600;

  // Degree bounds → radius mapping
  const degrees  = data.nodes.map(n => n.degree || 1);
  const minDeg   = d3.min(degrees);
  const maxDeg   = d3.max(degrees);
  const rScale   = d3.scaleLinear().domain([minDeg, maxDeg]).range([8, 24]).clamp(true);

  // Edge weight → opacity mapping
  const weights  = data.links.map(l => l.weight || 0.5);
  const minW     = d3.min(weights) ?? 0;
  const maxW     = d3.max(weights) ?? 1;
  const opScale  = d3.scaleLinear().domain([minW, maxW]).range([0.2, 0.8]).clamp(true);

  // ── Zoom container ──
  const root = svgSel.append('g').attr('class', 'zoom-root');

  const zoom = d3.zoom()
    .scaleExtent([0.3, 5])
    .on('zoom', (event) => root.attr('transform', event.transform));
  svgSel.call(zoom);

  // Center on load
  svgSel.call(zoom.transform, d3.zoomIdentity.translate(W / 2, H / 2).scale(1));

  // ── Force simulation ──
  if (d3Simulation) d3Simulation.stop();

  d3Simulation = d3.forceSimulation(data.nodes)
    .force('link',      d3.forceLink(data.links).id(d => d.id).distance(120))
    .force('charge',    d3.forceManyBody().strength(-500))
    .force('center',    d3.forceCenter(0, 0))
    .force('collision', d3.forceCollide(28));

  // ── Links ──
  const linkSel = root.append('g').attr('class', 'links')
    .selectAll('line')
    .data(data.links)
    .join('line')
    .attr('stroke', COLOR_LINK)
    .attr('stroke-width', 1.5)
    .attr('stroke-opacity', d => opScale(d.weight ?? 0.5));

  // ── Nodes ──
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

  // Glow ring (hidden by default, shown on select)
  nodeSel.append('circle')
    .attr('class', 'node-ring')
    .attr('r', d => rScale(d.degree || 1) + 6)
    .attr('fill', 'none')
    .attr('stroke', COLOR_HIGHLIGHT)
    .attr('stroke-width', 2)
    .attr('opacity', 0);

  // Main circle
  nodeSel.append('circle')
    .attr('class', 'node-circle')
    .attr('r', d => rScale(d.degree || 1))
    .attr('fill', COLOR_NORMAL)
    .attr('stroke', 'rgba(255,255,255,0.6)')
    .attr('stroke-width', 1.5);

  // Labels
  nodeSel.append('text')
    .attr('class', 'node-label')
    .text(d => d.label || d.id)
    .attr('dx', d => rScale(d.degree || 1) + 4)
    .attr('dy', '0.35em')
    .attr('font-size', '10px')
    .attr('fill', '#64748b')
    .attr('pointer-events', 'none');

  // ── Tooltip events ──
  nodeSel
    .on('mouseover', (event, d) => {
      tooltip.innerHTML = `<strong>${d.label || d.id}</strong><br/>Degree: ${d.degree ?? '—'}<br/>Score: ${d.disease_score?.toFixed ? d.disease_score.toFixed(3) : (d.score ?? '—')}`;
      tooltip.style.display = 'block';
    })
    .on('mousemove', (event) => {
      tooltip.style.left = `${event.clientX + 14}px`;
      tooltip.style.top  = `${event.clientY - 10}px`;
    })
    .on('mouseleave', () => {
      tooltip.style.display = 'none';
    });

  // ── Click to select ──
  nodeSel.on('click', (event, d) => {
    event.stopPropagation();
    selectNode(d, nodeSel);
  });

  // Deselect on SVG background click
  svgSel.on('click', () => {
    selectedNode = null;
    nodeSel.selectAll('.node-ring').attr('opacity', 0);
    nodeSel.selectAll('.node-circle').attr('fill', COLOR_NORMAL);
    geneNameEl.textContent   = 'None';
    geneDegreeEl.textContent = '—';
    geneScoreEl.textContent  = '—';
  });

  // ── Tick ──
  d3Simulation.on('tick', () => {
    linkSel
      .attr('x1', d => d.source.x)
      .attr('y1', d => d.source.y)
      .attr('x2', d => d.target.x)
      .attr('y2', d => d.target.y);

    nodeSel.attr('transform', d => `translate(${d.x},${d.y})`);
  });

  // Store refs on window for post-simulation access
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

  // Reset all rings, then highlight selected
  nodeSel.selectAll('.node-ring').attr('opacity', 0);
  nodeSel.selectAll('.node-circle').attr('fill', COLOR_NORMAL);

  nodeSel
    .filter(n => n.id === d.id)
    .selectAll('.node-ring')
    .attr('opacity', 1);

  nodeSel
    .filter(n => n.id === d.id)
    .selectAll('.node-circle')
    .attr('fill', COLOR_HIGHLIGHT);
}

// ── Drag helpers ──────────────────────────────────────────────
function dragStart(event, d) {
  if (!event.active) d3Simulation.alphaTarget(0.3).restart();
  d.fx = d.x;
  d.fy = d.y;
}
function dragged(event, d) {
  d.fx = event.x;
  d.fy = event.y;
}
function dragEnd(event, d) {
  if (!event.active) d3Simulation.alphaTarget(0);
  d.fx = null;
  d.fy = null;
}

// ============================================================
// SIMULATE BUTTON
// ============================================================
btnSimulate.addEventListener('click', async () => {
  if (!selectedNode) {
    showError('Please select a gene by clicking a node');
    return;
  }
  hideError();

  const strength = parseFloat(slider.value) / 100;

  // Disable button & show loading state
  btnSimulate.disabled     = true;
  btnSimulate.textContent  = 'Computing…';

  const body = {
    gene_id:          selectedNode.id,
    perturbation_type: perturbationType,
    strength,
  };

  try {
    const res = await fetch(`${API_BASE}/api/perturb`, {
      method:  'POST',
      headers: { 'Content-Type': 'application/json' },
      body:    JSON.stringify(body),
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const result = await res.json();

    updateMetrics(result);
    applyVisualPerturbation(selectedNode, perturbationType, strength);
    renderImpactChart(result.top_impacted || []);
  } catch (err) {
    console.warn('Perturb API unavailable – using mock result.', err);
    const mockResult = generateMockPerturbResult(selectedNode, strength);
    updateMetrics(mockResult);
    applyVisualPerturbation(selectedNode, perturbationType, strength);
    renderImpactChart(mockResult.top_impacted || []);
  } finally {
    btnSimulate.disabled    = false;
    btnSimulate.textContent = 'Simulate Perturbation';
  }
});

// ── Update right-panel metrics ─────────────────────────────────
function updateMetrics(result) {
  const baseline = result.efficiency_baseline ?? 0;
  const current  = result.efficiency_perturbed ?? 0;
  const delta    = result.efficiency_delta_pct ?? 0;
  const impact   = result.impact_score ?? 0;

  metricBaselineEl.textContent = baseline.toFixed(3);
  metricCurrentEl.textContent  = current.toFixed(3);
  metricChangeEl.textContent   = `${delta >= 0 ? '+' : ''}${delta.toFixed(2)}%`;
  metricImpactEl.textContent   = typeof impact === 'number' ? impact.toFixed(3) : impact;

  // Color metric-change
  metricChangeEl.style.color = delta < 0 ? '#ef4444' : '#10b981';
}

// ── Apply visual perturbation to node ──────────────────────────
function applyVisualPerturbation(node, type, strength) {
  const rScale = window._rScale;
  if (!rScale) return;

  const baseR = rScale(node.degree || 1);

  if (type === 'knockout') {
    const nodeSel = window._nodeSel;
    if (!nodeSel) return;

    nodeSel
      .filter(n => n.id === node.id)
      .transition().duration(600)
      .style('opacity', 0)
      .remove();

    window._linkSel
      .filter(l => l.source.id === node.id || l.target.id === node.id)
      .transition().duration(600)
      .style('opacity', 0)
      .remove();

    graphData.nodes = graphData.nodes.filter(n => n.id !== node.id);
    graphData.links = graphData.links.filter(
      l => (l.source.id ?? l.source) !== node.id &&
           (l.target.id ?? l.target) !== node.id
    );
    d3Simulation.nodes(graphData.nodes);
    d3Simulation.force('link').links(graphData.links);
    d3Simulation.alpha(0.3).restart();

    selectedNode = null;
    geneNameEl.textContent   = 'None';
    geneDegreeEl.textContent = '—';
    geneScoreEl.textContent  = '—';

  } else if (type === 'inhibition') {
    // Live DOM query — avoids stale cached selection with D3 "removed" state
    const newR = Math.max(baseR * 0.4, 5);
    const liveGroup = d3.select('#network-svg').selectAll('.node-group')
      .filter(n => n.id === node.id);

    liveGroup.style('opacity', 1);  // reset opacity in case prior knockout faded it

    liveGroup.selectAll('.node-circle')
      .transition().duration(500)
      .attr('r', newR)
      .attr('fill', COLOR_PERTURBED);

    liveGroup.selectAll('.node-ring')
      .transition().duration(500)
      .attr('r', newR + 6);

  } else if (type === 'overexpression') {
    const newR = baseR * (1 + strength * 0.5);
    const liveGroup = d3.select('#network-svg').selectAll('.node-group')
      .filter(n => n.id === node.id);

    liveGroup.style('opacity', 1);

    liveGroup.selectAll('.node-circle')
      .transition().duration(500)
      .attr('r', newR)
      .attr('fill', COLOR_PERTURBED);

    liveGroup.selectAll('.node-ring')
      .transition().duration(500)
      .attr('r', newR + 6);
  }
}

// ============================================================
// IMPACT BAR CHART (D3, horizontal bars)
// ============================================================
function renderImpactChart(topImpacted) {
  // Clear old content (including placeholder)
  impactChartEl.innerHTML = '';

  if (!topImpacted || topImpacted.length === 0) {
    impactChartEl.innerHTML = '<p style="color:#94a3b8;font-size:0.8rem;padding:8px 0;">No impact data available.</p>';
    return;
  }

  const containerWidth = impactChartEl.getBoundingClientRect().width || 260;
  const labelW   = 60;
  const barW     = containerWidth - labelW - 24;
  const barH     = 18;
  const gap      = 10;
  const svgH     = topImpacted.length * (barH + gap);

  const maxDelta = d3.max(topImpacted, d => Math.abs(d.delta_bc || 0)) || 1;
  const xScale   = d3.scaleLinear().domain([0, maxDelta]).range([0, barW]).clamp(true);

  const chartSvg = d3.select(impactChartEl)
    .append('svg')
    .attr('width', '100%')
    .attr('height', svgH)
    .attr('viewBox', `0 0 ${containerWidth} ${svgH}`);

  const row = chartSvg.selectAll('g.bar-entry')
    .data(topImpacted)
    .join('g')
    .attr('class', 'bar-entry')
    .attr('transform', (d, i) => `translate(0, ${i * (barH + gap)})`);

  // Gene label
  row.append('text')
    .text(d => d.gene || d.gene_id || '')
    .attr('x', labelW - 6)
    .attr('y', barH / 2 + 1)
    .attr('text-anchor', 'end')
    .attr('dominant-baseline', 'middle')
    .attr('font-size', '11px')
    .attr('fill', '#475569')
    .attr('font-family', 'Inter, sans-serif');

  // Background track
  row.append('rect')
    .attr('x', labelW)
    .attr('y', 0)
    .attr('width', barW)
    .attr('height', barH)
    .attr('rx', 4)
    .attr('fill', '#f1f5f9');

  // Filled bar (animated)
  row.append('rect')
    .attr('class', 'impact-bar')
    .attr('x', labelW)
    .attr('y', 0)
    .attr('width', 0)
    .attr('height', barH)
    .attr('rx', 4)
    .attr('fill', COLOR_PERTURBED)
    .transition()
    .duration(700)
    .delay((d, i) => i * 60)
    .attr('width', d => xScale(Math.abs(d.delta_bc || 0)));

  // Delta value label
  row.append('text')
    .text(d => {
      const v = d.delta_bc || 0;
      return v.toFixed ? v.toFixed(3) : v;
    })
    .attr('x', labelW + barW + 4)
    .attr('y', barH / 2 + 1)
    .attr('dominant-baseline', 'middle')
    .attr('font-size', '10px')
    .attr('fill', '#94a3b8')
    .attr('font-family', 'Inter, sans-serif');
}

// ============================================================
// RESET BUTTON
// ============================================================
btnReset.addEventListener('click', () => {
  selectedNode   = null;
  slider.value   = 50;
  sliderValueEl.textContent = '50%';
  hideError();
  resetMetrics();

  // Re-enable slider group (in case knockout was selected)
  document.getElementById('slider-group').style.display = 'block';

  const cancerType = cancerTypeSelect.value;
  loadGraph(`${cancerType}_cancer`);
});

// ── Reset metric display ───────────────────────────────────────
function resetMetrics() {
  metricBaselineEl.textContent = '0.842';
  metricCurrentEl.textContent  = '0.842';
  metricChangeEl.textContent   = '0.00%';
  metricChangeEl.style.color   = '';
  metricImpactEl.textContent   = '0';
  geneNameEl.textContent       = 'None';
  geneDegreeEl.textContent     = '—';
  geneScoreEl.textContent      = '—';

  // Restore placeholder chart
  impactChartEl.innerHTML = `
    <div class="chart-placeholder">
      <div class="bar-row"><span>BRCA1</span><div class="bar" style="width:80%"></div></div>
      <div class="bar-row"><span>TP53</span><div class="bar" style="width:65%"></div></div>
      <div class="bar-row"><span>EGFR</span><div class="bar" style="width:50%"></div></div>
      <div class="bar-row"><span>MYC</span><div class="bar" style="width:40%"></div></div>
    </div>`;
}

// ── Clear graph ───────────────────────────────────────────────
function clearGraph() {
  if (d3Simulation) d3Simulation.stop();
  d3.select('#network-svg').selectAll('*').remove();
  window._nodeSel = null;
  window._linkSel = null;
}

// ── Error helpers ─────────────────────────────────────────────
function showError(msg) {
  errorMsg.textContent   = msg;
  errorMsg.style.display = 'block';
}
function hideError() {
  errorMsg.style.display = 'none';
}

// ============================================================
// MOCK DATA (fallback when API is unavailable)
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
    source: s,
    target: t,
    weight: parseFloat((Math.random() * 0.8 + 0.2).toFixed(2)),
  }));

  return { nodes: genes, links };
}

function generateMockPerturbResult(node, strength) {
  const baseline = 0.842;
  const drop     = strength * 0.18 * (node.disease_score || 0.5);
  const current  = Math.max(0, baseline - drop);
  const delta    = ((current - baseline) / baseline) * 100;
  const topGenes = [
    { gene: 'BRCA1', delta_bc: 0.312 * strength },
    { gene: 'TP53',  delta_bc: 0.258 * strength },
    { gene: 'EGFR',  delta_bc: 0.194 * strength },
    { gene: 'MYC',   delta_bc: 0.147 * strength },
    { gene: 'PTEN',  delta_bc: 0.098 * strength },
  ].filter(g => g.gene !== node.id);

  return {
    efficiency_baseline:  baseline,
    efficiency_perturbed: parseFloat(current.toFixed(4)),
    efficiency_delta_pct: parseFloat(delta.toFixed(3)),
    impact_score:         parseFloat((strength * (node.disease_score || 0.5)).toFixed(3)),
    top_impacted:         topGenes,
  };
}

// ============================================================
// INIT
// ============================================================
document.addEventListener('DOMContentLoaded', () => {
  // Wire error message element after button group
  const buttonGroup = document.querySelector('.button-group');
  buttonGroup.after(errorMsg);

  // Set default slider display
  slider.value = 50;
  sliderValueEl.textContent = '50%';

  // Setup perturbation toggle
  setupPerturbationToggle();

  // Load default graph
  loadGraph('breast_cancer');
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
