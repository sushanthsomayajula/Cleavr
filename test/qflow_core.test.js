/* Rigorous test harness for the qFlow compute core.
 *
 * 1. Extracts the exact QFLOW-CORE block from website/qflow.html (single source
 *    of truth -- no reimplementation, so tests can't drift from shipped code).
 * 2. Compares it number-for-number against the real Python functions in
 *    code/receptor_quantification.py and code/flow_cytometry_validation.py.
 * 3. Exercises every other module (units, stats, indices, compensation, FCS).
 *
 * Run:  node test/qflow_core.test.js       (from the repo root)
 * Needs node and python3 on PATH.
 */
const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');

const HTML = process.env.QFLOW_HTML || path.join(__dirname, '..', 'website', 'qflow.html');
const DATA = JSON.parse(fs.readFileSync(path.join(__dirname, 'test_data.json'), 'utf8'));

/* ---- extract the core block and load it ---- */
function loadCore(){
  const html = fs.readFileSync(HTML, 'utf8');
  const startTag = html.indexOf('QFLOW-CORE-START');
  const endTag = html.indexOf('QFLOW-CORE-END');
  if(startTag < 0 || endTag < 0) throw new Error('core markers not found in qflow.html');
  const startIdx = html.indexOf('(function(root, factory){', startTag);
  let chunk = html.slice(startIdx, endTag);
  chunk = chunk.slice(0, chunk.lastIndexOf('});') + 3);
  const fn = new Function('module', 'exports', chunk);
  const m = { exports: {} };
  fn(m, m.exports);
  return m.exports;
}
const CORE = loadCore();

/* ---- tiny assertion framework ---- */
let passed = 0, failed = 0;
const fails = [];
function ok(name, cond, detail){
  if(cond){ passed++; }
  else { failed++; fails.push(name + (detail ? '  ['+detail+']' : '')); }
}
function approx(a, b, tol){ tol = tol===undefined ? 1e-9 : tol; return Math.abs(a-b) <= tol * (1 + Math.abs(b)); }
function okApprox(name, a, b, tol){ ok(name, approx(a,b,tol), `got ${a}, expected ${b}`); }

/* ===================================================================
   1. CROSS-LANGUAGE PARITY WITH THE PYTHON SCRIPTS
   =================================================================== */
const py = JSON.parse(execFileSync('python3', [path.join(__dirname, 'py_parity.py')], {encoding:'utf8'}));

const beads = DATA.beads;
const linBeads = DATA.linear_beads;

const jsCurve = CORE.fitCalibrationCurve(beads, 'auto');
const jsLinCurve = CORE.fitCalibrationCurve(linBeads, 'auto');

ok('parity: auto curve method matches Python', jsCurve.method === py.curve.method, `js=${jsCurve.method} py=${py.curve.method}`);
okApprox('parity: curve slope', jsCurve.slope, py.curve.slope, 1e-9);
okApprox('parity: curve intercept', jsCurve.intercept, py.curve.intercept, 1e-9);
okApprox('parity: curve R^2', jsCurve.r_squared, py.curve.r_squared, 1e-9);

ok('parity: linear curve method matches Python', jsLinCurve.method === py.linear_curve.method, `js=${jsLinCurve.method} py=${py.linear_curve.method}`);
okApprox('parity: linear slope', jsLinCurve.slope, py.linear_curve.slope, 1e-9);
okApprox('parity: linear intercept', jsLinCurve.intercept, py.linear_curve.intercept, 1e-9);

DATA.measurements.forEach((m, i) => {
  const at = CORE.applyCurve(jsCurve, m.mfi);
  const ai = CORE.applyCurve(jsCurve, m.isotype);
  okApprox(`parity: apply_curve target row ${i}`, at, py.quant[i].abc_target, 1e-7);
  okApprox(`parity: apply_curve isotype row ${i}`, ai, py.quant[i].abc_isotype, 1e-7);
  okApprox(`parity: net receptors row ${i}`, at-ai, py.quant[i].net, 1e-7);
});

ok('parity: normalizeCellLine MDA-MB-453', CORE.normalizeCellLine('MDA-MB-453') === py.normalize['MDA-MB-453']);
ok('parity: normalizeCellLine spaces', CORE.normalizeCellLine('mda mb 453') === py.normalize['mda mb 453']);
Object.keys(py.subtype_lookup).forEach(k => {
  const j = CORE.CELL_LINE_SUBTYPE[k];
  ok(`parity: subtype lookup ${k} subtype`, j && j.subtype === py.subtype_lookup[k].subtype);
  ok(`parity: subtype lookup ${k} stable`, j && j.stable === py.subtype_lookup[k].stable);
});

/* ===================================================================
   2. CALIBRATION EDGE CASES
   =================================================================== */
ok('insufficient beads -> null', CORE.fitCalibrationCurve([[1,2],[3,4]], 'auto') === null);
ok('forced linear ignores log-log', CORE.fitCalibrationCurve(beads, 'linear').method === 'linear');
ok('forced log-log ignores linear', CORE.fitCalibrationCurve(beads, 'log-log').method === 'log-log');
ok('applyCurve rejects MFI<=0', CORE.applyCurve(jsCurve, 0) === null && CORE.applyCurve(jsCurve, -5) === null);
// log-log unavailable when a value is non-positive: linear still returned under 'auto'
const negPts = [[10, -5],[20, 100],[30, 400]];
const negCurve = CORE.fitCalibrationCurve(negPts, 'auto');
ok('auto falls back to linear when log-log impossible', negCurve && negCurve.method === 'linear');

// negative ABC extrapolation: a linear curve with a negative x-intercept goes < 0
// for MFI below that intercept -- the "isotype below curve range" case the tool
// flags as an unreliable, non-physical receptor count.
const lc = CORE.fitCalibrationCurve([[100,5000],[200,25000],[300,45000]], 'linear');
okApprox('linear curve has negative intercept', lc.intercept, -15000, 1e-6);
ok('linear curve extrapolates below zero (negative ABC detectable)', CORE.applyCurve(lc, 10) < 0);

/* ===================================================================
   3. UNIT CONVERSION ROUND-TRIP (applyCurve <-> invertCurve)
   =================================================================== */
[jsCurve, jsLinCurve].forEach((c, ci) => {
  [150, 1200, 8400].forEach(mfi => {
    const u = CORE.applyCurve(c, mfi);
    const back = CORE.invertCurve(c, u);
    okApprox(`unit round-trip curve${ci} mfi=${mfi}`, back, mfi, 1e-6);
  });
});
const logCurveForInvert = CORE.fitCalibrationCurve(beads, 'log-log');
ok('invertCurve rejects non-positive units in log-log', CORE.invertCurve(logCurveForInvert, -10) === null);

/* ===================================================================
   4. DESCRIPTIVE STATISTICS
   =================================================================== */
okApprox('median odd', CORE.median([3,1,2]), 2);
okApprox('median even', CORE.median([1,2,3,4]), 2.5);
okApprox('mean', CORE.mean([2,4,6]), 4);
okApprox('percentile p0', CORE.percentile([1,2,3,4,5], 0), 1);
okApprox('percentile p100', CORE.percentile([1,2,3,4,5], 100), 5);
okApprox('percentile p50 == median', CORE.percentile([1,2,3,4,5], 50), 3);
// sample SD of [2,4,4,4,5,5,7,9] is 2.138... (n-1)
okApprox('sample stddev', CORE.stddev([2,4,4,4,5,5,7,9]), 2.13808993529939, 1e-9);
okApprox('geometric mean of [1,10,100]', CORE.geometricMean([1,10,100]), Math.pow(1000, 1/3), 1e-9);
ok('geomean ignores non-positive', !isNaN(CORE.geometricMean([-1,0,4,4])));
okApprox('pctPositive gate', CORE.pctPositive([1,2,3,4,5,6,7,8,9,10], 5), 50);
{
  const arr = []; for(let i=1;i<=1000;i++) arr.push(i);
  const rsd = CORE.robustStd(arr);
  ok('robustStd positive and sensible', rsd > 0 && rsd < 1000);
}

/* ===================================================================
   5. STAIN / SEPARATION INDEX
   =================================================================== */
okApprox('stainIndex formula', CORE.stainIndex(1000, 100, 50), (1000-100)/(2*50));
okApprox('separationIndex uses rSD', CORE.separationIndex(1000, 100, 25), (1000-100)/(2*25));
ok('stainIndex NaN on zero SD', isNaN(CORE.stainIndex(1000,100,0)));
{
  const pos = [980,1000,1020,1010,990];
  const neg = [95,100,105,98,102];
  const r = CORE.indicesFromEvents(pos, neg);
  okApprox('indicesFromEvents posMFI', r.posMFI, CORE.median(pos));
  okApprox('indicesFromEvents negMFI', r.negMFI, CORE.median(neg));
  okApprox('indicesFromEvents stainIndex consistency', r.stainIndex, (r.posMFI-r.negMFI)/(2*r.negSD), 1e-9);
}

/* ===================================================================
   6. COMPENSATION (round-trip: true -> measured -> compensated)
   =================================================================== */
{
  const pctS = [[100,12,3],[8,100,5],[2,6,100]];
  const S = CORE.spilloverFromPercent(pctS);
  ok('spillover diagonal is 1', S[0][0]===1 && S[1][1]===1 && S[2][2]===1);
  okApprox('spillover off-diagonal fraction', S[0][1], 0.12);
  const trueVec = [1000, 500, 200];
  const measured = CORE.matVecRow(trueVec, S);       // measured = true * S
  const comp = CORE.compensate(S, measured);          // should recover true
  trueVec.forEach((t,i)=> okApprox(`compensation recovers true[${i}]`, comp[i], t, 1e-7));
  const I = CORE.spilloverFromPercent([[100,0,0],[0,100,0],[0,0,100]]);
  const same = CORE.compensate(I, [7,8,9]);
  ok('identity compensation is a no-op', approx(same[0],7)&&approx(same[1],8)&&approx(same[2],9));
  const inv = CORE.matInvert(S);
  const prod = S.map((row,i)=> inv[0].map((_,j)=> row.reduce((s,_,k)=> s + S[i][k]*inv[k][j], 0)));
  for(let i=0;i<3;i++) for(let j=0;j<3;j++) okApprox(`S*S^-1 identity [${i}][${j}]`, prod[i][j], i===j?1:0, 1e-7);
  ok('singular matrix -> null inverse', CORE.matInvert([[1,2],[2,4]]) === null);
}

/* ===================================================================
   7. FCS PARSER (build a synthetic FCS 3.1 file, parse it back)
   =================================================================== */
function buildSyntheticFCS(){
  const params = [
    {name:'FSC-A', label:'', bits:32, range:262144},
    {name:'PE-A',  label:'CD44', bits:32, range:262144},
    {name:'APC-A', label:'CD24', bits:32, range:262144},
  ];
  const events = [
    [1000.5, 880.0, 45.0],
    [ 950.0, 910.0, 50.0],
    [1100.0, 120.0, 600.0],
    [1020.0, 890.0, 47.0],
  ];
  const nPar = params.length, nTot = events.length;
  const delim = '/';
  const kv = {
    '$BEGINANALYSIS':'0','$ENDANALYSIS':'0','$BEGINSTEXT':'0','$ENDSTEXT':'0',
    '$BYTEORD':'1,2,3,4','$DATATYPE':'F','$MODE':'L','$NEXTDATA':'0',
    '$PAR':String(nPar),'$TOT':String(nTot),
  };
  params.forEach((p,i)=>{ const n=i+1; kv['$P'+n+'N']=p.name; kv['$P'+n+'B']=String(p.bits); kv['$P'+n+'R']=String(p.range); kv['$P'+n+'E']='0,0'; if(p.label) kv['$P'+n+'S']=p.label; });
  let body = '';
  for(const k in kv){ body += k + delim + kv[k] + delim; }
  const HEADER = 58;
  const textStart = HEADER;
  let dataStartField = '00000000', dataEndField='00000000';
  let textWithOffsets, textEnd, dataStart, dataEnd;
  for(let iter=0; iter<5; iter++){
    textWithOffsets = delim + '$BEGINDATA' + delim + dataStartField + delim + '$ENDDATA' + delim + dataEndField + delim + body.slice(1);
    textEnd = textStart + Buffer.byteLength(textWithOffsets, 'latin1') - 1;
    dataStart = textEnd + 1;
    dataEnd = dataStart + nPar*nTot*4 - 1;
    dataStartField = String(dataStart).padStart(8,'0');
    dataEndField = String(dataEnd).padStart(8,'0');
  }
  const textBuf = Buffer.from(textWithOffsets, 'latin1');
  const dataBuf = Buffer.alloc(nPar*nTot*4);
  let off=0;
  events.forEach(ev=> ev.forEach(v=>{ dataBuf.writeFloatLE(v, off); off+=4; }));
  const header = Buffer.alloc(HEADER, 0x20);
  header.write('FCS3.1', 0, 'latin1');
  const w8 = (val, pos)=> header.write(String(val).padStart(8,' '), pos, 'latin1');
  w8(textStart, 10); w8(textEnd, 18); w8(dataStart, 26); w8(dataEnd, 34); w8(0, 42); w8(0, 50);
  return { buffer: Buffer.concat([header, textBuf, dataBuf]), params, events };
}
{
  const { buffer, params, events } = buildSyntheticFCS();
  const parsed = CORE.parseFCS(buffer);
  ok('FCS: version', parsed.version.slice(0,6) === 'FCS3.1', parsed.version);
  ok('FCS: event count', parsed.nEvents === events.length, `got ${parsed.nEvents}`);
  ok('FCS: param count', parsed.params.length === params.length);
  ok('FCS: param names', parsed.params.map(p=>p.name).join(',') === 'FSC-A,PE-A,APC-A', parsed.params.map(p=>p.name).join(','));
  ok('FCS: label read from $PnS', parsed.params[1].label === 'CD44');
  let allMatch = true;
  events.forEach((ev,e)=> ev.forEach((v,p)=>{ if(!approx(parsed.matrix[e][p], v, 1e-4)) allMatch=false; }));
  ok('FCS: event matrix values match', allMatch);
  ok('FCS: byParam column extraction', approx(CORE.median(parsed.byParam['PE-A']), CORE.median(events.map(r=>r[1])), 1e-4));
  let threw=false; try{ CORE.parseFCS(Buffer.from('NOTFCSxxxx')); }catch(e){ threw=true; }
  ok('FCS: rejects non-FCS input', threw);
}

/* ===================================================================
   8. channelStats integration
   =================================================================== */
{
  const vals = [10,20,30,40,50,60,70,80,90,100];
  const st = CORE.channelStats(vals, 55);
  okApprox('channelStats %positive', st.pctPositive, 50);
  okApprox('channelStats median', st.median, 55);
  ok('channelStats null threshold -> null pctPositive', CORE.channelStats(vals, null).pctPositive === null);
}

/* ---- report ---- */
console.log(`\n${passed} passed, ${failed} failed  (${passed+failed} assertions)\n`);
if(failed){ console.log('FAILURES:'); fails.forEach(f=>console.log('  x '+f)); process.exit(1); }
else { console.log('All qFlow core tests passed.'); }
