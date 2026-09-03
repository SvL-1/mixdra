# Markup and styling for validation_report_html(). Kept apart from the
# generator so the R stays readable. Typeset in IBM Plex (serif for headings,
# sans for prose, mono for every figure) -- a family with the right technical
# provenance for a measurement document, and monospaced figures so a workbook
# value and a fitted value line up digit for digit down the page.

.vh_template <- function() r"---(<title>mixdra validation report</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=IBM+Plex+Sans:wght@400;450;600&family=IBM+Plex+Serif:wght@500;600&display=swap">
<style>
:root{
  --paper:#fbfbfd; --panel:#ffffff; --ink:#15181d; --muted:#5c6572;
  --rule:#e2e6ec; --rule-soft:#eef1f5;
  --accent:#1c5f6b; --accent-soft:#dcebee;
  --ok:#2e6b4f; --diverge:#8a5a20; --fail:#a3352c; --quiet:#7a828e;
  --serif:"IBM Plex Serif",Georgia,"Times New Roman",serif;
  --sans:"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
  --mono:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;
}
@media (prefers-color-scheme:dark){
  :root:not([data-theme="light"]){
    --paper:#121419; --panel:#181b21; --ink:#e6e9ee; --muted:#98a1ad;
    --rule:#272c34; --rule-soft:#1f242b;
    --accent:#6fbac7; --accent-soft:#1d3339;
    --ok:#74c099; --diverge:#d8a85f; --fail:#e07a70; --quiet:#8b939f;
  }
}
:root[data-theme="dark"]{
  --paper:#121419; --panel:#181b21; --ink:#e6e9ee; --muted:#98a1ad;
  --rule:#272c34; --rule-soft:#1f242b;
  --accent:#6fbac7; --accent-soft:#1d3339;
  --ok:#74c099; --diverge:#d8a85f; --fail:#e07a70; --quiet:#8b939f;
}

*{box-sizing:border-box}
body{margin:0;background:var(--paper);color:var(--ink);
  font-family:var(--sans);font-size:15px;line-height:1.6;
  -webkit-font-smoothing:antialiased}
.wrap{max-width:56rem;margin:0 auto;padding:3.5rem 1.5rem 5rem;
  display:flex;flex-direction:column;gap:3rem}

/* ---- header ------------------------------------------------------------- */
.head{display:flex;flex-direction:column;gap:1rem}
.eyebrow{margin:0;font-family:var(--mono);font-size:.72rem;letter-spacing:.14em;
  text-transform:uppercase;color:var(--accent)}
h1{margin:0;font-family:var(--serif);font-weight:600;line-height:1.15;
  font-size:clamp(1.8rem,1.2rem + 2.2vw,2.55rem);text-wrap:balance;
  letter-spacing:-.01em}
.lede{margin:0;max-width:62ch;color:var(--muted);font-size:1.02rem}
.lede strong{color:var(--ink);font-weight:600}

.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(9rem,1fr));
  gap:1px;margin:.5rem 0 0;padding:0;background:var(--rule);
  border:1px solid var(--rule);border-radius:3px;overflow:hidden}
.summary>div{background:var(--panel);padding:.85rem 1rem;display:flex;
  flex-direction:column;gap:.15rem}
.summary dt{font-family:var(--mono);font-size:.68rem;letter-spacing:.1em;
  text-transform:uppercase;color:var(--muted)}
.summary dd{margin:0;font-family:var(--mono);font-size:1.5rem;font-weight:500;
  font-variant-numeric:tabular-nums;line-height:1.1}
.summary dd small{font-size:.85rem;color:var(--quiet)}
.summary .is-ok dd{color:var(--ok)}
.summary .is-div dd{color:var(--diverge)}
.summary .is-fail dd{color:var(--fail)}
.stamp{margin:0;font-family:var(--mono);font-size:.75rem;color:var(--quiet)}

/* ---- legend ------------------------------------------------------------- */
.legend{border-top:1px solid var(--rule);border-bottom:1px solid var(--rule);
  padding:1.25rem 0;display:grid;gap:.7rem}
.legend h2{margin:0 0 .2rem;font-family:var(--sans);font-size:.72rem;
  letter-spacing:.12em;text-transform:uppercase;color:var(--muted);
  font-weight:600}
.legend p{margin:0;font-size:.9rem;color:var(--muted);max-width:70ch}
.legend b{font-family:var(--mono);font-size:.82rem;font-weight:500;
  color:var(--ink)}

/* ---- panels ------------------------------------------------------------- */
.panel{display:flex;flex-direction:column}
.panel h2{margin:0;font-family:var(--serif);font-weight:600;font-size:1.3rem;
  line-height:1.3;text-wrap:balance}
.panel__id{margin:.15rem 0 1.25rem;font-family:var(--mono);font-size:.75rem;
  color:var(--quiet)}

.row{display:grid;align-items:baseline;gap:.35rem 1.25rem;
  grid-template-columns:minmax(10rem,1.5fr) 7.5rem 7.5rem minmax(6rem,1fr) 8.5rem;
  padding:.85rem 0 .85rem .9rem;border-top:1px solid var(--rule-soft);
  border-left:2px solid var(--rule)}
.row--match{border-left-color:var(--accent)}
.row--divergent{border-left-color:var(--diverge)}
.row--pinned,.row--unidentified,.row--bounded{border-left-color:var(--rule)}
.row__q{font-weight:450;display:flex;flex-direction:column;gap:.1rem}
.row__param{font-family:var(--mono);font-size:.73rem;color:var(--quiet)}
.row__v{font-family:var(--mono);font-size:.95rem;font-variant-numeric:tabular-nums;
  display:flex;flex-direction:column;gap:.1rem;overflow-wrap:anywhere}
.row__cap{font-size:.63rem;letter-spacing:.1em;text-transform:uppercase;
  color:var(--quiet);font-family:var(--sans);font-weight:600}
.row__d{display:flex;align-items:center;gap:.55rem}
.row__k{font-size:.75rem;color:var(--muted);text-align:right}
.row--divergent .row__k{color:var(--diverge)}

.bar{position:relative;flex:1;min-width:2.5rem;height:4px;border-radius:2px;
  background:var(--rule)}
.bar--na{background:repeating-linear-gradient(90deg,var(--rule) 0 3px,
  transparent 3px 6px)}
.bar__fill{position:absolute;top:0;bottom:0;left:0;border-radius:2px;
  background:var(--accent)}
.bar--over .bar__fill{background:var(--fail)}
.bar__pct{font-family:var(--mono);font-size:.68rem;color:var(--quiet);
  white-space:nowrap;font-variant-numeric:tabular-nums}
.bar--over + .bar__pct{color:var(--fail)}

.note{margin:0;padding:.1rem 0 .9rem .9rem;border-left:2px solid var(--rule);
  font-size:.84rem;color:var(--muted);max-width:78ch}
.note span{font-family:var(--mono);font-size:.74rem;color:var(--quiet);
  display:block}

/* ---- footer ------------------------------------------------------------- */
.foot{border-top:1px solid var(--rule);padding-top:1.5rem;font-size:.88rem;
  color:var(--muted);display:grid;gap:.8rem;max-width:70ch}
.foot h2{margin:0;font-family:var(--serif);font-size:1.05rem;color:var(--ink);
  font-weight:600}
code{font-family:var(--mono);font-size:.85em;background:var(--accent-soft);
  color:var(--ink);padding:.1em .35em;border-radius:2px}

@media (max-width:760px){
  .row{grid-template-columns:1fr 1fr;gap:.6rem 1rem}
  .row__q{grid-column:1/-1}
  .row__d{grid-column:1/-1}
  .row__k{grid-column:1/-1;text-align:left}
}
@media (prefers-reduced-motion:reduce){*{animation:none!important;transition:none!important}}
</style>

<!--BODY-->
<div class="wrap">
<header class="head">
  <p class="eyebrow">Generated validation report</p>
  <h1>Does mixdra reproduce the published workbooks?</h1>
  <p class="lede">Every figure below was produced by re-fitting the reference
  datasets with the current version of the package. The <strong>workbook</strong>
  column comes from <code>inst/validation/oracles.csv</code>, which records each
  published value together with the sheet and column it was read from &mdash; and
  which is also what the automated test suite asserts against. A number cannot
  appear on this page without being tested, or be tested without appearing here.</p>
  <dl class="summary">
    <div class="is-ok"><dt>Workbook values reproduced</dt>
      <dd>{{PASS}}<small>&thinsp;/&thinsp;{{MATCH}}</small></dd></div>
    <div class="is-div"><dt>Deliberate differences</dt><dd>{{DIV}}</dd></div>
    <div class="is-fail"><dt>Failures</dt><dd>{{FAIL}}</dd></div>
  </dl>
  <p class="stamp">Generated {{GENERATED}}</p>
</header>

<section class="legend">
  <h2>What each row claims</h2>
  <p><b>reproduces workbook</b> &mdash; the engine must land on the published
  value within the stated tolerance. The bar shows how much of that tolerance the
  result actually used: a short bar is a strong agreement, a full bar only just
  qualifies.</p>
  <p><b>deliberate difference</b> &mdash; the engine departs from the workbook on
  purpose, and the note says why. These are not failures, but they are the rows
  worth reading rather than skimming. The workbook fits the dose&ndash;response
  curves jointly on the mixture data; this engine fits them from the
  single-compound data alone and then holds them fixed, so interaction shows up
  as interaction instead of being absorbed into the curves.</p>
  <p><b>fixed in both</b> &mdash; the parameter is pinned in the workbook and
  pinned here, so agreement is arithmetic rather than evidence.
  <b>not identified by the data</b> and <b>inequality only</b> mark quantities the
  data cannot pin to a point value; only a bound is claimed, and that bound is
  asserted in the test suite.</p>
</section>

{{PANELS}}

<footer class="foot">
  <h2>What this report does not show</h2>
  <p>Two implementations agreeing shows that they agree. It does not show that
  either is right. That is what the synthetic recovery tests are for: data is
  generated forward from known parameters through the same model the fitter uses,
  and the engine has to invert it back to those parameters. Those live in
  <code>tests/testthat/test-recovery.R</code> and run in the same CI job that
  produced this page.</p>
  <p>Reference data are <em>Folsomia candida</em> soil reproduction assays
  (chlorpyrifos, microplastics, imidacloprid and FBSA, in mg&thinsp;kg&#8315;&sup1;),
  fitted with the Jonker et&nbsp;al. (2005) mixture dose&ndash;response method.</p>
</footer>
</div>
)---"
