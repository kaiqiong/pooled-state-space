#!/usr/bin/env python3
"""Render manuscript displays, explicitly choosing reference summaries or recomputed outputs."""
import argparse
import csv
import hashlib
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--mode', '--source', dest='mode', choices=['reference', 'recomputed'], default='reference')
p.add_argument('--only', choices=['all', 'schematic', 'simulation', 'references', 'application', 'window', 'tables'], default='all')
p.add_argument('--results-dir', default=str(ROOT / 'results'), help='Root of recomputed simulation/application results and rendered outputs')
p.add_argument('--include-application', action='store_true', help='Also render Figure 3 and S2 from local application results; no participant traces are shipped')
p.add_argument('--rscript', default='Rscript')
p.add_argument('--python', default=sys.executable, help='Python interpreter with numpy and matplotlib for the local-data window illustration')
a = p.parse_args()
results = Path(a.results_dir).resolve()
figures = results / 'figures'
tables = results / 'tables'
figures.mkdir(parents=True, exist_ok=True)
tables.mkdir(parents=True, exist_ok=True)

def run_r(script, *args):
    subprocess.check_call([a.rscript, str(ROOT/'scripts/figures'/script)] + [str(x) for x in args], cwd=str(ROOT))

def needed(part):
    return a.only in ['all', part]

if needed('schematic'):
    run_r('render_schematic.R', figures)

if a.mode == 'reference':
    sim = ROOT/'reference/figures/simulation'
    boundary = ROOT/'reference/tables/boundary_coverage.csv'
    supplement = ROOT/'reference/simulation/tables'
    app = ROOT/'reference/figures/application'
else:
    sim = figures/'recomputed_summary_inputs'
    sim.mkdir(parents=True, exist_ok=True)
    source = results/'simulation'
    if needed('simulation') or needed('tables'):
        sources = {
            'benchmark.csv': source/'tables/part_g1_M60_beta_benchmark.csv',
            'information_ablation.csv': source/'free_variance_ablation/results/summary.csv',
            'clustered_free_variance.csv': source/'matched/paired_posvar200_free_variance_comparison.csv',
            'episode_size.csv': source/'tables/part_j_NT_sweep_summary.csv',
        }
        for name, src in sources.items():
            if not src.exists():
                raise SystemExit('Missing recomputed summary: {}. Run the corresponding simulation job first; see docs/OUTPUT_MAP.csv.'.format(src))
            shutil.copyfile(str(src), str(sim/name))
    boundary = source/'tables/part_m_boundary_v2_working_summary.csv'
    supplement = figures/'recomputed_supplement_inputs'
    supplement.mkdir(parents=True, exist_ok=True)
    if needed('tables'):
        summary_names = ['part_a_se_analytic_vs_numeric.csv', 'part_e1_mv_se_analytic_vs_numeric.csv',
          'part_h1_interaction_se_analytic_vs_numeric.csv', 'part_e2_mv_coverage_summary.csv',
          'part_m_boundary_v2_working_summary.csv', 'part_p_mixed_comparator_v2_summary.csv',
          'part_n_final_summary.csv', 'part_l_heavytail_summary.csv', 'part_o_direct_optimization_v2.csv',
          'cap_sensitivity_application_ms.csv', 'cap_sensitivity_simulation_ms.csv']
        for name in summary_names:
            shutil.copyfile(str(source/'tables'/name), str(supplement/name))
        for name in ['paired_posvar200_free_variance_comparison.csv', 'free_fixed_matched200_comparison.csv']:
            shutil.copyfile(str(source/'matched'/name), str(supplement/name))
        fixed_check = source/'tables/fixed_variance_se_analytic_vs_numeric.csv'
        if fixed_check.exists():
            shutil.copyfile(str(fixed_check), str(supplement/fixed_check.name))
    app = results/'application' 

if needed('simulation'):
    run_r('render_simulations.R', sim, figures)
if needed('references'):
    run_r('render_references.R', app, figures)
if needed('tables'):
    run_r('render_table_1.R', sim, boundary, tables)
    run_r('render_supplement_tables.R', supplement, tables, a.mode)
    app_tables = ROOT/'reference/application' if a.mode == 'reference' else results/'application'
    for stem in ['table_S11_cohort_flow', 'table_S12_variance_profile', 'table_S13_primary_regression', 'table_S14_reference_sensitivity']:
        src = app_tables/(stem+'.csv')
        if not src.exists():
            raise SystemExit('Missing application table summary: {}. Run scripts/application/run.R tables first.'.format(src))
        shutil.copyfile(str(src), str(tables/src.name))
        with src.open(newline='') as f:
            rows = list(csv.reader(f))
        def tex_escape(value):
            return value.replace('_', r'\_').replace('%', r'\%').replace('&', r'\&')
        lines = [r'\begin{tabular}{'+'l'*len(rows[0])+'}', r'\toprule']
        for idx, row in enumerate(rows):
            lines.append(' & '.join(tex_escape(x) for x in row)+r' \\')
            if idx == 0:
                lines.append(r'\midrule')
        lines.extend([r'\bottomrule', r'\end{tabular}'])
        (tables/(stem+'.tex')).write_text('\n'.join(lines)+'\n')

if a.mode == 'recomputed' or a.include_application:
    local_app = results/'application'
    if needed('application'):
        run_r('render_application.R', local_app, figures)
    if needed('window'):
        subprocess.check_call([a.python, str(ROOT/'scripts/figures/render_window.py'), str(local_app), str(figures)], cwd=str(ROOT))
elif a.only in ['application', 'window']:
    raise SystemExit('Figures 3A and S2 require local CGMacros episode data; use --source recomputed after the application pipeline.')
elif a.only == 'all':
    print('Reference mode does not create Figure 3 or S2: individual glucose traces are not redistributed. Run the application pipeline and --source recomputed for those figures.')

print('Rendered outputs from {} inputs; no model fitting performed.'.format(a.mode))
