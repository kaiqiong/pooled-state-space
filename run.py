#!/usr/bin/env python3
"""Entry points for the manuscript replication package (Python standard library)."""
import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent

def call(command):
    print('+ '+' '.join(str(x) for x in command),flush=True)
    env = os.environ.copy()
    for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
        env.setdefault(key,'1')
    subprocess.run([str(x) for x in command],cwd=str(ROOT),env=env,check=True)

def rscript(script,*args):
    exe = shutil.which('Rscript')
    if exe is None:
        raise RuntimeError('Rscript is not on PATH. Install R before running this command.')
    call([exe,'--vanilla',ROOT/script]+list(args))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='command')
    check=sub.add_parser('check',help='Check runtime dependencies and package file integrity.')
    check.add_argument('--install',action='store_true',help='Install missing R packages from CRAN.')
    sub.add_parser('smoke',help='Run small synthetic numerical and bootstrap checks; no CGMacros download needed.')
    data=sub.add_parser('data',help='Download and verify CGMacros v1.0.0, or extract a local ZIP.')
    data.add_argument('--archive',type=Path)
    data.add_argument('--output',type=Path)
    app=sub.add_parser('application',help='Run a CGMacros analysis stage.')
    app.add_argument('stage',choices=['preprocess','prepare','primary','profile','reference','diagnostics','tables','all'])
    app.add_argument('--data',type=Path)
    app.add_argument('--out',type=Path)
    app.add_argument('--cores',type=int,default=1)
    sim=sub.add_parser('simulation',help='Run a simulation stage with published replication counts.')
    sim.add_argument('stage',choices=['help','smoke','independent','targets','clustered','summaries','robustness','numerical','diagnostics','caps','all'])
    sim.add_argument('--cores',type=int,default=1)
    fig=sub.add_parser('figures',help='Render outputs from shipped reference summaries or newly computed results.')
    fig.add_argument('--source',choices=['reference','recomputed'],default='reference')
    fig.add_argument('--include-application',action='store_true',help='Also render raw-window/application figures from local analysis outputs.')
    fig.add_argument('--python',type=Path,help='Python interpreter with numpy/matplotlib for the illustrated window.')
    fig.add_argument('--only',choices=['all','schematic','simulation','references','application','window','tables'],default='all')
    args=parser.parse_args()
    if args.command is None:
        parser.print_help();return
    if args.command=='check':
        rscript('scripts/check_environment.R',*(['--install'] if args.install else []))
        call([sys.executable,ROOT/'scripts/validate_release.py'])
    elif args.command=='smoke':
        for script in ['check_fixed_information.R','check_nuisance_inference.R','check_free_information.R']:
            rscript('tests/'+script)
    elif args.command=='data':
        command=[sys.executable,ROOT/'scripts/download_cgmacros.py']
        for key in ['archive','output']:
            value=getattr(args,key)
            if value is not None:command+=['--'+key,value.resolve()]
        call(command)
    elif args.command=='application':
        if args.cores<1:parser.error('--cores must be positive')
        os.environ['PSS_CORES']=str(args.cores)
        options=[]
        for key in ['data','out']:
            value=getattr(args,key)
            if value is not None:options+=['--'+key,str(value.resolve())]
        rscript('scripts/application/run.R',args.stage,*options)
    elif args.command=='simulation':
        if args.cores<1:parser.error('--cores must be positive')
        rscript('scripts/simulation/run.R',args.stage,str(args.cores))
    elif args.command=='figures':
        command=[sys.executable,ROOT/'scripts/figures/run.py','--source',args.source,'--only',args.only]
        if args.include_application:command.append('--include-application')
        if args.python:command+=['--python',args.python.resolve()]
        call(command)

if __name__=='__main__':
    try:
        main()
    except (RuntimeError,OSError,subprocess.CalledProcessError) as exc:
        sys.exit(str(exc))
