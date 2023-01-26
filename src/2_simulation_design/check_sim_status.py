import glob
from termcolor import colored

outdir = '/lustre/scratch117/cellgen/team205/ed6/PBMC_CZI_integration_filtered/tmp/'
diff_methods = ['meld', 'milo', 'cna']
emb_methods = ['scArches', 'scVI']
n_cells = ['500', '1000']

def _check_output(diff_method, emb_method, n='500'):
    print(f'Subsample {n} cells | {emb_method} embedding | {diff_method}')
    for d in ['ACR', 'AR', 'CR']:
        match_dirs = glob.glob(f'{outdir}/*{n}cells*/')
        match_files = glob.glob(f'{outdir}/*{n}cells*/{d}_design.{emb_method}_{diff_method}.h5ad')
        if len(match_files) == len(match_dirs):
            print(colored(f'{d} design | {len(match_files)}/{len(match_dirs)} completed', 'green'))
        elif len(match_files) == 0:
            print(colored(f'{d} design | {len(match_files)}/{len(match_dirs)} completed', 'red'))
        else:
            print(colored(f'{d} design | {len(match_files)}/{len(match_dirs)} completed', 'yellow'))


n=500
for diff_method in diff_methods:
    _check_output(diff_method, emb_method='scArches', n = n)
for m in emb_methods:
    _check_output(diff_method = 'milo', emb_method=m, n = n)

n=1000
# for diff_method in diff_methods:
#     _check_output(diff_method, emb_method='scArches', n)
for m in emb_methods:
    _check_output(diff_method = 'milo', emb_method=m, n = n)
