import os
import pandas as pd
from pathlib import Path


OUTPUT_CLUSTER = 'result.csv'
OUTPUT_CLUSTER_WTIH_GENES = 'result.with_genes.tsv'
OUTPUT_GENES = 'result_genes.csv'
OUTPUT_GENES_JSON = 'result_genes.json'
OUTPUT_RUNTIME = 'runtime.txt'
OUTPUT_SAMPLES = 'samples.txt'

"""Something went wrong with some of the result files, maybe it was due to the scikit learn update."""
def _fix_result(df_result_broken):
    assert all(df_result_broken['samples'].map(lambda x: '{' in x or ('(' in x and 's' in x and 'e' in x and 't' in x and ')' in x)))
    df_result_broken['samples'] = df_result_broken['samples'].map(lambda x: set(x.replace(',', '').replace('{', '').replace('}', '').replace("'", '').split(' ')))
    df_result_broken['samples'] = df_result_broken['samples'].map(lambda x: set() if len(x) and list(x)[0] == 'set()' else x)
    return df_result_broken

def write_runtime(output_path, runtime):
    with open(os.path.join(output_path, OUTPUT_RUNTIME), "w") as text_file:
        text_file.write(f"{runtime}")

def write_samples(output_path, samples):
    samples = ','.join(list(samples))
    with open(os.path.join(output_path, OUTPUT_SAMPLES), "w") as text_file:
        text_file.write(f"{samples}")
        
def write_output(output_path, df_data):
    # df_data['samples'] = df_data['samples'].map(lambda x: ','.join(str(x)))
    df_data.to_csv(os.path.join(output_path, OUTPUT_CLUSTER))

def save(df_result, runtime, output_path):
    write_output(output_path, df_result)
    write_runtime(output_path, runtime)
    print(f'Saved {output_path}.')

def create_or_get_result_folder(output):
    Path(output).mkdir(parents=True, exist_ok=True)
    return os.path.isfile(os.path.join(output, OUTPUT_CLUSTER))

def read_result(output_path):
    file_path = None
    if os.path.isfile(os.path.join(output_path, OUTPUT_CLUSTER_WTIH_GENES)):
        file_path = os.path.join(output_path, OUTPUT_CLUSTER_WTIH_GENES)
        df_data = pd.read_csv(file_path, index_col=0, sep='\t')
        df_data = df_data.reset_index(drop=True)
        df_data['genes'] = df_data['genes'].map(lambda x: set(x.split(' ')) if isinstance(x, str) else set()) # is sometimes float nan
        df_data['genes_up'] = df_data['genes_up'].map(lambda x: set(x.split(' ')) if isinstance(x, str) else set()) # is sometimes float nan
        df_data['genes_down'] = df_data['genes_down'].map(lambda x: set(x.split(' ')) if isinstance(x, str) else set()) # is sometimes float nan
        df_data['samples'] = df_data['samples'].map(lambda x: set(x.split(' ')))
        return df_data
    else:
        file_path = os.path.join(output_path, OUTPUT_CLUSTER)
        df_data = pd.read_csv(file_path, index_col=0)

    # fucked up saving data, some contain the string 'set()' instead of NaN or an actual set {...} 
    try:
       df_data = _fix_result(df_data)
    except Exception:
        df_data['samples'] = df_data['samples'].map(lambda x: set(x.split(',')) if not isinstance(x, float) else set())
        
    if OUTPUT_CLUSTER_WTIH_GENES in file_path:
        # genes in first result file
        return df_data
        
    try:
        df_data_genes = pd.read_csv(os.path.join(output_path, OUTPUT_GENES), index_col=0)
        # fucked up saving data, some contain the string 'set()' instead of NaN or an actual set {...} 
        try:
            df_data_genes = _fix_result(df_data_genes)
        except Exception:
            df_data_genes['samples'] = df_data_genes['samples'].map(lambda x: set(x.split(',')) if not isinstance(x, float) else set())
        df_data['genes'] = df_data_genes['samples']
    except:
        # methods do not return genes, will throw file not found exception
        pass
    
    return df_data

def read_runtime(output_path):
    with open(os.path.join(output_path, OUTPUT_RUNTIME), 'r') as f:
        return float(f.read())

def read_samples(output_path):
    with open(os.path.join(output_path, OUTPUT_SAMPLES), 'r') as f:
        samples = f.read().split(',')
        return set(samples)
    