# This script computes correaltion between clusters and each feature in metadata.
import pandas as pd
from scipy.stats import pearsonr
import argparse
import os

# check data and count NA
def column_stats(df):

    results = []
    
    for col in df.columns:
        na_count = df[col].isna().sum()
        dtype = df[col].dtype
        if dtype in ['int64', 'float64']:
            unique_values = df[col].unique()
            if len(unique_values) <= 10:
                data_type = 'Categorical'
            else:
                data_type = 'Numerical'
        else:
            data_type = 'Categorical'

        results.append({
            'Column': col,
            'Data Type': data_type,
            'NA Count': na_count
        })
    
    return pd.DataFrame(results)

def process_metadata_files(df, outpath):
    stats = column_stats(df)
    stats = stats[stats['NA Count'] != len(df)]
    stats.to_csv(outpath, sep='\t', index=False)
        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare eigensp networks')
    parser.add_argument('--prefix', type=str, required=True, help='Output file prefix')
    parser.add_argument('--idir', type=str, required=True, help='input dir')
    parser.add_argument('--odir', type=str, required=True, help='output dir')
    parser.add_argument('--pheno', type=str, required=True, help='Phenotype')
    args = parser.parse_args()
    prefix = args.prefix
    pheno = args.pheno
    idir = args.idir
    odir = args.odir
    # expr_df = pd.read_csv(f'data_large_scale/{pheno}/{prefix}/abd.tsv', sep='\t', index_col=0)
    # meta_df = pd.read_csv(f'data_large_scale/{pheno}/{prefix}/metadata.tsv', sep='\t')
    # read file
    eigensp_df = pd.read_csv(os.path.join(odir, pheno, '{}.together.eigensp.csv'.format(prefix)),index_col=0,header=0,sep="\t")
    metadata_df = pd.read_csv(os.path.join(idir, pheno, prefix,'metadata.tsv'), sep='\t', header=0)
    metadata_stat_path = os.path.join(odir,pheno,'{}.metadata.stat.tsv'.format(prefix))
    process_metadata_files(metadata_df, metadata_stat_path)
    metadata_stat_df = pd.read_csv(metadata_stat_path, sep='\t',header=0)


# select datatype == numerical columns from metadata
selected_col_list = metadata_stat_df[metadata_stat_df['Data Type'] == 'Numerical']['Column'].tolist()


# select sample id and other numerical columns from metadata
selected_metadata_df = metadata_df[['sample_id'] + selected_col_list]

# compute module and features' pearson and p_value for each column
results = []
for module in eigensp_df['module'].unique():
    module_df = eigensp_df[eigensp_df['module'] == module]
    
    # merge eigenspand metadata
    merged_df = module_df.merge(selected_metadata_df, left_on='sample', right_on='sample_id', how='inner')
    
    for col in selected_col_list:
        valid_rows = merged_df[[col, 'eigensp']].dropna()

        corr, p_value = pearsonr(valid_rows[col], valid_rows['eigensp'])

        results.append({
            'module': module,
            'column': col,
            'project': prefix,
            'correlation': corr,
            'p_value': p_value
        })

result_df = pd.DataFrame(results)
result_df.to_csv(os.path.join(odir, pheno, '{}.module_correlation.tsv'.format(prefix)), index=False,sep="\t")
