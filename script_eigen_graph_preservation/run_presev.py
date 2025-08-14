import pandas as pd
import argparse
import os

# compute preservation of two network
def preserv_matrix(network1, network2):
    common_modules = list(set(network1.index) & set(network2.index))
    preserv_matrix = pd.DataFrame(index=common_modules, columns=common_modules)
    
    for i, module1 in enumerate(common_modules):
        for j, module2 in enumerate(common_modules):
            #preserv_matrix.iloc[i, j] = 1 - max(abs(network1.loc[module1, module2] - network2.loc[module1, module2]),
                                              #abs(network1.loc[module1, module2] - network2.loc[module2, module1]))
            preserv_matrix.iloc[i, j] = 1 - (max(network1.loc[module1, module2], network2.loc[module1, module2]) - min(network1.loc[module1, module2], network2.loc[module1, module2]))
    
    return preserv_matrix

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compare eigensp networks')
    parser.add_argument('--n1', type=str, required=True)
    parser.add_argument('--n2', type=str, required=True)
    parser.add_argument('--odir', type=str, required=True, help='output dir')
    args = parser.parse_args()
    n1 = args.n1
    n2 = args.n2
    odir = args.odir
    
    network1 = pd.read_csv("{}.eigensp_cor.tsv".format(n1), sep='\t',index_col=0, header=0)
    network2 = pd.read_csv("{}.eigensp_cor.tsv".format(n2), sep='\t',index_col=0, header=0)
    preserv_matrix = preserv_matrix(network1, network2)
    preserv_matrix = preserv_matrix.astype(float)

    n2 = os.path.split(n2)[-1]
    n1 = os.path.split(n1)[-1]
    prefix = "{}.{}".format(n1,n2)
    preserv_matrix.to_csv(os.path.join(odir, "{}.preserv_matrix.tsv".format(prefix)),sep='\t')



