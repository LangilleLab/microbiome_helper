# Create a simple dataframe 
  
# importing pandas as pd 
import pandas as pd 
import argparse, sys, textwrap

parser=argparse.ArgumentParser()
parser.add_argument('--filename', help = 'picrust2 output')
parser.add_argument('--outfilename', help = 'Output filename to store stratified output')

def main():
    # creating a dataframe
    args = parser.parse_args()
    fileN = args.filename
    outputfileN = args.outfilename
    
    #df = pd.read_table(file, sep='\t', header=0, index_col=0, lineterminator='\n')
    df = pd.read_csv(fileN, sep = '\t')
    print (df.head())
    #df[['function','sequence']] = df.function.str.split("|",expand=True)
    #df = df[ ['taxon'] + [ col for col in df.columns if col != 'taxon' ] ]
    df['function'] = df['pathway'].str.cat(df['sequence'], sep ="|")
    df.drop(['pathway','sequence'], axis = 1, inplace = True)
    cols = df.columns.tolist()
    cols.insert(0, cols.pop(cols.index('function')))
    df = df.reindex(columns= cols)
    print (df.head())
    
    pd.DataFrame.to_csv(df, path_or_buf=outputfileN, sep='\t', na_rep='', header=True, index=False, mode='w', line_terminator='\n', escapechar=None, decimal='.')
    
if __name__ == "__main__":
        main();
