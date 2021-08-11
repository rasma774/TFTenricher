
import argparse

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'


def parse():
    DESC = """.

    Dependancies:
        - Python 3.9
        - NumPy
        - Scipy
        - Matplotlib
        - Pandas

    Author:
        - Rasmus Magnusson

    COPYRIGHT:
        - Rasmus Magnusson Linköping 2021

    LICENCE:
        - - GNU Affero General Public License v3.0

    Further reference:
        - https://github.com/rasma774/TFTenricher


    """
    # Create the parser instance
    parser = argparse.ArgumentParser(description=DESC)
    parser.add_argument('--tfs',
                        type=str,
                        nargs='*',
                        help='Names of the TFs that are to be mapped to target genes. Input should be a textfile of one row.')
    parser.add_argument('--sep',
                        default=' ',
                        type=str,
                        nargs=1,
                        help='Separator of the tf text file. Default is one blank space')
    parser.add_argument('--multiple_test_corr',
                        default='BenjaminiHochberg',
                        type=str,
                        nargs=1,
                        help='Multiple testing correction function. Either BenjaminiHochberg or Bonferroni. Default is BH. For user-defined multiple testing corrections, please use the TFTenricher toolbox.')

    parser.add_argument('--db',
                        type=str,
                        nargs=1,
                        default='GO',
                        help='What to compare the putative target genes to. Can be\n\
                        of types {"KEGG", "REACTOME", "GO", "GWAS"}')
    parser.add_argument('--results_savename',
                        type=str,
                        default=['results.csv'],
                        nargs=1,
                        help='The name of the file to which the enrichment results are saved')
    parser.add_argument('--ngenes',
                        default=None,
                        type=int,
                        nargs=1,
                        help='Number of top target genes to be analysed. Default is\n\
                            to use built-in monte carlo estimation')
    parser.add_argument('--FDR',
                        default=.05,
                        type=float,
                        nargs=1,
                        help='Set the FDR for enrichment analysis')
    parser.add_argument('--silent',
                        default=0,
                        type=int,
                        nargs=1,
                        help='Run TRenricher in silent mode')
    parser.add_argument('--plotname',
                        type=str,
                        default='-1',
                        nargs=1,
                        help='Set a name to save the plot')
    parser.add_argument('--plot_n_top',
                        type=int,
                        default=[15],
                        nargs=1,
                        help='The maximum number of top enrichment to plot')

    return parser.parse_args()


if __name__ == '__main__':
    parse()
