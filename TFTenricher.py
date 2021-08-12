from src import enrich_utils
from src import map2trgt_utils
from src import stat_utils
from src import plot_utils
from src import parse_utils
from src import citation_handler
from src.build_corrmat import check_corrmat

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2021 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'
__LICENSE__ = 'GNU Affero General Public License v3.0'
__version__ = '0.02'
__cite__ = ''


# TODO: put the TFTenricher ref into citation when possible

class TFTenricher:
    def __init__(self,
                 TFs,
                 mapmethod='corr',
                 silent=False,
                 top_n_genes=None):
        """
        The transcription factor downstream annotation enricher (TFTenricher)
        package is a bioinformatics tool to enable users to do an enrichment
        analysis on lists of transcription factors (TFs).


        Parameters
        ----------
        TFs : list
            List of transcription factors (TFs) in SYMBOL annotations.

        mapmethod : str or function, optional
            The function to map TFs to target genes. If an independent function
            is given, the function should take a list of TFs, and return a list
            of target genes. The default is 'corr', a built-in multiple testing
            corrected correlation cut-off.

        silent : bool, optional
            Specify whether to print the output. The default is False.

        top_n_genes : int or None, optional
            In the case of an execive number of output genes from 'mapmethod',
            TFTenricher can be set to only consider the top N genes in
            top_n_genes. The default is None, equaling to include all genes
            returned from the mapmethod.


        Attributes
        -------
        target_genes : The estimated target genes.


        How to cite
        -------
        <put cite here>
        Moreover, additional citation from the third-party materials can be
        accessed using the TFTenricher.cite method.

        Thank you for using TFTenricher!
        """

        assert len(TFs) > 0

        self.TFs = TFs.copy()
        self.silent = silent

        # For the references
        self.used_methods = []

        # Attributes that will be filled in other methods
        self.enrichments = None
        self.multtest_fun = None


        if mapmethod == 'corr':
            # If corr matrix does not exist, build it
            check_corrmat()

            self.mapmethod = map2trgt_utils.correlation_genes
            self.used_methods.append('corrs')
        else:
            self.mapmethod = mapmethod

        self.target_genes = self.mapmethod(TFs,
                                           silent=silent,
                                           top_n_genes=top_n_genes
                                           )


    def downstream_enrich(self,
                          db='GO',
                          FDR=0.05,
                          multiple_testing_correction='BenjaminiHochberg',
                          ):
        """


        Parameters
        ----------
        db : str or dict, optional
            If a string, 'db' should be in ['GO', REACTOME', 'KEGG', GWAS']
            Else, 'db' can be a dict with annotations as keys, and associated
            genes as elements, e.g. {'annot' : ['gene1', 'gene2'}. Here, the
            gene names are in the SYMBOL id. Note that, in the case of 'db'
            being a dict, TFTenricher.downstream_enrich by default removes
            annotations with fewer target genes than 10.
            The default is 'GO'.

        FDR : float, optional
            The accepted false discovery rate of associated annotations. The
            default is 0.05.

        multiple_testing_correction : str or function, optional
            The name of the multiple testing correction method. The default is
            'BenjaminiHochberg'. For examples of how to pass own funciton here,
            please see the github page or
           stat_utils.benjaminihochberg_correction, which is also the default

        Attributes
        -------
        enrichments : pandas DataFrame of the results, with annotation, P-value
        and a bool value of whether the test passed a multiple testing
        correction

        """

        self.db = db
        if multiple_testing_correction == 'BenjaminiHochberg': # use this unless otherwise told
            self.multtest_fun = stat_utils.benjaminihochberg_correction
        elif multiple_testing_correction == 'Bonferroni':
            self.multtest_fun = stat_utils.bonferroni_correction
        else:
            # Here, we let the user define the testing function
            self.multtest_fun = multiple_testing_correction

        if isinstance(db, str):
            db = db.upper()
            self.used_methods.append(db)

        res = enrich_utils.set_enrichments(self.target_genes,
                                           mult_test_corr=self.multtest_fun,
                                           db=db,
                                           FDR=FDR,
                                           )
        self.enrichments = res

    def plot(self,
             savename=None,
             plot_Ntop='all',
             textlength=None,
             sorton='OR',
             remove_non_FDR=True,
             padding=0.7,
             tick_font_size=14):
        """


        Parameters
        ----------
        savename : str, optional
            IF specified, save the fig to this path. The default is None.

        plot_Ntop : str or int, optional
            If int, plot the N top associations. If 'all' plot all
            associations that passed mult. testing correction. The default is
            'all'.

        textlength : int, optional
            At what position to wrap the annotations headings in the figure.
            TFTenricher.plot searches for the most suitable place closest to the
            textlengh int and inserts a new line. The default is None.

        sorton : {'OR', 'p'}, optional
            Which variable of the results DataFrame to sort on. The default is
            'OR'.
        remove_non_FDR : bool, optional
            If set to False, also annotations that do not pass an FDR can be
            plotted. The default is True.

        padding : float, optional
            Adds some padding between the bars in the plot. The default is 0.7.

        tick_font_size : int, optional
            The size of the plot ticks. The default is 14.

        Returns
        -------
        f : matplotlib figure.

        ax : Axes of a matplotlib figure.

        """
        if plot_Ntop == 'all':
            plot_Ntop = self.enrichments.shape[0]

        f, ax = plot_utils.plot_res(self.enrichments.copy(),
                                    savename,
                                    plot_Ntop,
                                    textlength=textlength,
                                    sorton=sorton,
                                    remove_non_FDR=remove_non_FDR,
                                    padding=padding,
                                    tick_font_size=tick_font_size,
                                    )
        return f, ax

    def cite(self):
        """
        Write the references to the databases used, and the TFTenricher
        reference.
        """
        citation_handler.write_citations(self.used_methods)



if __name__ == '__main__':
    args = parse_utils.parse()

    # Unpack the TF names
    with open(args.tfs[0], 'r') as f:
        TFs = f.read().strip('\n').split(args.sep[0])

    # Map TFs to targets
    enr = TFTenricher(TFs, silent=args.silent, top_n_genes=args.ngenes)

    # Calculate the overlaps between putative downstream genes and gene sets
    enr.downstream_enrich(
        db=args.db,
        FDR=args.FDR,
        multiple_testing_correction=args.multiple_test_corr
        )

    if args.plotname != '-1':
        enr.plot(savename=args.plotname[0], plot_Ntop=args.plot_n_top[0])

    enr.enrichments.to_csv(args.results_savename[0])

