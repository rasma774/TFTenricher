import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

__author__ = 'Rasmus Magnusson'
__COPYRIGHT__ = 'Copyright (C) 2020 Rasmus Magnusson'
__contact__ = 'rasma774@gmail.com'

_golden = (1 + (5**(.5)))/2

# TODO: clean this horrible code up...

def _split_lines(lines, textlength):
    # Dont think we should break lines more than 1 time (?)
    for i, lbl in enumerate(lines):
        if len(lbl) > textlength:
            # get where whitespace
            lbl = list(lbl)
            ws = np.where(np.array(lbl) == ' ')[0]
            if len(ws) < 2:
                continue
            if not ws[-2] <= textlength:
                # First place that breaks textlength
                place = ws[np.where(ws > textlength)][0]
                lbl = list(lbl)
                lbl[place] = '\n'
                lbl = ''.join(lbl)
                lines[i] = lbl
    return lines

def _get_template():
    f, ax = plt.subplots(
        1,
        1,
        figsize=(4, 4*_golden))
        
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    
    ax.spines['bottom'].set_position(('outward', 20))
        

    ax.tick_params(axis='both', 
                   width=2, 
                   length=4, 
                   labelsize=14)
    return f, ax
    
def plot_res(enrichments, 
               savename=None, 
               plot_Ntop=25, 
               color='#850000ff',
               textlength=30,
               cmap=mpl.cm.OrRd,
               sorton='OR',
               remove_non_FDR=True,
               padding=0.7,
               tick_font_size=14,
               ):
    
    if remove_non_FDR:
        enrichments = enrichments[enrichments.FDR]
        
    # Handle input variable                             
    if plot_Ntop is None:
        # TODO: what it there are too many here?
        nplot = enrichments.FDR.sum()
    else:
        nplot = np.min((plot_Ntop, enrichments.shape[0]))
        
    # Clean input data
    enrichments = enrichments.iloc[:, :2]
    enrichments.p = -np.log10(enrichments.p)
    enrichments = enrichments.sort_values(sorton).iloc[::-1, :]
    enrichments = enrichments.iloc[:nplot, :]
    

    
    # For the colorbar
    norm = mpl.colors.Normalize(vmin=np.min((enrichments.p.min(), -np.log10(0.05))), vmax=enrichments.p.max())

    if textlength is not None:
        enrichments.index = _split_lines(enrichments.index.values, textlength)
        
    f, ax = _get_template()
    ax.tick_params(axis='both', labelsize=tick_font_size)
    
    
    for i in range(nplot):
        ax.barh(padding*(nplot - i), width=enrichments.OR[i], height=0.4, color=cmap(norm(enrichments.p[i])))    
    
    ax.set_xlim([0, enrichments.OR[:i].max()*1.2])
    ax.set_yticks(padding*np.array(range(1, 1 + nplot)))
    ax.set_yticklabels(enrichments.index[:plot_Ntop][::-1])
    ax.set_xlabel('Odds Ratio', fontsize=17)
    colorbar = f.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax, shrink=0.4, pad = 0.15)
    colorbar.set_label(r'-log$_{10}$ P', fontsize=15, labelpad=-50)
    
    if  savename is not None:
        f.savefig(savename, bbox_inches='tight')
        
    return f, ax