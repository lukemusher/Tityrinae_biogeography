################################################################################
#
# bg_parse.py
# 
# This file is used to extract information from the posterior from the
# 		mnCharacterHistoryNewick monitor in RevBayes for style="events".
#
# authors: Michael Landis
#
################################################################################


import re

def get_events(fn='../output/bg_2rate.events.txt'):
    """ Extract events from mnCharacterHistoryNewick(...,style="events") file
   
        Returns a dictionary of character histories. The dictionary keys match
            branch indexes you get from the RevBayes function. In RevBayes, to
            get the branch index for the MRCA of e.g. "wolf" and "coyote"

                RevBayes > names <- my_tree.names()
                RevBayes > mrcaIndex(tree=my_tree, clade=clade("wolf", "coyote"))
                    16 
        
        Then, access the character history for the branch from Python using

            > my_histories = get_events(fn='my_events.txt')
            > my_histories[16]      # Warning, lots of output

        
        Keywords
            fn : the full filepath to your events file

    """
    f = open(fn,'r')
    lines = f.readlines()
    f.close()
    events = {}
    for l in lines:
        y = l.split('\t')
        x = re.findall('\\[([^]]+)\\]', y[-1])
        n = len(x)
        for i in range(1,n):
            s = x[i][1:-1]
            toks = s.split(';')
            taxon_index = int(toks[0].split('=')[1])
            if not events.has_key( taxon_index ):
                events[ taxon_index ] = {}
                events[ taxon_index ]['iteration'] = []
                events[ taxon_index ]['posterior'] = []
                events[ taxon_index ]['likelihood'] = []
                events[ taxon_index ]['prior'] = []
            events[taxon_index]['iteration'].append(int(y[0]))
            events[taxon_index]['posterior'].append(float(y[1]))
            events[taxon_index]['likelihood'].append(float(y[2]))
            events[taxon_index]['prior'].append(float(y[3]))
            for j in range(1,len(toks)):
                k,v = toks[j].split('=')
                if not events[ taxon_index ].has_key( k ):
                    events[ taxon_index ][ k ] = []
                if k == 'nd' or k == 'ch0' or k == 'ch1' or k == 'pa':
                    v = [ int(b) for b in v ]
                elif k == 'cs':
                    if v == 's':
                        v = 'subset_sympatry'
                    elif v == 'n':
                        v = 'narrow_sympatry'
                    elif v == 'w':
                        v = 'widespread_sympatry'
                    elif v == 'a':
                        v = 'allopatry'
                    v = v
                elif k == 'bn':
                    v = int(v)
                elif k == 'ev':
                    v = re.sub( '[{|}]', '', v[1:-1] ).split(',')
                    tmp = []
                    for vi in range(len(v)/4):
                        time=float(v[vi*4+0][2:])
                        age=float(v[vi*4+1][2:])
                        state=int(v[vi*4+2][2:])
                        idx=int(v[vi*4+3][2:])
                        tmp.append( {"time":time,"age":age,"state":state,"idx":idx} )
                    v = tmp

                events[ taxon_index ][ k ].append(v)
    return events

def get_best(d,n=5,f=None,p='posterior'):
    """ Get the best sampled histories from a branch history dictionary. 
            e.g. for my_data[16]

                > get_best(d=my_data[16])

        Returns the n best samples according the criterion p (posterior, prior, likelihood).
            Or, best proportion of samples by frequency f, if f is defined.

        Keywords
            d : The branch history dictionary for a single branch.
            n : The count of best samples to find.
            f : The proportion of the best samples to find. (Overrides n.)
            p : The criteria used to define "best", which can be
                    "posterior", "prior", or "likelihood"
    """
    if p not in ['posterior','prior','likelihood']:
        print('WARNING: p=\'' + p + '\' invalid, set p=\'posterior\'')
        p = 'posterior'
    K = d[p]
    if n < 0:
        print('ERROR: n < 0')
    if f is not None:
        if f > 1.:
            f = 1.
        if f < 0.:
            print('ERROR: f < 0.')
        n = f * len(K)
    best = sorted(range(len(K)), key=lambda x: K[x])[-n:][::-1]
    ret = {}
    for k in d.keys():
        ret[k] = []
    for i in best:
        for k in d.keys():
            ret[k].append(d[k][i])
    return(ret)

def get_gain_loss(d, freqs=True):
    """ Get the mean posterior number of area gain and loss events per area.

        Returns a 2d-vector, where the first index is what type of state was
            acquired (e.g. 0: loss, 1: gain) and the second index is what
            character (e.g. area) underwent that change.

        Keywords
            d     : The branch history dictionary for a single branch.
            freqs : If False, the numbers are raw counts
                        (i.e. not divided by the number of samples)
    """
    num_char = len(d['nd'][0])
    v = []
    for i in range(2):
        v.append([0.]*num_char)
    n = 0
    for ev in d['ev']:
        n += 1
        if len(ev) is 0:
            next
        else:
            for e in ev:
                v[e['state']][e['idx']] += 1.
    if freqs and n == 0:
        return([0.]*num_char)
    if freqs:
        for i in range(2):
            for j in range(num_char):
                v[i][j] = v[i][j] / n
    return(v)


def get_area_pair(d,k='nd',freqs=True):
    """ Get the posterior that two areas were occupied at a given node.

        Returns a 2d-vector, indexed by the areas' indices. The diagonal
            elements are the marginal probability of that area being
            in the node's range. The off-diagonal elements give the
            marginal probability the pair of areas we co-occupied by
            that node's range.
        
        Keywords
            d     : The branch history dictionary for a single branch.
            freqs : If False, the numbers are raw counts
                        (i.e. not divided by the number of samples)
    """
    num_char = len(d[k][0])
    v = []
    for i in range(num_char):
        v.append([0.]*num_char)
   
    if not d.has_key(k) or k not in ['nd','ch0','ch1']:
        print('ERROR: \'' + k + '\' invalid key')
        return

    n = 0
    for x in d[k]:
        n += 1
        for i in range(num_char):
            for j in range(i,num_char):
                if x[i] is 1 and x[j] is 1:
                    v[i][j] += 1.
                    v[j][i] = v[i][j]
    if freqs and n == 0:
        return([0.]*num_char)
    if freqs:
        for i in range(num_char):
            for j in range(num_char):
                v[i][j] = v[i][j] / n
    return(v)

def get_clado_state(d,minSize=1,freqs=True):
    """ Get the posterior for cladogenic state.

        Returns a dictionary for whose keys are the possible cladogenic states. States
            are defined as:
            'narrow_sympatry'       10000 -> 10000 | 10000 -> 'narrow_sympatry'
            'subset_sympatry'       11110 -> 10000 | 11110 -> 'subset_sympatry'
            'allopatry'             11110 -> 11000 | 00110 -> 'allopatry'
            'widespread_sympatry'   11110 -> 11110 | 11110 -> 'widespread_sympatry'

            Dictionary values give the frequency that cladogenic state for the branch
                was found in the posterior distribution. Use get_clado_prob()
                for more detailed information.
        
        Keywords
            d     : The branch history dictionary for a single branch.
            freqs : If False, the numbers are raw counts
                        (i.e. not divided by the number of samples)
    """
    if not d.has_key('ch0') or not d.has_key('ch1'):
        print('ERROR: no cladogenic state recorded')
        return

    num_char = len(d['nd'][0])
    v = {'subset_sympatry':0.,'narrow_sympatry':0.,'allopatry':0.,'widespread_sympatry':0.}
   
    n = 0
    for i,x in enumerate(d['cs']):

        if sum(d['nd'][i]) >= minSize:
            n += 1
            v[x] += 1

    if freqs:
        if n == 0:
            return([0.]*num_char)
        for k in v.keys():
            v[k] = v[k]/n

    return(v)


def get_clado_prob(d,freqs=True):
    """ Get the posterior for cladogenic configuration probs (incomplete)

        Returns a dictionary for whose keys are the possible cladogenic states. States
            are defined as:
            'narrow_sympatry'       10000 -> 10000|10000 -> 20000
            'subset_sympatry'       11110 -> 10000|21110 -> 21110
            'allopatry'             11110 -> 11000|00110 -> 11330
            'widespread_sympatry'   11110 -> 11110|11110 -> 22220

            1 means only one descendant lineage retained the area
            2 means both descendant lineages retained thea area
            3 means only the second descedant lineage retained the area

            Dictionary values give the frequency that each particular cladogenic state for the branch
                was found in the posterior distribution.

        Example:
            > get_clado_prob(dd[23]).keys()
                ['0110', '0010', '0100', '1110', '1100', '1010', '1000']
            > get_clado_prob(dd[23])['0110']
                {'0120': 0.1236,
                 '0130': 0.0200,
                 '0220': 0.2788}
        
        Keywords
            d     : The branch history dictionary for a single branch.
            freqs : If False, the numbers are raw counts
                        (i.e. not divided by the number of samples)
    """
    if not d.has_key('ch0') or not d.has_key('ch1'):
        print('ERROR: no cladogenic state recorded')
        return

    num_char = len(d['nd'][0])

    v = {}
    n = 0
    for k in range(len(d['nd'])):
        if d['cs'][k] == 'allopatry':
            n += 1
            nd_s = "".join([str(x) for x in d['nd'][k]])
            if not v.has_key(nd_s):
                v[nd_s] = {}
            y = d['ch0'][k]
            z = d['ch1'][k]
            ch0_s = "".join([str(x) for x in y])
            ch1_s = "".join([str(x) for x in z])
            y_tmp = y
            z_tmp = z
            if ch0_s < ch1_s:
                y_tmp = [ a*3 for a in y_tmp ]
            else:
                z_tmp = [ a*3 for a in z_tmp ]
            b = [ sum(a) for a in zip(y_tmp,z_tmp) ]
            b_s = "".join([str(x) for x in b])
            if not v[nd_s].has_key(b_s):
                v[nd_s][b_s] = 0.
            v[nd_s][b_s] += 1.
        else:
            n += 1
            nd_s = "".join([str(x) for x in d['nd'][k]])
            if not v.has_key(nd_s):
                v[nd_s] = {}
            y = d['ch0'][k]
            z = d['ch1'][k]
            ch0_s = "".join([str(x) for x in y])
            ch1_s = "".join([str(x) for x in z])
            b = [ sum(a) for a in zip(y,z) ]
            b_s = "".join([str(x) for x in b])
            if not v[nd_s].has_key(b_s):
                v[nd_s][b_s] = 0.
            v[nd_s][b_s] += 1.

    if freqs and n == 0:
        return([0.]*num_char)
    if freqs:
        for ki in v.keys() :
            for kj in v[ki].keys():
                v[ki][kj] = v[ki][kj] / n

    return(v)

def area_heatmap(x,title='',fn='out.pdf'):
    try:
        from matplotlib import pyplot as plt
        from pylab import colorbar
        import numpy as np
    except ImportError:
        'ERROR: cannot import matplotlib.pyplot'

    m = np.array(np.matrix(x))
    fig,ax = plt.subplots()
    hm = plt.pcolor(m,
            cmap=plt.get_cmap('Blues',20))
    plt.clim(0.,1.)
    cbar = plt.colorbar(hm)
    fig.suptitle(title,
            fontsize=14,fontweight='bold')
    n_areas = m.shape[0]

    tick_txt = [ str(s) if (s)%5==0 else '' for s in range(n_areas) ]
    ax.set_xticks(np.arange(n_areas)+0.5,minor=False)
    ax.set_yticks(np.arange(n_areas)+0.5,minor=False)
    ax.set_xticklabels(tick_txt)
    ax.set_yticklabels(tick_txt)

    plt.show()
    #plt.savefig(fn)
    return

