# try:
#     import SSRCalc
# except ImportError:
#     SSRCalc = None
# try:
#     from pyteomics import biolccc
# except ImportError:
#     pass
import sys
import re
import string
from pyteomics import parser, mass, auxiliary as aux, achrom, pepxml
import numpy as np
from scipy.stats import scoreatpercentile
from copy import copy
from scipy.spatial import cKDTree
from collections import Counter, defaultdict
from operator import itemgetter
from itertools import izip, izip_longest
import logging
logger = logging.getLogger(__name__)

try:
    from configparser import RawConfigParser
except ImportError:
    from ConfigParser import RawConfigParser

infiles_dict = dict()


def get_modif_name(modif, aa=''):
    return str(round(modif['mass'], 4)) + aa

class CustomRawConfigParser(RawConfigParser):
    def get(self, section, option):
        val = RawConfigParser.get(self, section, option)
        if isinstance(val, basestring):
            return val.split('|')[0]
        return val

    def get_choices(self, section, option):
        val = RawConfigParser.get(self, section, option)
        if isinstance(val, basestring) and len(val.split('|')) > 1:
            return val.split('|')[1]
        else:
            return ''

def get_dbname(prot):
    return prot['protein']

_modchars = set(string.ascii_lowercase + string.digits)
def custom_split_label(mod):
    j = 0
    while mod[j] in _modchars:
        j += 1
    if j == 0:
        return mod[1:], '-', ']'
    if len(mod[j:]) > 1 and '[' in mod:
        return mod[:j], mod[j:].replace('[', ''), '['
    elif len(mod[j:]) > 1 and ']' in mod:
        return mod[:j], mod[j:].replace(']', ''), ']'
    elif len(mod[j:]) == 1:
        if mod[0] == '-' or ']' in mod:
            return mod[:j], '-', ']'
        elif mod[-1] == '-' or '[' in mod:
            return mod[:j], '-', '['
        else:
            return mod[:j], mod[j:], ''

def get_aa_mass(settings):
    aa_mass = mass.std_aa_mass.copy()
    aa_mass['-'] = 0.0
    fmods = settings.get('modifications', 'fixed')
    # if fmods:
    #     for mod in re.split(r'[,;]\s*', fmods):
    #         m, aa = parser._split_label(mod)
    #         aa_mass[aa] += settings.getfloat('modifications', m)
    #         aa_mass[mod] = aa_mass[aa] + settings.getfloat('modifications', m)
    vmods = settings.get('modifications', 'variable')
    if fmods:
        if vmods:
            vmods = ','.join([vmods, fmods])
        else:
            vmods = fmods
    if vmods:
        mods = [custom_split_label(mod) for mod in re.split(r',\s*', vmods)]#[(l[:-1], l[-1]) for l in re.split(r',\s*', vmods)]
        if settings.getboolean('advanced options', 'snp'):
            for k, v in mass.std_aa_mass.items():
                for kk in mass.std_aa_mass.keys():
                    aa_mass['snp' + kk.lower()] = v
                # mods.append(('snp' + k.lower(), round(v, 4), 'snp'))
        # for (mod, aa, term), char in zip(mods, string.punctuation):
        for (mod, aa, term), char in zip(mods, string.punctuation):
            if term == '[' and aa == '-':
                aa_mass[mod + '-'] = settings.getfloat('modifications', mod.replace('-', '')) + settings.getfloat('modifications', 'protein nterm cleavage')
            elif term == ']' and aa == '-':
                aa_mass['-' + mod] = settings.getfloat('modifications', mod.replace('-', '')) + settings.getfloat('modifications', 'protein cterm cleavage')
            else:
                aa_mass[mod + aa] = aa_mass[aa] + settings.getfloat('modifications', mod)
    return aa_mass

def filter_evalue_prots(prots, FDR=1.0, remove_decoy=True, dec_prefix='DECOY_'):

    proteins = prots.items()

    isdecoy = lambda x: x[0].startswith(dec_prefix)
    escore = lambda x: float(x[1]['expect'])
    filtered_proteins = aux.filter(proteins, fdr=float(FDR) / 100, key=escore, is_decoy=isdecoy,
                                   remove_decoy=False, formula=1, full_output=True)
    qvals_e = aux.qvalues(filtered_proteins, key=escore, is_decoy=isdecoy, reverse=False, remove_decoy=False, formula=1,
                          full_output=True)
    new_prots = {}
    for val in qvals_e:
        val[-1][1]['qval'] = val[-2]
        if (not remove_decoy or not val[-1][0].startswith(dec_prefix)):
            new_prots[val[-1][0]] = val[-1][1]
    logger.info('Actual protein-level FDR = %.2f%%', aux.fdr(filtered_proteins, is_decoy=isdecoy) * 100)
    return new_prots

def get_settings(fname=None, default_name='default.cfg'):
    """Read a configuration file and return a :py:class:`RawConfigParser` object.
    """
    kw = {'inline_comment_prefixes': ('#', ';')
            } if sys.version_info.major == 3 else {}
    kw['dict_type'] = dict
    kw['allow_no_value'] = True

    raw_config = CustomRawConfigParser(**kw)
    if default_name:
        raw_config.read(default_name)
    if fname:
        raw_config.read(fname)
    return raw_config


def FDbinSize(X):
    """Calculates the Freedman-Diaconis bin size for
    a data set for use in making a histogram
    Arguments:
    X:  1D Data set
    Returns:
    h:  F-D bin size
    """
    X = np.sort(X)
    upperQuartile = scoreatpercentile(X, 75)
    lowerQuartile = scoreatpercentile(X, 25)
    IQR = upperQuartile - lowerQuartile
    h = 2. * IQR / len(X) ** (1. / 3.)
    return h

class Descriptor():
    def __init__(self, name, single_formula=False, massive_formula=False, group='A', binsize='auto'):
        self.name = name
        self.single_formula = single_formula
        self.massive_formula = massive_formula
        self.group = group
        self.binsize = binsize
        self.normK = None
        self.array = None

    def formula(self, peptides):
        if self.single_formula:
            return [self.single_formula(peptide) for peptide in peptides.peptideslist]
        else:
            return self.massive_formula(peptides)

    def get_array(self, peptides):
        # tmp = np.array([self.formula(peptide) for peptide in peptides.peptideslist])
        tmp = np.array(self.formula(peptides))
        tmp.sort()
        self.array = tmp
        return tmp

    def get_binsize(self, peptides=None):
        if self.binsize == 'auto':
            if self.group == 'B':
                return 0.25
            else:
                return FDbinSize(self.get_array(peptides))
        else:
            return self.binsize


class PeptideList:
    def __init__(self, settings=None, mods=None):
        self.listing = ['peptideslist', 'spectrumlist', 'RT_exp', 'RT_predicted']
        self.listing_nparrays = ['spectrumlist', 'RT_exp', 'RT_predicted']
        self.peptideslist = []
        self.spectrumlist = []
        self.RT_exp = []
        self.RT_predicted = []
        self.calibrate_coeff = None
        self.RC = False
        self.infiles = set()
        self.settings = settings
        if not mods:
            self.modification_list = {}
        else:
            self.modification_list = mods
        self.total_number_of_PSMs = 0
        self.total_number_of_PSMs_decoy = 0
        self.total_number_of_peptides_in_searchspace = 0
        self.total_number_of_proteins_in_searchspace = 0
        self.total_number_of_spectra = 0
        self.nterm_mass = self.settings.getfloat('modifications', 'protein nterm cleavage')
        self.cterm_mass = self.settings.getfloat('modifications', 'protein cterm cleavage')
        self.aa_list = get_aa_mass(settings)
        self.proteins_dict = defaultdict(list)

        fmods = self.settings.get('modifications', 'fixed')
        # if fmods:
        #     for mod in re.split(r'[,;]\s*', fmods):
        #         m, aa = parser._split_label(mod)
        #         self.modification_list[str(round(mass.std_aa_mass[aa] + settings.getfloat('modifications', m), 4)) + aa] = m
        
        vmods = settings.get('modifications', 'variable')
        if fmods:
            if vmods:
                vmods = ','.join([vmods, fmods])
            else:
                vmods = fmods
        if vmods:
            mods = [custom_split_label(mod) for mod in re.split(r',\s*', vmods)]#[(l[:-1], l[-1]) for l in re.split(r',\s*', vmods)]
            if settings.getboolean('advanced options', 'snp'):
                for k, v in mass.std_aa_mass.items():
                    mods.append(('snp' + k.lower(), round(v, 4), 'snp'))
            for (mod, aa, term), char in zip(mods, string.punctuation):
                if term == 'snp':
                    self.modification_list[aa] = mod
                elif term == '[' and aa == '-':
                    self.modification_list[str(round(self.nterm_mass + settings.getfloat('modifications', mod), 4))+'n'] = mod + '-'
                elif term == ']' and aa == '-':
                    self.modification_list[str(round(self.cterm_mass + settings.getfloat('modifications', mod), 4))+'c'] = '-' + mod
                else:
                    self.modification_list[str(round(self.aa_list[aa] + settings.getfloat('modifications', mod), 4))+aa] = mod

    def __len__(self):
        return len(self.peptideslist)

    def get_number_of_peptides(self):
        return len(set(p.sequence for p in self.peptideslist))

    def get_infiles(self):
        if not self.infiles:
            self.infiles.update(p.infile for p in self.peptideslist)
        return self.infiles

    def is_decoy_array(self):
        return np.array([peptide.note2 == 'wr' for peptide in self.peptideslist])


    def rem_elements(self, js):
        js.sort(reverse=True)
        js_set = set(js)
        js_stay = [x for x in range(self.spectrumlist.size) if x not in js_set]
        for l in self.listing:
            tmp = getattr(self, l)
            if len(tmp):
                if isinstance(tmp, list):
                    for j in js:
                        tmp.pop(j)
                    setattr(self, l, tmp)
                elif isinstance(tmp, np.ndarray):
                    tmp = tmp[js_stay]
                    setattr(self, l, tmp)
                else:
                    logger.critical('Unknown type of PeptidesList attribute. Smth wrong in the Code!')

    def get_number_of_spectra(self):
        """Returns the number of MS/MS spectra used for the search. If mgf file is not available,
         returns number of identified PSMs as approximation.
        """
        if self.total_number_of_spectra:
            return self.total_number_of_spectra
        else:
            return self.total_number_of_PSMs

    def get_from_pepxmlfile(self, pepxmlfile, min_charge=1, max_charge=0, allowed_peptides=False, prefix='DECOY_', FDR_type=None, termini=set([2,1,0])):
        if allowed_peptides:
            allowed_peptides_set = set([x.strip() for x in open(allowed_peptides)])

        try:
            pepxml_params = {k: v for d in pepxml.iterfind(pepxmlfile, 'parameter name', read_schema=False) for k, v in d.items()}
            self.total_number_of_peptides_in_searchspace = int(pepxml_params.get('modelling, total peptides used', self.total_number_of_peptides_in_searchspace))
            self.total_number_of_proteins_in_searchspace = int(pepxml_params.get('modelling, total proteins used', self.total_number_of_proteins_in_searchspace))
            self.total_number_of_spectra = int(pepxml_params.get('modelling, total spectra used', self.total_number_of_spectra))
        except Exception as e:
            logger.critical('Error reading pepXML file: %s, %s ', e, e.args)
            return 0

        best_scores = {}
        standard_aminoacids = set(k for k in mass.std_aa_comp if '-' not in k)
        first_psm = True
        for record in pepxml.read(pepxmlfile, read_schema=False):
            if 'search_hit' in record:
                if int(min_charge) <= int(record['assumed_charge']) and (int(record['assumed_charge']) <= int(max_charge) or not max_charge):
                    if first_psm:
                        if 'num_missed_cleavages' not in record['search_hit'][0]:
                            logger.warning('Missed cleavages are missing in pepXML file, using 0 for all peptides')
                        try:
                            float(record['retention_time_sec'])
                        except:
                            try:
                                float(record['spectrum'].split(',')[2].split()[0])
                            except Exception:
                                logger.warning('RT experimental is missing in pepXML file, using 0 value for all peptides')
                        first_psm = False
                    if 'peptide' in record['search_hit'][0]:
                        sequence = record['search_hit'][0]['peptide']
                        num_tol_term_tmp = record['search_hit'][0]['proteins'][0]['num_tol_term']
                        if num_tol_term_tmp in termini and (not allowed_peptides or sequence in allowed_peptides_set):
                            try:
                                evalue = record['search_hit'][0]['search_score']['expect']
                                # evalue = 1/record['search_hit'][0]['search_score']['hyperscore']
                            except:
                                try:
                                    evalue = 1.0 / float(record['search_hit'][0]['search_score']['ionscore'])
                                except IOError:
                                    'Cannot read e-value!'
                            tags = {}
                            for k in record['search_hit'][0]['search_score'].keys():
                                if k.startswith('tmt'):
                                    tags[k] = float(record['search_hit'][0]['search_score'][k])
                            if not (FDR_type.startswith('peptide') and best_scores.get(sequence, 1e6) < evalue) and not set(sequence).difference(standard_aminoacids):
                                if FDR_type.startswith('peptide'):
                                    best_scores[sequence] = evalue
                                mc = record['search_hit'][0].get('num_missed_cleavages', 0)
                                modifications = record['search_hit'][0]['modifications']
                                try:
                                    sumI = 10 ** float(record['search_hit'][0]['search_score']['sumI'])
                                except:
                                    sumI = 0
                                try:
                                    frag_mt = float(record['search_hit'][0]['search_score']['fragmentMT'])
                                except:
                                    frag_mt = None
                                spectrum = record['spectrum']
                                pcharge = record['assumed_charge']
                                mass_exp = record['precursor_neutral_mass']

                                if pepxmlfile not in infiles_dict:
                                    infiles_dict[pepxmlfile] = len(infiles_dict)
                                infile_current = infiles_dict[pepxmlfile]
                                pept = Peptide(sequence=sequence, settings=self.settings, evalue=evalue, pcharge=pcharge, mass_exp=mass_exp, modifications=modifications, modification_list=self.modification_list, custom_aa_mass=self.aa_list, sumI=sumI, mc=mc, infile=infile_current, frag_mt=frag_mt, tags=tags)
                                try:
                                    RT_exp = float(record['retention_time_sec']) / 60
                                except:
                                    try:
                                        RT_exp = float(spectrum.split(',')[2].split()[0])
                                    except:
                                        RT_exp = 0

                                if not all(protein['protein'].startswith(prefix) for protein in record['search_hit'][0]['proteins']):
                                    pept.note = 'target'
                                else:
                                    pept.note = 'decoy'
                                pept.num_tol_term = num_tol_term_tmp
                                pept.next_aa = record['search_hit'][0]['proteins'][0]['peptide_next_aa']
                                pept.prev_aa = record['search_hit'][0]['proteins'][0]['peptide_prev_aa']

                                if pept.sequence not in self.proteins_dict:
                                    for prot in record['search_hit'][0]['proteins']:
                                        prot_name = get_dbname(prot)
                                        if prot_name not in [protein.dbname for protein in self.proteins_dict[pept.sequence]]:
                                            #pept.num_tol_term = prot['num_tol_term']
                                            self.proteins_dict[pept.sequence].append(Protein(dbname=get_dbname(prot), description=prot.get('protein_descr', None)))
                                            #pept.parentproteins.append(Protein(dbname=get_dbname(prot), description=prot.get('protein_descr', None)))

                                if len(self.proteins_dict[pept.sequence]) and (not modifications or Counter(v['position'] for v in modifications).most_common(1)[0][1] <= 1):
                                    self.add_elem((pept, spectrum, RT_exp))
                                    # self.peptideslist.append(pept)
                                    # self.spectrumlist.append(spectrum)

        self.spectrumlist = np.array(self.spectrumlist)

        if FDR_type.startswith('peptide'):
            js = []
            j = len(self.peptideslist) - 1
            while j >= 0:
                if self.peptideslist[j].evalue > best_scores.get(self.peptideslist[j].sequence, 1e6):
                    js.append(j)
                j -= 1
            self.rem_elements(js)

    def get_RC(self):
        try:
            seqs = [pept.modified_sequence for pept in self.peptideslist]
            RTexp = self.RT_exp#[pept.RT_exp for pept in self.peptideslist]
            RC_def = achrom.RCs_gilar_rp
            RC_def['aa'].setdefault('U', RC_def['aa'].get('C', 0.0))
            RC_def['aa'].setdefault('O', RC_def['aa'].get('K', 0.0))
            aa_labels = set(RC_def['aa'].keys())
            for pept in self.peptideslist:
                for v in pept.modification_list.itervalues():
                    aa_labels.add(v)
            xdict = {}
            for key, val in RC_def['aa'].items():
                xdict[key] = [val, None]
            RC_dict = achrom.get_RCs_vary_lcp(seqs, RTexp, labels=aa_labels)
            for key, val in RC_dict['aa'].items():
                try:
                    xdict[key][1] = val
                except:
                    xdict[key] = [None, val]
            a, b, _, _ = aux.linear_regression([x[0] for x in xdict.values() if all(v != None for v in x)], [x[1] for x in xdict.values() if all(v != None for v in x)])
            for key, x in xdict.items():
                if x[1] == None:
                    x[1] = x[0] * a + b
                RC_dict['aa'][key] = x[1]
            if 'C' not in RC_dict['aa']:
                RC_dict['aa']['C'] = RC_dict['aa']['C*']
        except:
            logger.error('Error in get_RC for achrom model. Using RCs_gilar_rp')
            RC_dict = achrom.RCs_gilar_rp
        self.RC = RC_dict

    def calc_RT(self, calibrate_coeff=(1, 0, 0, 0), RTtype='achrom'):
        self.RT_predicted = []
        if RTtype == 'ssrcalc':
            if SSRCalc is None:
                logger.critical('SSRCalc not available. Make sure that mechanize is installed.')
                sys.exit(1)
            ps = list(set([peptide.sequence for peptide in self.peptideslist]))
            SSRCalc_RTs = SSRCalc.calculate_RH(ps[:], pore_size=100, ion_pairing_agent='FA')

        for peptide in self.peptideslist:
            if RTtype == 'achrom':
                RT_predicted = achrom.calculate_RT(peptide.modified_sequence, self.RC, raise_no_mod=False)
                if np.isinf(RT_predicted):
                    elems = peptide.modified_sequence.split('-')
                    if len(elems) > 1:
                        if not all(el in parser.std_amino_acids for el in elems[0]) and str(elems[0] + '-') not in self.RC['aa']:
                            self.RC['aa'][str(elems[0] + '-')] = self.RC['aa']['H-']
                        elif not all(el in parser.std_amino_acids for el in elems[-1]) and str('-' + elems[-1]) not in self.RC['aa']:
                            self.RC['aa'][str('-' + elems[-1])] = self.RC['aa']['-OH']
                    RT_predicted = achrom.calculate_RT(peptide.modified_sequence, self.RC, raise_no_mod=False)

            elif RTtype == 'ssrcalc':
                SSRCalc_RT = SSRCalc_RTs[peptide.sequence]
                if SSRCalc_RT is not None:
                    RT_predicted = float(SSRCalc_RT) * calibrate_coeff[0] + calibrate_coeff[1]
                else:
                    RT_predicted = 0
                    logger.error('SSRCalc error')
            elif RTtype == 'biolccc':
                RT_predicted = biolccc.calculateRT(peptide.sequence, biolccc.rpAcnTfaChain, biolccc.standardChromoConditions)
            else:
                logger.error('RT_type error')
            self.RT_predicted.append(RT_predicted)
        self.check_arrays()

    def filter_RT(self, RT_tolerance):
        js = []
        j = len(self.peptideslist) - 1
        while j >= 0:
            if abs(float(self.RT_predicted[j]) - float(self.RT_exp[j])) > float(RT_tolerance):
                js.append(j)
            j -= 1
        self.rem_elements(js)

    def get_calibrate_coeff(self):
        peptides = []
        peptides_added = {}
        for peptide, RT_exp, RT_predicted in izip(self.peptideslist, self.RT_exp, self.RT_predicted):
            if peptide.sequence not in peptides_added:
                peptides_added[peptide.sequence] = [RT_exp, ]
                peptides.append([RT_predicted, RT_exp])
            else:
                if any(abs(RT_exp - v) < 2 for v in peptides_added[peptide.sequence]):
                    pass
                else:
                    peptides_added[peptide.sequence].append(RT_exp)
                    peptides.append([RT_predicted, RT_exp])
        aux_RT = aux.linear_regression([val[0] for val in peptides], [val[1] for val in peptides])
        return aux_RT

    def filter_modifications(self, RT_type=None):
        #TODO
        pass

    def filter_decoy(self):
        js = []
        j = len(self.peptideslist) - 1
        while j >= 0:
            if self.peptideslist[j].note == 'decoy':
                js.append(j)
                #self.peptideslist.pop(j)
            j -= 1
        self.rem_elements(js)

    def filter_evalue_new(self, FDR=1, FDR2=1, useMP=True, drop_decoy=True, toprint=True):
        "A function for filtering PSMs by e-value and MP-score with some FDR"
        isdecoy = lambda x: x[0].note == 'decoy'
        escore = lambda x: float(x[0].evalue)
        mscore = lambda x: -float(x[0].peptscore)

        new_peptides = self.copy_empty()
        for infile in self.get_infiles():
            infile_peptides = []
            for val in self.get_izip_full():
            # for peptide, spectrum in izip(self.peptideslist, self.spectrumlist):
            #     if peptide.infile == infile:
                if val[0].infile == infile:
                    infile_peptides.append(val)
            filtered_peptides = aux.filter(infile_peptides, fdr=float(FDR)/100, key=escore, is_decoy=isdecoy, remove_decoy=False, formula=1, full_output=True)
            qvals_e = aux.qvalues(filtered_peptides, key=escore, is_decoy=isdecoy, reverse=False, remove_decoy=False, formula=1, full_output=True)
            try:
                best_cut_evalue = max(escore(p) for p in filtered_peptides)
                real_FDR = round(aux.fdr(filtered_peptides, is_decoy=isdecoy) * 100, 1)
            except:
                best_cut_evalue = 0
                real_FDR = 0
            if toprint:
                logger.info('%s %s e-value', real_FDR, best_cut_evalue)
            best_cut_peptscore = 1.1
            if useMP:
                tmp_peptides = []
                for p in infile_peptides:
                    if escore(p) > best_cut_evalue:
                        tmp_peptides.append(p)
                filtered_peptides = aux.filter(tmp_peptides, fdr=float(FDR2)/100, key=mscore, is_decoy=isdecoy, remove_decoy=False, formula=1, full_output=True)
                qvals_m = aux.qvalues(filtered_peptides, key=mscore, is_decoy=isdecoy, reverse=False, remove_decoy=False, formula=1, full_output=True)
                try:
                    best_cut_peptscore = min(float(p[0].peptscore) for p in filtered_peptides)
                    real_FDR = round(aux.fdr(filtered_peptides, is_decoy=isdecoy) * 100, 1)
                except:
                    best_cut_peptscore = 1.1
                    real_FDR = 0
                if toprint:
                    logger.info('%s %s MP score', real_FDR, best_cut_peptscore)
            for val in qvals_e:
                val[-1][0].qval = val[-2]
                new_peptides.add_elem(val[-1])
                # new_peptides.peptideslist.append(val[-1][0])
                # new_peptides.peptideslist[-1].qval = val[-2]
                # new_peptides.spectrumlist.append(val[-1][1])
            if useMP:
                for val in qvals_m:
                    val[-1][0].qval = val[-2]
                    new_peptides.add_elem(val[-1])
                    # new_peptides.peptideslist.append(val[-1][0])
                    # new_peptides.peptideslist[-1].qval = val[-2]
                    # new_peptides.spectrumlist.append(val[-1][1])
        # new_peptides.spectrumlist = np.array(new_peptides.spectrumlist)
        new_peptides.check_arrays()
        if drop_decoy:
            new_peptides.filter_decoy()
        return (new_peptides, best_cut_evalue, best_cut_peptscore)

    def get_izip_full(self):
        #TODO check the memory usage for this place
        # for val in izip(*(getattr(self, l) for l in self.listing)):
        for val in izip_longest(*(getattr(self, l) for l in self.listing)):
            yield val

    def copy_empty(self):
        new_peptides = PeptideList(self.settings)
        new_peptides.total_number_of_spectra = self.total_number_of_spectra
        new_peptides.total_number_of_PSMs = self.total_number_of_PSMs
        new_peptides.total_number_of_PSMs_decoy = self.total_number_of_PSMs_decoy
        new_peptides.total_number_of_proteins_in_searchspace = self.total_number_of_proteins_in_searchspace
        new_peptides.total_number_of_peptides_in_searchspace = self.total_number_of_peptides_in_searchspace
        new_peptides.calibrate_coeff = self.calibrate_coeff
        new_peptides.RC = self.RC
        new_peptides.modification_list = self.modification_list
        new_peptides.infiles = self.infiles
        new_peptides.proteins_dict = self.proteins_dict
        return new_peptides

    def update(self, new_peptides):
        self.settings = new_peptides.settings
        self.total_number_of_spectra = new_peptides.total_number_of_spectra
        self.total_number_of_PSMs = new_peptides.total_number_of_PSMs
        self.total_number_of_PSMs_decoy = new_peptides.total_number_of_PSMs_decoy
        self.total_number_of_proteins_in_searchspace = max(new_peptides.total_number_of_proteins_in_searchspace, self.total_number_of_proteins_in_searchspace)
        self.total_number_of_peptides_in_searchspace = max(new_peptides.total_number_of_peptides_in_searchspace, self.total_number_of_peptides_in_searchspace)
        self.calibrate_coeff = new_peptides.calibrate_coeff
        self.RC = new_peptides.RC
        self.modification_list = new_peptides.modification_list

        for l in self.listing:
            new_peptides_attr = getattr(new_peptides, l)
            if isinstance(new_peptides_attr, list):
                tmp = getattr(self, l)
                tmp.extend(new_peptides_attr)
                setattr(self, l, tmp)
            elif isinstance(new_peptides_attr, np.ndarray):
                setattr(self, l, np.append(getattr(self, l), new_peptides_attr))

        self.proteins_dict.update(new_peptides.proteins_dict)

    def remove_duplicate_spectra(self):
        sdict = dict()
        for peptide, spectrum in izip(self.peptideslist, self.spectrumlist):
            if spectrum not in sdict or peptide.evalue < sdict[spectrum]:
                sdict[spectrum] = peptide
        js = []
        for idx, spectrum in enumerate(self.spectrumlist):
            if spectrum not in sdict:
                js.append(idx)
        self.rem_elements(js)

    def add_elem(self, val):
        for l, v in zip(self.listing, val):
            getattr(self, l).append(v)

    def check_arrays(self):
        for l in self.listing_nparrays:
            tmp = getattr(self, l)
            if not isinstance(tmp, np.ndarray):
                setattr(self, l, np.array(tmp))

    def remove_duplicate_sequences(self):
        edict = dict()
        for peptide in self.peptideslist:
            edict[peptide.sequence] = min(float(peptide.evalue), edict.get(peptide.sequence, np.inf))
        new_peptides = self.copy_empty()
        for val in self.get_izip_full():
            peptide_ind = self.listing.index('peptideslist')
            if val[peptide_ind].evalue == edict.get(val[peptide_ind].sequence, None):
                new_peptides.add_elem(val)
                del edict[val[peptide_ind].sequence]
        new_peptides.check_arrays()
        return new_peptides

    def cut_left(self, msize):
        for l in self.listing:
            setattr(self, l, getattr(self,l)[msize:])

    def get_right(self, extpeptides, msize):
        for l in self.listing:
            setattr(self, l, getattr(extpeptides,l)[:msize])

class Protein:
    __slots__ = ['dbname', 'description']

    def __init__(self, dbname, description='Unknown'):
        self.dbname = dbname
        self.description = description

class Peptide:
    __slots__ = ['sequence', 'modified_sequence', 'modification_list', 'pcharge',
                 'aa_mass', 'pmass', 'mz', 'evalue', 'massdiff',
                 'num_missed_cleavages', 'mc', 'note', 'note2', 'note3', 'protscore2', 'peptscore',
                 'peptscore2', 'spectrum_mz', 'fragment_mt', 'sumI', 'it',
                 'infile', 'fragments', 'valid_sequence', 'prev_aa', 'next_aa', 'num_tol_term']
    def __init__(self, sequence, settings, pcharge=0, evalue=0, note='unknown', mass_exp=0, modifications=[], modification_list={}, custom_aa_mass=None, sumI=0, mc=None, infile='unknown', frag_mt=None, tags=None):
        self.sequence = sequence
        self.modified_sequence = sequence
        self.modification_out_str = ''
        self.modification_list = modification_list
        self.pcharge = int(pcharge)
        self.aa_mass = custom_aa_mass
        self.pmass = float(mass.fast_mass(sequence=self.sequence, charge=0)) - 18.0105646837 + settings.getfloat('modifications', 'protein nterm cleavage') + settings.getfloat('modifications', 'protein cterm cleavage')
        for modif in modifications:
            self.pmass += modif['mass']
            if modif['position'] not in [0, len(self.sequence) + 1]:
                aminoacid = self.sequence[modif['position'] - 1]
                self.pmass -= mass.std_aa_mass[aminoacid]
            else:
                if modif['position'] == 0:
                    self.pmass -= settings.getfloat('modifications', 'protein nterm cleavage')
                else:
                    self.pmass -= settings.getfloat('modifications', 'protein cterm cleavage')

        self.mz = (mass_exp + pcharge * 1.007276) / pcharge
        self.modified_peptide(modifications)
        # self.RT_exp = RT_exp
        # self.RT_predicted = False
        self.evalue = float(evalue)
        #self.parentproteins = []
        self.massdiff = float(mass_exp) - float(self.pmass)
        self.num_missed_cleavages = dict()
        self.mc = mc
        self.note = note
        self.note2 = ''
        self.note3 = ''
        self.protscore2 = 1
        self.peptscore = 1
        self.peptscore2 = 1
        self.spectrum_mz = None
        self.fragment_mt = frag_mt
        self.sumI = sumI# / self.pcharge
        self.it = 1.0
        self.infile = infile
        self.fragments = defaultdict(dict)
        self.valid_sequence = dict()
        self.tags = tags if len(tags) else None

    def theor_spectrum(self, types=('b', 'y'), maxcharge=None, **kwargs):
        peaks = {}
        if not maxcharge:
            maxcharge = max(self.pcharge, 1)
        for ion_type in types:
            ms = []
            for i in range(1, len(self.modified_sequence)):
                if self.modified_sequence[i - 1] in parser.std_amino_acids and self.modified_sequence[i] != '-':
                    for charge in range(1, maxcharge + 1):
                        if ion_type[0] in 'abc':
                            ms.append(mass.fast_mass2(
                                str(self.modified_sequence)[:i], ion_type=ion_type, charge=charge,
                                **kwargs))
                        else:
                            ms.append(mass.fast_mass2(
                                str(self.modified_sequence)[i:], ion_type=ion_type, charge=charge,
                                **kwargs))
            marr = np.array(ms)
            marr.sort()
            peaks[ion_type] = marr
        return peaks

    def get_missed_cleavages(self, protease='trypsin'):
        if protease not in self.num_missed_cleavages:
            self.num_missed_cleavages[protease] = parser.num_sites(self.sequence, protease)
        return self.num_missed_cleavages[protease]

    def count_modifications(self, label):
        if '-' in label:
            naa = 1
            nmods = self.modified_sequence.count(label)
        else:
            naa = self.modified_sequence.count(label[-1])
            nmods = self.modified_sequence.count(label)
        if naa:
            if nmods:
                return 1
            else:
                return 0
        else:
            return -1.0

    def get_median_fragment_mt(self, settings=None):
        if not self.fragment_mt and len(self.spectrum_mz) > 1:
            temp = settings.get('fragment mass', 'ion types')
            ion_types = (x.strip() for x in temp.split(','))
            acc = settings.getfloat('fragment mass', 'mass accuracy')
            spectrum_mz = copy(self.spectrum_mz)
            int_array = copy(self.spectrum_i)
            i = int_array > int_array.max() / 100
            spectrum_mz = spectrum_mz[i]
            int_array = int_array[i]
            theor = self.theor_spectrum(types=ion_types, aa_mass=self.aa_mass, maxcharge=1)
            spectrum_KDTree = cKDTree(spectrum_mz.reshape((spectrum_mz.size, 1)))
            dist_total = np.array([])
            for itype, fragments in theor.iteritems():
                self.fragments[itype]['m/z'] = fragments
                n = fragments.size
                dist, ind = spectrum_KDTree.query(fragments.reshape((n, 1)),
                    distance_upper_bound=acc)
                dist_total = np.append(dist_total, dist[dist != np.inf])
                self.fragments[itype]['intensity'] = np.zeros(len(fragments))
                for idx in range(len(dist)):
                    if dist[idx] != np.inf:
                        self.fragments[itype]['intensity'][idx] += int_array[ind[idx]]

            for itype, val in self.fragments.iteritems():
                if itype == 'b':
                    try:
                        flag = max(idx for idx, intensity in enumerate(val['intensity']) if intensity) + 1
                    except:
                        flag = 0
                    self.valid_sequence[itype] = self.sequence[:flag] + self.sequence[flag:].lower()
                elif itype == 'y':
                    try:
                        flag = len(self.sequence) - (max(idx for idx, intensity in enumerate(val['intensity']) if intensity) + 1)
                    except:
                        flag = len(self.sequence)
                    self.valid_sequence[itype] = self.sequence[:flag].lower() + self.sequence[flag:]

            if dist_total.size:
                self.fragment_mt = np.median(dist_total)
            else:
                self.fragment_mt = acc
        return self.fragment_mt

    def mass_diff(self):
        """Calculates a difference between theoretical and experimental masses. Takes into account an isotope mass difference error"""
        return (self.massdiff - round(self.massdiff, 0) * 1.0033548378) / (self.pmass - round(self.massdiff, 0) * 1.0033548378) * 1e6

    def modified_peptide(self, modifications):
        def add_modification(arg, term=None):
            i = ''
            done_flag = 0
            while 1:
                for x in string.ascii_lowercase:
                    if i + x not in self.modification_list.values():
                        self.modification_list[arg] = i + x
                        if term and term == 'c':
                             self.modification_list[arg] = '-' + self.modification_list[arg]
                        elif term and term == 'n':
                             self.modification_list[arg] += '-'
                        else:
                            logger.warning('label for %s modification is missing in parameters, using %s label', arg, self.modification_list[arg])
                        done_flag = 1
                        break
                if not done_flag:
                    if i and i[-1] != string.ascii_lowercase[-1]:
                        i = i[:-1] + string.ascii_lowercase.index(i[-1] + 1)
                    else:
                        i += string.ascii_lowercase[0]
                else:
                    break

        self.modified_sequence = str(self.sequence)
        for modif in sorted(modifications, key=itemgetter('position'), reverse=True):
            if modif['position'] == 0:
                modname = get_modif_name(modif, 'n')
                try:
                    self.modified_sequence = self.modification_list[modname] + self.modified_sequence
                except:
                    add_modification(modname, term='n')
                    logger.warning('label for %s nterm modification is missing in parameters, using %s label', modname, self.modification_list[modname])
                    self.aa_mass[self.modification_list[modname]] = float(modif['mass'])
                    self.modified_sequence = self.modification_list[modname] + self.modified_sequence
                self.modification_out_str += '[n-term] %s,' % (self.modification_list[modname].replace('-', ''))
            elif modif['position'] == len(self.sequence) + 1:
                modname = get_modif_name(modif, 'c')
                try:
                    self.modified_sequence = self.modified_sequence + self.modification_list[modname]
                except:
                    add_modification(modname, term='c')
                    logger.warning('label for %s cterm modification is missing in parameters, using label %s', modname, self.modification_list[modname])
                    self.aa_mass[self.modification_list[modname]] = float(modif['mass'])
                    self.modified_sequence = self.modified_sequence + self.modification_list[modname]
                self.modification_out_str += '[c-term] %s,' % (self.modification_list[modname].replace('-', ''))
            else:
                pos = modif['position']
                modname = get_modif_name(modif, self.modified_sequence[pos - 1])
                try:
                    self.modified_sequence = self.modified_sequence[:pos - 1] + self.modification_list[modname] + self.modified_sequence[pos - 1:]
                except:
                    add_modification(modname)
                    self.aa_mass[self.modification_list[modname]] = float(modif['mass'])
                    self.modified_sequence = self.modified_sequence[:pos - 1] + self.modification_list[modname] + self.modified_sequence[pos - 1:]
                self.modification_out_str += '[%d] %s (%s),' % (modif['position'], self.modification_list[modname], self.sequence[pos - 1])
        if self.modification_out_str:
            self.modification_out_str = self.modification_out_str[:-1]
            self.modification_out_str = ', '.join(self.modification_out_str.split(',')[::-1])
        else:
            self.modification_out_str = ' '
            
