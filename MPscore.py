from MPlib import PeptideList, Descriptor, get_settings, filter_evalue_prots, FDbinSize
import sys
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
from os import path, listdir
from pyteomics import mzml, fasta, auxiliary, mgf, parser, pepxmltk, mass
from Queue import Empty
import multiprocessing
import shutil
from time import sleep, time
import pickle
from copy import copy, deepcopy
from collections import defaultdict, Counter
from itertools import izip
import logging
import logging.config
logger = logging.getLogger(__name__)

try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#ffffff'})
    seaborn.set_style('whitegrid')
except ImportError:
    pass

manager = multiprocessing.Manager()
protsC = {}
protsL = dict()#manager.dict()
protsS = dict()#manager.dict()
protsN = dict()#manager.dict()
stime = time()
redcolor='#FC6264'
bluecolor='#70aed1'
greencolor='#8AA413'
def get_output_string(obj, spectrum, RT_exp, type, fragments_info=False, fragments_info_zeros=False, peptide_count=False, proteins_dict={}, tags=None):
    if type == 'psm':
        out = '%s\t' % (obj.sequence, )
        if peptide_count:
            out += '%s\t' % (peptide_count, )
        out += '%s\t%s\t%0.3f\t%d\t%0.1f\t%d\t%d\t%s\t%s\t%0.2E\t%0.2E\t%0.2f\t%s\t' % (obj.modified_sequence, obj.modification_out_str, obj.mz, obj.pcharge, obj.mass_diff(), obj.mc,
                                                             obj.num_tol_term, obj.prev_aa, obj.next_aa, obj.evalue, obj.peptscore, RT_exp, spectrum)
        for protein in proteins_dict[obj.sequence]:
            out += '%s;' % (protein.dbname,)
        out += '\t'
        for protein in proteins_dict[obj.sequence]:
            out += '%s;' % (protein.description,)
        out += '\t%s\t%0.3f\t%s' % (obj.sumI, obj.massdiff, obj.note)

        if tags:
            out += '\t'
            out += '\t'.join(str(obj.tags[t]) for t in tags)
        if fragments_info:
            for itype, val in obj.fragments.iteritems():
                out += '\t'
                for idx, mz in enumerate(val['m/z']):
                    if fragments_info_zeros or int(val['intensity'][idx]):
                        out += '%s:%s;' % (round(mz, 3), int(val['intensity'][idx]))
            out += '\t%s\t%s' % (obj.valid_sequence['b'], obj.valid_sequence['y'])
        out += '\n'
    return out

def calc_sq(protein, peptides):
    if not protein:
        return 0
    psq = [False for x in protein]
    plen = len(protein)
    for pep in peptides:
        csize = len(pep)
        for j in range(plen):
            if protein[j:j+csize] == pep:
                for y in range(csize):
                    psq[j + y] = True
    return float(sum(psq)) / len(psq) * 100


def handle(q, q_output, settings, protsL):
    while 1:
        try:
            filenames = q.get(timeout=1)
        except Empty:
            q_output.put('1')
            break

        FDR_type = settings.get('options', 'FDR_type')
        if FDR_type.startswith('protein'):
            FDR = 1.0
        else:
            FDR = settings.getfloat('options', 'FDR')

        RT_type = settings.get('retention time', 'model')

        min_charge = settings.getint('charges', 'min charge')
        max_charge = settings.getint('charges', 'max charge')

        valid_proteins = []
        valid_proteins_input = settings.get('options', 'valid proteins')
        if path.isfile(valid_proteins_input):
            for line in open(valid_proteins_input, 'r'):
                if len(line.strip().split()) > 1:
                    dbname, conc = line.strip().split()
                    valid_proteins.append(dbname)
                    protsC[dbname] = float(conc)
                else:
                    dbname = line.strip().split()[0]
                    valid_proteins.append(dbname)
        else:
            valid_proteins = []

        spectra_dict = dict()
        spectra_dict_intensities = dict()

        peptides = PeptideList(settings)
        mods = manager.dict()
        iprocs = []
        inprocs = peptides.settings.getint('options', 'threads')
        iq = multiprocessing.Queue()
        iq_output = multiprocessing.Queue()

        def getpepxml(iq, iq_output, settings, mods=False):
            allowed_termini = set(int(z.strip()) for z in settings.get('missed cleavages', 'number_of_enzyme_termini').split(','))
            for curfile in iter(iq.get, None):
                qpeptides = PeptideList(settings, mods)
                qpeptides.get_from_pepxmlfile(curfile['.pep'], min_charge=min_charge, max_charge=max_charge,
                allowed_peptides=settings.get('advanced options', 'allowed peptides'),
                prefix=settings.get('input', 'decoy prefix'), FDR_type=settings.get('options', 'FDR_type'),
                termini=allowed_termini)

                if len(qpeptides.peptideslist):
                    mzmlfile = curfile.get('.mzML', None)
                    if mzmlfile:
                        logger.info('Processing mzML')
                        # isolation_window = settings.getfloat('precursor ion fraction', 'isolation window')
                        # mass_acc = settings.getfloat('precursor ion fraction', 'mass accuracy')
                        # spectra_ms1 = []
                        # spectra_ms2 = []
                        for spectrum in mzml.read(mzmlfile):
                            if spectrum['ms level'] == 2:
                                spectra_dict[spectrum['id'].strip()] = spectrum['m/z array']
                                spectra_dict_intensities[spectrum['id'].strip()] = spectrum['intensity array']
                                #spectra_ms2.append(x)
                            # elif x['ms level'] == 1:
                            #     spectra_ms1.append(x)
                        # for psm in qpeptides.peptideslist:
                        #     j = len(spectra_ms1) - 1
                        #     while spectra_ms1[j]['scanList']['scan'][0]['scan start time'] > psm.RT_exp:
                        #         j -= 1
                        #     basemz = spectra_ms1[j]
                        #     I = []
                        #     Ip = 0
                        #     for idx, mz in enumerate(spectra_ms1[j]['m/z array']):
                        #         if abs(mz - (float(psm.mass_exp + 1.007825 * psm.pcharge) / psm.pcharge)) <= isolation_window:
                        #             if any(abs(mz - (float(psm.mass_exp + k + 1.007825 * psm.pcharge) / psm.pcharge)) <= mz * mass_acc * 1e-6 for k in [-2, -1, 0, 1, 2]):
                        #                 Ip += float(spectra_ms1[j]['intensity array'][idx])
                        #             I.append(float(spectra_ms1[j]['intensity array'][idx]))
                        #     PIF = Ip / sum(I) * 100
                        #     psm.I = Ip
                        #     psm.PIF = PIF

                        # itimes = dict()
                        # for x in spectra_ms2:
                        #     sc_n = x['id'].split('scan=')[-1]
                        #     itimes[sc_n] = float(x['scanList']['scan'][0]['ion injection time'])
                        # for peptide in qpeptides.peptideslist:
                        #     try:
                        #         peptide.it = itimes[peptide.spectrum.split('.')[1]]
                        #     except:
                        #         print 'Smth wrong with mzML indexes'
                        #         print peptide.spectrum.split('.')[1], itimes[-1]

                    if qpeptides.settings.getboolean('descriptors', 'fragment mass tolerance, Da') and not qpeptides.peptideslist[0].fragment_mt:
                        mgffile = curfile.get('.mgf', None)
                        if mgffile:
                            logger.info('Processing MGF')
                            spectra = mgf.read(mgffile)
                            for spectrum in spectra:
                                spectra_dict[spectrum['params']['title'].strip()] = spectrum['m/z array']
                                spectra_dict_intensities[spectrum['params']['title'].strip()] = spectrum['intensity array']
                            protsL['total spectra'] = len(spectra_dict)
                            if not qpeptides.total_number_of_spectra:
                                qpeptides.total_number_of_spectra = len(spectra_dict)
                        if not spectra_dict and qpeptides.settings.getboolean('descriptors', 'fragment mass tolerance, Da'):
                            qpeptides.settings.set('descriptors', 'fragment mass tolerance, Da', '0')
                            logger.warning('Fragment mass tolerance was turned off due to missing MGF file')
                        if qpeptides.settings.getboolean('descriptors', 'fragment mass tolerance, Da'):
                            for peptide, spectrum in izip(qpeptides.peptideslist, qpeptides.spectrumlist):
                                try:
                                    peptide.spectrum_mz = spectra_dict[spectrum.split(' RTINSECONDS=')[0].strip()]
                                    peptide.spectrum_i = spectra_dict_intensities[spectrum.split(' RTINSECONDS=')[0].strip()]
                                except:
                                    try:
                                        peptide.spectrum_mz = spectra_dict[spectrum.strip() + ' min']
                                        peptide.spectrum_i = spectra_dict_intensities[spectrum.strip() + ' min']
                                    except:
                                        peptide.spectrum_mz = spectra_dict[spectrum.strip()]
                                        peptide.spectrum_i = spectra_dict_intensities[spectrum.strip()]
                            for peptide in qpeptides.peptideslist:
                                peptide.get_median_fragment_mt(qpeptides.settings)
                                peptide.spectrum_mz = None
                                peptide.spectrum_i = None

                    tmp_peptides = qpeptides.copy_empty()
                    msize = 10000
                    while len(qpeptides.peptideslist):
                        tmp_peptides.get_right(qpeptides, msize)
                        # tmp_peptides.peptideslist = qpeptides.peptideslist[:msize]
                        # tmp_peptides.spectrumlist = qpeptides.spectrumlist[:msize]
                        iq_output.put(copy(tmp_peptides))
                        qpeptides.cut_left(msize)
                        # qpeptides.peptideslist = qpeptides.peptideslist[msize:]
                        # qpeptides.spectrumlist = qpeptides.spectrumlist[msize:]
                iq_output.put(None)

        for filename in filenames:
            iq.put(filename)
        for i in range(inprocs):
            iq.put(None)

        for i in range(inprocs):
            p = multiprocessing.Process(target=getpepxml, args=(iq, iq_output, settings, mods))
            iprocs.append(p)
            p.start()

        j = 0
        while j < len(filenames):
            for res_peptides in iter(iq_output.get, None):
                peptides.update(res_peptides)
            j += 1
        if settings.getboolean('advanced options', 'choose_best_spectra_results'):
            for peptide in peptides.peptideslist:
                peptide.infile = 0
            peptides.remove_duplicate_spectra()
        peptides.total_number_of_PSMs = len(peptides.peptideslist)

        for p in iprocs:
            p.terminate()

        logger.info('Total number of PSMs = %s', len(peptides.peptideslist))
        logger.info('Total number of peptides: %s', len(set(pept.sequence for pept in peptides.peptideslist)))
        if len(peptides.peptideslist):
            prots_dict = defaultdict(int)
            pepts_dict = defaultdict(int)
            for peptide in peptides.peptideslist:
                pepts_dict[peptide.sequence] += 1
                for protein in peptides.proteins_dict[peptide.sequence]:#peptide.parentproteins:
                    prots_dict[protein.dbname] += 1
                    if peptide.note == 'decoy':
                        protein.note = 'W'
                        if peptide.note2 != 'tr' or peptide.note == 'decoy':
                            peptide.note2 = 'wr'
                    else:
                        if protein.dbname in valid_proteins:
                            peptide.note3 = 'valid'
                        protein.note = 'Valid'
                        peptide.note2 = 'tr'
            peptides.total_number_of_PSMs_decoy = sum(1 for pept in peptides.peptideslist if pept.note2 == 'wr')

            if FDR_type.startswith('peptide'):
                logger.debug('Choosing best PSM per peptide sequence')
                peptidesdict = dict()
                for peptide, spectrum in izip(peptides.peptideslist, peptides.spectrumlist):
                    if peptide.sequence not in peptidesdict:
                        peptidesdict[peptide.sequence] = [spectrum, peptide.evalue]
                    elif peptide.evalue < peptidesdict[peptide.sequence][1]:
                        peptidesdict[peptide.sequence] = [spectrum, peptide.evalue]
                passed = set([x[0] for x in peptidesdict.values()])

                js = []
                j = len(peptides.peptideslist) - 1
                while j >= 0:
                    if peptides.spectrumlist[j] not in passed:
                        js.append(j)
                    j -= 1
                peptides.rem_elements(js)
                del js

            logger.debug('Calculating number of peptides per protein')
            misprotflag = 0
            for peptide in peptides.peptideslist:
                peptide.peptscore2 = pepts_dict[peptide.sequence]
                for protein in peptides.proteins_dict[peptide.sequence]:
                    if protein.dbname not in protsL:
                        misprotflag = +1
                        protsL[protein.dbname] = 5000
                        protsN[protein.dbname] = 50
                        protsL['total proteins'] += 1
                        protsL['total peptides'] += 50
                    if peptide.protscore2 < float(prots_dict[protein.dbname]) / protsL[protein.dbname] * 500:
                        peptide.protscore2 = float(prots_dict[protein.dbname]) / protsL[protein.dbname] * 500

            if misprotflag:
                logger.warning('Some proteins are missing in FASTA. Using length 5000 and 50 theoretical peptides for normalization and emPAI calculation\n'
                      'It is highly recommended to check the FASTA file for correct work of MPscore')
            pepts_dict = None
            prots_dict = None
            logger.info('Starting first FDR filtering...')
            copy_peptides, threshold0, _ = peptides.filter_evalue_new(FDR=FDR, useMP=False)

            logger.info('Default filtering:')
            numPSMs, numpeptides_true, numprots_true = PSMs_info(copy_peptides, valid_proteins, settings)
            if numPSMs > 1:
                descriptors = []
                if numPSMs >= 50:
                    dname = 'RT difference, min'
                    if peptides.settings.getboolean('descriptors', dname) and not np.allclose(peptides.RT_exp, 0):
                        descriptors.append(Descriptor(name=dname, massive_formula=lambda peptides: peptides.RT_exp - peptides.RT_predicted, group='A', binsize='auto'))
                    dname = 'precursor mass difference, ppm'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.mass_diff(), group='A', binsize='auto'))
                    dname = 'missed cleavages, protease 1'
                    if peptides.settings.getboolean('descriptors', dname.split(',')[0]):
                        protease1 = [x.strip() for x in settings.get('missed cleavages', 'protease1').split(',')]
                        expasy1 = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules else protease for protease in protease1))
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.get_missed_cleavages(expasy1), group='A', binsize=1))
                    dname = 'missed cleavages, protease 2'
                    if peptides.settings.getboolean('descriptors', dname.split(',')[0]):
                        protease2 = [x.strip() for x in settings.get('missed cleavages', 'protease2').split(',')]
                        expasy2 = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules else protease for protease in protease2))
                        if protease2[0]:
                            descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.get_missed_cleavages(expasy2), group='A', binsize=1))

                    dname = 'charge states'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.pcharge, group='A', binsize=1))
                    dname = 'potential modifications'
                    if peptides.settings.getboolean('descriptors', dname):
                        temp = settings.get('modifications', 'variable')
                        if temp:
                            for mod in temp.replace(' ', '').split(','):
                                descriptors.append(Descriptor(name='%s, %s' % (dname, mod), single_formula=lambda peptide, mod=mod: peptide.count_modifications(mod), group='A', binsize=1))
                    dname = 'isotopes mass difference, Da'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: round(peptide.massdiff, 0), group='A', binsize=1))
                    dname = 'PSMs per protein'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.protscore2, group='B'))
                    dname = 'PSM count'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.peptscore2, group='B'))
                    dname = 'fragment mass tolerance, Da'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.get_median_fragment_mt(), group='A', binsize = 'auto'))
                    dname = 'PIF'
                    if peptides.settings.getboolean('descriptors', dname):
                        descriptors.append(Descriptor(name=dname, single_formula=lambda peptide: peptide.PIF, group='B'))

                if 'RT difference, min' in [d.name for d in descriptors] and numPSMs >= 50:
                    if RT_type == 'achrom':
                        copy_peptides.get_RC()
                        peptides.RC = copy_peptides.RC
                        if peptides.settings.getint('advanced options', 'saveRC'):
                            pickle.dump(peptides.RC, open('RC.pickle', 'w'))
                        copy_peptides.calc_RT(RTtype=RT_type)
                        logger.debug('Calibration coefficients: %s', copy_peptides.get_calibrate_coeff())
                        peptides.calc_RT(RTtype=RT_type)
                    else:
                        copy_peptides.filter_modifications(RT_type=RT_type)
                        copy_peptides.calc_RT(RTtype=RT_type)
                        calibrate_coeff = copy_peptides.get_calibrate_coeff()
                        copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
                        copy_peptides.filter_RT(RT_tolerance=3 * calibrate_coeff[3])
                        copy_peptides.calc_RT(RTtype=RT_type)
                        calibrate_coeff = copy_peptides.get_calibrate_coeff()
                        logger.debug('Calibration coefficients: %s', calibrate_coeff)
                        peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
                        copy_peptides, _, _ = peptides.filter_evalue_new(FDR=FDR, useMP=False)
                        copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)

                curfile = filenames[-1]['.pep']
                if descriptors:
                    descriptors = prepare_hist(descriptors, copy_peptides, first=False)
                    fig = plot_histograms(descriptors, peptides, FDR, curfile, savesvg=settings.getboolean('advanced options', 'saveSVG'), sepfigures=settings.getboolean('advanced options', 'separatefigures'))

                    for idx, cpeptscore in calc_peptscore(peptides, descriptors):
                        peptides.peptideslist[idx].peptscore = cpeptscore

                    js = []
                    j = len(peptides.peptideslist) - 1
                    while j >= 0:
                        if peptides.peptideslist[j].peptscore == 0:
                            js.append(j)
                        j -= 1
                    peptides.rem_elements(js)
                    del js

                    k_temp = []
                    while len(k_temp) < 3 or k_temp[-1] != k_temp[-3]:
                        if not k_temp:
                            copy_peptides = peptides.filter_evalue_new(FDR=FDR, useMP=False, toprint=False)[0]
                        else:
                            FDR_new = (((FDR / 100 - float(peptides.total_number_of_PSMs_decoy) / numPSMs_true * k)) / (1 - k)) * 100
                            copy_peptides = peptides.filter_evalue_new(FDR=FDR, FDR2=FDR_new, useMP=True, toprint=False)[0]

                        numPSMs_true = len(copy_peptides)
                        try:
                            k = float(copy_peptides.get_number_of_peptides()) / peptides.total_number_of_peptides_in_searchspace
                        except:
                            k = 0
                        FDR_new = (((FDR / 100 - float(peptides.total_number_of_PSMs_decoy) / numPSMs_true * k)) / (1 - k)) * 100
                        k_temp.append(float(k))

                    plot_MP(descriptors, peptides, fig, FDR, FDR_new, valid_proteins, settings, threshold0, curfile)
                else:
                    fig = plt.figure(figsize=(16, 12))
                    DPI = fig.get_dpi()
                    fig.set_size_inches(2000.0/float(DPI), 2000.0/float(DPI))
                    plot_MP(descriptors, peptides, fig, FDR, 0, valid_proteins, settings, threshold0, curfile)


def find_optimal_xy(descriptors):
    x, y = 1, 1
    while x * y < len(descriptors) + 11:
        if x > y:
            y += 1
        else:
            x += 1
    return x, y


def nsaf(prots, norm=True):
    sumI = 0
    for protein in prots:
        sumI += float(prots[protein]['PSMs'])/protsL[protein]
    for protein in prots:
        prots[protein]['NSAF'] = float(prots[protein]['PSMs']) / (1 if not norm else sumI) / protsL[protein]
    return prots, sumI


def calc_emPAI(prots, protsN, norm=True):
    try:
        for dbname in prots:
            PAI = float(prots[dbname]['PSMs']) / max(protsN[dbname], 1)
            emPAI = 10 ** PAI - 1
            prots[dbname]['emPAI'] = emPAI

        sum_emPAI = sum(val['emPAI'] for val in prots.itervalues())
        if norm:
            for dbname in prots.keys():
                prots[dbname]['emPAI'] = prots[dbname]['emPAI'] / sum_emPAI
    except OverflowError:
        logger.warning('OverflowError encountered in emPAI calculation. emPAI calculation disabled.')
        for dbname, prot in prots.iteritems():
            prot['emPAI'] = 1
        sum_emPAI = 1
    return prots, sum_emPAI


def PSMs_info(peptides, valid_proteins, settings, fig=False, printresults=True, tofile=False, curfile=False, loop=True, ox=False, oy=False):
    def keywithmaxval(d):
        #this method is much faster than using max(prots.iterkeys(), key=(lambda key: prots[key]))
        # v=list(d.values())
        # k=list(d.keys())
        maxv = max(d.values())
        tmp = [k for k, v in d.iteritems() if v == maxv]
        return sorted(tmp)[0]
    tostay = set()
    prots = defaultdict(int)
    prots_pep = defaultdict(set)
    peptides_added = defaultdict(set)
    for peptide in peptides.peptideslist:
        if peptide.sequence not in peptides_added:
            if peptide.note2 == 'wr':
                add_label = ''#settings.get('input', 'decoy prefix')
            else:
                add_label = ''
            for protein in peptides.proteins_dict[peptide.sequence]:#peptide.parentproteins:
                tmp_dbname = add_label + protein.dbname
                prots[tmp_dbname] += 1
                prots_pep[tmp_dbname].add(peptide.sequence)
                peptides_added[peptide.sequence].add(tmp_dbname)
    while peptides_added and loop:
        bestprot = keywithmaxval(prots)
        tostay.add(bestprot)
        for pep in prots_pep[bestprot]:
            for k in peptides_added[pep]:
                prots[k] -= 1
            del peptides_added[pep]

    prots = dict()
    peptides_added = set()
    true_prots = set()
    Total_prots = set()

    for peptide in peptides.peptideslist:
        # Peptide sumI normalization!
        if tofile:
            peptide.sumI = round(peptide.sumI / peptide.it, 2)
        if peptide.note2 == 'wr':
            add_label = ''#settings.get('input', 'decoy prefix')
        else:
            add_label = ''

        for protein in peptides.proteins_dict[peptide.sequence]:
            tmp_dbname = add_label + protein.dbname
            Total_prots.add(tmp_dbname)
            try:
                prots[tmp_dbname]['PSMs'] += 1
                prots[tmp_dbname]['sumI'] += peptide.sumI
                prots[tmp_dbname]['pept'].add(peptide.sequence)
            except:
                prots[tmp_dbname] = dict()
                prots[tmp_dbname]['pept'] = set([peptide.sequence, ])
                prots[tmp_dbname]['PSMs'] = 1
                prots[tmp_dbname]['sumI'] = peptide.sumI
                prots[tmp_dbname]['evalues'] = []
                prots[tmp_dbname]['expect'] = 1
                prots[tmp_dbname]['description'] = protein.description
            if tmp_dbname in valid_proteins and peptide.note != 'decoy':
                true_prots.add(tmp_dbname)
        if peptide.sequence not in peptides_added:
            for protein in peptides.proteins_dict[peptide.sequence]:
                tmp_dbname = add_label + protein.dbname
                try:
                    prots[tmp_dbname]['Peptides'] += 1
                except:
                    prots[tmp_dbname]['Peptides'] = 1
        peptides_added.add(peptide.sequence)

    def calc_expect_log(es):
        n = len(es)
        es_new = []
        for x in es:
            if x > 1:
                x = 1
            elif x == 0:
                x = 1e-15
            es_new.append(x)
        es = list(es_new)
        if n == 1:
            return np.log10(es[0])
        expect = sum([np.log10(x) for x in es])
        return expect

    for dbname in prots.keys():
        if 'Peptides' not in prots[dbname]:
            del prots[dbname]

    for k in prots:
        prots[k]['fullgroup'] = set()

    for pep in peptides.peptideslist:
        tprots = set([pr.dbname for pr in peptides.proteins_dict[pep.sequence]])
        for k in tostay:
            if k in tprots:
                prots[k]['fullgroup'].update(tprots)
    for k in tostay:
        prots[k]['fullgroup'] = ';'.join(prots[k]['fullgroup'])

    prots_full = deepcopy(prots)

    for dbname in prots.keys():
        if loop and dbname not in tostay:
            del prots[dbname]

    if tofile:
        new_peptides = peptides.remove_duplicate_sequences()
        dec_label = settings.get('input', 'decoy prefix')
        for peptide in new_peptides.peptideslist:
            if peptide.note2 == 'wr':
                add_label = ''#dec_label
            else:
                add_label = ''
            for protein in new_peptides.proteins_dict[peptide.sequence]:
                tmp_dbname = add_label + protein.dbname
                if tmp_dbname in prots: # <---- WTF??? not works with tmp_dbname in tostay
                    prots[tmp_dbname]['evalues'].append(peptide.qval)
                if tmp_dbname in prots_full:
                    prots_full[tmp_dbname]['evalues'].append(peptide.qval)
        for k in prots:
            prots[k]['expect'] = calc_expect_log(prots[k]['evalues'])
        for k in prots_full:
            prots_full[k]['expect'] = calc_expect_log(prots_full[k]['evalues'])

        FDR_type = settings.get('options', 'FDR_type')
        remove_decoy = settings.getint('advanced options', 'remove_decoy')
        if FDR_type.startswith('protein'):
            protFDR = settings.getfloat('options', 'FDR')
            prots = filter_evalue_prots(prots, FDR=protFDR, remove_decoy=remove_decoy, dec_prefix=settings.get('input', 'decoy prefix'))
            all_prots = set()
            for v in prots.itervalues():
                all_prots.update(v['fullgroup'].split(';'))
            for k in prots_full.keys():
                if k not in all_prots:
                    del prots_full[k]
        else:
            prots = filter_evalue_prots(prots, FDR=100.0, remove_decoy=remove_decoy, dec_prefix=settings.get('input', 'decoy prefix'))
            if remove_decoy:
                peptides.filter_decoy()
                for k in prots.keys():
                    if k.startswith(settings.get('input', 'decoy prefix')):
                        del prots[k]
                for k in prots_full.keys():
                    if k.startswith(settings.get('input', 'decoy prefix')):
                        del prots_full[k]

        prots_full, _ = nsaf(prots_full, norm=False)
        prots_full, _ = calc_emPAI(prots_full, protsN, norm=False)
        sumI_norm = sum(x['sumI'] for x in prots.itervalues())
        if not sumI_norm:
            sumI_norm = 1.0
        for k in prots_full.keys():
            prots_full[k]['sumI'] = prots_full[k]['sumI'] / sumI_norm / protsL[k]

        prots, nsaf_norm = nsaf(prots)
        prots, emPAI_norm = calc_emPAI(prots, protsN)
        sumI_norm = sum(x['sumI'] for x in prots.itervalues())
        if not sumI_norm:
            sumI_norm = 1.0
        for k in prots.keys():
            prots[k]['sumI'] = prots[k]['sumI'] / sumI_norm / protsL[k]

        for k in prots_full:
            prots_full[k]['NSAF'] /= nsaf_norm
            prots_full[k]['emPAI'] /= emPAI_norm

        ffolder = path.dirname(path.realpath(curfile))
        if peptides.settings.get('options', 'files') == 'union':
            fname = 'union'
        else:
            fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]

        try:
            fig = plot_quantiation(prots, curfile, peptides.settings, fig, separatefigs=settings.getboolean('advanced options', 'separatefigures'), ox=ox, oy=oy)
            fig = plot_useful_histograms(peptides, curfile, fig, separatefigs=settings.getboolean('advanced options', 'separatefigures'), savesvg=settings.getboolean('advanced options', 'saveSVG'), ox=ox, oy=oy)
        except Exception as e:
            logger.error('Cannot plot quantitation figures: %s', e)
        output_proteins = open('%s/%s_proteins.tsv' % (ffolder, fname), 'w')
        output_proteins.write('dbname\tdescription\tPSMs\tpeptides\tsequence coverage\tLFQ(SIn)\tLFQ(NSAF)\tLFQ(emPAI)\tprotein LN(e-value)\tq-value\tall proteins\n')
        output_proteins_full = open('%s/%s_proteins_full.tsv' % (ffolder, fname), 'w')
        output_proteins_full.write('dbname\tdescription\tPSMs\tpeptides\tsequence coverage\tLFQ(SIn)\tLFQ(NSAF)\tLFQ(emPAI)\tprotein LN(e-value)\tall proteins\n')
        output_PSMs = open('%s/%s_PSMs.tsv' % (ffolder, fname), 'w')
        output_PSMs_pepxml = '%s/%s_PSMs.pep.xml' % (ffolder, fname)

        tags = peptides.peptideslist[0].tags.keys() if peptides.peptideslist[0].tags else None
        tags_header = ('\t' + '\t'.join(tags)) if tags else ''
        output_PSMs.write('sequence\tmodified_sequence\tmodifications\tm/z exp\tcharge\tm/z error in ppm\tmissed cleavages\tnum tol term\tprev_aa\tnext_aa\te-value\tMPscore\tRT exp\tspectrum\tproteins\tproteins description\tSIn\tmassdiff\tis decoy%s' % (tags_header, ))
        output_peptides_detailed = open('%s/%s_peptides.tsv' % (ffolder, fname), 'w')
        output_peptides_detailed.write('sequence\tPSM count\tmodified_sequence\tmodifications\tm/z exp\tcharge\tm/z error in ppm\tmissed cleavages\tnum tol term\tprev_aa\tnext_aa\te-value\tMPscore\tRT exp\tspectrum\tproteins\tproteins description\tSIn\tmassdiff\tis decoy%s' % (tags_header, ))
        framents_info = settings.getboolean('advanced options', 'fragments_info')
        framents_info_zeroes = settings.getboolean('advanced options', 'fragments_info_zeros')
        if framents_info:
            for itype in peptides.peptideslist[0].fragments:
                output_peptides_detailed.write('\t%s_ions' % (itype))
                output_PSMs.write('\t%s_ions' % (itype))
            output_peptides_detailed.write('\tb_sequence\ty_sequence')
            output_PSMs.write('\tb_sequence\ty_sequence')
        output_PSMs.write('\n')
        output_peptides_detailed.write('\n')
        pickle.dump(peptides.RC, open('%s/%s_RC.pickle' % (ffolder, fname), 'w'))
        if protsC:
            output_proteins_valid = open('%s/%s_proteins_valid.tsv' % (ffolder, fname), 'w')
            temp_data = []

        for k, v in prots_full.iteritems():
            sq_tmp = calc_sq(protsS.get(k, []), v['pept'])
            prots_full[k]['sq'] = sq_tmp
            if k in prots:
                prots[k]['sq'] = sq_tmp

        for k, v in prots.items():
            if protsC and k in valid_proteins:
                output_proteins_valid.write('%s,%s,%s,%s,%s\n' % (k, v['PSMs'], v['Peptides'], v['sumI'], protsC[k]))
                temp_data.append([float(v['sumI']), protsC[k]])
            if int(v['Peptides']) > 0:
                sqc = v['sq']#calc_sq(protsS.get(k, []), v['pept'])
                output_proteins.write('%s\t%s\t%s\t%s\t%0.1f\t%0.2E\t%0.2E\t%0.2E\t%0.2E\t%0.4f\t%s\n' % (k, v['description'], v['PSMs'], v['Peptides'], sqc, v['sumI'],v['NSAF'], v['emPAI'], v['expect'], v['qval'], v['fullgroup']))

        for k, v in prots_full.items():
            if int(v['Peptides']) > 0:
                sqc = v['sq']#calc_sq(protsS.get(k, []), v['pept'])
                output_proteins_full.write('%s\t%s\t%s\t%s\t%0.1f\t%0.2E\t%0.2E\t%0.2E\t%0.2E\t%s\n' % (k, v['description'], v['PSMs'], v['Peptides'], sqc, v['sumI'],v['NSAF'], v['emPAI'], v['expect'], v['fullgroup']))

        for val in peptides.get_izip_full():
            if any(protein.dbname in prots for protein in peptides.proteins_dict[val[0].sequence]):#peptide.parentproteins):
                output_PSMs.write(get_output_string(val[0], val[1], val[2], type='psm', fragments_info=framents_info, fragments_info_zeros=framents_info_zeroes, proteins_dict=peptides.proteins_dict, tags=tags))
        peptides_best = dict()
        peptides_best_sp = dict()
        peptides_count = Counter()
        for peptide, spectrum in izip(peptides.peptideslist, peptides.spectrumlist):
            if peptide.evalue < peptides_best.get(peptide.sequence, np.inf):
                peptides_best[peptide.sequence] = peptide.evalue
                peptides_best_sp[peptide.sequence] = spectrum
            peptides_count[peptide.sequence] += 1
        for val in peptides.get_izip_full():
            if val[1] == peptides_best_sp.get(val[0].sequence, None):
                if any(protein.dbname in prots for protein in peptides.proteins_dict[val[0].sequence]):#peptide.parentproteins):
                    output_peptides_detailed.write(get_output_string(val[0], val[1], val[2], type='psm', fragments_info=framents_info, fragments_info_zeros=framents_info_zeroes, peptide_count=peptides_count[val[0].sequence], proteins_dict=peptides.proteins_dict, tags=tags))
                    del peptides_best_sp[val[0].sequence]
        output_peptides_detailed.close()
        output_PSMs.close()

        try:
            filt_pepxml = settings.getint('advanced options', 'add_filtered_pepxml')
        except:
            filt_pepxml = 1
        if filt_pepxml:
            pepxmltk.easy_write_pepxml([curfile], output_PSMs_pepxml, {val[1] for val in peptides.get_izip_full()})
        if protsC:
            temp_sum = sum([x[0] for x in temp_data])
            temp_data = [[x[0] / temp_sum, x[1]] for x in temp_data]
            logger.debug('conc: %s', auxiliary.linear_regression([x[0] for x in temp_data], [x[1] for x in temp_data]))

        output_proteins.close()
        output_proteins_full.close()
    if printresults:
        logger.info('PSMs: %s', sum(x.note2 == 'tr' for x in peptides.peptideslist))
        logger.info('Peptides: %s', len(set(p.sequence for p in peptides.peptideslist if p.note2 == 'tr')))
        logger.info('Protein groups: %s', sum(1 for k in prots if not k.startswith(settings.get('input', 'decoy prefix'))))
        logger.info('Protein groups with >= 2 peptides: %s', sum(1 for k, v in prots.iteritems() if v['Peptides'] >= 2 and not k.startswith(settings.get('input', 'decoy prefix'))))
        if valid_proteins:
            logger.info('PSMs_true: %s', sum(1 for x in peptides.peptideslist if x.note3))
            logger.info('Peptides_true: %s', len(set(x.sequence for x in peptides.peptideslist if x.note3)))
            logger.info('Protein groups_true: %s', len(true_prots))
            logger.info('Real FDR = %s', (100 * float(sum(1 for x in peptides.peptideslist if not x.note3)) / len(peptides.peptideslist)))
    return sum(x.note2 == 'tr' for x in peptides.peptideslist), len(set(p.sequence for p in peptides.peptideslist)), sum(v['Peptides'] > 1 for v in prots.values())

def plot_useful_histograms(peptides, curfile, fig, separatefigs=False, savesvg=False, ox=False, oy=False):
    formulas = [
        (lambda peptide: peptide.RT_exp, 'RT experimental', 'PSMs, RT experimental, min'),
        (lambda peptide: peptide.mz, 'precursor mass', 'PSMs, precursor m/z'),
        (lambda peptide: len(peptide.sequence), 'peptide length', 'PSMs, peptide length'),
        (lambda peptide: peptide.RT_exp, 'RT experimental, peptides', 'peptides, RT experimental, min'),
        (lambda peptide: peptide.mz, 'precursor mass, peptides', 'peptides, precursor m/z'),
        (lambda peptide: len(peptide.sequence), 'peptide length, peptides', 'peptides, peptide length')
    ]
    tmp_dict = {}
    for peptide in peptides.peptideslist:
        if peptide.note2 == 'tr':
            tmp_dict[peptide.sequence] = min(tmp_dict.get(peptide.sequence, 1e6), peptide.evalue)

    for idx, form in enumerate(formulas):
        if form[1].startswith('RT experimental'):
            if form[1].endswith(', peptides'):
                decarray = np.array([peptide.evalue != tmp_dict.get(peptide.sequence, None) for peptide in peptides.peptideslist])
                array_valid = form[0](peptides)[~decarray]
            else:
                decarray = peptides.is_decoy_array()
                array_valid = form[0](peptides)[~decarray]
        else:
            if form[1].endswith(', peptides'):
                array_valid = [form[0](peptide) for peptide in peptides.peptideslist if peptide.evalue == tmp_dict.get(peptide.sequence, None)]
            else:
                array_valid = [form[0](peptide) for peptide in peptides.peptideslist if peptide.note2 == 'tr']
        if separatefigs:
            plt.clf()
            fig = plt.figure()
            DPI = fig.get_dpi()
            fig.set_size_inches(500.0/float(DPI), 400.0/float(DPI))
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.add_subplot(ox, oy, idx+1)
        lbin = min(array_valid)
        rbin = max(array_valid)
        if form[1].startswith('peptide length'):
            binsize = 1
            lbin = lbin-0.5
        else: binsize = FDbinSize(array_valid)
        if lbin != rbin:
            H1, _ = np.histogram(array_valid, bins=np.arange(lbin, rbin+binsize, binsize))
            ind = np.arange(lbin, rbin, binsize)
            width = binsize
            ax.bar(ind, H1, width, align='edge', color=redcolor, alpha=0.8,edgecolor='#EEEEEE')
            ax.set_ylabel('# of identifications')
            ax.set_xlabel(form[2])
            plt.grid(color='#EEEEEE')

            from matplotlib.ticker import MaxNLocator
            ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6))


            if separatefigs:
                plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)

                if peptides.settings.get('options', 'files') == 'union':
                    fname = 'union'
                else:
                    fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
                # tmpfname = '%s/%s_%s' % (path.dirname(path.realpath(curfile)), fname, form[1])
                tmpfname = '%s_%s' % (fname, form[1])
                tmpfname = tmpfname.replace(' ', '_').replace(',', '')
                tmpfname = path.join(path.dirname(path.realpath(curfile)), tmpfname)
                plt.grid(color='#EEEEEE')
        #        seaborn.despine()
                plt.savefig(tmpfname + '.png')
                if savesvg:
                    plt.savefig(tmpfname + '.svg')
                plt.close()
    del tmp_dict
    return fig


def plot_histograms(descriptors, peptides, FDR, curfile, savesvg=False, sepfigures=False):
    if not sepfigures:
        fig = plt.figure(figsize=(16, 12))
        ox, oy = find_optimal_xy(descriptors)
        DPI = fig.get_dpi()
        fig.set_size_inches(2000.0/float(DPI), 2000.0/float(DPI))
    copy_peptides, _, _ = peptides.filter_evalue_new(FDR=FDR, useMP=False)

    for idx, descriptor in enumerate(descriptors):
        if sepfigures:
            plt.clf()
            fig = plt.figure()
            DPI = fig.get_dpi()
            fig.set_size_inches(500.0/float(DPI), 400.0/float(DPI))
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.add_subplot(ox, oy, idx + 10)
        darray = np.array(descriptor.formula(peptides))
        dnotes = peptides.is_decoy_array()
        array_wrong = darray[dnotes==True]
        array_valid = darray[dnotes==False]

        if descriptor.group == 'B':
            array_wrong = np.log10(array_wrong)
            array_valid = np.log10(array_valid)
        binsize = float(descriptor.get_binsize(copy_peptides))
        if binsize < float(max(np.append(array_wrong, array_valid)) - min(np.append(array_wrong, array_valid))) / 300:
            binsize = float(max(np.append(array_wrong, array_valid)) - min(np.append(array_wrong, array_valid))) / 300
        lbin_s = scoreatpercentile(np.append(array_wrong, array_valid), 1.0)
        lbin = min(np.append(array_wrong, array_valid))
        if lbin_s and abs((lbin - lbin_s) / lbin_s) > 1.0:
            lbin = lbin_s * 1.05
        rbin_s = scoreatpercentile(np.append(array_wrong, array_valid), 99.0)
        rbin = max(np.append(array_wrong, array_valid))
        if rbin_s and abs((rbin - rbin_s) / rbin_s) > 1.0:
            rbin = rbin_s * 1.05
        rbin += 1.5 * binsize
        if descriptor.name.startswith('potential modifications'):
            lbin = -1
            rbin = 2.5
            binsize = 1
        H1, _ = np.histogram(array_wrong, bins=np.arange(lbin, rbin+binsize, binsize))
        H2, bins_valid = np.histogram(array_valid, bins=np.arange(lbin, rbin+binsize, binsize))
        if descriptor.group == 'B':
            H3, _ = np.histogram(np.log10(descriptor.formula(copy_peptides)), bins=np.arange(lbin, rbin+binsize, binsize))
        else:
            H3, _ = np.histogram(np.array(descriptor.formula(copy_peptides)), bins=np.arange(lbin, rbin+binsize, binsize))

        if descriptor.name == 'precursor mass difference, ppm':
            mass_arr = np.array(descriptor.formula(copy_peptides))
            median_mass = np.median(mass_arr)
            logger.info('MAX BIN precursor mass difference of top PSMs = %s ppm', bins_valid[:-1][H2 == H2.max()])
            logger.info('STD precursor mass difference of top PSMs = %s ppm', np.std(mass_arr - median_mass))
        if descriptor.group == 'B':
            H1 = H1.clip(1)
            H2 = H2.clip(1)
            H3 = H3.clip(1)
            H1 = np.log10(H1)
            H2 = np.log10(H2)
            H3 = np.log10(H3)
        ind = np.arange(lbin, rbin+binsize, binsize)
        width = binsize
        if descriptor.name.startswith('potential modifications'):
            ind=np.append(-2,ind)
            H1=np.append([H1[0],0],H1[1:])
            H2=np.append([H2[0],0],H2[1:])
            H3=np.append([H3[0],0],H3[1:])
      #  step_ind = ind[:-1]+float(width)/2
        
        if len(ind)>50: 
            ax.bar(ind[:-1], H1, width, align='edge',color=redcolor, alpha=0.4, edgecolor=redcolor)
            ax.bar(ind[:-1], H2, width, align='edge',color=bluecolor, alpha=0.4, edgecolor=bluecolor)
            ax.bar(ind[:-1], H3, width, align='edge',color=greencolor, alpha=1, edgecolor=greencolor)
            ind=np.append(ind[0],ind)
            H1=np.append(np.append(0,H1),0)
            H2=np.append(np.append(0,H2),0)
            
            ax.step(ind, H2, where='post', color=bluecolor,alpha=0.8)
            ax.step(ind, H1, where='post', color=redcolor,alpha=0.8)
        else:
            ax.bar(ind[:-1], H1, width, align='edge',color=redcolor, alpha=0.4,edgecolor='#EEEEEE')
            ax.bar(ind[:-1], H2, width, align='edge',color=bluecolor, alpha=0.4,edgecolor='#EEEEEE')
            ax.bar(ind[:-1], H3, width, align='edge',color=greencolor, alpha=1,edgecolor='#EEEEEE')
            ind=np.append(ind[0],ind)
            H1=np.append(np.append(0,H1),0)
            H2=np.append(np.append(0,H2),0)
            ax.step(ind, H2, where='post', color=bluecolor,alpha=0.8)
            ax.step(ind, H1, where='post', color=redcolor,alpha=0.8)
        if any(descriptor.name.startswith(clabel) for clabel in ['missed cleavages', 'charge states', 'isotopes mass difference, Da']):
            logger.debug('%s %s', descriptor.name, ind)
            ax.set_xticks(np.arange(0.5, 5.5, 1.0))
            fig.canvas.draw()
            labels = [item.get_text() for item in ax.get_xticklabels()]
            ax.set_xticklabels([int(float(l)) for l in labels])
        elif descriptor.name.startswith('potential modifications'):
            logger.debug('%s %s', descriptor.name, ind)
            ax.set_xticks(np.arange(-1.5,2.5,1))
            ax.set_xticklabels(['NA','']+list(np.arange(0, 3, 1)))
            ax.set_xlim(-2, 2.5)
       # else:
        #    from matplotlib.ticker import MaxNLocator
         #   ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6))
        
            

        if descriptor.name.startswith('PSM count'): ax.set_xlabel('log(PSMs per peptide)')
        elif descriptor.name.startswith('PSMs per protein'): ax.set_xlabel('log(PSMs per protein)')
        else: ax.set_xlabel(descriptor.name)
        if peptides.settings.get('options', 'files') == 'union':
            fname = 'union'
        else:
            fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
        if descriptor.group == 'A':
            ax.set_ylabel('# of identifications')
        else:
            ax.set_ylabel('Log(# of identifications)')
        if sepfigures:
            plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)
            # tmpfname = '%s/%s_%s' % (path.dirname(path.realpath(curfile)), fname, descriptor.name)
            tmpfname = '%s_%s' % (fname, descriptor.name)
            tmpfname = tmpfname.replace(' ', '_').replace(',', '')
            tmpfname = path.join(path.dirname(path.realpath(curfile)), tmpfname)
            plt.grid(color="#EEEEEE")
         #   seaborn.despine()
            plt.savefig(tmpfname + '.png' )
            if savesvg:
                plt.savefig(tmpfname + '.svg')
            plt.close()
    return fig


def plot_quantiation(prots, curfile, settings, fig, separatefigs, ox, oy):
    for i, idx in enumerate(['sumI', 'NSAF', 'emPAI']):
        dots = [np.log10(v[idx]) for v in prots.itervalues()]
        if separatefigs:
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            DPI = fig.get_dpi()
            fig.set_size_inches(500.0/float(DPI), 400.0/float(DPI))
        else:
            ax = fig.add_subplot(ox, oy, i+7)
        ax.hist(dots, bins = 10, alpha=0.8, color=bluecolor,edgecolor='#EEEEEE')
        plt.grid(color="#EEEEEE")
        #seaborn.despine()
        ax.set_xlabel('LOG10(%s)' % (idx, ))
        ax.set_ylabel('Number of proteins')
        ax.locator_params(axis='x', nbins=4)
        if settings.get('options', 'files') == 'union':
            fname = 'union'
        else:
            fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
        if separatefigs:
            plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)
            # tmpfname = '%s/%s_%s' % (path.dirname(path.realpath(curfile)), fname, idx)
            tmpfname = '%s_%s' % (fname, idx)
            tmpfname = tmpfname.replace(' ', '_').replace(',', '')
            tmpfname = path.join(path.dirname(path.realpath(curfile)), tmpfname)
            plt.savefig(tmpfname + '.png')
            if settings.getboolean('advanced options', 'saveSVG'):
                plt.savefig(tmpfname + '.svg')
            plt.close()
    return fig


def plot_MP(descriptors, peptides, fig, FDR, FDR2, valid_proteins, settings, threshold0=False, curfile=False):
    ox, oy = find_optimal_xy(descriptors)
    sepfigures = settings.getboolean('advanced options', 'separatefigures')
    copy_peptides, threshold1, threshold2 = peptides.filter_evalue_new(FDR=FDR, FDR2=FDR2, useMP=True, drop_decoy=False)

    threshold1 = -np.log(threshold1)
    try:
        threshold2 = np.log(threshold2)
    except:
        pass

    zero_peptscore = np.log(min([p.peptscore for p in peptides.peptideslist if p.peptscore != 0])) - 1
    zero_evalue = -np.log(min([p.evalue for p in peptides.peptideslist if p.evalue != 0])) + 1
    PSMs_wrong = [[(-np.log(pept.evalue) if pept.evalue != 0 else zero_evalue), (np.log(pept.peptscore) if pept.peptscore != 0 else zero_peptscore)] for pept in peptides.peptideslist if pept.note2 == 'wr']
    PSMs_true = [[(-np.log(pept.evalue) if pept.evalue != 0 else zero_evalue), (np.log(pept.peptscore) if pept.peptscore != 0 else zero_peptscore)] for pept in peptides.peptideslist if pept.note2 == 'tr']

    logger.info('MP filtering:')
    PSMs_info(copy_peptides, valid_proteins, settings, fig=fig, tofile=True, curfile=curfile, ox=ox, oy=oy)
    ffolder = path.dirname(path.realpath(curfile))
    if peptides.settings.get('options', 'files') == 'union':
        fname = 'union'
    else:
        fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
    output_PSMs_full = open('%s/%s_PSMs_full.tsv' % (ffolder, fname), 'w')
    tags = peptides.peptideslist[0].tags.keys() if peptides.peptideslist[0].tags else None
    tags_header = ('\t' + '\t'.join(tags)) if tags else ''
    output_PSMs_full.write('sequence\tmodified_sequence\tmodifications\tm/z exp\tcharge\tm/z error in ppm\tmissed cleavages\tnum tol term\tprev_aa\tnext_aa\te-value\tMPscore\tRT exp\tspectrum\tproteins\tproteins description\tSIn\tmassdiff\tis decoy%s\n' % (tags_header, ))
    for val in peptides.get_izip_full():
        output_PSMs_full.write(get_output_string(val[0], val[1], val[2], type='psm', fragments_info=False, fragments_info_zeros=False, proteins_dict=peptides.proteins_dict, tags=tags))
    logger.info('Without filtering, after removing outliers:')
    PSMs_info(peptides, valid_proteins, settings, loop=False)

    all_seqs = set([p.sequence for p in copy_peptides.peptideslist])
    std_aa_list = list(mass.std_aa_mass.keys())
    std_aa_list.remove('O')
    std_aa_list.remove('U')
    aa_exp = Counter()
    for pep in all_seqs:
        for aa in pep:
            aa_exp[aa] += 1

    aa_theor = Counter()
    for pep in all_seqs:
        for prot_obj in peptides.proteins_dict[pep]:
            prot_name = prot_obj.dbname
            prot_seq = protsS[prot_name]
            for aa in prot_seq:
                aa_theor[aa] += 1
            
    aa_exp_sum = sum(aa_exp.values())
    aa_theor_sum = sum(aa_theor.values())
    lbls, vals = [], []
    for aa in sorted(std_aa_list):
        if aa in aa_theor:
            lbls.append(aa)
            vals.append((float(aa_exp.get(aa, 0))/aa_exp_sum)/(float(aa_theor.get(aa, 0))/aa_theor_sum))
        else:
            logger.warning('Amino acid %s is missing in theor sequences' % (aa, ))
    std_val = np.std(vals)
    clrs = [greencolor if abs(x-1)<=2*std_val else redcolor for x in vals]

    if sepfigures:
        plt.clf()
        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(500.0/float(DPI), 400.0/float(DPI))
        ax = fig.add_subplot(1, 1, 1)
    else:
        ax = fig.add_subplot(ox, oy, ox*oy-1)

    ax.bar(range(len(vals)), vals, color=clrs)
    ax.set_xticks(range(len(lbls)))
    ax.set_xticklabels(lbls)
    ax.hlines(1.0, range(len(vals))[0]-1, range(len(vals))[-1]+1)
    ax.set_ylabel('amino acid ID rate')

    if sepfigures:
        plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)
        # tmpfname = '%s/%s_%s' % (path.dirname(path.realpath(curfile)), fname, 'scores')
        tmpfname = '%s_%s' % (fname, 'aa_stats')
        tmpfname = tmpfname.replace(' ', '_').replace(',', '')
        tmpfname = path.join(path.dirname(path.realpath(curfile)), tmpfname)
        plt.savefig(tmpfname + '.png')
        if settings.getboolean('advanced options', 'saveSVG'):
            plt.savefig(tmpfname + '.svg')
        plt.close()


    if sepfigures:
        plt.clf()
        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(500.0/float(DPI), 400.0/float(DPI))
        ax = fig.add_subplot(1, 1, 1)
    else:
        ax = fig.add_subplot(ox, oy, ox*oy)
    col=[redcolor,]*len(PSMs_wrong)+[bluecolor,]*len(PSMs_true)
    PSMs_all = PSMs_wrong+PSMs_true
    idx_range = range(len(col))
    np.random.shuffle(idx_range)
    x=[]
    y=[]
    ncol=[]
    for a in idx_range:
        x.append(PSMs_all[a][0])
        y.append(PSMs_all[a][1])
        ncol.append(col[a])
    col = ncol
    plt.scatter(x,y, s=2, color=col)
    ax.axvline(threshold1, color=greencolor)
    if threshold2:
        ax.axhline(threshold2, color=greencolor)
    if threshold0:
        threshold0 = -np.log(threshold0)
        ax.axvline(threshold0, color=redcolor)
    ax.set_ylim(min(y) - 1, max(y) + 1)
    ax.set_xlim(min(x) - 1, max(x) + 1)
    if peptides.settings.get('options', 'files') == 'union':
        fname = 'union'
    else:
        fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
    ax.set_xlabel('-LOG(evalue)')
    ax.set_ylabel('LOG(MPscore)')
    plt.grid(color='#EEEEEE')
  #  seaborn.despine()
    if sepfigures:
        plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)
        # tmpfname = '%s/%s_%s' % (path.dirname(path.realpath(curfile)), fname, 'scores')
        tmpfname = '%s_%s' % (fname, 'scores')
        tmpfname = tmpfname.replace(' ', '_').replace(',', '')
        tmpfname = path.join(path.dirname(path.realpath(curfile)), tmpfname)
        plt.savefig(tmpfname + '.png')
        if settings.getboolean('advanced options', 'saveSVG'):
            plt.savefig(tmpfname + '.svg')
        plt.close()
    else:
        plt.gcf().subplots_adjust(bottom=0.05, left=0.05, top=0.95, right=0.95)
        # tmpfname = '%s/%s' % (path.dirname(path.realpath(curfile)), fname)
        tmpfname = str(fname)
        tmpfname = tmpfname.replace(' ', '_').replace(',', '')
        tmpfname = path.join(path.dirname(path.realpath(curfile)), tmpfname)
        plt.tight_layout()
        plt.savefig(tmpfname + '.png')
        if settings.getboolean('advanced options', 'saveSVG'):
            plt.savefig(tmpfname + '.svg')
        plt.close()


def prepare_hist(descriptors, copy_peptides, first=False):
    for descriptor in descriptors:
        descriptor.get_array(copy_peptides)
        if descriptor.group == 'A':
            binsize = float(descriptor.get_binsize(copy_peptides))
            if binsize < float(max(descriptor.array) - min(descriptor.array)) / 400:
                binsize = float(max(descriptor.array) - min(descriptor.array)) / 400

            lbin, rbin = min(descriptor.array), max(descriptor.array) + 1.5 * binsize
            if lbin == rbin:
                rbin = lbin + binsize
            descriptor.hist = np.histogram(descriptor.array, bins=np.arange(lbin, rbin + binsize, binsize))

            if descriptor.name.startswith('potential modifications'):
                descriptor.hist[0][0] = max(descriptor.hist[0][1:])
                if not descriptor.hist[0][0]:
                    descriptor.hist[0][0] = 1
            del descriptor.array
            descriptor.hsum = sum(descriptor.hist[0])
    return descriptors


def calc_peptscore(cq, descriptors):
    desc_dict = {}
    for descriptor in descriptors:
        desc_dict[descriptor.name] = np.array(descriptor.formula(cq))
    for idx, peptide in enumerate(cq.peptideslist):
        tmp_peptscore = peptide.peptscore
        for descriptor in descriptors:
            # descriptor_value = descriptor.single_formula(peptide)
            descriptor_value = desc_dict[descriptor.name][idx]
            if descriptor.group == 'A':
                if descriptor_value < descriptor.hist[1][0] or descriptor_value >= descriptor.hist[1][-1]:
                    tmp_peptscore = 0
                else:
                    j = descriptor.hist[1].searchsorted(descriptor_value)
                    if descriptor_value < descriptor.hist[1][j]:
                        j -= 1
                    try:
                        tmp_peptscore *= float(descriptor.hist[0][j]) / descriptor.hsum#sum(descriptor.hist[0])
                    except:
                        logger.error('Error when multiplying by descriptor %s', descriptor.name)
                        logger.error('%s', descriptor.hist[0])
                        logger.error('%s', descriptor.hist[1])
                        logger.error('%s %s %s %s', descriptor_value, descriptor.hist[1][j], descriptor.hist[1][descriptor.hist[1].searchsorted(descriptor_value)], descriptor.hist[1].searchsorted(descriptor_value))
            elif descriptor.group == 'B':
                tmp_peptscore *= float(descriptor.array.searchsorted(descriptor_value, side='right')) / descriptor.array.size

        yield (idx, tmp_peptscore)

def main(argv_in, union_custom=False):
    # inputfile = argv_in[1]
    files = {}
    fastafile = None
    configfile = None

    def update_dict(inputdict, path_to_file=None):
        if path_to_file:
            extension = path.splitext(path_to_file)[-1]
            if extension == '.xml':
                extension = path.splitext(path.splitext(path_to_file)[0])[-1]
                filename = path.basename(path.splitext(path.splitext(path_to_file)[0])[0])
            else:
                filename = path.basename(path.splitext(path_to_file)[0])
            inputdict.setdefault(filename, {})[extension] = path_to_file
        else:
            for k, v in inputdict.items():
                if '.pep' not in v:
                    del inputdict[k]
                else:
                    for ext in ('.mgf', '.mzML'):
                        if ext not in v:
                            path_to_file = path.join(path.dirname(v['.pep']), k) + (ext if not ext == '.t' else ext + '.xml')
                            if path.isfile(path_to_file):
                                inputdict[k][ext] = path_to_file
        return inputdict

    for arg in argv_in:
        if path.splitext(arg)[-1].lower() in {'.fasta', '.faa'}:
            fastafile = arg
        elif path.splitext(arg)[-1] == '.cfg':
            configfile = arg
        elif path.isdir(arg):
            for filename in listdir(arg):
                files = update_dict(files, path.join(arg, filename))
        else:
            files = update_dict(files, arg)
    files = update_dict(files)

    if configfile:
        defpath = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'default.cfg')
        logger.debug('Loading defaults from %s ...', defpath)
        settings = get_settings(defpath)
        logger.debug('Loading config from %s ...', configfile)
        settings.read(configfile)
        if union_custom:
            settings.set('options', 'files', 'union')

        if settings.get('advanced options', 'copy_params_to_output_folder') == 'yes':
            try:
                shutil.copy(configfile, os.path.join(os.path.dirname(os.path.abspath(files.values()[0]['.pep'])), os.path.basename(configfile)))
            except:
                pass

        proteases = [x.strip() for x in settings.get('missed cleavages', 'protease1').split(',')]
        proteases.extend([x.strip() for x in settings.get('missed cleavages', 'protease2').split(',') if x.strip()])
        expasy = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules else protease for protease in proteases))
        try:
            mc = settings.getint('missed cleavages', 'number of missed cleavages')
        except Exception as e:
            logger.warning('Number of missed cleavages is missing in parameters, using 2 for normalization and emPAI calculation')
            logger.debug(e)
            mc = 2
        # fprocs = []
        # fnprocs = 12
        # fq = multiprocessing.Queue()
        # fq_output = multiprocessing.Queue()

        protsL['total proteins'] = 0
        protsL['total peptides'] = 0
        def get_number_of_peptides(prot, expasy, mc, minl, maxl):
            return sum(minl <= len(x) <= maxl for x in parser.cleave(prot, expasy, mc))
        def protein_handle(fq, fq_output, protsL, protsN, protsS, expasy, mc, minl, maxl):
            while 1:
                try:
                    x = fq.get(timeout=1)
                except Empty:
                    fq_output.put('1')
                    break

                try:
                    dbname = x[0].split(' ')[0]
                    protsS[dbname] = x[1]
                    protsL[dbname] = len(x[1])
                    protsN[dbname] = get_number_of_peptides(x[1], expasy, mc, minl, maxl)
                    protsL['total proteins'] += 1
                    protsL['total peptides'] += protsN.get(dbname, 0)
                except Exception as e:
                    logger.error('Error reading of FASTA file: %s', e)

        if fastafile:
            if settings.get('input', 'add decoy') == 'yes':
                decoy_method = settings.get('input', 'decoy method')
                decoy_prefix = settings.get('input', 'decoy prefix')
                fastagen = fasta.decoy_db(fastafile, mode=decoy_method, prefix=decoy_prefix)
                # for x in fasta.decoy_db(fastafile, mode=decoy_method, prefix=decoy_prefix):
                #     fq.put(x)
            else:
                fastagen = fasta.read(fastafile)
                # for x in fasta.read(fastafile):
                #     fq.put(x)
        else:
            logger.critical('FASTA file not specified.')
            return 1

        minl = settings.getint('search', 'peptide minimum length')
        maxl = settings.getint('search', 'peptide maximum length')

        for x in fastagen:
            dbname = x[0].split(' ')[0]
            protsS[dbname] = x[1]
            protsL[dbname] = len(x[1])
            protsN[dbname] = get_number_of_peptides(x[1], expasy, mc, minl, maxl)
            protsL['total proteins'] += 1
            protsL['total peptides'] += protsN.get(dbname, 0)

        # for i in range(fnprocs):
        #     p = multiprocessing.Process(target=protein_handle, args=(fq, fq_output, protsL, protsN, protsS, expasy, mc, minl, maxl))
        #     fprocs.append(p)
        #     p.start()

        # while fq_output.qsize() != fnprocs:
        #     sleep(10)
        # for p in fprocs:
        #     p.terminate()

        procs = []
        nprocs = 1
        q = multiprocessing.Queue()
        q_output = multiprocessing.Queue()

        files_processing = settings.get('options', 'files')
        if files_processing == 'union':
            q.put(sorted(files.values()))
        else:
            for filename in sorted(files.values()):
                q.put([filename, ])
        for i in range(nprocs):
            p = multiprocessing.Process(target=handle, args=(q, q_output, settings, protsL))
            procs.append(p)
            p.start()

        while q_output.qsize() != nprocs:
            sleep(10)
        for p in procs:
            p.terminate()
    else:
        logger.critical('.cfg file with parameters is missing')

if __name__ == '__main__':
    LOGGING = {
        'version': 1,
        'disable_existing_loggers': True,
        'formatters': {
            'simple': {
                'format': '%(levelname)8s: %(asctime)s %(message)s',
                'datefmt': '[%H:%M:%S]',
            },
        },
        'handlers': {
            'console': {
                'level': 'DEBUG',
                'class': 'logging.StreamHandler',
                'formatter': 'simple',
            },
        },
        'loggers': {
            '__main__': {
                'handlers': ['console'],
                'level': 'DEBUG',
            },
            'MPlib': {
                'handlers': ['console'],
                'level': 'DEBUG',
            }
        }
    }

    logging.config.dictConfig(LOGGING)
    main(sys.argv)
