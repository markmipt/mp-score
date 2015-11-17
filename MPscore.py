from MPlib import PeptideList, Descriptor, get_settings, filter_evalue_prots, FDbinSize
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
from os import path, listdir
from pyteomics import mzml, fasta, auxiliary, mgf, parser
from Queue import Empty
import multiprocessing
from time import sleep, time
import pickle
from copy import copy, deepcopy
from collections import defaultdict
try:
    import seaborn
    seaborn.set(rc={'axes.facecolor':'#eeeeff'})
except:
    pass

manager = multiprocessing.Manager()
protsC = {}
protsL = manager.dict()
protsS = manager.dict()
protsN = manager.dict()
stime = time()

def calc_sq(protein, peptides):
    if not protein:
        return 0
    psq = [False for x in protein]
    for pep in peptides:
        csize = len(pep)
        for j in range(len(protein)):
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
        print 'inputfile = %s' % (','.join(f['.pep'] for f in filenames), )
        FDR = settings.getfloat('options', 'FDR')
        FDR_type = settings.get('options', 'FDR_type')

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

        Fragment_intensities = {}
        spectra_dict = dict()
        spectra_dict_intensities = dict()

        peptides = PeptideList(settings)
        mods = manager.dict()
        iprocs = []
        inprocs = peptides.settings.getint('options', 'threads')
        iq = multiprocessing.Queue()
        iq_output = multiprocessing.Queue()

        def getpepxml(iq, iq_output, settings, mods=False):
            for curfile in iter(iq.get, None):
                qpeptides = PeptideList(settings, mods)
                qpeptides.get_from_pepxmlfile(curfile['.pep'], min_charge=min_charge, max_charge=max_charge, allowed_peptides=settings.get('advanced options', 'allowed peptides'))

                if len(qpeptides.peptideslist):
                    mzmlfile = curfile.get('.mzML', None)
                    if mzmlfile:
                        print 'mzml is processing'
                        # isolation_window = settings.getfloat('precursor ion fraction', 'isolation window')
                        # mass_acc = settings.getfloat('precursor ion fraction', 'mass accuracy')
                        # spectra_ms1 = []
                        spectra_ms2 = []
                        for spectrum in mzml.read(mzmlfile):
                            if spectrum['ms level'] == 2:
                                spectra_dict[spectrum['spectrum title'].strip()] = spectrum['m/z array']
                                spectra_dict_intensities[spectrum['spectrum title'].strip()] = spectrum['intensity array']
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

                    mgffile = curfile.get('.mgf', None)
                    if mgffile:
                        print 'mgf is processing'
                        spectra = mgf.read(mgffile)
                        for spectrum in spectra:
                            spectra_dict[spectrum['params']['title'].strip()] = spectrum['m/z array']
                            spectra_dict_intensities[spectrum['params']['title'].strip()] = spectrum['intensity array']
                        protsL['total spectra'] = len(spectra_dict)
                        if not qpeptides.total_number_of_spectra:
                            qpeptides.total_number_of_spectra = len(spectra_dict)
                    if not spectra_dict:
                        qpeptides.settings.set('descriptors', 'fragment mass tolerance, Da', '0')
                        print 'fragment mass tolerance was turned off due to missed mgf file'
                    if qpeptides.settings.getboolean('descriptors', 'fragment mass tolerance, Da'):
                        for peptide in qpeptides.peptideslist:
                            try:
                                peptide.spectrum_mz = spectra_dict[peptide.spectrum.split(' RTINSECONDS=')[0].strip()]
                                peptide.spectrum_i = spectra_dict_intensities[peptide.spectrum.split(' RTINSECONDS=')[0].strip()]
                            except:
                                try:
                                    peptide.spectrum_mz = spectra_dict[peptide.spectrum.strip() + ' min']
                                    peptide.spectrum_i = spectra_dict_intensities[peptide.spectrum.strip() + ' min']
                                except:
                                    peptide.spectrum_mz = spectra_dict[peptide.spectrum.strip()]
                                    peptide.spectrum_i = spectra_dict_intensities[peptide.spectrum.strip()]
                        for peptide in qpeptides.peptideslist:
                            peptide.get_median_fragment_mt(qpeptides.settings)
                            peptide.spectrum_mz = None
                            peptide.spectrum_i = None

                    tmp_peptides = qpeptides.copy_empty()
                    msize = 10000
                    while len(qpeptides.peptideslist):
                        tmp_peptides.peptideslist = qpeptides.peptideslist[:msize]
                        iq_output.put(copy(tmp_peptides))
                        qpeptides.peptideslist = qpeptides.peptideslist[msize:]
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
        peptides.total_number_of_PSMs = len(peptides.peptideslist)

        for p in iprocs:
            p.terminate()

        print 'Total number of PSMs = %d' % (len(peptides.peptideslist),)
        print 'Total number of peptides: %s' % (len(set(pept.sequence for pept in peptides.peptideslist)), )
        if len(peptides.peptideslist):
            prots_dict = defaultdict(int)
            pepts_dict = defaultdict(int)
            for peptide in peptides.peptideslist:
                pepts_dict[peptide.sequence] += 1
                for protein in peptide.parentproteins:
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

            if FDR_type == 'peptide':
                peptidesdict = dict()
                for peptide in peptides.peptideslist:
                    if peptide.sequence not in peptidesdict:
                        peptidesdict[peptide.sequence] = [peptide.spectrum, peptide.evalue]
                    elif peptide.evalue < peptidesdict[peptide.sequence][1]:
                        peptidesdict[peptide.sequence] = [peptide.spectrum, peptide.evalue]
                passed = set([x[0] for x in peptidesdict.values()])
                j = len(peptides.peptideslist) - 1
                while j >= 0:
                    if peptides.peptideslist[j].spectrum not in passed:
                        peptides.peptideslist.pop(j)
                    j -= 1

            for peptide in peptides.peptideslist:
                try:
                    peptide.sumI = Fragment_intensities[peptide.start_scan]
                except:
                    pass
            Fragment_intensities = None

            for peptide in peptides.peptideslist:
                peptide.peptscore2 = pepts_dict[peptide.sequence]
                for protein in peptide.parentproteins:
                    if protein.dbname not in protsL:
                        print 'protein %s is missed in fasta, 5000 length and 50 theoretical peptides is used for normalization and emPAI calculation' % (protein.dbname, )
                        protsL[protein.dbname] = 5000
                        protsN[protein.dbname] = 50
                        protsL['total proteins'] += 1
                        protsL['total peptides'] += 50
                    if peptide.protscore2 < float(prots_dict[protein.dbname]) / protsL[protein.dbname] * 500:
                        peptide.protscore2 = float(prots_dict[protein.dbname]) / protsL[protein.dbname] * 500
            pepts_dict = None
            prots_dict = None
            copy_peptides, threshold0, _ = peptides.filter_evalue_new(FDR=FDR, useMP=False)

            print 'Default filtering:'
            numPSMs, numpeptides_true, numprots_true = PSMs_info(copy_peptides, valid_proteins, settings)
            if numPSMs > 1:
                descriptors = []
                dname = 'RT difference, min'
                if peptides.settings.getboolean('descriptors', dname) and not np.allclose([pep.RT_exp for pep in peptides.peptideslist], 0):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.RT_exp - peptide.RT_predicted, group='A', binsize='auto'))
                dname = 'precursor mass difference, ppm'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.mass_diff(), group='A', binsize='auto'))
                dname = 'missed cleavages, protease 1'
                if peptides.settings.getboolean('descriptors', dname.split(',')[0]):
                    protease1 = [x.strip() for x in settings.get('missed cleavages', 'protease1').split(',')]
                    expasy1 = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules else protease for protease in protease1))
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.get_missed_cleavages(expasy1), group='A', binsize=1))
                dname = 'missed cleavages, protease 2'
                if peptides.settings.getboolean('descriptors', dname.split(',')[0]):
                    protease2 = [x.strip() for x in settings.get('missed cleavages', 'protease2').split(',')]
                    expasy2 = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules else protease for protease in protease2))
                    if protease2[0]:
                        descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.get_missed_cleavages(expasy2), group='A', binsize=1))

                dname = 'charge states'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.pcharge, group='A', binsize='1'))
                dname = 'potential modifications'
                if peptides.settings.getboolean('descriptors', dname):
                    labeldict = dict()
                    temp = settings.get('modifications', 'variable')
                    if temp:
                        for mod in temp.replace(' ', '').split(','):
                            descriptors.append(Descriptor(name='%s, %s' % (dname, mod), formula=lambda peptide, mod=mod: peptide.count_modifications(mod), group='A', binsize='1'))
                dname = 'isotopes mass difference, Da'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: round(peptide.massdiff, 0), group='A', binsize='1'))
                dname = 'PSMs per protein'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.protscore2, group='B'))
                dname = 'PSM count'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.peptscore2, group='B'))
                dname = 'fragment mass tolerance, Da'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.get_median_fragment_mt(), group='A', binsize = 'auto'))
                dname = 'PIF'
                if peptides.settings.getboolean('descriptors', dname):
                    descriptors.append(Descriptor(name=dname, formula=lambda peptide: peptide.PIF, group='B'))

                if 'RT difference, min' in [d.name for d in descriptors]:
                    if RT_type == 'achrom':
                        copy_peptides.get_RC()
                        peptides.RC = copy_peptides.RC
                        if peptides.settings.getint('advanced options', 'saveRC'):
                            pickle.dump(peptides.RC, open('RC.pickle', 'w'))
                        copy_peptides.calc_RT(RTtype=RT_type)
                        print copy_peptides.get_calibrate_coeff()
                        peptides.calc_RT(RTtype=RT_type)
                    else:
                        copy_peptides.filter_modifications(RT_type=RT_type)
                        copy_peptides.calc_RT(RTtype=RT_type)
                        calibrate_coeff = copy_peptides.get_calibrate_coeff()
                        copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
                        copy_peptides.filter_RT(RT_tolerance=3 * calibrate_coeff[3])
                        copy_peptides.calc_RT(RTtype=RT_type)
                        calibrate_coeff = copy_peptides.get_calibrate_coeff()
                        print calibrate_coeff
                        peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)
                        copy_peptides, _, _ = peptides.filter_evalue_new(FDR=FDR, useMP=False)
                        copy_peptides.calc_RT(calibrate_coeff=calibrate_coeff, RTtype=RT_type)

                curfile = filenames[-1]['.pep']
                if descriptors:
                    descriptors = prepare_hist(descriptors, copy_peptides, first=False)
                    fig = plot_histograms(descriptors, peptides, FDR, curfile, savesvg=settings.getboolean('advanced options', 'saveSVG'), sepfigures=settings.getboolean('advanced options', 'separatefigures'))

                    if len(copy_peptides.peptideslist) > 100:
                        jk = manager.dict()
                        cprocs = []
                        cnprocs = peptides.settings.getint('options', 'threads')
                        cq = multiprocessing.Queue()
                        cq_output = multiprocessing.Queue()
                        cq_finish = multiprocessing.Queue()

                        for i in range(cnprocs):
                            p = multiprocessing.Process(target=calc_peptscore, args=(cq, cq_output, descriptors, jk, cq_finish))
                            cprocs.append(p)
                            p.start()

                        counter = 0
                        for idx, peptide in enumerate(peptides.peptideslist):
                            # if counter < 10000:
                            cq.put([idx, peptide])
                            counter += 1
                            if counter >= 10000:
                                while counter != 0:
                                    ind, pscore = cq_output.get()
                                    peptides.peptideslist[ind].peptscore = float(pscore)
                                    counter -= 1
                        while counter != 0:
                            ind, pscore = cq_output.get()
                            peptides.peptideslist[ind].peptscore = float(pscore)
                            counter -= 1

                        for i in range(cnprocs):
                            cq.put(None)
                        while cq_finish.qsize() != cnprocs:
                            sleep(5)

                        for p in cprocs:
                            p.terminate()

                        print jk
                        j = len(peptides.peptideslist) - 1
                        while j >= 0:
                            if peptides.peptideslist[j].peptscore == 0:
                                peptides.peptideslist.pop(j)
                            else:
                                peptides.peptideslist[j].peptscore = 1.0
                            j -= 1

                    copy_peptides, _, _ = peptides.filter_evalue_new(FDR=FDR, useMP=False)
                    descriptors = prepare_hist(descriptors, copy_peptides)

                    jk = manager.dict()
                    cprocs = []
                    cnprocs = peptides.settings.getint('options', 'threads')
                    cq = multiprocessing.Queue()
                    cq_output = multiprocessing.Queue()
                    cq_finish = multiprocessing.Queue()

                    for i in range(cnprocs):
                        p = multiprocessing.Process(target=calc_peptscore, args=(cq, cq_output, descriptors, jk, cq_finish))
                        cprocs.append(p)
                        p.start()

                    counter = 0
                    for idx, peptide in enumerate(peptides.peptideslist):
                        # if counter < 10000:
                        cq.put([idx, peptide])
                        counter += 1
                        if counter >= 10000:
                            while counter != 0:
                                ind, pscore = cq_output.get()
                                peptides.peptideslist[ind].peptscore = float(pscore)
                                counter -= 1
                    while counter != 0:
                        ind, pscore = cq_output.get()
                        peptides.peptideslist[ind].peptscore = float(pscore)
                        counter -= 1

                    for i in range(cnprocs):
                        cq.put(None)
                    while cq_finish.qsize() != cnprocs:
                        sleep(5)

                    for p in cprocs:
                        p.terminate()

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
                    print k, 'k factor'

                    plot_MP(descriptors, peptides, fig, FDR, FDR_new, valid_proteins, settings, threshold0, curfile)
                else:
                    fig = plt.figure(figsize=(16, 12))
                    DPI = fig.get_dpi()
                    fig.set_size_inches(2000.0/float(DPI), 2000.0/float(DPI))
                    plot_MP(descriptors, peptides, fig, FDR, 0, valid_proteins, settings, threshold0, curfile)


def find_optimal_xy(descriptors):
    x, y = 1, 1
    while x * y < len(descriptors) + 10:
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
    for dbname in prots:
        PAI = float(prots[dbname]['PSMs']) / protsN[dbname]
        emPAI = 10 ** PAI - 1
        prots[dbname]['emPAI'] = emPAI

    sum_emPAI = sum(val['emPAI'] for val in prots.itervalues())
    if norm:
        for dbname in prots.keys():
            prots[dbname]['emPAI'] = prots[dbname]['emPAI'] / sum_emPAI
    return prots, sum_emPAI


def PSMs_info(peptides, valid_proteins, settings, fig=False, printresults=True, tofile=False, curfile=False, loop=True, ox=False, oy=False):
    def keywithmaxval(d):
        #this method is much faster than using max(prots.iterkeys(), key=(lambda key: prots[key]))
        v=list(d.values())
        k=list(d.keys())
        return k[v.index(max(v))]
    # full_sequences = set((peptide.sequence for peptide in peptides.peptideslist))
    tostay = set()
    prots = defaultdict(int)
    prots_pep = defaultdict(set)
    peptides_added = defaultdict(set)
    for peptide in peptides.peptideslist:
        if peptide.sequence not in peptides_added:
            if peptide.note2 == 'wr':
                add_label = 'L'
            else:
                add_label = ''
            for protein in peptide.parentproteins:
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
            add_label = 'L'
        else:
            add_label = ''

        for protein in peptide.parentproteins:
            tmp_dbname = add_label + protein.dbname
            Total_prots.add(tmp_dbname)
            try:
                prots[tmp_dbname]['PSMs'] += 1
                prots[tmp_dbname]['sumI'] += peptide.sumI
                prots[tmp_dbname]['pept'].add(peptide.sequence)
            except:
                prots[tmp_dbname] = dict()
                prots[tmp_dbname]['sumI_mod'] = defaultdict(float)
                prots[tmp_dbname]['pept'] = set([peptide.sequence, ])
                prots[tmp_dbname]['PSMs'] = 1
                prots[tmp_dbname]['sumI'] = peptide.sumI
                prots[tmp_dbname]['evalues'] = []
                prots[tmp_dbname]['expect'] = 1
                prots[tmp_dbname]['description'] = protein.description
            prots[tmp_dbname]['sumI_mod'][peptide.sequence] += peptide.sumI
            if tmp_dbname in valid_proteins and peptide.note != 'decoy':
                true_prots.add(tmp_dbname)
        if peptide.sequence not in peptides_added:
            for protein in peptide.parentproteins:
                tmp_dbname = add_label + protein.dbname
                try:
                    prots[tmp_dbname]['Peptides'] += 1
                except:
                    prots[tmp_dbname]['Peptides'] = 1
        peptides_added.add(peptide.sequence)

    def calc_expect_log(es, s, N, T):
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
#        beta = float(N) / T
#        expect += sum([np.log10((s - i) / (n - i)) for i in range(n)])
#        expect += n * np.log10(beta) + (s - n) * np.log10(1 - beta) - np.log10(s) - (n - 1) * np.log10(N)
        return expect

    for dbname in prots.keys():
        if 'Peptides' not in prots[dbname]:
            del prots[dbname]

    for k in prots:
        prots[k]['fullgroup'] = set()

    for pep in peptides.peptideslist:
        tprots = set([pr.dbname for pr in pep.parentproteins])
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
        # Add normal s, T calculation
        N = len(new_peptides.peptideslist)
        s = new_peptides.get_number_of_spectra()
        T = protsL['total peptides']
        for peptide in new_peptides.peptideslist:
            if peptide.note2 == 'wr':
                add_label = 'L'
            else:
                add_label = ''
            for protein in peptide.parentproteins:
                tmp_dbname = add_label + protein.dbname
                if tmp_dbname in prots: # <---- WTF??? not works with tmp_dbname in tostay
                    prots[tmp_dbname]['evalues'].append(peptide.qval)
                if tmp_dbname in prots_full:
                    prots_full[tmp_dbname]['evalues'].append(peptide.qval)
        for k in prots:
            prots[k]['expect'] = calc_expect_log(prots[k]['evalues'], s, N, T)
        for k in prots_full:
            prots_full[k]['expect'] = calc_expect_log(prots_full[k]['evalues'], s, N, T)

        protFDR = peptides.settings.getfloat('options', 'protFDR')
        if protFDR >= 0:
            prots = filter_evalue_prots(prots, FDR=protFDR)
            all_prots = set()
            for v in prots.itervalues():
                all_prots.update(v['fullgroup'].split(';'))
            for k in prots_full.keys():
                if k not in all_prots:
                    del prots_full[k]
        else:
            peptides.filter_decoy()
            for k in prots.keys():
                if k.startswith('L'):
                    del prots[k]
            for k in prots_full.keys():
                if k.startswith('L'):
                    del prots_full[k]
        #protein new sumI
        for k in prots:
            prots[k]['sumI_mod'] = np.median(prots[k]['sumI_mod'].values())
        sumI_mod_norm = sum(x['sumI_mod'] for x in prots.itervalues())
        if not sumI_mod_norm:
            sumI_mod_norm = 1.0
        for k in prots.keys():
            prots[k]['sumI_mod'] = prots[k]['sumI_mod'] / sumI_mod_norm / protsL[k]
        #protein sumI normalization


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
        except:
            print 'Cannot plot quantitation figures'
        output_proteins = open('%s/%s_proteins.csv' % (ffolder, fname), 'w')
        output_proteins.write('dbname\tdescription\tPSMs\tpeptides\tsequence coverage\tLFQ(SIn)\tLFQ(NSAF)\tLFQ(emPAI)\tprotein LN(e-value)\tall proteins\n')
        output_proteins_full = open('%s/%s_proteins_full.csv' % (ffolder, fname), 'w')
        output_proteins_full.write('dbname\tdescription\tPSMs\tpeptides\tsequence coverage\tLFQ(SIn)\tLFQ(NSAF)\tLFQ(emPAI)\tprotein LN(e-value)\tall proteins\n')
        output_PSMs = open('%s/%s_PSMs.csv' % (ffolder, fname), 'w')
        output_PSMs.write('sequence\tmodified_sequence\tm/z exp\tcharge\tm/z error in ppm\tmissed cleavages\te-value\tMPscore\tRT exp\tspectrum\tproteins\tproteins description\tSIn\n')
        output_peptides_detailed = open('%s/%s_peptides.csv' % (ffolder, fname), 'w')
        output_peptides_detailed.write('sequence\tmodified_sequence\tm/z exp\tcharge\tmissed cleavages\te-value\tMPscore\tRT exp\tspectrum\tproteins\tproteins description\tSIn\n')
        pickle.dump(peptides.RC, open('%s/%s_RC.pickle' % (ffolder, fname), 'w'))
        if protsC:
            output_proteins_valid = open('%s/%s_proteins_valid.csv' % (ffolder, fname), 'w')
            temp_data = []
        for k, v in prots.items():
            if protsC and k in valid_proteins:
                output_proteins_valid.write('%s,%s,%s,%s,%s\n' % (k, v['PSMs'], v['Peptides'], v['sumI'], protsC[k]))
                temp_data.append([float(v['sumI']), protsC[k]])
            if int(v['Peptides']) > 0:
                sqc = calc_sq(protsS.get(k, []), v['pept'])
                output_proteins.write('%s\t%s\t%s\t%s\t%0.1f\t%s\t%s\t%s\t%s\t%s\n' % (k, v['description'], v['PSMs'], v['Peptides'], sqc, v['sumI'],v['NSAF'], v['emPAI'], v['expect'], v['fullgroup']))

        for k, v in prots_full.items():
            if int(v['Peptides']) > 0:
                sqc = calc_sq(protsS.get(k, []), v['pept'])
                output_proteins_full.write('%s\t%s\t%s\t%s\t%0.1f\t%s\t%s\t%s\t%s\t%s\n' % (k, v['description'], v['PSMs'], v['Peptides'], sqc, v['sumI'],v['NSAF'], v['emPAI'], v['expect'], v['fullgroup']))

        for peptide in peptides.peptideslist:
            if any(protein.dbname in prots for protein in peptide.parentproteins):
                output_PSMs.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (peptide.sequence, peptide.modified_sequence, peptide.mz, peptide.pcharge, peptide.mass_diff(), peptide.mc, peptide.evalue, peptide.peptscore, peptide.RT_exp, peptide.spectrum))
                for protein in peptide.parentproteins:
                    output_PSMs.write('%s;' % (protein.dbname, ))
                output_PSMs.write('\t')
                for protein in peptide.parentproteins:
                    output_PSMs.write('%s;' % (protein.description, ))
                output_PSMs.write('\t%s\n' % (peptide.sumI, ))

        peptides_best = dict()
        peptides_best_sp = dict()
        for peptide in peptides.peptideslist:
            if peptide.sequence in peptides_best:
                if peptide.evalue < peptides_best[peptide.sequence]:
                    peptides_best[peptide.sequence] = peptide.evalue
                    peptides_best_sp[peptide.sequence] = peptide.spectrum
            else:
                peptides_best[peptide.sequence] = peptide.evalue
                peptides_best_sp[peptide.sequence] = peptide.spectrum
        for peptide in peptides.peptideslist:
            if peptide.spectrum == peptides_best_sp[peptide.sequence]:
                if any(protein.dbname in prots for protein in peptide.parentproteins):
                    output_peptides_detailed.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (peptide.sequence, peptide.modified_sequence, peptide.mz, peptide.pcharge, peptide.mc, peptide.evalue, peptide.peptscore, peptide.RT_exp, peptide.spectrum))
                    for protein in peptide.parentproteins:
                        output_peptides_detailed.write('%s;' % (protein.dbname, ))
                    output_peptides_detailed.write('\t')
                    for protein in peptide.parentproteins:
                        output_peptides_detailed.write('%s;' % (protein.description, ))
                    output_peptides_detailed.write('\t%s\n' % (peptide.sumI, ))
        if protsC:
            temp_sum = sum([x[0] for x in temp_data])
            temp_data = [[x[0] / temp_sum, x[1]] for x in temp_data]
            print 'conc: ', auxiliary.linear_regression([x[0] for x in temp_data], [x[1] for x in temp_data])

        output_proteins.close()
        output_proteins_full.close()
    if printresults:
        print 'PSMs: %s' % (len([1 for x in peptides.peptideslist if x.note2 == 'tr']), )
        print 'Peptides: %d' % (len(set(p.sequence for p in peptides.peptideslist if p.note2 == 'tr')))
        print 'Protein groups: %s' % (sum(1 for k in prots if not k.startswith('L')))
        print 'Protein groups with >= 2 peptides: %s' % (sum([1 for k, v in prots.iteritems() if v['Peptides'] >= 2 and not k.startswith('L')]))
        if valid_proteins:
            print 'PSMs_true: %s' % (len([1 for x in peptides.peptideslist if x.note3]), )
            print 'Peptides_true: %d' % (len(set(x.sequence for x in peptides.peptideslist if x.note3)), )
            print 'Protein groups_true: %s' % (len(true_prots), )
            print 'Real FDR = %s' % (100 * float(len([1 for x in peptides.peptideslist if not x.note3])) / len(peptides.peptideslist) )
        print '\n'
    return (len([1 for x in peptides.peptideslist if x.note2 == 'tr']), len(set(p.sequence for p in peptides.peptideslist)), len([v for v in prots.values() if v['Peptides'] > 1]))

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
        print form[2]
        if form[1].endswith(', peptides'):
            array_valid = [form[0](peptide) for peptide in peptides.peptideslist if peptide.evalue == tmp_dict.get(peptide.sequence, None)]
        else:
            array_valid = [form[0](peptide) for peptide in peptides.peptideslist if peptide.note2 == 'tr']
        if separatefigs:
            plt.clf()
            fig = plt.figure()
            DPI = fig.get_dpi()
            fig.set_size_inches(300.0/float(DPI), 300.0/float(DPI))
            ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.add_subplot(ox, oy, idx+1)
        lbin = min(array_valid)
        rbin = max(array_valid)
        if form[1].startswith('peptide length'):
            binsize = 1
            lbin = lbin-0.5
        else: binsize = FDbinSize(array_valid)
        H1, _ = np.histogram(array_valid, bins=np.arange(lbin, rbin+binsize, binsize))
        ind = np.arange(lbin, rbin, binsize)
        width = binsize
        ax.bar(ind, H1, width, color='#AE0066', alpha=0.8)
        ax.set_ylabel('# of identifications')
        ax.set_xlabel(form[2])
        
        from matplotlib.ticker import MaxNLocator
        ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6))

        if separatefigs:
            plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)

            if peptides.settings.get('options', 'files') == 'union':
                fname = 'union'
            else:
                fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]

            plt.savefig('%s/%s_%s.png' % (path.dirname(path.realpath(curfile)), fname, form[1]))
            if savesvg:
                plt.savefig('%s/%s_%s.svg' % (path.dirname(path.realpath(curfile)), fname, form[1]))
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
            fig.set_size_inches(300.0/float(DPI), 300.0/float(DPI))
            ax = fig.add_subplot(1, 1, 1)
        else:
        # fig = plt.figure()
            ax = fig.add_subplot(ox, oy, idx + 10)
        array_wrong = [descriptor.formula(peptide) for peptide in peptides.peptideslist if peptide.note2 == 'wr']
        array_valid = [descriptor.formula(peptide) for peptide in peptides.peptideslist if peptide.note2 == 'tr']

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
        H1, _ = np.histogram(array_wrong, bins=np.arange(lbin, rbin+binsize, binsize))
        H2, _ = np.histogram(array_valid, bins=np.arange(lbin, rbin+binsize, binsize))
        if descriptor.group == 'B':
            H3, _ = np.histogram([np.log10(descriptor.formula(peptide)) for peptide in copy_peptides.peptideslist], bins=np.arange(lbin, rbin+binsize, binsize))
        else:
            H3, _ = np.histogram([descriptor.formula(peptide) for peptide in copy_peptides.peptideslist], bins=np.arange(lbin, rbin+binsize, binsize))

        if descriptor.name in ['precursor mass difference, ppm']:
            print 'MEDIAN precursor mass difference of top PSMs=%s ppm' % (np.median([descriptor.formula(peptide) for peptide in copy_peptides.peptideslist]), )
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
            ind=np.append(ind[0]-1,ind)
            H1=np.append(H1[0],np.append(0,H1[1:]))
            H2=np.append(H2[0],np.append(0,H2[1:]))
            H3=np.append(H3[0],np.append(0,H3[1:]))
        if len(ind)>50: 
            ax.bar(ind[:-1], H1, width, color='#AE0066', alpha=0.4, edgecolor='#AE0066')
            ax.bar(ind[:-1], H2, width, color='#000099', alpha=0.4, edgecolor='#000099')
            ax.bar(ind[:-1], H3, width, color='#007E08', alpha=1, edgecolor='#007E08')
            ax.step(ind, np.append(0,H2), color='#000099',alpha=0.8)
            ax.step(ind, np.append(0,H1), color='#AE0066',alpha=0.8)
        else:
            ax.bar(ind[:-1], H1, width, color='#AE0066', alpha=0.4)
            ax.bar(ind[:-1], H2, width, color='#000099', alpha=0.4)
            ax.bar(ind[:-1], H3, width, color='#007E08', alpha=1)
            ax.step(ind, np.append(0,H2), color='#000099',alpha=0.8)
            ax.step(ind, np.append(0,H1), color='#AE0066',alpha=0.8)
        if any(descriptor.name.startswith(clabel) for clabel in ['missed cleavages', 'charge states', 'isotopes mass difference, Da']):
            ax.set_xticks(np.arange(0.5, 5.5, 1.0))
            fig.canvas.draw()
            labels = [item.get_text() for item in ax.get_xticklabels()]
            ax.set_xticklabels([int(float(l)) for l in labels])
        elif descriptor.name.startswith('potential modifications'):
            ax.set_xticks(np.arange(-1.5,5.5,1))
            ax.set_xticklabels(['NA','']+list(np.arange(0, 6, 1)))
        else:
            from matplotlib.ticker import MaxNLocator
            ax.get_xaxis().set_major_locator(MaxNLocator(nbins=6))
        
            

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
            plt.savefig('%s/%s_%s.png' % (path.dirname(path.realpath(curfile)), fname, descriptor.name))
            if savesvg:
                plt.savefig('%s/%s_%s.svg' % (path.dirname(path.realpath(curfile)), fname, descriptor.name))
    return fig


def plot_quantiation(prots, curfile, settings, fig, separatefigs, ox, oy):
    for i, idx in enumerate(['sumI', 'NSAF', 'emPAI']):
        dots = [np.log10(v[idx]) for v in prots.itervalues()]
        if separatefigs:
            plt.clf()
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            DPI = fig.get_dpi()
            fig.set_size_inches(300.0/float(DPI), 300.0/float(DPI))
        else:
            ax = fig.add_subplot(ox, oy, i+7)
        ax.hist(dots, bins = 10, alpha=0.8, color='#000099')
        ax.set_xlabel('LOG10(%s)' % (idx, ))
        ax.set_ylabel('Number of proteins')
        ax.locator_params(axis='x', nbins=4)
        if settings.get('options', 'files') == 'union':
            fname = 'union'
        else:
            fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
        if separatefigs:
            plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)
            plt.savefig('%s/%s_%s.png' % (path.dirname(path.realpath(curfile)), fname, idx))
            if settings.getboolean('advanced options', 'saveSVG'):
                plt.savefig('%s/%s_%s.svg' % (path.dirname(path.realpath(curfile)), fname, idx))
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

    print 'MP filtering:'
    PSMs_info(copy_peptides, valid_proteins, settings, fig=fig, tofile=True, curfile=curfile, ox=ox, oy=oy)
    ffolder = path.dirname(path.realpath(curfile))
    if peptides.settings.get('options', 'files') == 'union':
        fname = 'union'
    else:
        fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
    output_PSMs_full = open('%s/%s_PSMs_full.csv' % (ffolder, fname), 'w')
    output_PSMs_full.write('sequence\tmodified_sequence\tm/z exp\tm/z error in ppm\tmissed cleavages\te-value\tMPscore\tRT exp\tspectrum\tproteins\tproteins description\tSIn\tdecoy\n')
    for peptide in peptides.peptideslist:
        output_PSMs_full.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t' % (peptide.sequence, peptide.modified_sequence, peptide.mz, peptide.mass_diff(), peptide.mc, peptide.evalue, peptide.peptscore, peptide.RT_exp, peptide.spectrum))
        for protein in peptide.parentproteins:
            output_PSMs_full.write('%s;' % (protein.dbname, ))
        output_PSMs_full.write('\t')
        for protein in peptide.parentproteins:
            output_PSMs_full.write('%s;' % (protein.description, ))
        output_PSMs_full.write('\t%s\t%s\n' % (peptide.sumI, peptide.note2))
    print 'Without filtering, after removing outliers:'
    PSMs_info(peptides, valid_proteins, settings, loop=False)

    if sepfigures:
        plt.clf()
        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(300.0/float(DPI), 300.0/float(DPI))
        ax = fig.add_subplot(1, 1, 1)
    else:
        ax = fig.add_subplot(ox, oy, ox*oy)
    ax.plot([x[0] for x in PSMs_wrong], [x[1] for x in PSMs_wrong], 'o', markersize=2, color='#AE0066')
    ax.plot([x[0] for x in PSMs_true], [x[1] for x in PSMs_true], 'o', markersize=2, color='#000099')
    ax.axvline(threshold1, color='g')
    if threshold2:
        ax.axhline(threshold2, color='g')
    if threshold0:
        threshold0 = -np.log(threshold0)
        ax.axvline(threshold0, color='#AE0066')
    ax.set_ylim(min(x[1] for x in PSMs_true) - 1, max(x[1] for x in PSMs_true) + 1)
    if peptides.settings.get('options', 'files') == 'union':
        fname = 'union'
    else:
        fname = path.splitext(path.splitext(path.basename(curfile))[0])[0]
    ax.set_xlabel('-LOG(evalue)')
    ax.set_ylabel('LOG(MPscore)')
    if sepfigures:
        plt.gcf().subplots_adjust(bottom=0.15, left=0.2, top=0.95, right=0.9)
        plt.savefig('%s/%s_%s.png' % (path.dirname(path.realpath(curfile)), fname, 'scores'))
        if settings.getboolean('advanced options', 'saveSVG'):
            plt.savefig('%s/%s_%s.svg' % (path.dirname(path.realpath(curfile)), fname, 'scores'))
    else:
        plt.gcf().subplots_adjust(bottom=0.05, left=0.05, top=0.95, right=0.95)
        plt.savefig('%s/%s.png' % (path.dirname(path.realpath(curfile)), fname))
        if settings.getboolean('advanced options', 'saveSVG'):
            plt.savefig('%s/%s.svg' % (path.dirname(path.realpath(curfile)), fname))


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
    return descriptors


def calc_peptscore(cq, cq_output, descriptors, jk, cq_finish):
    for idx, peptide in iter(cq.get, None):
        tmp_peptscore = peptide.peptscore
        for descriptor in descriptors:
            descriptor_value = descriptor.formula(peptide)
            if tmp_peptscore != 0:
                flag = 1
            else:
                flag = 0
            if descriptor.group == 'A':
                if descriptor_value < descriptor.hist[1][0] or descriptor_value >= descriptor.hist[1][-1]:
                    tmp_peptscore = 0
                else:
                    j = descriptor.hist[1].searchsorted(descriptor_value)
                    if descriptor_value < descriptor.hist[1][j]:
                        j -= 1
                    try:
                        tmp_peptscore *= float(descriptor.hist[0][j]) / sum(descriptor.hist[0])
                    except:
                        print descriptor.name
                        print descriptor.hist[0]
                        print descriptor.hist[1]
                        print descriptor_value, descriptor.hist[1][j], descriptor.hist[1][descriptor.hist[1].searchsorted(descriptor_value)], descriptor.hist[1].searchsorted(descriptor_value)
            elif descriptor.group == 'B':
                tmp_peptscore *= float(descriptor.array.searchsorted(descriptor_value, side='right')) / descriptor.array.size

            if flag and not tmp_peptscore:
                try:
                    jk[descriptor.name] += 1
                except:
                    jk[descriptor.name] = 1
        cq_output.put([idx, tmp_peptscore])
    cq_finish.put(True)


def main(argv_in, union_custom=False):
    inputfile = argv_in[1]
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
        if path.splitext(arg)[-1].lower() == '.fasta':
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
        settings = get_settings(configfile)
    else:
        settings = get_settings('default.cfg')
    try:
        settings.getboolean('advanced options', 'saveSVG')
    except:
        settings.set('advanced options', 'saveSVG', '0')
    try:
        settings.get('advanced options', 'allowed peptides')
    except:
        settings.set('advanced options', 'allowed peptides', '')
    try:
        settings.getboolean('advanced options', 'separatefigures')
    except:
        settings.set('advanced options', 'separatefigures', '0')
    try:
        settings.getboolean('input', 'add decoy')
    except:
        if 'input' not in settings.sections():
            settings.add_section('input')
        settings.set('input', 'add decoy', '0')
    if union_custom:
        settings.set('options', 'files', 'union')

    proteases = [x.strip() for x in settings.get('missed cleavages', 'protease1').split(',')]
    proteases.extend([x.strip() for x in settings.get('missed cleavages', 'protease2').split(',')])
    expasy = '|'.join((parser.expasy_rules[protease] if protease in parser.expasy_rules else protease for protease in proteases))
    try:
        mc = settings.getint('missed cleavages', 'number of missed cleavages')
    except:
        print 'Number of missed cleavages is missed in parameters, using 2 for normalization and emPAI calculation'
        mc = 2
    fprocs = []
    fnprocs = 12
    fq = multiprocessing.Queue()
    fq_output = multiprocessing.Queue()

    protsL['total proteins'] = 0
    protsL['total peptides'] = 0
    def protein_handle(fq, fq_output, protsL, protsN, protsS, expasy, mc):
        while 1:
            try:
                x = fq.get(timeout=1)
            except Empty:
                fq_output.put('1')
                break

            try:
                if not any(x[0].startswith(tag) for tag in ['sp', 'tr', 'DECOY_sp', 'DECOY_tr']):
                    if any(tag in x[0] for tag in ['SWISS-PROT:', 'TREMBL:']):
                        dbname = x[0].split(' ')[0]
                    else:
                        dbname = x[0]#.replace('>', ' ')
                    if 'DECOY_' not in x[0]:
                        protsS[dbname] = x[1]
                else:
                    dbname = x[0].split('|')[1]
                    if 'DECOY_' not in x[0]:
                        protsS[x[0].split('|')[1]] = x[1]
                if 'DECOY_' not in x[0]:
                    protsL[dbname] = len(x[1])
                    protsN[dbname] = len(parser.cleave(x[1], expasy, mc))
                protsL['total proteins'] += 1
                protsL['total peptides'] += protsN.get(dbname, 0)
            except:
                print 'Smth wrong with reading of FASTA file'

    if fastafile:
        if settings.getboolean('input', 'add decoy'):
            decoy_method = settings.get('input', 'decoy method')
            decoy_prefix = settings.get('input', 'decoy prefix')
            for x in fasta.decoy_db(fastafile, mode=decoy_method, prefix=decoy_prefix):
                fq.put(x)
        else:
            for x in fasta.read(fastafile):
                fq.put(x)

    for i in range(fnprocs):
        p = multiprocessing.Process(target=protein_handle, args=(fq, fq_output, protsL, protsN, protsS, expasy, mc))
        fprocs.append(p)
        p.start()

    while fq_output.qsize() != fnprocs:
        sleep(10)
    for p in fprocs:
        p.terminate()

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

if __name__ == '__main__':
    main(argv)
