
import sys
import HapFileIO
import PcaFileIO

from LineageTrack import get_matched_info, Trackbuild
from PcaPlot import PcaPlot

def run_hap():
    arguments = HapFileIO.hap_parser()
    filename, input_file, remove_indel = HapFileIO.judge_inputfile_type(arguments.input, HapFileIO.work_path)
    if arguments.output:
        output_file = set_outputfile(arguments.output, HapFileIO.work_path)
    else:
        if input_file.endswith(r'.vcf.gz|.inp.gz'):
            output_file = input_file.rstrip('.gz').rstrip(r'vcf|inp') + 'hap'
        else:
            output_file = input_file.rstrip(r'vcf|inp') + 'hap'
    HapFileIO.judge_missing(arguments.missing)
    HapFileIO.set_log(input_file, filename, remove_indel, arguments)

    (info_matched,
    ind_hap_dict,
    ind_pos_dict,
    ind_keyhap_dict,
    nohap_ind,
    key_hap_num) = get_matched_info(arguments, input_file)

    tracker = Trackbuild(key_hap_num, HapFileIO.data_path)
    freq_info_ind = tracker.count_freq(ind_hap_dict)
    freq_info_sum = tracker.count_freq(ind_pos_dict)
    freq_info_ind = tracker.num_inter(freq_info_ind, freq_info_sum)
    hap_info, key_hap_info = tracker.search_final_hap(arguments.filter, freq_info_ind, ind_keyhap_dict)
    hap_info, key_hap_info = tracker.search_mut(hap_info, key_hap_info, info_matched)
    key_hap_info = tracker.search_population(info_matched, key_hap_info)
    hap_info = tracker.tree_track(hap_info)
    info_gathered = tracker.tree_iter(hap_info, key_hap_info, freq_info_ind)
    hap_data = tracker.write_info(HapFileIO.start, info_gathered, output_file, nohap_ind)

    if arguments.pca:
        pca_output_file = output_file.rstrip('hap') + 'png'
        population_file, extra_file, classify, hap_type = HapFileIO.get_pca_file(arguments.pca ,HapFileIO.work_path)
        pca = PcaPlot(population_file)
        hap_df = pca.hap_read(hap_data, extra_file, hap_type, arguments.filter, classify)
        pca.pca_plot(hap_df, population_file, pca_output_file)


def run_pca():
    arguments = PcaFileIO.pca_parser()
    pca_hap_df, pca_inputfile, pca_inputfile_extra, filename = PcaFileIO.judge_pca_inputfile(arguments.input, HapFileIO.work_path)
    if arguments.output:
        output_file = PcaFileIO.set_outputfile(pca_inputfile, HapFileIO.work_path)
    else:
        output_file = pca_inputfile.rstrip('hap') + 'png'
    PcaFileIO.set_log(pca_inputfile, filename, arguments)
    population_file = PcaFileIO.judge_info_file(arguments.info, HapFileIO.work_path)

    pca = PcaPlot(population_file)
    hap_df = pca.hap_read(pca_hap_df, pca_inputfile_extra, arguments.type, arguments.filter, arguments.mode[0])
    pca.pca_plot(hap_df, output_file, arguments.freq)


def run_phy():
    pass

if __name__ == '__main__':
    if sys.argv[1] == 'hap':
        run_hap()
    elif sys.argv[1] == 'pca':
        run_pca()
    elif sys.argv[1] == 'phy':
        pass
