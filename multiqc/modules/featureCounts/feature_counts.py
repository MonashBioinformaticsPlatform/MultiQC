#!/usr/bin/env python

""" MultiQC module to parse output from featureCounts """

from __future__ import print_function
from collections import OrderedDict
import logging

from multiqc import config, BaseMultiqcModule, plots

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='featureCounts', 
        anchor='featurecounts', target='Subread featureCounts', 
        href='http://bioinf.wehi.edu.au/featureCounts/', 
        info="is a highly efficient general-purpose read summarization program"\
        " that counts mapped reads for genomic features such as genes, exons,"\
        " promoter, gene bodies, genomic bins and chromosomal locations.")

        # Find and load any featureCounts reports
        self.featurecounts_data = dict()
        self.featurecounts_keys = list()
        for f in self.find_log_files(config.sp['featurecounts']):
            self.parse_featurecounts_report(f)

        if len(self.featurecounts_data) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.featurecounts_data)))

        # Write parsed report data to a file
        self.write_data_file(self.featurecounts_data, 'multiqc_featureCounts')

        for f in self.featurecounts_data:
            # Basic Stats Table
            # Report table is immutable, so just updating it works
            self.featurecounts_stats_table(f)

            # Assignment bar plot
            # Only one section, so add to the intro
            self.intro += self.featureCounts_chart(f)


    def parse_featurecounts_report (self, f):
        """ Parse the featureCounts log file. """
        
        file_names = list()
        parsed_data = dict()
        for l in f['f'].splitlines():
            s = l.split("\t")
            if len(s) < 2:
                continue
            if s[0] == 'Status':
                for f_name in s[1:]:
                    file_names.append(f_name)
            else:
                k = s[0]
                parsed_data[k] = list()
                if k not in self.featurecounts_keys:
                    self.featurecounts_keys.append(k)
                for val in s[1:]:
                    parsed_data[k].append(int(val))
        
        for idx, f_name in enumerate(file_names):
            
            # Clean up sample name
            s_name = self.clean_s_name(f_name, f['root'])
            
            # Reorganised parsed data for this sample
            # Collect total count number
            data = dict()
            data['Total'] = 0
            for k in parsed_data:
                data[k] = parsed_data[k][idx]
                data['Total'] += parsed_data[k][idx]
            
            # Calculate the percent aligned if we can
            if 'Assigned' in data:
                data['percent_assigned'] = (float(data['Assigned'])/float(data['Total'])) * 100.0
            
            # Add to the main dictionary
            if len(data) > 1:
                file_stem = f['s_name']
                if s_name in self.featurecounts_data:
                    log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                self.add_data_source(f, s_name)
                if file_stem not in self.featurecounts_data:
                    self.featurecounts_data[file_stem] = {}
                self.featurecounts_data[file_stem][s_name] = data
        

    def featurecounts_stats_table(self,f):
        """ Take the parsed stats from the featureCounts report and add them to the
        basic stats table at the top of the report """
        
        headers = OrderedDict()
        headers['percent_assigned'] = {
            'title': '% Assigned : '+f,
            'description': '% Assigned reads',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn',
            'format': '{:.1f}%'
        }
        headers['Assigned'] = {
            'title': 'M Assigned : '+f,
            'description': 'Assigned reads (millions)',
            'min': 0,
            'scale': 'PuBu',
            'modify': lambda x: float(x) / 1000000,
            'shared_key': 'read_count'
        }
        self.general_stats_addcols(self.featurecounts_data[f], headers)


    def featureCounts_chart (self,f):
        """ Make the featureCounts assignment rates plot """
        
        # Config for the plot
        config = {
            'id': 'featureCounts_assignment_plot_'+f,
            'title': 'featureCounts Assignments : '+f,
            'ylab': '# Reads',
            'cpswitch_counts_label': 'Number of Reads'
        }
        
        return plots.bargraph.plot(self.featurecounts_data[f], self.featurecounts_keys, config)
