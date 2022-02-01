from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
from PyQt5.QtCore import *
from PyQt5 import uic

import sys
import os
import gzip
import pysam
import random

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib
from matplotlib.figure import Figure
from matplotlib import pyplot as plt


def init(inp_dict,key,value):
    if not key in inp_dict:     inp_dict[key] = value
    
def getRangeOvlp(range1,range2):
    num_bp_ovlp =-( ( max([range1[-1],range2[-1]]) - min(range1[0],range2[0]) ) - ( (range1[-1]-range1[0]) + (range2[-1]-range2[0]) ) )
    return num_bp_ovlp #Return: positive vals = numBPovlp. Neg vals: Distance ranges are apart

def computeOvlpRanges_wIdxs(hit_ranges,attribute=None):
    """
    inputs: array of ranges, input of form [start,identifier,end]
    outputs: array of ranges and the active hit_ranges index for that range
    """   
    events = {}
    for hit_range in hit_ranges:
        start = hit_range[0]
        end = hit_range[-1]
        query = hit_range[1]
        entry_id = hit_range[2]
        
        # Check if we have attribute we should parse from hit_range[1]
        if attribute != None: query = query.attribute
        
        # init pos
        if not start in events:     events[start] = {'starts':{},'ends':{}}
        if not end in events:     events[end] = {'starts':{},'ends':{}}

        # init query
        if not query in events[start]['starts']:     events[start]['starts'][query] = set()
        if not query in events[end]['ends']:         events[end]['ends'][query] = set()
        
        # add entry
        events[start]['starts'][query].add(entry_id)
        events[end]['ends'][query].add(entry_id)

    ranges = []
    active = {}
    old_active = {}
    for pos in sorted(events):
        old_active = dict(active)
        
        ## update active
        # add query+entry_id from START
        for query in events[pos]['starts']:
            # check if query it doesnt exist, then start
            if not query in active:      active[query] = set()
            # add all query entries
            active[query].update(events[pos]['starts'][query])
          
        # remove query+entry_id from END
        remove_list = set()
        for query in events[pos]['ends']:
            # remove end count
            active[query].difference_update(events[pos]['ends'][query])
            # check if empty, then add to remove list
            if not active[query]:
                remove_list.add(query)
                
        # remove empty entries
        for query in remove_list:
            del active[query]
        ##/
        
        ## handle range
        #check if start new range
        if not old_active:
            ranges.append([pos,{},None])
            for query in active:
                for entry_id in active[query]:
                    if not query in ranges[-1][1]:      ranges[-1][1][query] = set()
                    ranges[-1][1][query].add(entry_id)
            
        #check if close old range
        if (not active) and old_active:
            ranges[-1][-1] = pos

        #check if close old + start new (tts overlap change)          
        if active and old_active and set(old_active) != set(active):
            ranges[-1][-1] = pos #close old
            ranges.append([pos,{},None])
            for query in active:
                for entry_id in active[query]:
                    if not query in ranges[-1][1]:      ranges[-1][1][query] = set()
                    ranges[-1][1][query].add(entry_id)
        ##/
        
    return ranges

### RANGEOVERLAPS, 2 functions
def rangeOverlaps_makeMappingDict(input_dict,bin_step,coordsKey=None,startKey='rstart',endKey='rend',sortByKey=None,ignoreIfHasKeys=None):
    ## Function1: make map
    """"
    # Make functions for making overlap-queries to a mapped dict, to semi-replace computeOvlpRanges
    # For example, map GTF/SV, query alignments at them. Or, map alignments, query alignments, and so on.
    # Two functions: 1/ make map, 2/ look for overlap
    # It will procude bins of poses, per rname, containing which IDXs exist in each bin.
    # Then we can dump GTF_idx_map into it, recieve back GTF_idxs_to_check -> Bin_pos -> GTF_idxs.
    # Then have another function to look in "RANGES", take start+end coord, calc start/end bin, scan all bins in-between
    # take out all all "GTF idxs" and corresponding GTF ranges, then look for overlap between RANGEE and GTF ranges
    
    ## EXAMPLE USAGE, 1/ map GTF file with bin size 1000. 2/ Query a range on X.
    GTF_idx_map_ROMD = rangeOverlaps_makeMappingDict(GTF_idx_map,1000,coordsKey='rcoords',sortByKey='rname') #ROMD = "range ovlp mapping dict"
    test_ovlp = rangeOverlaps_lookup([132000,138000] , GTF_idx_map_ROMD , 1000 , map_lookup_key='X')

    """
    coordBins = {} # "sortByKey" -> coordBins -> Idxs present
    for idx in input_dict:
        entry = input_dict[idx]
        
        # Check if we specified to filter
        skipEntry = False
        if ignoreIfHasKeys != None:
            for key in ignoreIfHasKeys:
                if key in entry and entry[key]:
                    skipEntry = True
                    break
        if skipEntry: continue
        #/
        
        # Parse start/end
        if coordsKey != None:
            start,end = entry[coordsKey]
        else:
            start,end = entry[startKey],entry[endKey]
            
        sortBy = None
        if sortByKey != None:           sortBy = entry[sortByKey] #update sortBy if we provided key
        
        # Round up/down coords
        start_roundD = int(start/float(bin_step)) *bin_step
        end_roundU = int(end/float(bin_step)) *bin_step + bin_step
            
        # Save at map, add idx to all bins in-between start/end rounded coords
        init(coordBins,sortBy,{})
        for bin_to_save in range(start_roundD,end_roundU+bin_step,bin_step): #rangee is open end
            init(coordBins[sortBy],bin_to_save,[])
            coordBins[sortBy][bin_to_save].append([start,idx,end])
        #/Save
        
    # If no sortByKey is provided, return dict of bin->idxs
    if sortByKey == None:       return coordBins[sortByKey]
    # Else, return full ("default")
    return coordBins
    ##/Function1

def rangeOverlaps_lookup(inp_coords,mapping_dict,bin_step,map_lookup_key=None,haltOnNumBPOvlp=None,skipNumBPOvlpBelow=None):
    ## Function2: overlap scan, using map (Query-function for input vs mapping dict)
        
    # pre-flight, if we provided key that the mapping dict is sorted by, parse it out
    if map_lookup_key != None:
        if not map_lookup_key in mapping_dict: return []
        mapping_dict = mapping_dict[map_lookup_key]
    
    # round input coords
    start,end = inp_coords
    start_roundD = int(start/float(bin_step)) *bin_step
    end_roundU = int(end/float(bin_step)) *bin_step + bin_step
    
    ## Get ovlp ranges
    ovlps = {} #run per target idx
    haltOnNumBPOvlp_toggle = False #for haltOnNumBPOvlp input
    for bin_to_scan in range(start_roundD,end_roundU+bin_step,bin_step): #rangee is open end
        if not bin_to_scan in mapping_dict: continue #skip binn if it doesnt exist
        # Traverse bin entries
        for (target_start,target_idx,target_end) in mapping_dict[bin_to_scan]:
            # check if we did target_idx already
            if target_idx in ovlps: continue
        
            # check if ovlp at all
            if not getRangeOvlp(inp_coords,[target_start,target_end]) >= 0: continue #0 means border ovlp/coord ovlp
            
            # calc ovlp
            ovlp_start = max([start,target_start])
            ovlp_end = min([end,target_end])
            
            numBP_ovlp = (ovlp_end-ovlp_start)
            
            # Check if we have requirements for minimal numBP ovlp for save
            if skipNumBPOvlpBelow != None and numBP_ovlp <= skipNumBPOvlpBelow:
                continue
            
            # save ovlp
            ovlps[target_idx] = [ovlp_start,ovlp_end]
            
            # Check if we want to break on first found occurrance, for instance on repeatmasking.
            if haltOnNumBPOvlp != None and (numBP_ovlp >= haltOnNumBPOvlp):
                haltOnNumBPOvlp_toggle = True
                break
        
        # Check if break outer as well
        if haltOnNumBPOvlp_toggle:
            break
    ##/
        
    # Return on "RANGES" format
    ovlps_arr = []
    for target_idx,(ovlp_start,ovlp_end) in ovlps.items():
        ovlps_arr.append([ovlp_start,target_idx,ovlp_end])
    
    return ovlps_arr
    ##/Function2
###/RANGEOVERLAPS, 2 functions

def GTF_GFF_parser(input_file,gtf_attribute_split=' ',gff_attribute_split='=',gtf_gene_tag='gene_id',gff_gene_tag='ID',gene_column_identifiers={'gene','transcript'},skip_features=set()):
    isGFF = False
    isGTF = False
    if input_file.endswith('.gff') or input_file.endswith('.gff3') or input_file.endswith('.gff2'):     isGFF = True
    if input_file.endswith('.gtf'):     isGTF = True
    
    if not (isGFF or isGTF): return 'file_not_recognized'
    
    features_pointers = {} # feature tag -> location
    genes = {}
    with open(input_file,'r') as f:
        for line in f:
            if line[0] == '#': continue #skip comment line
            line = line.split('\t')
            line[-1] = line[-1].strip('\n')
            
            rname,_source,typee,rstart,rend,_score,strand,_frame,data = line
            
            if skip_features and typee in skip_features: continue #skip entry if feature is in "skip_features". For performance reasons, lets skip entires which should be covered by an overview entry (i.e. skip exon/CDS as it is covered in gene/transcript entry)
            
            try:
                rstart,rend = int(rstart),int(rend)
            except:
                continue
            
            if data.find('; ') != -1:
                data = data.split('; ')
                if data[-1][-1] == ';':     data[-1] = data[:-1] # correct for last entry
            else:
                data = data.split(';')
            
            for entry in data:
                if entry[0] == '#': continue #skip commentary entries
                if isGFF:
                    try:
                        tagType,tag = entry.split(gff_attribute_split)
                    except:
                        continue
                elif isGTF:
                    try:
                        tagType,tag = entry.split(gtf_attribute_split)
                        tag = tag.replace('"','')
                    except:
                        continue
                
                # Save to pointer
                tag = tag.lower()
                init(features_pointers,tag,[])
                features_pointers[tag].append([rname,[rstart,rend],strand,typee])
                #/
                
                # Save to gene
                if typee in gene_column_identifiers:
                    if tag in genes: continue # take first add of entry only. This is because i.e. human GTF file contain supplemental "fix" chromosomes with corresponding annotation appearing later in file.
                    tmp_save = {'rname':rname,'rcoords':[rstart,rend],'strand':strand,'idx':len(genes),'name':tag}
                    if isGFF and tagType == gff_gene_tag:
                        genes[tag] = tmp_save
                    elif isGTF and tagType == gtf_gene_tag:
                        genes[tag] = tmp_save
                #/
    return features_pointers,genes

def repeat_parser_repeatMasker(input_file):
    repeat_entries = {}
    with open(input_file,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split()
            
            # Identify entry line by having a number at first entry
            if not line: continue
            if not line[0].isdigit():
                continue
            
            # Example lines
            # OUT-file
            #   319    4.3  0.0  0.0  211000022278031                                   1      336      (685) + (TATAA)n          Simple_repeat           1    336     (0)      1  
            #   997    4.9  0.0  0.0  211000022278031                                 345      467      (554) C HETA              LINE/I-Jockey      (4357)   1724    1602      2  
            
            SW_score,percDic,percDel,percIns,rname,rstart,rend,rstartLeftMost,strand,qname,repeat_class,qstart_rightMost,qend,qstart_leftMost,ID = line[:15]
            rstart = int(rstart)
            rend = int(rend)
            rcoords = [min([rstart,rend]),max([rstart,rend])]
            
            tmp_data = {}
            tmp_data['repeat_id'] = qname
            tmp_data['rname'] = rname
            tmp_data['rcoords'] = rcoords
            tmp_data['strand'] = strand
            tmp_data['type'] = 'repeat'
            tmp_data['len'] = rcoords[-1]-rcoords[0]
            tmp_data['strand'] = strand
            tmp_data['class'] = repeat_class
            
            repeat_entries[len(repeat_entries)] = tmp_data
    
    repeats_mapping_dict = rangeOverlaps_makeMappingDict(repeat_entries,100,coordsKey='rcoords',sortByKey='rname')
    
    return repeat_entries,repeats_mapping_dict

def VCF_parser(vcf_file,rname2_key='CHR2',rpos2_key='END',SV_type_key='SVTYPE',read_supp_key='RNAMES',SV_len_key='SVLEN'):
    def is_bnd(alt):
        if (alt.find(']') != -1 or alt.find('[') != -1) and alt.find(':') != -1:
            return True
        else:
            return False
    def is_SV(alt):
        if (alt.find('<') != -1 and alt.find('>') != -1):
            return True
        else:
            return False
    
    rname_vcfs = {}
    read_vcfs = {}
    with open(vcf_file,'r') as f:
        for ln,line in enumerate(f):
            if line[0] == '#': continue #skip commentary lines
            line = line.split()
            line[-1] = line[-1].strip('\n')
            
            rname,rpos,ID,ref,alt,qual,filter_,info = line[:8]
            rpos = int(rpos)
            
            info = info.split(';')
            data = {}
            for entry in info:
                try:
                    key,val = entry.split('=')
                    data[key] = val
                except:
                    continue
            
            # Check if ins/del from ref+alt
            rname2 = None
            rpos2 = None
            SVtype = None
            SVlen = None
            read_supp = set()
            if not is_SV(alt) and not is_bnd(alt):
                deleted_bases = len(ref)
                inserted_bases = len(alt)
                rname2 = rname
                rpos2 = rpos + deleted_bases
                if deleted_bases == 1 and inserted_bases > deleted_bases:
                    SVtype = 'INS'
                    SVlen = inserted_bases
                elif inserted_bases == 1 and deleted_bases > inserted_bases:
                    SVtype = 'DEL'
            #/
            # Parse breakend
            elif is_bnd(alt):
                rname2,rpos2 = alt.split(':')
                if rname2.find('[') != -1:      rname2 = rname2.split('[')[1]
                if rname2.find(']') != -1:      rname2 = rname2.split(']')[1]
                rpos2 = int(rpos2.split('[')[0].split(']')[0])
                SVtype = 'BND'
            #/
            # Parse SV type
            elif is_SV(alt):
                try:
                    rname2 = data[rname2_key]
                except:
                    print('Line: '+str(ln+1)+' Could not parse chr2 from vcf! Skipping...')
                    continue
                try:
                    rpos2 = int(data[rpos2_key])
                except:
                    print('Line: '+str(ln+1)+' Could not parse pos2 from vcf! Skipping...')
                    continue
                try:
                    SVtype = data[SV_type_key]
                except:
                    SVtype = alt.replace('<','').replace('>','')
                if SVtype == 'INS':
                    try:
                        SVlen = data[SV_len_key]
                    except:
                        SVlen = 0
                        print('Line: '+str(ln+1)+' Could not parse SV insertion length from vcf! Setting to 0.')
            #/
            # Parse read support
            try:
                read_supp.update(data[read_supp_key].split(','))
            except:
                print('Could not parse read support from vcf! Will not be able to link SVs to readnames!')
                read_supp.add('None')
            #/
            
            # Save
            if SVlen == None and rpos2 != None:       SVlen = abs(rpos-rpos2)
            tmp_SV = {'rname':rname,'rpos':rpos,'rname2':rname2,'rpos2':rpos2,'len':SVlen,'type':SVtype,'read_supp':read_supp,'idx':ln}
            
            if not SVtype:
                print('Could not parse SV type! Setting as None')
            for read in read_supp:
                init(read_vcfs,read,[])
                read_vcfs[read].append(tmp_SV)
            #/
    return read_vcfs

def calc_aln_masking(aln,repeat_entries,repeats_mapping_dict):
    rangees = [ [aln['rcoords'][0],'base',aln['rcoords'][-1]] ]
    for oS,idx,oE in rangeOverlaps_lookup(aln['rcoords'],repeats_mapping_dict,100,map_lookup_key=aln['rname']):
        if oE == oS: continue #skip single-nucleotide ovlps
        rangees.append([oS,idx,oE])
    ovlps = computeOvlpRanges_wIdxs(rangees)
    base_numBP = 0
    masked_numBP = 0
    for rangee in ovlps:
        if not 'base'  in rangee[1]: continue
        base_numBP += rangee[-1]-rangee[0]
        if len(rangee[1]) > 1:
            masked_numBP += rangee[-1]-rangee[0]
    masked_covFrac = masked_numBP / base_numBP
    return base_numBP,masked_numBP,masked_covFrac


def calc_aChroms_from_alignments(entries,alns_chain_distance=10000,region_margin=0,lenReq=0,show_all_for_qnames={}):
    ### Define regions to parse
    # Setup rangees per rname with BASE and entry's
    rname_regions_to_parse = {}
    for enum,entry in enumerate(entries):
        if ((entry['qcoords'][-1]-entry['qcoords'][0] >= lenReq) or (show_all_for_qnames and entry['rname'] in show_all_for_qnames)):
            rname = entry['rname']
            rstart,rend = entry['rcoords']
            init(rname_regions_to_parse,rname,[])
            rstart_fattened = rstart-alns_chain_distance
            if rstart_fattened < 0: rstart_fattened = 0
            rend_fattened = rend+alns_chain_distance
            rname_regions_to_parse[rname].append([rstart_fattened,enum,rend_fattened])
    for rname,rangees in rname_regions_to_parse.items():
        rname_regions_to_parse[rname].append([0,'base',sorted(rangees,key=lambda x: x[-1])[-1][-1]])
    #/
    # Traverse rnames, compile ranges where entry's exist
    rname_regions_PRE = [] # [rname, [rposes], set(rsqIdxs)] to plot at this "artificial chrom"
    for rname, rangees in rname_regions_to_parse.items():
        ovlpRangees = computeOvlpRanges_wIdxs(rangees)
        for rangee in ovlpRangees:
            tmp_save = {'rname':rname,'rposes':[],'rposes_wFatten':[],'entry_idxs':set()}
            
            # check if base only, then init new artificialChrom and do nothing more...
            if len(rangee[1]) == 1 and 'base' in rangee[1]: 
                if not rname_regions_PRE or len(rname_regions_PRE[-1]) >= 1:
                    rname_regions_PRE.append(tmp_save)
            else:
                # Append coords, entry_idxs to existing artificialChrom
                if not rname_regions_PRE:                           rname_regions_PRE.append(tmp_save)
                if not rname_regions_PRE[-1]['rname'] == rname:     rname_regions_PRE.append(tmp_save) # check if rname was changed
                
                for member in rangee[1]:
                    if member == 'base': continue
                    entry = entries[member]
                    rname_regions_PRE[-1]['rposes_wFatten'].append(rangee[0])
                    rname_regions_PRE[-1]['rposes_wFatten'].append(rangee[-1])
                    rname_regions_PRE[-1]['rposes'] += entry['rcoords']
                    rname_regions_PRE[-1]['entry_idxs'].add(entry['idx'])
    #/
    # Extract "artifical chroms"
    aChroms = {}
    entry_aChromEnums = {}
    for enum,aChrom_PRE in enumerate(sorted(rname_regions_PRE,key=lambda x: x['rname'])):
        rstart,rend = max(0,min(aChrom_PRE['rposes']) - region_margin) , max(aChrom_PRE['rposes']) + region_margin # avoid rstart to be less < 0
        aChroms[enum] = {'rcoords':[rstart,rend],'rname':aChrom_PRE['rname'],'len':rend-rstart,'enum':enum}
        for entry_idx in aChrom_PRE['entry_idxs']:
            entry_aChromEnums[entry_idx] = enum
    #/
    # Return
    return aChroms,entry_aChromEnums
    #/
    
def parse_bam_ref_lens(bam_handle):
    ref_lens = {}
    for entry in bam_handle.header['SQ']:
        ref,length = entry['SN'],entry['LN']
        ref_lens[ref] = length
        if not length+1:
            print('Error: reference length was not an integer!!')
    return ref_lens
    
def parse_bam_entry(bam_entry):
    rname,rstart,rend = bam_entry.reference_name,bam_entry.reference_start,bam_entry.reference_end   
    qname,cigar = bam_entry.qname,bam_entry.cigar
    mapq = bam_entry.mapping_quality
    
    ## find qstart,qend from cigar
    tmp_qstart = 0
    qend = None
    # find start
    if cigar[0][0] in (4,5): #check if hard (5) or soft-clipped (4), else we start at coords 0
        tmp_qstart = cigar[0][1]
        
    # find end, by summing cigarstring - ignoring soft/hard-clips and DEL (2) made on query by reference
    cigar_sum = 0
    cigar_add = 0 #add soft/hard-clipping to get read length
    for i,j in cigar:
        if not i in (4,5,2):
            cigar_sum += j
        elif not i == 2:
            cigar_add += j

    # Compute read length from cigar. needed for fixing reverse mapped qstart
    read_len = cigar_sum + cigar_add
    
    # Check if entry is reverse and fix poses accordingly...
    strand = None
    if bam_entry.is_reverse:
        strand = '-'
        qend = read_len - tmp_qstart
        qstart = qend - cigar_sum
    else:
        strand = '+'
        qstart = tmp_qstart
        qend = tmp_qstart + cigar_sum
        
    data = {'qcoords':[qstart,qend],'rcoords':[rstart,rend],
            'rname':rname,'qname':qname,'read_len':read_len,
            'strand':strand,'mapq':mapq}
        
    return data
    ##/

def parse_reads_metadata(bam_handle,read_len_thresh=0):
    reads_metadata = {}
    aln_idx_map = {}
    for bam_entry in bam_handle.fetch():
        qname = bam_entry.qname
        if qname in reads_metadata: continue # skip qname if we already imported metadata for it
        if not bam_entry.has_tag('SA'): continue # skip if  no SA tag is available
        
        read_length = bam_entry.infer_read_length()
        if read_length == None: continue #skip if no cigar was output
        if read_length < read_len_thresh: continue #skip if read was too short
        
        aln_regions = []
        
        # Add primary alignment
        entry = parse_bam_entry(bam_entry)
        entry['idx'] = len(aln_idx_map)
        aln_regions.append(entry)
        aln_idx_map[entry['idx']] = entry
        #/
        
        # get alignment region and other regions from SA tag
        SA_tags = bam_entry.get_tag('SA').split(';')
        for SA_tag in SA_tags:
            if not SA_tag: continue
            SA_tag = SA_tag.split(',')
            
            # Parse SA tag
            rname = SA_tag[0]
            rstart = int(SA_tag[1])-1 # 1-based by SAM specification
            strand = SA_tag[2]
            mapq = int(SA_tag[4])
            
            letters_count = []
            number = ''
            for letter in SA_tag[3]:
                if letter.isdigit():
                    number += letter
                else:
                    letters_count.append([letter,int(number)])
                    number = ''
            
            clip1, clip2 = 0, 0
            if letters_count[0][0] == 'S':               clip1 = letters_count[0][1]
            if letters_count[-1][0] == 'S':              clip2 = letters_count[-1][1]
            
            if SA_tag[2] == '+':
                qstart = clip1
                qend = read_length - clip2
            else:
                qstart = clip2
                qend = read_length - clip1
                
            alen = 0
            for letter,count in letters_count:
                if not letter in ('S','I',):
                    alen += count
            rend = rstart + alen
            #/
            
            entry = {'qname':qname,'rname':rname,'rcoords':[rstart,rend],'qcoords':[qstart,qend],'strand':strand,'mapq':mapq}
            entry['idx'] = len(aln_idx_map)
            aln_regions.append(entry)
            aln_idx_map[entry['idx']] = entry
        #/
        # Save
        reads_metadata[qname] = [read_length,aln_regions]
        #/
        
    aln_idx_map_ROMD = rangeOverlaps_makeMappingDict(aln_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
    return reads_metadata,aln_idx_map,aln_idx_map_ROMD

def read_list_parser(input_file):
    reads = set()
    with open(input_file,'r') as f:
        for line in f:
            line = line.strip('\n')
            line = line.split()
            reads.add(line[0])
    return reads
          
def plot_on_ax(ax,fig=None,bam_handle='',read='',read_metadata='',alen_thresh=None,mapq_thresh=None,masking_thresh=None,mark_q_ovlps=True,SVs=None,
               GTFs=None,GTF_show_label=None,GTF_minLen=0,GTF_label_minLen=0,
               showScale=None,highlight_region=None,additional_reads_metadata={},
               alns_chain_distance=20000,region_margin=10000):
    """
    Hold plotting function separate, so we can use it in both GUI and command-line export
    """
    
    def filterAlns(read_alns):
        if alen_thresh != None or mapq_thresh != None or masking_thresh != None:
            read_alns_filt = []
            for aln in read_alns:
                if alen_thresh != None and aln['rcoords'][-1]-aln['rcoords'][0] < alen_thresh: continue
                if mapq_thresh != None and aln['mapq'] < mapq_thresh: continue
                if masking_thresh != None and masking_thresh >= 0 and 'masking' in aln and aln['masking'][2] > masking_thresh: continue
                read_alns_filt.append(aln)
        return read_alns_filt
    
    if bam_handle and read_metadata:
        read_alns = read_metadata[1]
        
        # Check if filt alns
        if alen_thresh != None or mapq_thresh != None or masking_thresh != None:
            """
            read_alns_filt = []
            for aln in read_alns:
                if alen_thresh != None and aln['rcoords'][-1]-aln['rcoords'][0] < alen_thresh: continue
                if mapq_thresh != None and aln['mapq'] < mapq_thresh: continue
                if masking_thresh != None and masking_thresh >= 0 and 'masking' in aln and aln['masking'][2] > masking_thresh: continue
                read_alns_filt.append(aln)
            read_alns = read_alns_filt
            """
            read_alns = filterAlns(read_alns)
        #/
        
        # Check if look for ovlps on query: For each alignment, mark qcoords which ovlp to other alignments.
        if mark_q_ovlps != None and mark_q_ovlps:
            qRanges = []
            for aln in read_alns:
                qRanges.append([aln['qcoords'][0],aln['idx'],aln['qcoords'][-1]])
            ovlps = computeOvlpRanges_wIdxs(qRanges)
            
            # 1. for each alignment, store qcoords which ovlp
            # 2. chain overlapping alignments (gonna plot chained alignments by individual colors)
            alns_qcoords_ovlp_raw = {}
            ovlps_chained = []
            for rangee in ovlps:
                # check if more than one alignment exist at current query coordinates
                if len(rangee[1]) > 1:
                    # 1.
                    # Check if chain to previous 
                    if ovlps_chained and ovlps_chained[-1][-1] == rangee[0]:
                        ovlps_chained[-1][-1] = rangee[-1] #update end qcoord
                        ovlps_chained[-1][1].update(rangee[1]) #update alignment enumerates at chain
                    # Else initiate new chain
                    else:
                        ovlps_chained.append([rangee[0],set(rangee[1]),rangee[-1]])
                    #
                    #/1
                    # 2
                    for idx in rangee[1]:
                        chain_enum = len(ovlps_chained)
                        init(alns_qcoords_ovlp_raw,idx,{})
                        init(alns_qcoords_ovlp_raw[idx],chain_enum,[])
                        alns_qcoords_ovlp_raw[idx][chain_enum].append(rangee[0])
                        alns_qcoords_ovlp_raw[idx][chain_enum].append(rangee[-1])
                    #/
            
            # Assign alignments with chains
            alns_ovlpChains = {}
            for chain_enum,chain in enumerate(ovlps_chained):
                for aln_idx in chain[1]:
                    alns_ovlpChains[aln_idx] = chain_enum
            #/
            # Calc min/max qcoords ovlp to each chain
            alns_qcoords_ovlp = {}
            for aln_idx in alns_qcoords_ovlp_raw:
                for chain_enum in alns_qcoords_ovlp_raw[aln_idx]:
                    qcoords_raw = alns_qcoords_ovlp_raw[aln_idx][chain_enum]
                    qcoords = [min(qcoords_raw),max(qcoords_raw)]
                    init(alns_qcoords_ovlp,aln_idx,{})
                    alns_qcoords_ovlp[aln_idx][chain_enum] = qcoords
            #/
        #/
        
        # Compute regions to plot
        aChroms,alns_aChromEnums = calc_aChroms_from_alignments(read_alns,alns_chain_distance=alns_chain_distance,region_margin=region_margin,lenReq=0)
        
        aChroms_totLen = 0
        for aChrom in aChroms.values():
            aChroms_totLen += aChrom['len']
        
        offset_x = 0
        for aChrom in aChroms.values():
            aChrom['x_offset'] = offset_x
            offset_x += aChrom['len'] + aChroms_totLen*0.1
        #/
        
        
        y_offset_scale = len(alns_aChromEnums)*0.05
        y_offset_scale_minor = max(1,y_offset_scale*0.25)
        y_offset = 0
        
        # Draw features
        gtf_offset = 0
        if GTFs != None:
            gtf_offset = y_offset_scale*0.5
            GTF_idx_map = GTFs[0]
            GTF_idx_map_ROMD = GTFs[1]
            
            for enum,aChrom in aChroms.items():
                GTFs_to_paint = []
                for oS,idx,oE in rangeOverlaps_lookup(aChrom['rcoords'],GTF_idx_map_ROMD,1000,map_lookup_key=aChrom['rname']):
                    GTF = GTF_idx_map[idx]
                    if GTF_minLen and (GTF['rcoords'][-1]-GTF['rcoords'][0]) < GTF_minLen: continue # skip GTF if too short
                    GTFs_to_paint.append(GTF)
                
                # paint
                for GTF in GTFs_to_paint:
                    GTF_left_shrink = aChrom['rcoords'][0] - GTF['rcoords'][0] # positive vals if rstart is left of aChrom
                    if GTF_left_shrink < 0:     GTF_left_shrink = 0
                    GTF_right_shrink = GTF['rcoords'][-1] - aChrom['rcoords'][-1] # positive vals if rend is right of aChrom
                    if GTF_right_shrink < 0:    GTF_right_shrink = 0
                    
                    if GTF['rcoords'][0] > GTF['rcoords'][-1]: print('GTF start coords larger than GTF end coords!')
                    
                    aChrom_rstart = aChrom['rcoords'][0]
                    rstart_offset = (GTF['rcoords'][0] + abs(GTF_left_shrink)) - aChrom_rstart  # we start plot at GTF_rcoords + GTF shrink
                    alen = (GTF['rcoords'][-1] - GTF_right_shrink) - (GTF['rcoords'][0] + abs(GTF_left_shrink))
                    x_start = aChrom['x_offset'] + rstart_offset
                    x_end = x_start + alen
                    plotColor = 'red'
                    ax.plot([x_start,x_end],[y_offset,y_offset],color=plotColor)
                    
                    if GTF_show_label:
                        if GTF_label_minLen and (GTF['rcoords'][-1]-GTF['rcoords'][0]) < GTF_label_minLen: continue
                        ax.text(x_start,y_offset,GTF['name'])
                #/
        #/  
        
        # Draw ref regions
        #@ init colormap for highlight (genome overview)
        if highlight_region != None:
            colEnums = []
            for i in range(len(highlight_region)):     colEnums.append(i/(len(highlight_region)-1))
            colormap = plt.cm.rainbow(colEnums)
        #@/
        y_offset -= gtf_offset
        for enum,aChrom in aChroms.items():
            x_start,x_end = aChrom['x_offset'],aChrom['x_offset'] + aChrom['len']

            # Check if we want to hightlight regions
            if highlight_region != None:
                for hl_enum,hl_reg in enumerate(highlight_region):
                    if aChrom['rname'] == hl_reg[0] and getRangeOvlp(aChrom['rcoords'],hl_reg[1]) >= 0:
                        # pre-correct coords if out of aChrom boundaries
                        reg_coords = [hl_reg[1][0],hl_reg[1][1]]
                        if reg_coords[0] < aChrom['rcoords'][0]:       reg_coords[0] = aChrom['rcoords'][0]
                        if reg_coords[-1] > aChrom['rcoords'][-1]:     reg_coords[-1] = aChrom['rcoords'][-1]
                        
                        # Plot reg
                        aChrom_rstart = aChrom['rcoords'][0]
                        rstart_offset = reg_coords[0] - aChrom_rstart
                        paint_len = reg_coords[-1]-reg_coords[0]
                        reg_start = aChrom['x_offset'] + rstart_offset
                        reg_end = reg_start + paint_len
                        
                        if hl_enum == 0: # for genome overview. Always have the clicked region appear in fixed color
                            color = 'brown'
                        else:
                            color = colormap[hl_enum]
                        ax.plot([reg_start,reg_end],[y_offset,y_offset],color=color,linewidth=5)
                        #/
            #/
            
            ax.plot([x_start,x_end],[y_offset,y_offset],color='blue')
            
            if 0 and 'static text':
                ax.text(x_start,y_offset+abs(gtf_offset)+y_offset_scale+y_offset_scale*enum,aChrom['rname']+':'+'-'.join(map(str,aChrom['rcoords'])))
            
            # write text as draggable annotation
            rcoords_formatted = [f'{aChrom["rcoords"][0]:,}',f'{aChrom["rcoords"][-1]:,}']
            ann = ax.annotate(aChrom['rname']+':'+'-'.join(rcoords_formatted) , xy = [x_start,y_offset],
                              #arrowprops=dict(facecolor='black',arrowstyle='->',relpos=(0,0),linestyle='--',linewidth=0.5),
                              xytext = [x_start,y_offset+abs(gtf_offset)+y_offset_scale+y_offset_scale*enum],
                              bbox=dict(fc='None',ec='None'))
            ann.draggable()
            #/
        #/
        
        # Draw scale
        if showScale != None:
            y_offset -= y_offset_scale
            scale_to_plot = 0
            if aChroms_totLen < 100:        scale_to_plot = 10
            elif aChroms_totLen < 1000:     scale_to_plot = 100
            elif aChroms_totLen < 10000:    scale_to_plot = 1000
            elif aChroms_totLen < 100000:   scale_to_plot = 10000
            elif aChroms_totLen < 1000000:  scale_to_plot = 100000
            else:                           scale_to_plot = 1000000
            x_start = 0
            x_end = scale_to_plot
            ax.plot([x_start,x_end],[y_offset*1.3,y_offset*1.3],color='red')
            ax.text(x_start,y_offset,'Scale: '+ f"{scale_to_plot:,}"+' bp')
        #/
        
        # Draw alignments
        if 0 and 'draw_vanilla':
            y_offset -= y_offset_scale
            prev_aln_point = None
            for aln in sorted(read_alns,key=lambda x: x['qcoords'][0]):
                # Plot aln
                aChrom = aChroms[alns_aChromEnums[aln['idx']]]
                aChrom_rstart = aChrom['rcoords'][0]
                rstart_offset = aln['rcoords'][0] - aChrom_rstart
                alen = aln['rcoords'][-1] - aln['rcoords'][0]
                x_start = aChrom['x_offset'] + rstart_offset
                x_end = x_start + alen
                ax.plot([x_start,x_end],[y_offset,y_offset],color='black')
                #/
                # Plot alignment direction
                arrow_y = y_offset_scale*0.1
                arrow_x = aChroms_totLen*0.01
                if aln['strand'] in ('+','fv','fw',):
                    ax.plot([x_end-arrow_x,x_end],[y_offset+arrow_y,y_offset],color='black')
                else:
                    ax.plot([x_start+arrow_x,x_start],[y_offset+arrow_y,y_offset],color='black')
                #/
                # Check if paint query alignment overlaps
                if mark_q_ovlps != None and mark_q_ovlps:
                    if aln['idx'] in alns_qcoords_ovlp:
                        for chain_enum,qcoords in alns_qcoords_ovlp[aln['idx']].items():
                            if qcoords[-1]-qcoords[0] >= 10:
                                #print('minimum qovlp thresh for plot: 10')
                                
                                # Get plotting ratio of rcoords vs. qcoords (this is not perfect. it should be done via CIGAR traversal, but we dont have CIGAR and it takes much computation.)
                                # Unless there are big insertions/deletions encoded into the CIGARS, this should have very low impact on plotting accuracy
                                q_to_r_ratio = (aln['rcoords'][-1]-aln['rcoords'][0]) / (aln['qcoords'][-1]-aln['qcoords'][0])
                                
                                if aln['strand'] in ('+','fv','fw',):
                                    rstart_mark = x_start + int( (qcoords[0]-aln['qcoords'][0])*q_to_r_ratio )
                                    rend_mark = rstart_mark + int( (qcoords[-1]-qcoords[0])*q_to_r_ratio )
                                else:
                                    rend_mark = x_end - int( (qcoords[0]-aln['qcoords'][0])*q_to_r_ratio )
                                    rstart_mark = rend_mark - int( (qcoords[-1]-qcoords[0])*q_to_r_ratio )
    
                                ax.plot([rstart_mark,rend_mark],[y_offset,y_offset],color='red')
                #/
                
                # Connect current alignment to previous point
                if prev_aln_point != None:
                    if aln['strand'] in ('+','fv','fw',):
                        x1,x2 = prev_aln_point[0],x_start
                    else:
                        x1,x2 = prev_aln_point[0],x_end
                    y1,y2 = prev_aln_point[1],y_offset
                    ax.plot([x1,x2],[y1,y2],color='black',linestyle='--',linewidth=0.5)
                #/
                # Update previous point
                if aln['strand'] in ('+','fv','fw',):
                    prev_aln_point = (x_end,y_offset)
                else:
                    prev_aln_point = (x_start,y_offset)
                #/
                
                y_offset -= y_offset_scale_minor
            #/
        
        # Draw alignments
        y_offset -= y_offset_scale
        mainRead_aChroms_yoffsets = {} # aChrom_enum -> "set of y_offsets" (used if plotting additional reads)
        prev_aln_point = None
        for aln in sorted(read_alns,key=lambda x: x['qcoords'][0]):
            # Plot aln
            aChrom = aChroms[alns_aChromEnums[aln['idx']]]
            aChrom_rstart = aChrom['rcoords'][0]
            rstart_offset = aln['rcoords'][0] - aChrom_rstart
            alen = aln['rcoords'][-1] - aln['rcoords'][0]
            x_start = aChrom['x_offset'] + rstart_offset
            x_end = x_start + alen
            ax.plot([x_start,x_end],[y_offset,y_offset],color='black')
            #/
            # Plot alignment direction
            arrow_y = y_offset_scale*0.1
            arrow_x = aChroms_totLen*0.01
            if aln['strand'] in ('+','fv','fw',):
                ax.plot([x_end-arrow_x,x_end],[y_offset+arrow_y,y_offset],color='black')
            else:
                ax.plot([x_start+arrow_x,x_start],[y_offset+arrow_y,y_offset],color='black')
            #/
            # Check if paint query alignment overlaps
            if mark_q_ovlps != None and mark_q_ovlps:
                if aln['idx'] in alns_qcoords_ovlp:
                    for chain_enum,qcoords in alns_qcoords_ovlp[aln['idx']].items():
                        if qcoords[-1]-qcoords[0] >= 10:
                            #print('minimum qovlp thresh for plot: 10')
                            
                            # Get plotting ratio of rcoords vs. qcoords (this is not perfect. it should be done via CIGAR traversal, but we dont have CIGAR and it takes much computation.)
                            # Unless there are big insertions/deletions encoded into the CIGARS, this should have very low impact on plotting accuracy
                            q_to_r_ratio = (aln['rcoords'][-1]-aln['rcoords'][0]) / (aln['qcoords'][-1]-aln['qcoords'][0])
                            
                            if aln['strand'] in ('+','fv','fw',):
                                rstart_mark = x_start + int( (qcoords[0]-aln['qcoords'][0])*q_to_r_ratio )
                                rend_mark = rstart_mark + int( (qcoords[-1]-qcoords[0])*q_to_r_ratio )
                            else:
                                rend_mark = x_end - int( (qcoords[0]-aln['qcoords'][0])*q_to_r_ratio )
                                rstart_mark = rend_mark - int( (qcoords[-1]-qcoords[0])*q_to_r_ratio )

                            ax.plot([rstart_mark,rend_mark],[y_offset,y_offset],color='red')
            #/
            # Connect current alignment to previous point
            if prev_aln_point != None:
                if aln['strand'] in ('+','fv','fw',):
                    x1,x2 = prev_aln_point[0],x_start
                else:
                    x1,x2 = prev_aln_point[0],x_end
                y1,y2 = prev_aln_point[1],y_offset
                ax.plot([x1,x2],[y1,y2],color='black',linestyle='--',linewidth=0.5)
            #/
            # Update previous point
            if aln['strand'] in ('+','fv','fw',):
                prev_aln_point = (x_end,y_offset)
            else:
                prev_aln_point = (x_start,y_offset)
            #/
            # Save yoffset for plot at aChrom
            init(mainRead_aChroms_yoffsets,aChrom['enum'],set())
            mainRead_aChroms_yoffsets[aChrom['enum']].add(y_offset)
            #/
            y_offset -= y_offset_scale_minor
        
        mainRead_yoffset_end = y_offset # Keep track of y_offsets for read + additional reads
        #/
        
        # Check if draw more reads
        if additional_reads_metadata:
            oreads_yoffsets = [] # keep track of yoffsets for additional reads
            for oread_enum,(oread,oread_metadata) in enumerate(additional_reads_metadata.items()):
                oread_alns = oread_metadata[1]
                plotColor = oread_metadata[2]
                # Check if filt alns
                if alen_thresh != None or mapq_thresh != None or masking_thresh != None:
                    oread_alns = filterAlns(oread_alns)
                #/
                # Find which y-offset to start additional read plotting
                aChroms_mainRead_yoffsets = set()
                for aln in sorted(oread_alns,key=lambda x: x['qcoords'][0]):
                    # check which aChrom aln overlaps
                    for enum,aChrom in aChroms.items():
                        if aChrom['rname'] == aln['rname'] and getRangeOvlp(aChrom['rcoords'],aln['rcoords']) >= 0:
                            aChroms_mainRead_yoffsets.update(mainRead_aChroms_yoffsets[aChrom['enum']])
                if not aChroms_mainRead_yoffsets: continue #skip if no hits to mainread aChroms were found (these reads do not have overlapping alignments)
                y_offset = max(aChroms_mainRead_yoffsets) - (y_offset_scale_minor*((oread_enum+1) / (len(additional_reads_metadata)+1)))  # reset y-offset. Keep "one additional reads' margin"
                #/
                # Draw alns
                prev_aln_point = None
                for aln in sorted(oread_alns,key=lambda x: x['qcoords'][0]):
                    # check which aChrom aln overlaps
                    aln_plotted = False
                    for enum,aChrom in aChroms.items():
                        if aChrom['rname'] == aln['rname'] and getRangeOvlp(aChrom['rcoords'],aln['rcoords']) >= 0:
                            # pre-correct coords if out of aChrom boundaries
                            aln_coords = [aln['rcoords'][0],aln['rcoords'][-1]]
                            if aln_coords[0] < aChrom['rcoords'][0]:       aln_coords[0] = aChrom['rcoords'][0]
                            if aln_coords[-1] > aChrom['rcoords'][-1]:     aln_coords[-1] = aChrom['rcoords'][-1]
                            
                            # Plot aln
                            aln_plotted = True
                            aChrom_rstart = aChrom['rcoords'][0]
                            rstart_offset = aln_coords[0] - aChrom_rstart
                            paint_len = aln_coords[-1]-aln_coords[0]
                            x_start = aChrom['x_offset'] + rstart_offset
                            x_end = x_start + paint_len
                            ax.plot([x_start,x_end],[y_offset,y_offset],color=plotColor)
                            #/
                            
                            # Plot alignment direction
                            arrow_y = y_offset_scale*0.1
                            arrow_x = aChroms_totLen*0.01
                            if aln['strand'] in ('+','fv','fw',):
                                ax.plot([x_end-arrow_x,x_end],[y_offset+arrow_y,y_offset],color=plotColor)
                            else:
                                ax.plot([x_start+arrow_x,x_start],[y_offset+arrow_y,y_offset],color=plotColor)
                            #/
                    
                    if aln_plotted:
                        # Connect current alignment to previous point
                        if prev_aln_point != None:
                            if aln['strand'] in ('+','fv','fw',):
                                x1,x2 = prev_aln_point[0],x_start
                            else:
                                x1,x2 = prev_aln_point[0],x_end
                            y1,y2 = prev_aln_point[1],y_offset
                            ax.plot([x1,x2],[y1,y2],color=plotColor,linestyle='--',linewidth=0.5)
                        #/
                        # Update previous point
                        if aln['strand'] in ('+','fv','fw',):
                            prev_aln_point = (x_end,y_offset)
                        else:
                            prev_aln_point = (x_start,y_offset)
                        #/
                        
                        y_offset -= y_offset_scale_minor
                #/
                # Store y-offsets for additional plots
                oreads_yoffsets.append(y_offset)
                #/
            # Re-set y_offset
            y_offset = min([mainRead_yoffset_end]+oreads_yoffsets) #minimum since we plot on negative axis
            #/
        #/
        
        # Draw SVs
        def drawTriangle(ax,topX,topY,width=1,height=1):
            botLeft = [topX-width/2,topY-height]
            botRight = [topX+width/2,topY-height]
            ax.add_patch(plt.Polygon([botLeft,botRight,[topX,topY]],closed=True,color='green'))
            ax.plot(0,0) # need to have a dummy-plot for ax patch to actually show up :O tried plt.show() and ax.autosclae_view() as stackexchange comments suggested
        
        if SVs != None and SVs:
            for SV in sorted(SVs,key=lambda x: x['type']):
                SV_plot_x = []
                SV_plot_y = []
                for aChrom in aChroms.values():
                    x_plot,y_plot = None,None
                    # Check if SV overlap full region
                    if SV['rname'] == aChrom['rname'] and SV['rname'] == SV['rname2'] and \
                     getRangeOvlp([min([SV['rpos'],SV['rpos2']]),max([SV['rpos'],SV['rpos2']])],aChrom['rcoords']) >= abs(SV['rpos']-SV['rpos2']):
                        aChrom_rstart = aChrom['rcoords'][0]
                        rstart_offset = SV['rpos'] - aChrom_rstart
                        paint_len = abs(SV['rpos2'] - SV['rpos'])
                        x_start = aChrom['x_offset'] + rstart_offset
                        x_end = x_start + paint_len
                        x_plot = [x_start,x_end]
                        y_plot = [y_offset,y_offset]

                    # else check overlap to each breakpoint
                    elif SV['rname'] == aChrom['rname'] and getRangeOvlp([SV['rpos']]*2,aChrom['rcoords']) >= 0:
                        # Paint from SV start and extend to right
                        aChrom_rstart = aChrom['rcoords'][0]
                        rstart_offset = SV['rpos'] - aChrom_rstart
                        paint_len = aChrom['rcoords'][-1] - SV['rpos']
                        x_start = aChrom['x_offset'] + rstart_offset
                        x_end = x_start + paint_len
                        x_plot = [x_start,x_end]
                        y_plot = [y_offset,y_offset-(y_offset_scale_minor*0.6)]
                        
                    elif SV['rname2'] == aChrom['rname'] and getRangeOvlp([SV['rpos2']]*2,aChrom['rcoords']) >= 0:
                        # Paint from SV start and extend to left (rstart will be aChrom start, rend will be SV start)
                        aChrom_rstart = aChrom['rcoords'][0]
                        rstart_offset = 0 # aChrom start
                        paint_len = SV['rpos2']-aChrom['rcoords'][0]
                        x_start = aChrom['x_offset'] + rstart_offset
                        x_end = x_start + paint_len
                        x_plot = [x_start,x_end]
                        y_plot = [y_offset-(y_offset_scale_minor*0.6),y_offset]
                                  
                    if x_plot != None and y_plot != None:
                        SV_plot_x += x_plot
                        SV_plot_y += y_plot
                
                if SV_plot_x and SV_plot_y:
                    y_offset -= y_offset_scale_minor
                    
                    if SV['type'].find('INS') != -1:
                        drawTriangle(ax,SV_plot_x[0],SV_plot_y[0],width=aChroms_totLen*0.02,height=0.2)
                        ax.text(SV_plot_x[0],SV_plot_y[0],SV['type']+':'+str(SV['len']))
                    elif SV['type'].find('DEL') != -1:
                        ax.text(SV_plot_x[0],SV_plot_y[0],SV['type']+':'+str(SV['len']))
                        ax.plot(SV_plot_x,SV_plot_y,color='red')
                    elif SV['type'] in ('TRA','BND',) and SV['rname'] != SV['rname2']:
                        ax.text(SV_plot_x[0],max(SV_plot_y),SV['type'])
                        ax.plot(SV_plot_x,SV_plot_y,color='black')
                    else:
                        ax.text(SV_plot_x[0],max(SV_plot_y),SV['type']+':'+str(SV['len']))
                        ax.plot(SV_plot_x,SV_plot_y,color='black')
                    
            #/
        #/
    
    alignments_out = sorted(read_alns,key=lambda x: x['qcoords'][0])
    
    # Reposition chromosome labels
    if fig != None:
        fig.canvas.draw()
        # Get reflabel positions
        ann_objs_data = []
        for child in ax.get_children():
            if isinstance(child,matplotlib.text.Annotation):
                #for attr in dir(child):
                #    print(attr,getattr(child,attr))
                bbox_patch = child.get_bbox_patch()
                extents = bbox_patch.get_extents().transformed(ax.transData.inverted())
                width = extents.x1-extents.x0
                height = extents.y1-extents.y0
                x,y = child.xyann
                
                ann_objs_data.append([child,[x,y],width,height])
        #/
        # Only continue if we actually have something plotted...
        if ann_objs_data:
            # Get x/y starts
            y_starts = []
            x_starts = []
            y_heights = []
            for ann_obj in ann_objs_data:
                y_starts.append(ann_obj[1][1])
                x_starts.append(ann_obj[1][0])
                y_heights.append(ann_obj[3])
            y_start = min(y_starts)
            x_start = min(x_starts)
            y_height = max(y_heights)
            #/
            # Set all y coordinates to lowest Y
            for ann_obj in ann_objs_data:
                ann_obj[1][1] = y_start
            #/
            ## Sort by x-coordinate, set yoffset to 0, traverse: check if overlap any previous label on current yoffset. If it does, increase yoffset and repeat until no ovlp. Move annotation.
            # Sort by x coord
            ann_objs_data.sort(key=lambda x: x[1][0])
            #/
            # Calculate new y coords for annotations
            y_label_poses = set([y_start]) # for ylim, keep track of highest added yval to set plot limits. sometimes text appears on border of plot
            for ann_enum,ann_obj in enumerate(ann_objs_data):
                if ann_enum == 0: continue #skip first entry, it is already placed correctly
                cur_y = y_start
                cur_x = [ann_obj[1][0],ann_obj[1][0]+ann_obj[2]] # x_start = xstart, x_end = xstart+width of label
                while True:
                    # check if ovlp to existing
                    ovlp_existing = False
                    for ann_enum2,exist_ann_obj in enumerate(ann_objs_data):
                        if ann_enum == ann_enum2: continue #skip self
                        exist_y = exist_ann_obj[1][1]
                        exist_x = [exist_ann_obj[1][0],exist_ann_obj[1][0]+exist_ann_obj[2]] # x_start = xstart, x_end = xstart+width of label
                        
                        if getRangeOvlp([exist_y]*2,[cur_y]*2) >= 0: # check only on same y pos
                            if getRangeOvlp(exist_x,cur_x) >= 0: # check if ovlp on x
                                ovlp_existing = True
                                break
                    #/
                    # If ovlped existing, increase ycoord and traverse again
                    if ovlp_existing:
                        cur_y += y_height
                    else: # else, move ann to current y
                        ann_obj[1][1] = cur_y
                        break
                    #/
                y_label_poses.add(cur_y)
            #/
            # Move annotation objects
            for exist_ann_obj in ann_objs_data:
                exist_ann_obj[0].set_y(exist_ann_obj[1][1])
            #/
            # Set ylim as highest coordinate. Sometimes text appears on border of plot
            ax.set_ylim(top=max(y_label_poses))
            #/
            ##/
        #/
    #/
    return aChroms, alignments_out

def exportImage(fig_handle,output_loc):
    fig_handle.savefig(output_loc)

### Window for settings. Window is called inside main window and its input-form is set via the main window
class inputSettings_window(QWidget):
    window_closed = pyqtSignal()
    
    def __init__(self):
        super().__init__()
        self.input_form = []
        self.output_form = []
    
    def executeWindow(self):
        
        self.setFixedWidth(500)
        
        self.mainLayout = QVBoxLayout()
        
        
        self.setLayout(self.mainLayout)
        
        # header
        self.header = QLabel('Header')
        #self.header.setAlignment(Qt.AlignCenter)
        self.mainLayout.addWidget(self.header)
        #/
        # form
        self.formLayout = QFormLayout()
        
        self.output_form = []
        for enum,(label,edit) in enumerate(self.input_form.items()):
            self.label = QLabel(label)
            self.lineEdit = QLineEdit(str(edit))
            self.formLayout.setWidget(enum,QFormLayout.LabelRole, self.label)
            self.formLayout.setWidget(enum,QFormLayout.FieldRole, self.lineEdit)            
            self.output_form.append([self.label,self.lineEdit])
        
        self.mainLayout.addLayout(self.formLayout)
        #/
        # save
        self.buttonLayout = QHBoxLayout()
        self.mainLayout.addLayout(self.buttonLayout)
        
        self.buttonSave = QPushButton('Save')
        self.buttonSave.clicked.connect(self.close)
        self.buttonLayout.addWidget(self.buttonSave)
        
        self.buttonRestore = QPushButton('Restore defaults')
        self.buttonRestore.clicked.connect(self.restore_defaults)
        self.buttonLayout.addWidget(self.buttonRestore)
        #/
        
        self.show()
        
    def restore_defaults(self):
        if self.output_form:
            self.output_form[0][1].setText('RESTORE_DEFAULTS')
        self.close()
    
    def closeEvent(self, event):
        self.window_closed.emit()
        event.accept()
###/
        
### WINDOW FOR GENOME OVERVIEW
class genomeOverview_window(QMainWindow):
    bin_clicked_signal = pyqtSignal()
    def __init__(self):
        super(genomeOverview_window, self).__init__()
        uic.loadUi('genome_overview.ui',self)
        
        ## Set default input variables (these must be specified from main window code)
        self.bam_handle = None
        self.maskings = None
        self.reads_metadata = None
        ##/
        
    def executeWindow(self):
        self.setWindowTitle('Genome Overview')
        
        # Init ref lens (import them if we have bam handle)
        self.ref_lens = None
        self.chroms_show_default = '' # init chromosome select default
        if self.bam_handle != None:
            self.ref_lens = parse_bam_ref_lens(self.bam_handle)
            
            for rname,rlen in sorted(self.ref_lens.items(),key=lambda x: x[0]):
                self.chroms_show_default += rname+'\n'
        
        self.select_chromosomes_reset() # set current chromosome select to the default
        #/
        
        ### Actions
        ## Menu
        # File
        self.actionExit.triggered.connect(self.hide)
        #/
        # Export
        self.actionExport_image.triggered.connect(self.export_image)    
        #/
        # View
        self.actionReset_zoom.triggered.connect(lambda: self.zoom('reset'))
        self.actionZoom_in_both.triggered.connect(lambda: self.zoom('in',zoom_size=100))
        self.actionZoom_in_X_axis.triggered.connect(lambda: self.zoom('in',zoom_axis='x',zoom_size=100))
        self.actionZoom_in_Y_axis.triggered.connect(lambda: self.zoom('in',zoom_axis='y',zoom_size=100))
        self.actionZoom_out_both.triggered.connect(lambda: self.zoom('out',zoom_size=100))
        self.actionZoom_out_X_axis.triggered.connect(lambda: self.zoom('out',zoom_axis='x',zoom_size=100))
        self.actionZoom_out_Y_axis.triggered.connect(lambda: self.zoom('out',zoom_axis='y',zoom_size=100))
        #/
        ##/
        
        ## Buttons
        self.button_selectChroms.clicked.connect(self.select_chromosomes)
        self.button_selectChroms_reset.clicked.connect(self.select_chromosomes_reset)
        self.button_run.clicked.connect(self.run_genomeOverview)
        ##/
        ###/
        
        ### Add plot
        self.updatePlot = False
        self.qscrollLayout = QVBoxLayout(self.scrollContents)
        self.initPlot()
        ###/
        
        # Init variable
        self.genome_bins = None
        self.read_alns_to_bins = None
        self.prev_artist = None
        #/
        self.show()
        
    def select_chromosomes(self):
        if self.ref_lens != None:
            # Filt chroms_show based on user input chrom min len
            min_chrom_len = 0
            if self.input_minChromLen.text() and self.input_minChromLen.text().isdigit():               min_chrom_len = max(0,int(self.input_minChromLen.text())) # never allow this to be < 0
            chroms_show = ''
            for rname in self.chroms_show.split():
                if rname in self.ref_lens and self.ref_lens[rname] >= min_chrom_len:
                    chroms_show += rname+'\n'
            #/
            
            box_input_str,_ = QInputDialog.getMultiLineText(self,'Chromosome selection','Chromosomes:', chroms_show)
            
            if box_input_str:
                self.chroms_show = box_input_str
    
    def select_chromosomes_reset(self):
        # Get input minimum chromosome length
        min_chrom_len = 0
        if self.input_minChromLen.text() and self.input_minChromLen.text().isdigit():               min_chrom_len = max(0,int(self.input_minChromLen.text())) # never allow this to be < 0
        #/
        chroms_show = ''
        for rname in self.chroms_show_default.split():
            if rname in self.ref_lens and self.ref_lens[rname] >= min_chrom_len:
                chroms_show += rname+'\n'
        self.chroms_show = chroms_show
        
    def initPlot(self):
        if not self.updatePlot:
            self.fig = Figure(tight_layout=True, dpi=100)
            self.canvas = FigureCanvas(self.fig)
            self.canvas_connect = None
            self.ax = self.canvas.figure.add_subplot(111)
            self.ax.set_axis_off()
    
            container = QWidget()
            lay = QVBoxLayout(container)
            lay.addWidget(self.canvas)
    
            self.qscrollLayout.addWidget(container)
            
            container.setMinimumHeight(700) # sets plot size
        
            self.updatePlot = True
            QScroller.grabGesture(self.scrollArea.viewport(), QScroller.LeftMouseButtonGesture)
            
    def zoom(self,zoom_direction,zoom_axis=None,zoom_size=100):
        if self.updatePlot:
            plot_widget = self.qscrollLayout.itemAt(0).widget()
            height = plot_widget.height()
            width = plot_widget.width()
            
            v_pos_ratio = self.scrollArea.verticalScrollBar().value() / max(1,self.scrollArea.verticalScrollBar().maximum())
            h_pos_ratio = self.scrollArea.horizontalScrollBar().value() / max(1,self.scrollArea.horizontalScrollBar().maximum())
            
            if zoom_direction == 'in':
                if zoom_axis in ('y','Y',):
                    plot_widget.setMinimumHeight(height + zoom_size)
                elif zoom_axis in ('x','X',):
                    plot_widget.setMinimumWidth(width + zoom_size)
                else:
                    plot_widget.setMinimumHeight(height + zoom_size)
                    plot_widget.setMinimumWidth(width + zoom_size)
                    
            elif zoom_direction == 'out':
                if zoom_axis in ('y','Y',):
                    plot_widget.setMinimumHeight(height - zoom_size)
                elif zoom_axis in ('x','X',):
                    plot_widget.setMinimumWidth(width - zoom_size)
                else:
                    plot_widget.setMinimumHeight(height - zoom_size)
                    plot_widget.setMinimumWidth(width - zoom_size)
                    
            elif zoom_direction == 'reset':
                if zoom_axis in ('y','Y',):
                    plot_widget.setMinimumHeight(700)
                elif zoom_axis in ('x','X',):
                    plot_widget.setMinimumWidth(700)
                else:
                    plot_widget.setMinimumHeight(700)
                    plot_widget.setMinimumWidth(700)
            
            self.scrollArea.verticalScrollBar().setValue( round(v_pos_ratio*self.scrollArea.verticalScrollBar().maximum()) )
            self.scrollArea.horizontalScrollBar().setValue( round(h_pos_ratio*self.scrollArea.horizontalScrollBar().maximum()) )
    
    def keyPressEvent(self, event):
        # Zoom plot
        if self.updatePlot and event.key() == Qt.Key_Plus or event.key() == Qt.Key_Minus:
            if event.key() == Qt.Key_Plus:
                if (QApplication.keyboardModifiers() & Qt.ShiftModifier):
                    self.zoom('in',zoom_axis='x',zoom_size=100)
                elif (QApplication.keyboardModifiers() & Qt.AltModifier):
                    self.zoom('in',zoom_axis='y',zoom_size=100)
                else:
                    self.zoom('in',zoom_axis='x',zoom_size=100)
                    self.zoom('in',zoom_axis='y',zoom_size=100)
                    
            elif event.key() == Qt.Key_Minus:
                if (QApplication.keyboardModifiers() & Qt.ShiftModifier):
                    self.zoom('out',zoom_axis='x',zoom_size=100)
                elif (QApplication.keyboardModifiers() & Qt.AltModifier):
                    self.zoom('out',zoom_axis='y',zoom_size=100)
                else:
                    self.zoom('out',zoom_axis='x',zoom_size=100)
                    self.zoom('out',zoom_axis='y',zoom_size=100)
        #/
    
    ## Functions to call freeze to objects. These will be used to block user inputs and/or show that the software is loading
    def freezeUI(self,obj):
        obj.setEnabled(False)
        QApplication.processEvents()
        
    def unfreezeUI(self,obj):
        QApplication.processEvents()
        obj.setEnabled(True)
    ##/
    
    def run_genomeOverview(self):
        if not self.updatePlot:     self.initPlot()
        ax = self.ax
        ax.cla()
        ax.set_axis_off()
        
        ## Parse input filters
        min_chrom_len = 0
        bin_size = 100000 # Set an arbitrary bin size so a blank field wont cause lag out
        read_len_req = 0
        read_supp_req = 1
        alen_req = 0
        mapq_req = 0
        masking_thresh = 0
        if self.input_minChromLen.text() and self.input_minChromLen.text().isdigit():               min_chrom_len = max(0,int(self.input_minChromLen.text())) # never allow this to be < 0
        if self.input_binSize.text() and self.input_binSize.text().isdigit():                       bin_size = max(1,int(self.input_binSize.text())) # never allow this to be < 1
        if self.input_minReadLen.text() and self.input_minReadLen.text().isdigit():                 read_len_req = max(0,int(self.input_minReadLen.text())) # never allow this to be < 0
        if self.input_maxAlnRepmaFrac.text():
            try:        masking_thresh = max(0,float(self.input_maxAlnRepmaFrac.text())) # never allow this to be < 0
            except:     masking_thresh = 0
        if self.input_minReadSupp.text() and self.input_minReadSupp.text().isdigit():               read_supp_req = max(0,int(self.input_minReadSupp.text())) # never allow this to be < 0
        if self.input_minAlen.text() and self.input_minAlen.text().isdigit():                       alen_req = max(0,int(self.input_minAlen.text())) # never allow this to be < 0
        if self.input_minMapq.text() and self.input_minMapq.text().isdigit():                       mapq_req = max(0,int(self.input_minMapq.text())) # never allow this to be < 0
        ##/
        
        # freeze UI while running calculations
        self.freezeUI(self.centralwidget)
        self.freezeUI(self.button_run)
        #/
        
        if self.bam_handle != None:
            # Parse reference lengths from bam
            ref_lens = self.ref_lens
            #/
            # Repeatmask alignments
            if self.maskings != None:
                for read,metadata in self.reads_metadata.items():
                    if metadata[0] < read_len_req: continue #skip short reads
                    for aln in metadata[1]:
                        if 'masking' in aln: continue
                        aln['masking'] = calc_aln_masking(aln,self.maskings[0],self.maskings[1])
            #/
            
            reads_metadata = self.reads_metadata

            ### Aim: find breakpoints reflective of rearrangements on chromosome-scale.
            ### Steps:
            ### -Initiate bins on chromosomes of size X bp.
            ### -Traverse + filter alignments (remove repetitive, etc.)
            ### -Sort alignments by query-coordinates. Save alignment order enumerate
            ### -Find overlaps to chromosome bins. Assign read alignment enumerates to bin
            ### -Traverse chromosome bins. Compare to next neighbor: Check if alignment orders of reads are without gaps (i.e. if a gap, means that read has an alignment in another bin)
            ###     Example: Bin1 has alignment orders 3,4,5 and bin2 has alignment orders 6,7 (No rearrangement)
            ###     Example: Bin1 has alignment orders 3,4,5 and bin2 has alignment orders 7,8 (Rearrangement, alignment with order 6 maps elsewhere)
            ## Initiate bins
            # store object 1. sorted by reference name and location, 2. in index-map dict
            update_all = False # keep track on which components we need to update
            if self.genome_bins == None or self.genome_bins[0] != bin_size:
                update_all = True # If we updated genome bins, then update everything else
                rname_bins = {}
                bin_idx_map = {}
                for rname,rlen in ref_lens.items():
                    init(rname_bins,rname,[])    
                    # Keep appending bins until reference length is hit
                    bin_enum = 0
                    while not rname_bins[rname] or rname_bins[rname][-1]['rcoords'][-1] < ref_lens[rname]:
                        tmp_bin = {'rname':rname,'rcoords':[(bin_enum*bin_size),(bin_enum*bin_size + bin_size)],
                                   'idx':len(bin_idx_map)}
                        rname_bins[rname].append(tmp_bin)
                        bin_idx_map[tmp_bin['idx']] = tmp_bin
                        bin_enum += 1
                #/
                # Make mapping dict for finding overlaps to read alignments
                bin_idx_map_mapping_dict = rangeOverlaps_makeMappingDict(bin_idx_map,1000,coordsKey='rcoords',sortByKey='rname')
                #/
                # Save
                self.genome_bins = [bin_size, [rname_bins,bin_idx_map,bin_idx_map_mapping_dict]]
                #/
            else:
                # Load previous values
                rname_bins,bin_idx_map,bin_idx_map_mapping_dict = self.genome_bins[1]
                #/
            ##/
            ## Traverse and filter which alignments to use. Assign selected alignments to bins. Find rearrangements between bins from alignments.
            if update_all or self.read_alns_to_bins == None or self.read_alns_to_bins[0] != [read_len_req,alen_req,mapq_req,masking_thresh]:
                # Init (or clear previous) alignment assignments to bins
                for binn in bin_idx_map.values():
                    binn['reads_aln_enums'] = {}
                #/
                
                # Find which alignments to use and assign to bins
                reads_alnEnums_bins = {}
                reads_alnEnums_entries = {}
                for read,metadata in reads_metadata.items():
                    # Filt alignments
                    alns_sorted_filt = []
                    for aln in sorted(metadata[1],key=lambda x: x['qcoords'][0]):
                        if (aln['rcoords'][-1]-aln['rcoords'][0]) < alen_req: continue
                        if aln['mapq'] < mapq_req: continue
                        if 'masking' in aln and aln['masking'][2] > masking_thresh: continue
                        alns_sorted_filt.append(aln)
                    #/
                    # Find overlapping bins
                    for aln_enum,aln in enumerate(alns_sorted_filt):
                        init(reads_alnEnums_entries,read,{})
                        reads_alnEnums_entries[read][aln_enum] = aln
                        for oS,bin_idx,oE in rangeOverlaps_lookup(aln['rcoords'],bin_idx_map_mapping_dict,1000,map_lookup_key=aln['rname']):
                            binE = bin_idx_map[bin_idx]
                            init(binE['reads_aln_enums'],read,[])
                            binE['reads_aln_enums'][read].append(aln_enum)
                            
                            init(reads_alnEnums_bins,read,{})
                            init(reads_alnEnums_bins[read],aln_enum,set())
                            reads_alnEnums_bins[read][aln_enum].add(binE['idx'])
                    #/
                #/
                # Traverse read alignments in-order and find alnEnums which do not have adjacent binEnums in pairwise alignment comparison.
                bins_rearr_reads = {} # bin comp_ID -> reads
                for read,read_alnEnums_bins in reads_alnEnums_bins.items():
                    if reads_metadata[read][0] < read_len_req: continue #skip short reads
                    if len(read_alnEnums_bins) < 2: continue #skip if not at least x2 alns
                    for comp_enum,alnEnum in enumerate(read_alnEnums_bins):
                        if comp_enum == 0: continue #backwards checking
                        prev_aln_binIdxs = read_alnEnums_bins[alnEnum-1]
                        cur_aln_binIdxs = read_alnEnums_bins[alnEnum]
                        
                        prev_aln_binIdxs_minMax = [min(prev_aln_binIdxs),max(prev_aln_binIdxs)]
                        cur_aln_binIdxs_minMax = [min(cur_aln_binIdxs),max(cur_aln_binIdxs)]
                        
                        if getRangeOvlp(prev_aln_binIdxs_minMax,cur_aln_binIdxs_minMax) < -1: # more than 1 bin apart
                            # get connection breakpoint
                            BP1_aln = reads_alnEnums_entries[read][alnEnum-1]
                            BP2_aln = reads_alnEnums_entries[read][alnEnum]
                            BP1_coord = None
                            BP2_coord = None
                            if BP1_aln['strand'] == '+':        BP1_coord = BP1_aln['rcoords'][0]
                            elif BP1_aln['strand'] == '-':      BP1_coord = BP1_aln['rcoords'][-1]
                            if BP2_aln['strand'] == '+':        BP2_coord = BP2_aln['rcoords'][0]
                            elif BP2_aln['strand'] == '-':      BP2_coord = BP2_aln['rcoords'][-1]
                            #/
                            # get connection bins
                            BP1_bin = None
                            BP2_bin = None
                            for oS,bin_idx,oE in rangeOverlaps_lookup([BP1_coord]*2,bin_idx_map_mapping_dict,1000,map_lookup_key=BP1_aln['rname']):
                                BP1_bin = bin_idx_map[bin_idx]
                                break
                            for oS,bin_idx,oE in rangeOverlaps_lookup([BP2_coord]*2,bin_idx_map_mapping_dict,1000,map_lookup_key=BP2_aln['rname']):
                                BP2_bin = bin_idx_map[bin_idx]
                                break    
                            #/
                            # save
                            bins_ID = '||'.join(map(str,sorted([BP1_bin['idx'],BP2_bin['idx']])))
                            
                            init(bins_rearr_reads,bins_ID,{})
                            init(bins_rearr_reads[bins_ID],read,[])
                            bins_rearr_reads[bins_ID][read].append([BP1_aln,BP2_aln])
                        
                        """
                        elif getRangeOvlp(prev_aln_binIdxs_minMax,cur_aln_binIdxs_minMax) == -1:
                            '1 bin apart'
                        elif getRangeOvlp(prev_aln_binIdxs_minMax,cur_aln_binIdxs_minMax) == 0:
                            '1 bin overlap'
                        elif getRangeOvlp(prev_aln_binIdxs_minMax,cur_aln_binIdxs_minMax) > 0:
                            'more than 1 bin overlap'
                        """
                #/
                # Save
                self.read_alns_to_bins = [[read_len_req,alen_req,mapq_req,masking_thresh], [reads_alnEnums_bins, reads_alnEnums_entries,bins_rearr_reads]]
                #/
            else:
                # Load previous values
                reads_alnEnums_bins, reads_alnEnums_entries, bins_rearr_reads = self.read_alns_to_bins[1]
                #/
            ##/
            
            ## Compile which bins to mark based on read support
            SK = 'rearr' # "save-key"
            # Init (or clear previous) rearrangement marker to bins
            for binn in bin_idx_map.values():
                binn[SK] = set()
            #/
            # Find bins with rearrangements supported by X amount of reads
            for bins_ID,reads in bins_rearr_reads.items():
                if len(reads) < read_supp_req: continue
                binA_idx,binB_idx = map(int,bins_ID.split('||'))
                binA,binB = bin_idx_map[binA_idx],bin_idx_map[binB_idx]

                init(binA,SK,set())
                init(binB,SK,set())
                binA[SK].add(bins_ID)
                binB[SK].add(bins_ID)
            #/
            ###/
            ### PLOTTING STAGE
            y_offset = 0
            yoffsets_rnames = {}
            #for rname,bins in rname_bins.items():
            for rname in self.chroms_show.split(): # plot chromosomes in order of appearance
                if not rname in rname_bins: continue
                bins = rname_bins[rname]
                #if not rname in ('2L','2R','3L','3R','X',): continue
                # Skip chromosome if too short
                if min_chrom_len and ref_lens[rname] < min_chrom_len: continue
                #if not rname == '2L': continue
            
                ax.text(0,y_offset-0.5,rname)
                yoffsets_rnames[y_offset] = rname

                # Paint chromosome backbone
                x_start = 0
                x_end = (ref_lens[rname] / bin_size)+1
                ax.plot([x_start,x_end],[y_offset,y_offset],alpha=0)
                ax.plot([x_start,x_end],[y_offset+1,y_offset+1],alpha=0)
                ax.fill_between([x_start,x_end],[y_offset,y_offset],[y_offset+1,y_offset+1],facecolor='blue')
                #/
                # Paint rearranged bins

                for binn in bins:
                    if 'rearr' in binn and binn['rearr']: # paint rearr bin with clickable
                        x_start,x_end = binn['rcoords'][0] / bin_size , binn['rcoords'][-1] / bin_size
                        
                        # paint between previous and current line (using fill_between)    
                        ax.plot([x_start,x_end],[y_offset,y_offset],alpha=0)
                        ax.plot([x_start,x_end],[y_offset+1,y_offset+1],alpha=0)
                        ax.fill_between([x_start,x_end],[y_offset,y_offset],[y_offset+1,y_offset+1],facecolor='red',picker=True,pickradius=0)

                y_offset -= 2
                
            def mpl_click(event):
                if event.mouseevent.button == 1: # respond to left-clicks only
                    #print(event.artist.__dict__)
                    #print(type(event.artist))
                    fill_color = event.artist._original_facecolor
                    if fill_color != 'blue': # blue color will be "non-rearranged" bins. make them unclickable
                        fill_paths = event.artist._paths[0]
                        
                        x_coords = []
                        y_coords = []
                        for x,y in fill_paths.vertices:
                            x_coords.append(x)
                            y_coords.append(y)
                        
                        x_start,x_end = int(min(x_coords)),int(max(x_coords))
                        y_start,y_end = int(min(y_coords)),int(max(y_coords))
                        
                        # Highlight selected rectangle
                        if self.prev_artist != None: # reset previous selection
                            self.prev_artist.set_facecolor('red')
                        
                        event.artist.set_facecolor('brown')
                        self.prev_artist = event.artist
                        self.fig.canvas.draw_idle()
                        #/

                        rname = yoffsets_rnames[y_start]    
                        binn = rname_bins[rname][x_start]
                        
                        binn_rearr_reads = set()
                        rearrBins_reads = {}
                        obins_coords = []
                        if 'rearr' in binn:
                            for bins_ID in binn['rearr']:
                                binn_rearr_reads.update(bins_rearr_reads[bins_ID])
                                
                                init(rearrBins_reads,bins_ID,set())
                                rearrBins_reads[bins_ID].update(bins_rearr_reads[bins_ID])
                                
                                # Get coords of other opposite bin
                                for binn_idx in list(map(int,bins_ID.split('||'))):
                                    if not binn_idx == binn['idx']:
                                        obins_coords.append([bin_idx_map[binn_idx]['rname'],bin_idx_map[binn_idx]['rcoords']])
                                
        
                        self.read_list = binn_rearr_reads
                        self.read_list_grouped = rearrBins_reads
                        self.highlight_region = [ [binn['rname'],binn['rcoords']] ]
                        if obins_coords:        self.highlight_region += obins_coords # if we had obins, add them as highlight entries
        
                        # Trigger data transfer to main window
                        self.trigger_bin_click()
                        #/
            
            if self.canvas_connect != None: self.canvas_connect = self.canvas.mpl_disconnect(self.canvas_connect) # Check if we have previous connect to the canvas for clicks. Then remove it
            self.canvas_connect = self.canvas.mpl_connect('pick_event', mpl_click)
            self.fig.canvas.draw_idle()
            ###/
            #/
        # unfreeze UI when calculations are done
        self.unfreezeUI(self.centralwidget)
        self.unfreezeUI(self.button_run)
        #/
        
    def trigger_bin_click(self):
        self.bin_clicked_signal.emit()
        return
    
    def export_image(self):
        export_path = QFileDialog.getSaveFileName(self, 'Save file', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not export_path: return
        export_path = export_path[0]
        try:
            if self.fig:
                exportImage(self.fig,export_path)
        except:
            self.messageBox('Could not write output file! Make sure you included file extension and have permissions to write in output folder.\nSupported extensions include *.png, *.svg, *.pdf')
            return
###/

class ApplicationWindow(QMainWindow):
    def messageBox(self,msg):
        # Make a simple message box
        message = QMessageBox()
        message.setText(msg)
        message.exec_()
        return True
    
    def __init__(self):
        super(ApplicationWindow, self).__init__()
        uic.loadUi('read_painter.ui',self)
        self.show()
        self.activateWindow()
        
        ### Add actions
        ## Menu
        # File
        self.actionImport_BAM.triggered.connect(self.import_bam)
        self.bam_handle = None
        self.actionImport_GTF.triggered.connect(self.import_gtf)
        self.features = None
        self.features_idxDict_ROMD = [None,None]
        self.actionImport_Repeats.triggered.connect(self.import_repeats)
        self.maskings = None
        self.actionImport_VCF.triggered.connect(self.import_vcf)
        self.SVs = None
        self.actionExit.triggered.connect(self.close)
        #/
        # Select
        self.actionInput_read.triggered.connect(self.import_read_name)
        self.actionImport_read_list.triggered.connect(self.import_read_list)
        self.read_list = None
        self.read_list_grouped = None # set in genomeOverview atm
        self.actionClear_selection.triggered.connect(self.clear_selection)
        self.additional_reads_metadata = {}
        #/
        # Export
        self.actionExport_current_read.triggered.connect(self.export_current_read)
        self.actionExport_read_highlight.triggered.connect(self.export_read_highlight)
        self.actionExport_read_selection.triggered.connect(self.export_read_selection)
        self.aChroms = None
        self.actionExport_regions.triggered.connect(self.export_regions)
        self.alignments = None
        self.actionExport_alignments.triggered.connect(self.export_alignments)
        
        self.actionExport_image.triggered.connect(self.export_image)
        #/
        # View
        self.actionReset_zoom.triggered.connect(lambda: self.zoom('reset'))
        self.actionZoom_in_both.triggered.connect(lambda: self.zoom('in',zoom_size=100))
        self.actionZoom_in_X_axis.triggered.connect(lambda: self.zoom('in',zoom_axis='x',zoom_size=100))
        self.actionZoom_in_Y_axis.triggered.connect(lambda: self.zoom('in',zoom_axis='y',zoom_size=100))
        self.actionZoom_out_both.triggered.connect(lambda: self.zoom('out',zoom_size=100))
        self.actionZoom_out_X_axis.triggered.connect(lambda: self.zoom('out',zoom_axis='x',zoom_size=100))
        self.actionZoom_out_Y_axis.triggered.connect(lambda: self.zoom('out',zoom_axis='y',zoom_size=100))
        #/
        # Settings (default values)
        self.actionShow_settings.triggered.connect(self.show_settings)
        self.show_settings = {'Alignment chain distance':20000,'Region margin':10000}
        self.show_settings_default = {'Alignment chain distance':20000,'Region margin':10000}
        self.actionBAM_import_settings.triggered.connect(self.bam_import_settings)
        self.bam_settings = {'Minimum read length':1000}
        self.bam_settings_default = {'Minimum read length':1000}
        self.actionVCF_import_settings.triggered.connect(self.vcf_import_settings)
        self.vcf_settings = {'Chr2':'CHR2','Pos2':'END','SV type':'SVTYPE','Supportive reads':'RNAMES','SV length':'SVLEN'}
        self.vcf_settings_default = {'Chr2':'CHR2','Pos2':'END','SV type':'SVTYPE','Supportive reads':'RNAMES','SV length':'SVLEN'}
        self.actionGTF_GFF_import_settings.triggered.connect(self.gtf_gff_import_settings)
        self.gtf_settings = {'GTF Attribute delimiter':' ','GTF Gene tag':'gene_id','GFF Attribute delimiter':'=','GFF Gene tag':'ID','Gene identifier (column 3)':'gene,transcript','Ignore types':'exon,CDS'}
        self.gtf_settings_default = {'GTF Attribute delimiter':' ','GTF Gene tag':'gene_id','GFF Attribute delimiter':'=','GFF Gene tag':'ID','Gene identifier (column 3)':'gene,transcript','Ignore types':'exon,CDS'}
        #/
        # Genome overview
        self.actionRun_genomeOverview.triggered.connect(self.genome_overview)
        self.genomeOverview_window = None
        self.highlight_region = None
        #/
        ##/
        
        # Update plot
        self.button_updatePlot.clicked.connect(self.update_plot)
        #/
        
        # Read search
        self.button_findReads.clicked.connect(self.search_reads)
        #/
        
        # Read list
        self.reads_metadata = False
        self.aln_idx_map = None
        self.aln_idx_map_ROMD = None
        self.listWidget.itemClicked.connect(self.listwidgetclicked)
        self.listWidget.installEventFilter(self)
        self.current_list_item = None
        self.listwidget_defaultColor = None
        #/
        ###/
        
        ### Add plot
        self.updatePlot = False
        self.qscrollLayout = QVBoxLayout(self.scrollContents)
        self.plot()
        ###/
    
    def listwidgetclicked(self,item):
        try:        read = item.text().split(') ')[1] # split away read length
        except:     return
        
        if self.current_list_item != None and self.listwidget_defaultColor != None:          self.current_list_item.setBackground(self.listwidget_defaultColor) # Remove previous highlight
        
        if self.listWidget.count() < 1000: # Only attempt to change color if below X reads are plotted. else UI becomes slow.
            item.setBackground( QColor('#D2B48C') ) # Set highlight
        # If clicked read is in highlight also, then remove it: else coloring becomes messy when playing around in the GUI
        if self.additional_reads_metadata != None and read in self.additional_reads_metadata:  del self.additional_reads_metadata[read]
        #/
        self.current_list_item = item
        if read in self.reads_metadata:
            self.plot(item.text())
        
    def listwidgetclicked_right(self,item):
        try:        read = item.text().split(') ')[1] # split away read length
        except:     return
        #/
        # Skip if read is the focus read (if this read is already in focus, then remove it. else coloring becomes messy when playing around in the GUI)
        if self.current_list_item != None and self.current_list_item.text().split(') ')[1] == read:
            if self.additional_reads_metadata != None and read in self.additional_reads_metadata:  del self.additional_reads_metadata[read]
            return
        #/
        # Add or remove read from additional-reads-to-plot
        if self.additional_reads_metadata != None and read in self.additional_reads_metadata: # check if remove
            del self.additional_reads_metadata[read]
            item.setBackground(self.listwidget_defaultColor) # remove highlight
                                       
        else: # else add
            if read in self.reads_metadata:
                if self.additional_reads_metadata == None:      self.additional_reads_metadata = {} # check if initiate dictionary
                
                # set color for plot
                numCols = 7 # Cycle between X colors
                colEnums = []
                for i in range(numCols):     colEnums.append(i/(numCols-1))
                colormap = plt.cm.rainbow(colEnums)
                #/
                #plotColor = colormap[ len(self.additional_reads_metadata) % (numCols-1) ]
                plotColor = colormap[random.randint(0,6)]
                #/
                self.additional_reads_metadata[read] = [self.reads_metadata[read][0],self.reads_metadata[read][1],plotColor]
                item.setBackground( QColor('#ABF1BC') ) # set highlight
        #/
        # Update plot
        self.update_plot()
        #/
    
    def zoom(self,zoom_direction,zoom_axis=None,zoom_size=100):
        if self.updatePlot:
            plot_widget = self.qscrollLayout.itemAt(0).widget()
            height = plot_widget.height()
            width = plot_widget.width()
            
            v_pos_ratio = self.scrollArea.verticalScrollBar().value() / max(1,self.scrollArea.verticalScrollBar().maximum())
            h_pos_ratio = self.scrollArea.horizontalScrollBar().value() / max(1,self.scrollArea.horizontalScrollBar().maximum())
            
            if zoom_direction == 'in':
                if zoom_axis in ('y','Y',):
                    plot_widget.setMinimumHeight(height + zoom_size)
                elif zoom_axis in ('x','X',):
                    plot_widget.setMinimumWidth(width + zoom_size)
                else:
                    plot_widget.setMinimumHeight(height + zoom_size)
                    plot_widget.setMinimumWidth(width + zoom_size)
                    
            elif zoom_direction == 'out':
                if zoom_axis in ('y','Y',):
                    plot_widget.setMinimumHeight(height - zoom_size)
                elif zoom_axis in ('x','X',):
                    plot_widget.setMinimumWidth(width - zoom_size)
                else:
                    plot_widget.setMinimumHeight(height - zoom_size)
                    plot_widget.setMinimumWidth(width - zoom_size)
                    
            elif zoom_direction == 'reset':
                if zoom_axis in ('y','Y',):
                    plot_widget.setMinimumHeight(700)
                elif zoom_axis in ('x','X',):
                    plot_widget.setMinimumWidth(700)
                else:
                    plot_widget.setMinimumHeight(700)
                    plot_widget.setMinimumWidth(700)
            
            self.scrollArea.verticalScrollBar().setValue( round(v_pos_ratio*self.scrollArea.verticalScrollBar().maximum()) )
            self.scrollArea.horizontalScrollBar().setValue( round(h_pos_ratio*self.scrollArea.horizontalScrollBar().maximum()) )
    
    def keyPressEvent(self, event):
        # Zoom plot
        if self.updatePlot and event.key() == Qt.Key_Plus or event.key() == Qt.Key_Minus:
            if event.key() == Qt.Key_Plus:
                if (QApplication.keyboardModifiers() & Qt.ShiftModifier):
                    self.zoom('in',zoom_axis='x',zoom_size=100)
                elif (QApplication.keyboardModifiers() & Qt.AltModifier):
                    self.zoom('in',zoom_axis='y',zoom_size=100)
                else:
                    self.zoom('in',zoom_axis='x',zoom_size=100)
                    self.zoom('in',zoom_axis='y',zoom_size=100)
                    
            elif event.key() == Qt.Key_Minus:
                if (QApplication.keyboardModifiers() & Qt.ShiftModifier):
                    self.zoom('out',zoom_axis='x',zoom_size=100)
                elif (QApplication.keyboardModifiers() & Qt.AltModifier):
                    self.zoom('out',zoom_axis='y',zoom_size=100)
                else:
                    self.zoom('out',zoom_axis='x',zoom_size=100)
                    self.zoom('out',zoom_axis='y',zoom_size=100)
        #/
        
    def eventFilter(self,source,event):
        # Check if right-click on item in readlist
        if source is self.listWidget and event.type() == QEvent.ContextMenu:
            try:        self.listwidgetclicked_right(source.itemAt(event.pos()))
            except:     pass
            return True
        #/
        return super(ApplicationWindow,self).eventFilter(source,event)
            

    def initPlot(self):
        if not self.updatePlot:
            self.fig = Figure(tight_layout=True, dpi=100)
            self.canvas = FigureCanvas(self.fig)
            self.ax = self.canvas.figure.add_subplot(111)
            self.ax.set_axis_off()

            container = QWidget()
            lay = QVBoxLayout(container)
            lay.addWidget(self.canvas)
    
            self.qscrollLayout.addWidget(container)

            container.setMinimumHeight(700) # sets plot size
        
            self.updatePlot = True
            QScroller.grabGesture(self.scrollArea.viewport(), QScroller.LeftMouseButtonGesture)
    
    def update_plot(self):
        try:
            self.plot( self.current_list_item.text() )
        except:
            pass

    def plot(self,read_clicked=''):
        if not self.updatePlot:     self.initPlot()
        ax = self.ax
        ax.cla()
        ax.set_axis_off()
        if self.bam_handle:
            # get alignment filters
            SVlen_thresh = None
            alen_thresh = None
            mapq_thresh = None
            masking_thresh = None
            if self.input_minSVlen_show.text() and self.input_minSVlen_show.text().isdigit():       SVlen_thresh = int(self.input_minSVlen_show.text())
            if self.input_minAlen_show.text() and self.input_minAlen_show.text().isdigit():         alen_thresh = int(self.input_minAlen_show.text())
            if self.input_minMapq_show.text() and self.input_minMapq_show.text().isdigit():         mapq_thresh = int(self.input_minMapq_show.text())
            if self.input_maxMaskFrac_show.text():
                try:        masking_thresh = max(0,float(self.input_maxMaskFrac_show.text()))
                except:     masking_thresh = 0
            if masking_thresh != None and masking_thresh > 1: # correct masking threshold if a percentage was input
                masking_thresh = masking_thresh / 100
            #/
            
            # Get read clicked
            read = read_clicked.split(') ')[1]
            #/
            # Check if we have masking file imported and if we have applied masking information to alignments
            if self.maskings != None and masking_thresh != None:
                for aln in self.reads_metadata[read][1]:
                    if 'masking' in aln: continue
                    aln['masking'] = calc_aln_masking(aln,self.maskings[0],self.maskings[1])
            #/
            # Check if we want to plot SVs
            SVs = None
            if self.SVs and read in self.SVs and self.checkbox_read_SVs.isChecked(): # read SVs
                SVs = self.SVs[read]
            if self.SVs and self.checkbox_SVs.isChecked(): # all SVs
                SVs = []
                SV_enums_added = set()
                for SV_read in self.SVs:
                    for SV in self.SVs[SV_read]:
                        if not SV['idx'] in SV_enums_added:
                            SVs.append(SV)
                            SV_enums_added.add(SV['idx'])
            #/
            # Check filter SVs by length
            if SVlen_thresh != None and SVs != None:
                SVs_filt = []
                for SV in SVs:
                    if SV['len'] >= SVlen_thresh or SV['type'] in ('BND','TRA',):
                        SVs_filt.append(SV)
                SVs = SVs_filt
            #/
            # Check filter SVs by BND
            if self.checkbox_ignoreTra.isChecked():
                SVs_filt = []
                for SV in SVs:
                    if not SV['type'] in ('BND','TRA',):
                        SVs_filt.append(SV)
                SVs = SVs_filt
            #/
            # Check if paint GTFs (check length thresh and label thresh)
            GTFs = None
            GTF_show_label = None
            GTF_minLen = 0
            GTF_label_minLen = 0
            if self.checkbox_geneTrack.isChecked() and self.features_idxDict_ROMD[0] != None:       GTFs = self.features_idxDict_ROMD
            if self.checkbox_geneLabels.isChecked():        GTF_show_label = True
            if self.input_minGeneLen_show.text() and self.input_minGeneLen_show.text().isdigit():                   GTF_minLen = max(0,int(self.input_minGeneLen_show.text()))
            if self.input_minGeneLen_label_show.text() and self.input_minGeneLen_label_show.text().isdigit():       GTF_label_minLen = max(0,int(self.input_minGeneLen_label_show.text()))
            #/
            # Check if paint scale
            showScale = None
            if self.checkbox_scale.isChecked():       showScale = True
            #/
            # Check if highlight region
            highlight_region = None
            if self.highlight_region != None:         highlight_region = self.highlight_region
            #/
            # Check if plot additional reads
            additional_reads_metadata = None
            if self.additional_reads_metadata != None: additional_reads_metadata = self.additional_reads_metadata
            #/
            # Get settings about region margin and alignment chaining distance (resulting in the same aChrom)
            alns_chain_distance = 20000
            try:        alns_chain_distance = int(self.show_settings['Alignment chain distance'])
            except:     pass
            region_margin = 10000
            try:        region_margin = int(self.show_settings['Region margin'])
            except:     pass
            #/
            aChroms,alignments = plot_on_ax(ax,fig=self.fig,bam_handle=self.bam_handle,read=read,read_metadata=self.reads_metadata[read],
                       alen_thresh=alen_thresh,mapq_thresh=mapq_thresh,masking_thresh=masking_thresh,SVs=SVs,
                       GTFs=GTFs,GTF_show_label=GTF_show_label,GTF_minLen=GTF_minLen,GTF_label_minLen=GTF_label_minLen,
                       showScale=showScale,highlight_region=highlight_region,additional_reads_metadata=additional_reads_metadata,
                       alns_chain_distance=alns_chain_distance,region_margin=region_margin)
            
            self.aChroms = aChroms
            self.alignments = alignments
        
        self.fig.canvas.draw_idle()
        #/
        
    ## MENU: FILE
    def import_bam(self):
        bam_path = QFileDialog.getOpenFileName(self, 'Alignment file import', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not bam_path: return
        bam_path = bam_path[0]
        # check so index file exist
        self.bam_path = bam_path
        
        # parse reads + lengths + one alignment (for quick lookup later)
        try:
            if bam_path:
                if not os.path.exists(bam_path+'.bai'): raise Exception
                if bam_path.endswith('.bam'):
                    bam_fo = pysam.AlignmentFile(bam_path,'rb')
                else:
                    bam_fo = pysam.AlignmentFile(bam_path,'r')
            else: return
        except:
            self.messageBox('Could not recognize input file.\nMake sure input file is a sorted BAM/SAM file and with available index-file.')
            return

        self.bam_handle = bam_fo
        self.setWindowTitle(bam_path)
        
        # If we had a genomeOverview window open, close it
        try:        self.genomeOverview_window.close()
        except:     pass
        #/
        
        # Check if we have specified read length for import
        read_len_thresh = 0
        try:    read_len_thresh = max(0,int(self.bam_settings['Minimum read length'])) #dont allow negative numbers
        except: pass
        #/
        
        # Parse read metadata (0=read_lenth + 1=locations from first alignment) + SA tag & save to outer
        try:
            self.freezeUI(self.centralwidget)
            self.reads_metadata,self.aln_idx_map,self.aln_idx_map_ROMD = parse_reads_metadata(self.bam_handle,read_len_thresh=read_len_thresh)
            self.reads_metadata_sorted = sorted(self.reads_metadata.items(),key=lambda x: x[1][0], reverse=True) # sort it by read length on import
            self.unfreezeUI(self.centralwidget)
        except:
            self.messageBox('Fatal error: Could not traverse alignment file. Is this a BAM/SAM file?')
            self.unfreezeUI(self.centralwidget)
            return
        #/
    
    def import_gtf(self):
        features_path = QFileDialog.getOpenFileName(self, 'Annotation file import', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not features_path: return
        features_path = features_path[0]
        
        try:
            if features_path:
                gene_column_identifiers = None
                try:                gene_column_identifiers = set(self.gtf_settings['Gene identifier (column 3)'].split(','))
                except:             gene_column_identifiers = set(['Failed_to_parse'])
                
                skip_features = set()
                try:
                    for entry in self.gtf_settings['Ignore types'].split(','):
                        skip_features.add(entry.strip())
                except:
                    skip_features = set(['Failed_to_parse'])
                
                try:
                    self.freezeUI(self.centralwidget)
                    self.features = GTF_GFF_parser(features_path,gtf_attribute_split=self.gtf_settings['GTF Attribute delimiter'],gff_attribute_split=self.gtf_settings['GFF Attribute delimiter'],
                                                   gtf_gene_tag=self.gtf_settings['GTF Gene tag'],gff_gene_tag=self.gtf_settings['GFF Gene tag'],
                                                   gene_column_identifiers=gene_column_identifiers,skip_features=skip_features)
            
                    self.features_idxDict_ROMD[0] = self.features[1]
                    self.features_idxDict_ROMD[1] = rangeOverlaps_makeMappingDict(self.features[1],1000,coordsKey='rcoords',sortByKey='rname')
                    self.unfreezeUI(self.centralwidget)
                except:
                    self.unfreezeUI(self.centralwidget)
                    self.messageBox('Fatal error: Could not traverse annotation file. Is this a GTF/GFF file?')
                    return
                
            else: return
        except:
            self.messageBox('Could not recognize input file.\nMake sure input file is in GTF/GFF format with corresponding file-extension.')
            return
        
    def import_repeats(self):
        masking_path = QFileDialog.getOpenFileName(self, 'Repeat annotation import', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not masking_path: return
        masking_path = masking_path[0]
        
        try:
            if masking_path:
                self.freezeUI(self.centralwidget)
                self.maskings = repeat_parser_repeatMasker(masking_path)
                self.unfreezeUI(self.centralwidget)
            else: return
        except:
            self.messageBox('Could not recognize input file.\nMake sure input file is repeatmasking .out file.')
            self.unfreezeUI(self.centralwidget)
            return
    
        # If we had a genomeOverview window open, close it
        try:        self.genomeOverview_window.close()
        except:     pass
        #/
    
    def import_vcf(self):
        vcf_path = QFileDialog.getOpenFileName(self, 'Structural variant import', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not vcf_path: return
        vcf_path = vcf_path[0]
        
        try:
            if vcf_path:
                self.freezeUI(self.centralwidget)
                self.SVs = VCF_parser(vcf_path,rname2_key=self.vcf_settings['Chr2'],rpos2_key=self.vcf_settings['Pos2'],
                                      SV_type_key=self.vcf_settings['SV type'],read_supp_key=self.vcf_settings['Supportive reads'],
                                      SV_len_key=self.vcf_settings['SV length'])
                self.unfreezeUI(self.centralwidget)
            else: return
        except:
            self.messageBox('Could not recognize input file.\nMake sure input file is in VCF format.')
            self.unfreezeUI(self.centralwidget)
            return
    ##/
    ## MENU:SELECT
    def import_read_name(self):
        inputbox,_ = QInputDialog.getMultiLineText(self,'Import read name(s)','Read name import (new-line or comma-separated):', '')
        if inputbox:
            read_list = []
            for line in inputbox.split('\n'):
                line = line.split(',')
                for entry in line:
                    read_list.append(entry)
            self.read_list = read_list
            self.search_reads()
            
    
    def import_read_list(self):
        read_list_path = QFileDialog.getOpenFileName(self, 'Open file', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not read_list_path: return
        read_list_path = read_list_path[0]
        
        try:
            if read_list_path:
                self.read_list = read_list_parser(read_list_path)
                self.search_reads()
        except:
            self.messageBox('Could not recognize input file.\nMake sure input file is in non-binary format and tab-separated with readname in first column.')
            return
        
    def clear_selection(self):
        self.read_list = None
        self.read_list_grouped = None
        self.highlight_region = None
        self.current_list_item = None
        self.additional_reads_metadata = None
        self.search_reads()
    ##/
    ## MENU:EXPORT
    def export_current_read(self):
        try:
            inputbox,_ = QInputDialog.getMultiLineText(self,'Read name export','Selected read name:', self.current_list_item.text().split(') ')[1])
        except:
            pass
    
    def export_read_highlight(self):
        try:
            reads_highlight = []
            if self.current_list_item != None:              reads_highlight.append(self.current_list_item.text().split(') ')[1])
            if self.additional_reads_metadata != None:      reads_highlight += list(self.additional_reads_metadata.keys())
            if reads_highlight:
                inputbox,_ = QInputDialog.getMultiLineText(self,'Read highlight export','Selected read name(s):', '\n'.join(reads_highlight))
        except:
            pass

    def export_read_selection(self):
        read_list_path = QFileDialog.getSaveFileName(self, 'Save file', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not read_list_path: return
        read_list_path = read_list_path[0]
        
        try:
            if read_list_path:
                reads = []
                for enum in range(self.listWidget.count()):
                    reads.append(self.listWidget.item(enum).text().split(') ')[1]) #split away the read length
                with open(read_list_path,'w') as nf:
                    nf.write('\n'.join(reads)+'\n')
        except:
            self.messageBox('Could not write output file! Make sure you have permissions to write in output folder.')
            return
        
    def export_regions(self):
        if self.aChroms:
            regions_str = ''
            for aChrom_enum,aChrom in self.aChroms.items():
                regions_str += aChrom['rname'] + ':' + '-'.join(map(str,aChrom['rcoords'])) + '\n'
                
            inputbox,_ = QInputDialog.getMultiLineText(self,'Region export','Regions:', regions_str)
    
    def export_alignments(self):
        if self.alignments:
            alignments_str = ''
            for aln in self.alignments:
                alignments_str += aln['rname'] + ':' + '-'.join(map(str,aln['rcoords'])) + '\t' + aln['strand'] + '\t' + '-'.join(map(str,aln['qcoords'])) + '\n'
                
            inputbox,_ = QInputDialog.getMultiLineText(self,'Alignment export','Reference position, strand, query position:', alignments_str)
        
    def export_image(self):
        export_path = QFileDialog.getSaveFileName(self, 'Save file', '/Datadrives/bulk1/jacob/vc_vsRef/SAKAMOTO/read_selections/test_data')
        if not export_path: return
        export_path = export_path[0]
        try:
            if self.fig and export_path:
                # check if no extension was given, then append default PNG
                if export_path.find('.') == -1:     export_path += '.png'
                #/
                # Check if file exists
                fileExists = False
                if os.path.exists(export_path):
                    fileExists = True
                    query_prompt = QMessageBox.question(self,'Overwrite file?','File already exists. Overwrite?', QMessageBox.Yes, QMessageBox.No)
                    if query_prompt == QMessageBox.Yes:
                        fileExists = False
                #/
                if not fileExists:
                    exportImage(self.fig,export_path)
        except:
            self.messageBox('Could not write output file! Make sure you included file extension and have permissions to write in output folder.\nSupported extensions include *.png, *.svg, *.pdf')
            return
        
    ##/
    ## MENU:SETTINGS
    def handle_settings_import(self,settings,default_settings=None):
        output_data = []
        for form_enum in range(self.settings_window.formLayout.count()):
            widget = self.settings_window.formLayout.itemAt(form_enum).widget()
            text = widget.text()
            if not output_data:                 output_data.append([])
            if len(output_data[-1]) == 2:      output_data.append([])
            
            output_data[-1].append(text)
            
        for key,val in output_data:
            settings[key] = val
            
        # Check if we applied default settings
        if default_settings and output_data[0][1] == 'RESTORE_DEFAULTS':
            for key,val in default_settings.items():
                settings[key] = val
        #/
    
    def show_settings(self):
        settings = self.show_settings
        default_settings = self.show_settings_default
        self.settings_window = inputSettings_window()
        self.settings_window.window_closed.connect(lambda: self.handle_settings_import(settings,default_settings))
        self.settings_window.input_form = settings
        self.settings_window.executeWindow()
        self.settings_window.setWindowTitle('BAM import settings')
        self.settings_window.header.setText('Enter settings for regions:')
        
    def bam_import_settings(self):
        settings = self.bam_settings
        default_settings = self.bam_settings_default
        self.settings_window = inputSettings_window()
        self.settings_window.window_closed.connect(lambda: self.handle_settings_import(settings,default_settings))
        self.settings_window.input_form = settings
        self.settings_window.executeWindow()
        self.settings_window.setWindowTitle('BAM import settings')
        self.settings_window.header.setText('Enter BAM import thresholds:')
    
    def vcf_import_settings(self):
        settings = self.vcf_settings
        default_settings = self.vcf_settings_default
        self.settings_window = inputSettings_window()
        self.settings_window.window_closed.connect(lambda: self.handle_settings_import(settings,default_settings))
        self.settings_window.input_form = settings
        self.settings_window.executeWindow()
        self.settings_window.setWindowTitle('VCF import settings')
        self.settings_window.header.setText('Fill in tags for keys:')
    
    def gtf_gff_import_settings(self):
        settings = self.gtf_settings
        default_settings = self.gtf_settings_default
        self.settings_window = inputSettings_window()
        self.settings_window.window_closed.connect(lambda: self.handle_settings_import(settings,default_settings))
        self.settings_window.input_form = settings
        self.settings_window.executeWindow()
        self.settings_window.setWindowTitle('GTF/GFF import settings')
        self.settings_window.header.setText('Fill in tags for keys:')
    ##/
    
    ## Functions to call freeze to objects. These will be used to block user inputs and/or show that the software is loading
    def freezeUI(self,obj):
        obj.setEnabled(False)
        QApplication.processEvents()
        
    def unfreezeUI(self,obj):
        QApplication.processEvents()
        obj.setEnabled(True)
    ##/
    
    def search_reads(self):
        self.listWidget.clear()
        self.current_list_item = None
        self.additional_reads_metadata = None
        if not self.reads_metadata:
            self.messageBox('No alignment file has been imported!')
            return
        
        self.freezeUI(self.listWidget) # freeze the list so user cannot spam-click it and lag it
        self.freezeUI(self.button_findReads) # freeze the list so user cannot spam-click it and lag it
        
        ## Check filters
        regs_req = []
        aln_reqs = []
        read_len_req = 0 # default is no length requirement
        numAlns_req = 1 # default is one
        read_SV_reqs = set()
        read_SV_len_req = 0
        # Region/gene/SV
        if self.input_region.text():
            for inp_text in self.input_region.text().split(';'):
                # Check if entry is genome overview, then skip
                if inp_text.find('<GenomeOverview>') != -1: continue
                #/
                
                # parse region
                try:
                    rname = inp_text.split(':')[0]
                    rcoords = list(map(int,inp_text.split(':')[1].replace(',','').split('-')))
                    regs_req.append([rname,rcoords])
                except:
                    # parse gene
                    try:
                        if self.features and inp_text.lower() in self.features[0]:
                            gtf_regs = [[],inp_text.lower(),'gtf']
                            for feature in self.features[0][inp_text.lower()]:
                                gtf_regs[0].append(feature[:2]) #import the rname,rcoords
                            regs_req.append(gtf_regs)
                        else:
                            regs_req.append(['notFound',[0,0]]) # dummy-add region if we didnt find
                    except:
                        pass
                    #/
                #/
        #/
        # SV-type / SV length
        if self.input_partOfSVtype.text():
            for entry in self.input_partOfSVtype.text().split(';'):
                read_SV_reqs.add(entry)
        
        if self.input_minSVlen.text() and self.input_minSVlen.text().isdigit():             read_SV_len_req = max(0,int(self.input_minSVlen.text())) # never allow this to be < 0
        #/
        # Read length / alignment length / MAPQ / masked fraction / Min. num. of alignments
        if self.input_minReadLen.text() and self.input_minReadLen.text().isdigit():         read_len_req = max(0,int(self.input_minReadLen.text())) # never allow this to be < 0
        if self.input_minNumAlns.text() and self.input_minNumAlns.text().isdigit():         numAlns_req = max(1,int(self.input_minNumAlns.text())) # never allow this to be set < 1
        if self.input_minAlen.text() and self.input_minAlen.text().isdigit():               aln_reqs.append(['alen',int(self.input_minAlen.text())])
        if self.input_minMapq.text() and self.input_minMapq.text().isdigit():               aln_reqs.append(['mapq',int(self.input_minMapq.text())])
        if self.input_maxMaskFrac.text(): # Try to parse fraction/percentage of alignment masked overlap. If not able to parse it, just pass.
            try:
                aln_reqs.append(['mask',max(0,float(self.input_maxMaskFrac.text()))])
            except:     pass
        #/
        ##/
        
        ## If we have regions specified, then look-up these reads. Else, traverse all reads
        reads_metadata = []
        if regs_req:
            reads_metadata_selection = set()
            for reg_req in regs_req:
                # Check if reg was parsed from GTF entry, then check GTF entries
                if reg_req[-1] == 'gtf':
                    for gtf_reg in reg_req[0]:
                        for oS,idx,oE in rangeOverlaps_lookup(gtf_reg[1],self.aln_idx_map_ROMD,1000,map_lookup_key=gtf_reg[0]):
                            read = self.aln_idx_map[idx]['qname']
                            reads_metadata_selection.add(read)
                #/
                else:
                    for oS,idx,oE in rangeOverlaps_lookup(reg_req[1],self.aln_idx_map_ROMD,1000,map_lookup_key=reg_req[0]):
                        read = self.aln_idx_map[idx]['qname']
                        reads_metadata_selection.add(read)
                        
            for read,data in self.reads_metadata_sorted:
                if read in reads_metadata_selection:
                    reads_metadata.append([read,data])
        else:
            reads_metadata = self.reads_metadata_sorted
        ##/        
        
        ## Traverse reads, filter, and populate list
        for read,data in reads_metadata:
            ## Run filters
            # Check if we have a read selection and skip if current read is not in selection
            if self.read_list and not read in self.read_list: continue
            #/
            # Skip read if it does not pass read length requirement (break, we have sorted reads by length.)
            if read_len_req and data[0] < read_len_req:
                break
            #/
            # Check if we have SV requirements
            if self.SVs and (read_SV_reqs or read_SV_len_req):
                if not read in self.SVs:        continue
                
                # Check if pass SV len req
                if read_SV_len_req > 0:
                    read_pass_SV_len_req = False
                    for SV in self.SVs[read]:
                        if SV['len'] != None and SV['len'] >= read_SV_len_req:
                            read_pass_SV_len_req = True
                    if not read_pass_SV_len_req: continue
                #/
                # Check if pass SV type req
                if read_SV_reqs:
                    pass_SV_reqs = []
                    for SV in self.SVs[read]:
                        if SV['type'] in read_SV_reqs:
                            # check if skip SV; if its too short
                            if read_SV_len_req > 0 and SV['len'] != None and SV['len'] < read_SV_len_req: continue 
                            #/
                            pass_SV_reqs.append(SV['type'])
                    if 'reqAny' and not pass_SV_reqs: continue
                    if 'reqAll' and not len(pass_SV_reqs) == len(read_SV_reqs): continue
                #/
            #/
            # Traverse alignments, check (1) region reqs and (2) alignment reqs
            pass_reg_reqs = set()
            pass_aln_reqs = []
            for aln in data[1]:
                # Check if pass regions
                for enum,reg_req in enumerate(regs_req):
                    # Check if reg was parsed from GTF entry, then check GTF entries
                    if reg_req[-1] == 'gtf':
                        for gtf_reg in reg_req[0]:
                            if gtf_reg[0] == aln['rname'] and getRangeOvlp(gtf_reg[1],aln['rcoords']) >= 0:
                                pass_reg_reqs.add(enum)
                                break # stop on first match
                    #/
                    else:
                        if reg_req[0] == aln['rname'] and getRangeOvlp(reg_req[1],aln['rcoords']) >= 0:
                            pass_reg_reqs.add(enum)
                #/
                # Check if pass alignment len
                aln_passes = []
                for aln_req in aln_reqs:
                    if aln_req[0] == 'alen' and (aln['rcoords'][-1]-aln['rcoords'][0] >= aln_req[1]):       aln_passes.append(aln_req)
                    if aln_req[0] == 'mapq' and aln['mapq'] >= aln_req[1]:                                  aln_passes.append(aln_req)
                    if aln_req[0] == 'mask':
                        if self.maskings != None: # Only run if we have a masking file imported
                            if not 'masking' in aln:        aln['masking'] = calc_aln_masking(aln,self.maskings[0],self.maskings[1])
                            fractional_cutoff = 0
                            if aln_req[1] > 1:     fractional_cutoff = aln_req[1]/100 # If percentage input, convert to fractional
                            else:                  fractional_cutoff = aln_req[1]
                            if aln['masking'][2] <= fractional_cutoff:
                                aln_passes.append(aln_req)
                                
                if len(aln_passes) == len(aln_reqs):
                    pass_aln_reqs.append(aln)
                #/
            #/
            # Skip if does not pass regions requirement (Only if we are not in genome_overview mode)
            if regs_req and 'passAny' and not pass_reg_reqs: continue
            if regs_req and 'passAll' and not len(pass_reg_reqs) == len(regs_req): continue
            #/

            # Skip if does not pass number of required alignments
            if numAlns_req and len(pass_aln_reqs) < numAlns_req: continue
            #/
            ##/
            self.listWidget.addItem('('+str(data[0])+') ' + read)
            
        # Check if we wanted to group reads by region
        if self.read_list_grouped != None:
            # Get list contents and add regional marker
            readlist_new_formatted = []
            for listEnum in range(self.listWidget.count()):
                read = self.listWidget.item(listEnum).text().split(') ')[1]
                read_regEnums = []
                for regEnum,(rearr_bins,readGroup) in enumerate(self.read_list_grouped.items()):
                    if read in readGroup:
                        read_regEnums.append(regEnum+1) # append it 1-based and not 0-based
                        
                readlist_new_formatted.append('R'+'+'.join(sorted(map(str,read_regEnums)))+' ('+str(self.reads_metadata[read][0])+') ' + read)
            #/
            # clear + re-populate list
            readlist_new_formatted = sorted(readlist_new_formatted,key=lambda x: (x.split(' ')[0],-int(x.split('(')[1].split(')')[0]) )) # sort by (1) region, (2) read length [use negative to have it sorted big>small]
            self.listWidget.clear()
            for entry in readlist_new_formatted:
                self.listWidget.addItem(entry)
            #/
        #/
        ##/
        
        self.label_reads_found.setText( self.label_reads_found.text().split(':')[0] + ': '+str(self.listWidget.count()) )
        self.listWidget.setCurrentRow(0)
        # Try get default background color (for left/right-click on reads)
        try:    self.listwidget_defaultColor = self.listWidget.item(0).background()
        except: pass
        #/
        # enable GUI again for user input
        self.unfreezeUI(self.listWidget)
        self.unfreezeUI(self.button_findReads)
        #/
    
    def handle_genomeOverview_trigger(self):
        self.read_list = self.genomeOverview_window.read_list
        self.read_list_grouped = self.genomeOverview_window.read_list_grouped
        self.highlight_region = self.genomeOverview_window.highlight_region
        self.input_region.setText('<GenomeOverview>')
        self.search_reads()
        return
        
    def genome_overview(self):
        if self.bam_handle == None:
            self.messageBox('You need to import an alignment file first!')
            return
        
        if self.maskings == None:
            self.messageBox('No masking file has been imported!\nWhile this module can run without masking data, it is highly recommended to include this data!')
        
        
        self.genomeOverview_window = genomeOverview_window()
        self.genomeOverview_window.bin_clicked_signal.connect(self.handle_genomeOverview_trigger)
        self.genomeOverview_window.bam_handle = self.bam_handle
        self.genomeOverview_window.maskings = self.maskings
        self.genomeOverview_window.reads_metadata = self.reads_metadata
        self.genomeOverview_window.executeWindow()
        
    def closeEvent(self, event):
        # if close mainwindow, then try to close other child windows too
        try:        self.genomeOverview_window.close()
        except:     pass
        try:        self.settings_window.close()
        except:     pass
        #/  
        event.accept()
    

def main():
    app = QApplication(sys.argv)
    application = ApplicationWindow()
    application.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()