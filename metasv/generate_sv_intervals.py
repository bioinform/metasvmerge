import sys
import os
import argparse
import multiprocessing
import logging
import collections
import itertools
import traceback
from functools import partial
import json
import base64
import time

import pysam
import pybedtools

from defaults import *
from sv_interval import *


def concatenate_files(files, output):
    with open(output, 'w') as outfile:
        for fname in files:
            if not os.path.isfile(fname): continue
            with open(fname) as infile:
                outfile.write(infile.read())


def find_softclip(aln):
    if aln.cigar is None:
        return None
    soft_clips = [(i, length) for i,(op, length) in enumerate(aln.cigar) if op == 4]
    if len(soft_clips) != 1:
        return None
    
    i, soft_clip = soft_clips[0]
    dist_L_end = sum(map(lambda x:x[1] if x[0] in [0,1,4] else 0, aln.cigar[0:i]))
    dist_R_end = sum(map(lambda x:x[1] if x[0] in [0,1,4] else 0, aln.cigar[i+1:]))

    return soft_clip, dist_L_end, dist_R_end

    


def is_good_candidate(aln, min_avg_base_qual=20, min_mapq=5, min_soft_clip=20, max_nm=10,
                      min_matches=50, skip_soft_clip=False, good_neigh_check= False ):
    if aln.is_duplicate:
        return False
    if aln.is_unmapped:
        return False
    if aln.mapq < min_mapq:
        return False
    if aln.cigar is None:
        return False

    if not good_neigh_check:
        soft_clip_tuple = find_softclip(aln)
        if not soft_clip_tuple:
            return False
        else:
            soft_clip, dist_L_end, dist_R_end = soft_clip_tuple
            if not (min_soft_clip <= soft_clip):
                return False
    else:
        if skip_soft_clip:
            soft_clip_tuple = find_softclip(aln)
            if soft_clip_tuple:
                if soft_clip_tuple[0] > min_soft_clip:
                    return False
    
            
    ins_lengths = sum([0] + [length for (op, length) in aln.cigar if op == 1])
    mismatches = int(aln.opt("XM")) if "XM" in aln.tags else 0
    matches = aln.alen - ins_lengths - mismatches
    nm = int(aln.opt("NM"))
    if nm > max_nm or matches < min_matches:
        return False

    if not good_neigh_check:
        if aln.cigar[0][0] == 4:
            avg_base_quality = float(sum(map(ord, aln.qual[:soft_clip]))) / soft_clip
        else:
            avg_base_quality = float(sum(map(ord, aln.qual[-soft_clip:]))) / soft_clip

        return avg_base_quality - 33 >= min_avg_base_qual
    else:
        return True

def get_interval(aln, pad=500):
    start = aln.pos
    end = aln.aend    

    if aln.cigar[0][0] == 4:
        return max(0, start - pad), start + pad
    return max(0, end - pad), end + pad


def merged_interval_features(feature, bam_handle):
    support_list = feature.name.split(",")
    locations = sorted(map(int, support_list[0:-1:3]))
    other_bp_ends = support_list[-1]
    num_unique_locations = len(set(locations))
    count_str = ",".join(["%s,%s" % (i, c) for (i, c) in collections.Counter(locations).items()])
    plus_support = len([i for i in support_list[2:-1:3] if i == "+"])
    minus_support = len(locations) - plus_support
    locations_span = max(locations) - min(locations)
    interval_readcount = bam_handle.count(reference=feature.chrom, start=feature.start, end=feature.end)
    info = {"SC_PLUS_SUPPORT":plus_support, "SC_MINUS_SUPPORT":minus_support, 
            "SC_LOCATIONS_SPAN":locations_span, "SC_NUM_UNIQUE_LOCATIONS":num_unique_locations,
            "SC_COUNT_STR": count_str, "SC_COVERAGE":interval_readcount, "SC_OTHER_BP_ENDS": other_bp_ends, 
            "SC_SC_BP_ENDS": "%s-%s"%(feature.start, feature.end)}
    name = "%s,%s,0,SC" % (
        base64.b64encode(json.dumps(info)), feature.fields[6].split(',')[0])

    return pybedtools.Interval(feature.chrom, feature.start, feature.end, name=name, score=feature.score,
                               otherfields=[str(interval_readcount)]+feature.fields[6:])



def generate_other_bp_interval(feature,pad):
    other_bp=int(feature.name.split(",")[1])
    return pybedtools.Interval(feature.chrom, max(other_bp-pad,0),other_bp+pad, name=feature.name, score=feature.score, strand=feature.strand, otherfields=[feature.fields[6]])


def add_other_bp_fields(feature,start,end):
    return pybedtools.Interval(feature.chrom, feature.start, feature.end, name='%s,%d-%d'%(feature.name,start,end), score=feature.score,
                               otherfields=feature.fields[6:])


def get_full_interval(feature,pad):
    name_fields = feature.name.split(",")
    info = json.loads(base64.b64decode(name_fields[0]))
    other_bp_ends=info["SC_OTHER_BP_ENDS"]
    start = feature.start
    end = feature.end
    sv_type = name_fields[1]
    if "-" in other_bp_ends:
        other_bp_start,other_bp_end=map(lambda x:int(x),other_bp_ends.split("-"))
        if other_bp_start != 0 or other_bp_end != 0:
            start = min((feature.start+feature.end)/2,(other_bp_start+other_bp_end)/2)
            end = max((feature.start+feature.end)/2,(other_bp_start+other_bp_end)/2)

    sv_len = 0 if sv_type == "INS" else max(end-start,0)
    info["SOURCES"] = "%s-%d-%s-%d-%d-SoftClip" % (feature.chrom, start, feature.chrom, end, sv_len)
    name = "%s,%s,%d,%s"%(base64.b64encode(json.dumps(info)),sv_type,sv_len,'SC')
    
    
    return pybedtools.Interval(feature.chrom, start, end, name=name, score=feature.score,
                                   otherfields=feature.fields[6:])

    

def generate_sc_intervals_callback(result, result_list):
    if result is not None:
        result_list.append(result)

def infer_svtype(aln, min_isize, max_isize):
    if aln.mate_is_unmapped:
        return "INS"
    if aln.tid != aln.rnext:
        return "CTX;INS"
    if (aln.is_reverse and aln.mate_is_reverse) or (not aln.is_reverse and not aln.mate_is_reverse):
        return "INV"
    if (aln.pos < aln.pnext and aln.is_reverse) or (aln.pos > aln.pnext and not aln.is_reverse):
        return "DUP;ITX"
    if abs(aln.tlen) > max_isize:
        return "DEL"
    if abs(aln.tlen) < min_isize:
        return "INS"
    return "NONE"

def find_other_bp(aln, isize_mean, svtype, soft_clip, dist_L_end, dist_R_end, 
                  soft_clip_location, skip_soft_clip =False, skip_neigh=True,  
                  min_dist_end=2, wiggle= 20):
    if (skip_soft_clip and soft_clip>0) or (skip_neigh and soft_clip<=0):
        return None
        
    other_bp = None      
    if soft_clip>0:
        is_left = dist_L_end <= min_dist_end
        is_right = dist_R_end <= min_dist_end
    else:
        is_left = aln.pos > (soft_clip_location - wiggle)
        is_right = aln.aend < (soft_clip_location + wiggle)
    
    is_first_pair = aln.pos <= aln.pnext
    isize_diff = abs(aln.tlen) - isize_mean
    isize_sum = abs(aln.tlen) + isize_mean
    
         
    if svtype == "INS":
        if soft_clip>0:
            if (is_left and aln.is_reverse) or (is_right and not aln.is_reverse): 
                return 1
        else:
            if not aln.is_reverse and is_right:
                if aln.mate_is_unmapped or aln.tid != aln.rnext or (aln.pnext > (soft_clip_location - wiggle)):
                    return 1
            elif aln.is_reverse and is_left:
                if aln.mate_is_unmapped or aln.tid != aln.rnext or ((aln.pnext+aln.rlen) < (soft_clip_location + wiggle)):
                    return 1
    elif svtype == "INV":
        if soft_clip>0:
            if is_right and not aln.is_reverse and not is_first_pair:
                other_bp = soft_clip_location - (isize_diff + 2*dist_L_end)
            elif is_left and not aln.is_reverse and not is_first_pair:
                other_bp = soft_clip_location + (-isize_diff - 2*soft_clip)
            elif is_left and aln.is_reverse and is_first_pair:
                other_bp = soft_clip_location + (isize_diff + 2*dist_R_end)
            elif is_right and aln.is_reverse and is_first_pair:
                other_bp = soft_clip_location - (-isize_diff - 2*soft_clip)
            elif is_left and aln.is_reverse and not is_first_pair:
                other_bp = soft_clip_location - (isize_sum - 2*dist_R_end)
            elif is_right and not aln.is_reverse and is_first_pair:
                other_bp = soft_clip_location + (isize_sum - 2*dist_L_end)
        else:
            if not aln.is_reverse and not is_first_pair and is_right:
                other_bp = soft_clip_location - (isize_diff + 2*(soft_clip_location-aln.pos))
            elif not aln.is_reverse and not is_first_pair and is_left:
                other_bp = soft_clip_location + (-isize_diff + 2*(aln.pos-soft_clip_location))
            elif aln.is_reverse and is_first_pair and is_left:
                other_bp = soft_clip_location + (isize_diff + 2*(aln.aend-soft_clip_location))
            elif aln.is_reverse and is_first_pair and is_right:
                other_bp = soft_clip_location - (-isize_diff + 2*(soft_clip_location-aln.aend))
            elif aln.is_reverse and not is_first_pair and is_left:
                other_bp = soft_clip_location - (isize_sum - 2*(aln.aend-soft_clip_location))
            elif not aln.is_reverse and is_first_pair and is_right:
                other_bp = soft_clip_location + (isize_sum - 2*(soft_clip_location-aln.pos))
    elif svtype == "DEL":
        if soft_clip>0:
            if is_left and aln.is_reverse and not is_first_pair:
                other_bp = soft_clip_location - isize_diff
            elif is_right and not aln.is_reverse and is_first_pair:
                other_bp = soft_clip_location + isize_diff
        else:
            if aln.is_reverse and is_left:
                other_bp = soft_clip_location - isize_diff
            elif not aln.is_reverse and is_right:
                other_bp = soft_clip_location + isize_diff
    elif svtype == "DUP":
        if soft_clip>0:
            if is_left and aln.is_reverse and is_first_pair:
                other_bp = soft_clip_location + isize_sum
            elif is_right and not aln.is_reverse and not is_first_pair:
                other_bp = soft_clip_location - isize_sum
        else:
            if aln.is_reverse and is_first_pair and is_left:
                other_bp = soft_clip_location + isize_sum
            elif not aln.is_reverse and not is_first_pair and aln.pos < (soft_clip_location + wiggle):
                other_bp = soft_clip_location - isize_sum       
    if not other_bp is None:
        return max(other_bp,0)
    else:
        return None


def check_overlap(start,end,chrom,interval,overlap_ratio=OVERLAP_RATIO):
    if interval.chrom != chrom:    
        return False,start,end
    if max(interval.start, start) \
            >= min(interval.end, end):
        return False,start,end
    self_length = float(end - start)
    other_length = float(interval.end - interval.start)
    overlap_length = min(interval.end , end ) - max(interval.start,start)
    if float(overlap_length) >= max(overlap_ratio * self_length,overlap_ratio * other_length):
        return True,  min(interval.start, start), max(interval.end, end)
    else:
        return False, start,end

def blind_merge(intervals,cols,ops):
    func_logger = logging.getLogger("%s-%s" % (blind_merge.__name__, multiprocessing.current_process()))
    try:
        columns=map(lambda x:int(x)+4,range(len(cols.split(","))))
        operations=ops.split(',')    
        if columns and not operations:
            func_logger.error("Aborting!")
            raise Exception("Bad column and operation combinations for merge: %s, %s\n" % (cols, ops))    
        if len(columns)>1 and len(operations)==1:
            operations=[ops for c in columns]

        if len(columns)!=len(operations):
            func_logger.error("Aborting!")
            raise Exception("Bad column and operation combinations for merge: %s, %s\n" % (cols, ops))    
    
        column_operations={c:operations[i] for i,c in enumerate(columns)}    
    
        if not intervals:
            return intervals
    
        start = intervals[0].start
        end = intervals[0].end
        chrom = intervals[0].chrom
        other_fields = {c:[intervals[0].fields[c-1]] if c<=len(intervals[0].fields) else [""]  for c,o in column_operations.iteritems()}
        for interval in intervals[1:]:
            start = min(start, int(interval.start))
            end = max(end, int(interval.end))
            for c in column_operations: 
                other_fields[c].append(interval.fields[c-1])
        operation_function={"collapse":(lambda x: ",".join(x)), 
                            "sum": "%s"%(lambda x:sum(map(lambda y:int(y) if y else 0,x))),
                            "min": "%s"%(lambda x:min(map(lambda y:int(y) if y else 0,x))),
                            "max": "%s"%(lambda x:max(map(lambda y:int(y) if y else 0,x))),
                            "first": "%s"%(lambda x:x[0]),
                            "last": "%s"%(lambda x:x[-1]),
                            "distinct": (lambda x:",".join(set(x))),
                            "count": (lambda x:len(x))}
        for c,o in column_operations.iteritems():
            if o in operation_function:
                other_fields[c] = operation_function[o](other_fields[c])
            else:
                func_logger.error("Aborting!")
                raise Exception("Not supported operation: %s\n" % (o))    
        name = other_fields[columns[0]] if columns else ""
        score = other_fields[columns[1]] if len(columns)>=2 else ""
        strand = other_fields[columns[2]] if len(columns)>=3 else ""
        otherfields = [other_fields[c] for c in columns[3:]]
        return pybedtools.Interval(chrom=chrom,start=start,end=end,name=name,score=score,strand=strand,otherfields=otherfields)    
    except Exception as e:
        func_logger.error('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()
        print()
        raise e

    
def merge_intervals_bed(bedtool, overlap_ratio , c ,o):
    bedtool=bedtool.sort().cut([0,1,2]+(map(lambda x: int(x)-1,c.split(',')) if c else []))
    new_intervals = []

    if bedtool.count()==0:
        return bedtool

    intervals=[intv for intv in bedtool]

    current_merged_interval_list = [intervals[0]]
    start = intervals[0].start
    end = intervals[0].end   
    chrom = intervals[0].chrom
    for i in xrange(len(intervals) - 1):
        next_interval = intervals[i+1]
        is_overlap, tmp_start, tmp_end=check_overlap(start,end,chrom,next_interval,overlap_ratio=overlap_ratio)
        if is_overlap:
            current_merged_interval_list.append(next_interval)
            start = tmp_start
            end = tmp_end
        else:
            new_intervals.append(current_merged_interval_list)
            current_merged_interval_list = [next_interval]
            start = next_interval.start
            end = next_interval.end  
            chrom = next_interval.chrom

    new_intervals.append(current_merged_interval_list)
    merged_bed = pybedtools.BedTool([blind_merge(intervals,c,o) for intervals in new_intervals])
                                     
    return merged_bed.sort()

def merge_for_each_sv(bedtool,c,o,svs_to_softclip=SVS_SOFTCLIP_SUPPORTED,
                      overlap_ratio=OVERLAP_RATIO,d=0, reciprocal_for_2bp=True,
                      sv_type_field = [3,1], inter_tools = False):
    merged_bedtool = pybedtools.BedTool([])
    for svtype in svs_to_softclip:
        sv_bedtool = bedtool.filter(lambda x: svtype in x.fields[sv_type_field[0]].split(',')[sv_type_field[1]]).sort()
        if sv_bedtool.count()==0: continue
        if svtype == "INS" or not reciprocal_for_2bp:
            sv_bedtool=sv_bedtool.merge(c=c, o=o, d=d)
        else:
            sv_bedtool = merge_intervals_bed(sv_bedtool,overlap_ratio=overlap_ratio,
                                                  c=c,o=o)
        merged_bedtool=sv_bedtool.cat(merged_bedtool,postmerge=False)
    return merged_bedtool.sort()
    

def fix_merged_fields(feature,inter_tools=True):
    name_fields = feature.name.split(",")
    n = len(name_fields)/4
    info = {}    
    sv_type = name_fields[1]
    sv_length = 0
    sv_tools=set([])
    
    sc_sub_intervals_info=[]
    if not inter_tools:
        info["SC_SUBINTERVAL_INFOs"]=[]
    for i in range(n):    
        sub_interval=name_fields[i*4:(i+1)*4]
        sub_info=json.loads(base64.b64decode(sub_interval[0]))
        if not inter_tools:
            info["SC_SUBINTERVAL_INFOs"].append(sub_info)
        else:
            info.update({k:v for k,v in sub_info.iteritems() if k not in ["SOURCES","SC_SUBINTERVAL_INFOs"]})
            if "SC_SUBINTERVAL_INFOs" in sub_info:
                sc_sub_intervals_info.extend(sub_info["SC_SUBINTERVAL_INFOs"])                
        sv_tools.update(set(map(lambda x: x.split('-')[-1],sub_info["SOURCES"].split(','))))
        if i==0:
            info["SOURCES"] = sub_info["SOURCES"]
        else:
            info["SOURCES"] += ","+sub_info["SOURCES"]
        sv_length = max(sv_length,int(sub_interval[2]))

    if sc_sub_intervals_info:
        for k in ["SC_COVERAGE", "SC_NEIGH_SUPPORT"]:
            info[k]=",".join(map(lambda x:"%s"%x[k],sc_sub_intervals_info))
        info["SC_READ_SUPPORT"]=",".join(map(lambda x:"%s"%(x["SC_PLUS_SUPPORT"]+x["SC_MINUS_SUPPORT"]),sc_sub_intervals_info))

    sv_methods=sorted(list(reduce(operator.add, [sv_sources_to_type[tool] for tool in list(sv_tools)])))
    info["NUM_SVMETHODS"] = len(sv_methods)
    info["NUM_SVTOOLS"] = len(sv_tools)
    
    if not inter_tools:
        info['SCORE'] = feature.score
        
    return pybedtools.Interval(feature.chrom, feature.start, feature.end, name="%s,%s,%d,%s" % (
            base64.b64encode(json.dumps(info)), sv_type, sv_length,
            ";".join(sv_methods)),
            score = feature.score if not inter_tools else "%d"%len(sv_methods), otherfields=feature.fields[6:])


def fine_tune_bps(feature,pad):
    name_fields = feature.name.split(",")
    sv_type = name_fields[1]
    sv_methods = name_fields[3]
    sv_length = int(name_fields[2])
    if sv_type  == "INS":
        return feature
    info = json.loads(base64.b64decode(name_fields[0]))
    if "SC_SUBINTERVAL_INFOs" in info:
        L_bps=[]
        R_bps=[]
        for source in info["SC_SUBINTERVAL_INFOs"]:
            sc_bp=sum(map(lambda x:int(x),source['SC_SC_BP_ENDS'].split('-')))/2
            other_bp=sum(map(lambda x:int(x),source['SC_OTHER_BP_ENDS'].split('-')))/2
            if abs(sc_bp-feature.start) < (abs(sc_bp-feature.end)) and abs(other_bp-feature.start) >= (abs(other_bp-feature.end)):
                L_bps.append(sc_bp)
            elif abs(sc_bp-feature.start) >= (abs(sc_bp-feature.end)) and abs(other_bp-feature.start) < (abs(other_bp-feature.end)):
                R_bps.append(sc_bp)
    
        if L_bps and R_bps:
            L_bp=sum(L_bps)/len(L_bps)
            R_bp=sum(R_bps)/len(R_bps)
            sv_length = R_bp-L_bp
            return pybedtools.Interval(feature.chrom, L_bp, R_bp, name="%s,%s,%d,%s" % (
                    base64.b64encode(json.dumps(info)), sv_type, sv_length,sv_methods),    score = feature.score , otherfields=feature.fields[6:])
        else:
            #will be omitted
            return pybedtools.Interval(feature.chrom, feature.start, feature.end, name=feature.name,    score = "-1", otherfields=feature.fields[6:])
    else:
        return feature
    
        
def add_INS_padding(feature,pad):
    name_fields = feature.name.split(",")
    sv_type = name_fields[1]
    return pybedtools.Interval(feature.chrom, max(feature.start-pad,0),
         feature.end+pad, name=feature.name, score = feature.score) if sv_type  == "INS" else feature
    
def find_coverage_frac(score,coverage):
    scores = map(lambda x: float(x),score.split(","))
    coverages = map(lambda x: float(x),coverage.split(","))
    return sum(map(lambda k, v: k/v , scores, coverages))/len(scores)
    
    
        
def add_neighbour_support(feature,bam_handle, min_mapq=SC_MIN_MAPQ,
                          min_soft_clip=SC_MIN_SOFT_CLIP, max_nm=SC_MAX_NM, 
                          min_matches=SC_MIN_MATCHES,skip_soft_clip=False, 
                          isize_mean=ISIZE_MEAN,
                          min_isize=ISIZE_MEAN-2*ISIZE_SD, max_isize=ISIZE_MEAN+2*ISIZE_SD, 
                          max_dist_sc= 100, max_dist_other_bp = 500, wiggle = 20):
    name_fields = feature.name.split(",")
    svtype = name_fields[1]
    sv_methods = name_fields[3]
    sv_length = int(name_fields[2])
    info = json.loads(base64.b64decode(name_fields[0]))

    other_bp_start, other_bp_end = map(int,info["SC_OTHER_BP_ENDS"].split('-'))
    num_neigh_support = 0
    soft_clip_location = (feature.start+feature.end)/2
    for aln in bam_handle.fetch(reference=feature.chrom, start=feature.start, end=feature.end):            
        if not is_good_candidate(aln, min_mapq=min_mapq,
                                 min_soft_clip=min_soft_clip, max_nm=max_nm,
                                 min_matches=min_matches,skip_soft_clip=skip_soft_clip, 
                                 good_neigh_check= True): continue

        svtype_neigh = infer_svtype(aln, min_isize, max_isize)
        
        if svtype_neigh == "CTX;INS":
            # TODO : Should be fixed to handle CTX
            svtype_neigh = "INS"

        if svtype_neigh == "DUP;ITX":
            # TODO : Should be fixed to handle ITX
            svtype_neigh = "DUP"


        if svtype != svtype_neigh:
            continue
        
        soft_clip, dist_L_end, dist_R_end = [-1, -1, -1]
        if not skip_soft_clip:
            soft_clip_tuple = find_softclip(aln)
            if soft_clip_tuple:
                soft_clip, dist_L_end, dist_R_end = soft_clip_tuple         
                soft_clip_location = sum(get_interval(aln))/2
                if abs(soft_clip_location -(feature.start+feature.end)/2) > max_dist_sc:
                    continue
        other_bp_neigh = find_other_bp(aln,isize_mean, svtype_neigh,
                                             soft_clip, dist_L_end,
                                             dist_R_end, soft_clip_location, 
                                             skip_soft_clip=skip_soft_clip, 
                                             skip_neigh=False, wiggle = wiggle)
        if other_bp_neigh is None: continue        
        if not svtype == "INS":
            if abs(other_bp_neigh -(other_bp_start+other_bp_end)/2) > max_dist_other_bp:
                continue

        num_neigh_support +=1
        
    info.update({"SC_NEIGH_SUPPORT": num_neigh_support})
    name = "%s,%s,%d,%s"%(base64.b64encode(json.dumps(info)),svtype,sv_length,sv_methods)

    return pybedtools.Interval(feature.chrom, feature.start, feature.end, name=name, score=feature.score,
                               otherfields=feature.fields[6:]+[str(num_neigh_support)])        
        
def generate_sc_intervals(bam, chromosome, workdir, min_avg_base_qual=SC_MIN_AVG_BASE_QUAL, min_mapq=SC_MIN_MAPQ,
                          min_soft_clip=SC_MIN_SOFT_CLIP,
                          pad=SC_PAD, min_support=MIN_SUPPORT, max_considered_isize=1000000000, 
                          min_support_frac=MIN_SUPPORT_FRAC_INS, max_nm=SC_MAX_NM, min_matches=SC_MIN_MATCHES, 
                          isize_mean=ISIZE_MEAN, isize_sd=ISIZE_SD, svs_to_softclip=SVS_SOFTCLIP_SUPPORTED,
                          overlap_ratio=OVERLAP_RATIO,merge_max_dist=-int(1*SC_PAD), 
                          mean_read_length=MEAN_READ_LENGTH, mean_read_coverage=MEAN_READ_COVERAGE, 
                          min_ins_cov_frac=MIN_INS_COVERAGE_FRAC, max_ins_cov_frac=MAX_INS_COVERAGE_FRAC,
                          num_sd = 2):
    func_logger = logging.getLogger("%s-%s" % (generate_sc_intervals.__name__, multiprocessing.current_process()))

    if not os.path.isdir(workdir):
        func_logger.error("Working directory %s doesn't exist" % workdir)
        return None

    func_logger.info("Generating candidate intervals from %s for chromsome %s" % (bam, chromosome))
    pybedtools.set_tempdir(workdir)

    
    min_isize = isize_mean - num_sd * isize_sd
    max_isize = isize_mean + num_sd * isize_sd

    unmerged_intervals = []
    start_time = time.time()
    ignore_none = True
    try:
        sam_file = pysam.Samfile(bam, "rb")
        for aln in sam_file.fetch(reference=chromosome):
            if abs(aln.tlen) > max_considered_isize:
                continue
            if not is_good_candidate(aln, min_avg_base_qual=min_avg_base_qual, min_mapq=min_mapq,
                                     min_soft_clip=min_soft_clip, max_nm=max_nm,
                                     min_matches=min_matches): continue
            interval = get_interval(aln, pad=pad)
            soft_clip_location = sum(interval) / 2
            strand = "-" if aln.is_reverse else "+"
            svtype = infer_svtype(aln, min_isize, max_isize)
            
            if svtype == "CTX;INS":
                # TODO : Should be fixed to handle CTX
                svtype = "INS"
            
            if svtype == "DUP;ITX":
                # TODO : Should be fixed to handle ITX
                svtype = "DUP"

            soft_clip_tuple = find_softclip(aln)
            if not soft_clip_tuple:    
                continue
            soft_clip, dist_L_end, dist_R_end = soft_clip_tuple 
            other_bp = find_other_bp(aln,isize_mean, svtype, soft_clip, dist_L_end,
                                                 dist_R_end, soft_clip_location)
            if other_bp is None: continue

            name = "%d,%d,%s" % (soft_clip_location, other_bp, strand)
            if ignore_none and svtype == "NONE":
                continue
            if svtype not in svs_to_softclip:
                continue

            unmerged_intervals.append(
                pybedtools.Interval(chromosome, interval[0], interval[1], name=name, score="1", strand=strand, otherfields=[svtype]))

        if not unmerged_intervals:
            sam_file.close()
            func_logger.warn("No intervals generated")
            return None

        unmerged_bed = os.path.join(workdir, "unmerged.bed")
        bedtool = pybedtools.BedTool(unmerged_intervals).sort().moveto(unmerged_bed)
        func_logger.info("%d candidate reads" % (bedtool.count()))



        merged_bed = os.path.join(workdir, "merged.bed")
        m_bedtool=merge_for_each_sv(bedtool,c="4,5,6,7",o="collapse,sum,collapse,collapse",
                                    svs_to_softclip=svs_to_softclip,d=merge_max_dist,
                                    reciprocal_for_2bp=False, sv_type_field = [6,0])
        m_bedtool = m_bedtool.moveto(merged_bed)
        func_logger.info("%d merged intervals" % (m_bedtool.count()))



        # Check if the other break point also can be merged for the merged intervals (for 2bp SVs)
        bp_merged_intervals = []
        for interval in m_bedtool:
            sv_type = interval.fields[6].split(',')[0]
            if len(set(interval.fields[6].split(',')))!=1:
                func_logger.warn("More than one svtypes: %s",(str(interval)))
            if  sv_type == "INS":
                bp_merged_intervals.append(add_other_bp_fields(interval,0,0))
            else:
                other_bp_bedtool=bedtool.filter(lambda x: x.name in interval.name and x.fields[6]==sv_type).each(partial(generate_other_bp_interval,pad=pad)).sort().merge(c="4,5,6,7", o="collapse,sum,collapse,collapse", d=merge_max_dist)
                for intvl in other_bp_bedtool:
                    name_fields = intvl.name.split(',')
                    sv_length = [abs(int(name_fields[3*i])-int(name_fields[3*i+1])) for i in range(len(name_fields)/3)]
                    sv_length = sum(sv_length)/len(sv_length)
                    if sv_length < pad:
                        L_bp_intvls=[]
                        R_bp_intvls=[]
                        for i in range(len(name_fields)/3):
                            if abs(int(name_fields[3*i])-(intvl.start+pad)) < abs(int(name_fields[3*i])-(intvl.end-pad)):
                                L_bp_intvls.append(",".join(name_fields[3*i:3*(i+1)]))
                            else:
                                R_bp_intvls.append(",".join(name_fields[3*i:3*(i+1)]))
                        if L_bp_intvls and R_bp_intvls:
                            L_merged_other_bp=max(map(lambda x: int(x.split(',')[1]),L_bp_intvls))
                            R_merged_other_bp=max(map(lambda x: int(x.split(',')[1]),R_bp_intvls))
                            L_other_end=max(map(lambda x: x.split(',')[1],L_bp_intvls))
                            bp_merged_intervals.extend(bedtool.filter(lambda x: x.name in L_bp_intvls and x.fields[6]==sv_type).sort().merge(c="4,5,6,7", o="collapse,sum,collapse,collapse", d=merge_max_dist).each(partial(add_other_bp_fields, start=max(L_merged_other_bp-pad,0), end=L_merged_other_bp+pad)).intervals)
                            bp_merged_intervals.extend(bedtool.filter(lambda x: x.name in R_bp_intvls and x.fields[6]==sv_type).sort().merge(c="4,5,6,7", o="collapse,sum,collapse,collapse", d=merge_max_dist).each(partial(add_other_bp_fields, start=max(R_merged_other_bp-pad,0), end=R_merged_other_bp+pad)).intervals)
                        else:
                            bp_merged_intervals.extend(bedtool.filter(lambda x: x.name in intvl.name and x.fields[6]==sv_type).sort().merge(c="4,5,6,7", o="collapse,sum,collapse,collapse", d=merge_max_dist).each(partial(add_other_bp_fields, start=intvl.start, end=intvl.end)).intervals)
                    else: 
                        bp_merged_intervals.extend(bedtool.filter(lambda x: x.name in intvl.name and x.fields[6]==sv_type).sort().merge(c="4,5,6,7", o="collapse,sum,collapse,collapse", d=merge_max_dist).each(partial(add_other_bp_fields, start=intvl.start, end=intvl.end)).intervals)
        
        bp_merged_bed = os.path.join(workdir, "bp_merged.bed")
        bedtool=pybedtools.BedTool(bp_merged_intervals).sort().moveto(bp_merged_bed)       
        func_logger.info("%d BP merged intervals" % (bedtool.count()))

        filtered_bed = os.path.join(workdir, "filtered_bp_merged.bed")
        bedtool = bedtool.filter(lambda x: int(x.score) >= min_support).each(
            partial(merged_interval_features, bam_handle=sam_file)).moveto(
            filtered_bed)
        func_logger.info("%d filtered intervals" % (bedtool.count()))
        
        # Now filter based on coverage
        coverage_filtered_bed = os.path.join(workdir, "coverage_filtered_bp_merged.bed")
        bedtool = bedtool.filter(lambda x: (x.fields[3].split(",")[1]!="INS" or 
                                           ((min_ins_cov_frac*mean_read_coverage)<=(float(x.fields[6])/abs(x.start-x.end+1)*mean_read_length)<=(max_ins_cov_frac*mean_read_coverage)))).moveto(coverage_filtered_bed)
        func_logger.info("%d coverage filtered intervals" % (bedtool.count()))


        thr_sv={"INS":MIN_SUPPORT_FRAC_INS, "INV":MIN_SUPPORT_FRAC_INV, 
                "DEL":MIN_SUPPORT_FRAC_DEL, "DUP": MIN_SUPPORT_FRAC_DUP}

        # Add number of neighbouring reads that support SC
        bedtool=bedtool.each(partial(add_neighbour_support,bam_handle=sam_file, min_mapq=min_mapq, 
                                     min_soft_clip=min_soft_clip, max_nm=max_nm, min_matches=min_matches,
                                     skip_soft_clip=False, isize_mean=isize_mean, min_isize=min_isize, max_isize=max_isize)).sort().moveto(coverage_filtered_bed)

        neigh_coverage_filtered_bed = os.path.join(workdir, "neigh_filtered_bp_merged.bed")
        bedtool = bedtool.filter(lambda x: (float(x.fields[6]) * thr_sv[x.fields[3].split(",")[1]] <= float(x.fields[8]))).moveto(neigh_coverage_filtered_bed)
        func_logger.info("%d neighbour support filtered intervals" % (bedtool.count()))

        # For 2bp SVs, the interval will be the cover of two intervals on the BP
        full_filtered_bed = os.path.join(workdir, "full_filtered_bp_merged.bed")
        bedtool = bedtool.each(partial(get_full_interval,pad=pad)).sort().moveto(full_filtered_bed)
        func_logger.info("%d full filtered intervals" % (bedtool.count()))

        # Now merge on full intervals
        merged_full_filtered_bed = os.path.join(workdir, "merged_full_filtered_bp_merged.bed")
        if bedtool.count()>0:
            bedtool=merge_for_each_sv(bedtool,c="4,5,6,7,9",o="collapse,collapse,collapse,collapse,collapse",
                                      svs_to_softclip=svs_to_softclip,
                                      overlap_ratio=overlap_ratio,
                                      reciprocal_for_2bp=True, 
                                      sv_type_field = [3,1], d=merge_max_dist)
        bedtool = bedtool.each(partial(fix_merged_fields,inter_tools=False)).each(partial(fine_tune_bps,pad=pad))
        bedtool = bedtool.filter(lambda x: x.score != "-1").sort().moveto(merged_full_filtered_bed)
        func_logger.info("%d merged full intervals" % (bedtool.count()))

        sam_file.close()
    except Exception as e:
        func_logger.error('Caught exception in worker thread')

        # This prints the type, value, and stack trace of the
        # current exception being handled.
        traceback.print_exc()

        print()
        raise e

    pybedtools.cleanup(remove_all=True)
    func_logger.info("Generated intervals in %g seconds for region %s" % ((time.time() - start_time), chromosome))

    return merged_full_filtered_bed


def parallel_generate_sc_intervals(bams, chromosomes, skip_bed, workdir, num_threads=1,
                                   min_avg_base_qual=SC_MIN_AVG_BASE_QUAL,
                                   min_mapq=SC_MIN_MAPQ, min_soft_clip=SC_MIN_SOFT_CLIP,
                                   pad=SC_PAD,
                                   min_support=MIN_SUPPORT, min_support_frac=MIN_SUPPORT_FRAC_INS, 
                                   max_intervals=MAX_INTERVALS, max_nm=SC_MAX_NM, min_matches=SC_MIN_MATCHES, 
                                   isize_mean=ISIZE_MEAN, isize_sd=ISIZE_SD,
                                   svs_to_softclip=SVS_SOFTCLIP_SUPPORTED,
                                   overlap_ratio=OVERLAP_RATIO, mean_read_length=MEAN_READ_LENGTH,
                                   mean_read_coverage=MEAN_READ_COVERAGE, min_ins_cov_frac=MIN_INS_COVERAGE_FRAC,
                                   max_ins_cov_frac=MAX_INS_COVERAGE_FRAC):
    func_logger = logging.getLogger(
        "%s-%s" % (parallel_generate_sc_intervals.__name__, multiprocessing.current_process()))

    if not os.path.isdir(workdir):
        func_logger.info("Creating directory %s" % workdir)
        os.makedirs(workdir)

    if not chromosomes:
        func_logger.info("Chromosome list unspecified. Inferring from the BAMs")
        for bam in bams:
            bamfile = pysam.Samfile(bam, "rb")
            chromosomes += list(bamfile.references)
            bamfile.close()
        chromosomes = sorted(list(set(chromosomes)))
        func_logger.info("Chromosome list inferred as %s" % (str(chromosomes)))

    if not chromosomes:
        func_logger.error("Chromosome list empty")
        return None


    merge_max_dist = -int(1 * pad)


    func_logger.info("SVs to soft-clip: %s" % (svs_to_softclip))

    pool = multiprocessing.Pool(num_threads)

    bed_files = []
    for index, (bam, chromosome) in enumerate(itertools.product(bams, chromosomes)):
        process_workdir = os.path.join(workdir, str(index))
        if not os.path.isdir(process_workdir):
            os.makedirs(process_workdir)

        args_list = [bam, chromosome, process_workdir]
        kwargs_dict = {"min_avg_base_qual": min_avg_base_qual, "min_mapq": min_mapq, "min_soft_clip": min_soft_clip,
                       "pad": pad, "min_support": min_support,
                       "min_support_frac": min_support_frac, "max_nm": max_nm, "min_matches": min_matches, 
                       "isize_mean": isize_mean, "isize_sd": isize_sd, "svs_to_softclip": svs_to_softclip, 
                       "merge_max_dist": merge_max_dist, "mean_read_length": mean_read_length,
                       "mean_read_coverage": mean_read_coverage, "min_ins_cov_frac": min_ins_cov_frac,
                       "max_ins_cov_frac": max_ins_cov_frac}
        pool.apply_async(generate_sc_intervals, args=args_list, kwds=kwargs_dict,
                         callback=partial(generate_sc_intervals_callback, result_list=bed_files))

    pool.close()
    pool.join()

    func_logger.info("Following BED files will be merged: %s" % (str(bed_files)))

    if not bed_files:
        func_logger.warn("No intervals generated")
        return None

    pybedtools.set_tempdir(workdir)
    bedtool = pybedtools.BedTool(bed_files[0])

    for bed_file in bed_files[1:]:
        bedtool = bedtool.cat(pybedtools.BedTool(bed_file), postmerge=False)

    bedtool = bedtool.sort().moveto(os.path.join(workdir, "all_intervals.bed"))

    func_logger.info("Selecting the top %d intervals based on normalized read support" % max_intervals)
    top_intervals_all_cols_file = os.path.join(workdir, "top_intervals_all_cols.bed")
    if bedtool.count() <= max_intervals:
        bedtool = bedtool.saveas(top_intervals_all_cols_file)
    else:
        # Sample the top intervals
        top_fraction_cutoff = \
            sorted([find_coverage_frac(interval.fields[7], interval.fields[6]) for interval in bedtool], reverse=True)[
                max_intervals - 1]
        func_logger.info("Normalized read support threshold: %0.3f" % top_fraction_cutoff)
        bedtool = bedtool.filter(lambda x: find_coverage_frac(x.fields[7],x.fields[6]) >= top_fraction_cutoff).moveto(
            top_intervals_all_cols_file)

    # Filter out the extra column added to simplify life later on
    bedtool = bedtool.cut(xrange(6)).saveas(os.path.join(workdir, "top_intervals.bed"))

    interval_bed = os.path.join(workdir, "intervals.bed")
    if skip_bed:        
        skip_bedtool = pybedtools.BedTool(skip_bed)
        sc_skip_bed = os.path.join(workdir, "sc_metasv.bed")
        if "INS" in svs_to_softclip:
            skip_bedtool = skip_bedtool.each(partial(add_INS_padding,pad=pad)).saveas(sc_skip_bed)
        nonsc_skip_bed = os.path.join(workdir, "non_sc_metasv.bed")
        func_logger.info(
            "Merging %d features with %d features from %s" % (bedtool.count(), skip_bedtool.count(), skip_bed))
        nonsc_skip_bedtool = skip_bedtool.filter(lambda x: x.name.split(',')[1] not in svs_to_softclip).saveas(nonsc_skip_bed)
        sc_skip_bedtool = skip_bedtool.filter(lambda x: x.name.split(',')[1] in svs_to_softclip).saveas(interval_bed)
        bedtool = bedtool.cat(sc_skip_bedtool, postmerge=False)
        bedtool = bedtool.sort()
        bedtool = merge_for_each_sv(bedtool,c="4",o="collapse",svs_to_softclip=svs_to_softclip,
                                  overlap_ratio=overlap_ratio, reciprocal_for_2bp=True, d=merge_max_dist)
        bedtool = bedtool.each(partial(fix_merged_fields,inter_tools=True)).sort().moveto(interval_bed)
        bedtool = bedtool.cat(nonsc_skip_bedtool, postmerge=False).sort().moveto(interval_bed)
        func_logger.info("After merging with %s %d features" % (skip_bed, bedtool.count()))
    else:    
        bedtool = bedtool.saveas(interval_bed)

    pybedtools.cleanup(remove_all=True)

    return bedtool.fn


if __name__ == "__main__":
    FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
    logging.basicConfig(level=logging.INFO, format=FORMAT)
    logger = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(
        description="Generate BED intervals for insertion detection using soft-clipped reads",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bams", nargs="+", help="BAMs", required=True)
    parser.add_argument("--chromosomes", nargs="+", help="Chromosomes", default=[])
    parser.add_argument("--workdir", help="Working directory", default="work")
    parser.add_argument("--num_threads", help="Number of threads to use", default=1, type=int)
    parser.add_argument("--min_avg_base_qual", help="Minimum average base quality", default=SC_MIN_AVG_BASE_QUAL,
                        type=int)
    parser.add_argument("--min_mapq", help="Minimum MAPQ", default=SC_MIN_MAPQ, type=int)
    parser.add_argument("--min_soft_clip", help="Minimum soft-clip", default=SC_MIN_SOFT_CLIP, type=int)
    parser.add_argument("--max_nm", help="Maximum number of edits", default=SC_MAX_NM, type=int)
    parser.add_argument("--min_matches", help="Minimum number of matches", default=SC_MIN_MATCHES, type=int)
    parser.add_argument("--isize_mean", help="Insert-size mean", default=ISIZE_MEAN, type=float)
    parser.add_argument("--isize_sd", help="Insert-size s.d.", default=ISIZE_SD, type=float)
    parser.add_argument("--pad", help="Padding on both sides of the candidate locations", default=SC_PAD, type=int)
    parser.add_argument("--min_support", help="Minimum supporting reads", default=MIN_SUPPORT, type=int)
    parser.add_argument("--min_support_frac", help="Minimum fraction of total reads for interval",
                        default=MIN_SUPPORT_FRAC_INS, type=float)
    parser.add_argument("--skip_bed", help="BED regions with which no overlap should happen", type=file)
    parser.add_argument("--max_intervals",
                        help="Maximum number of intervals to process. Intervals are ranked by normalized read-support",
                        type=int, default=MAX_INTERVALS)
    parser.add_argument("--svs_to_softclip", nargs="+", help="SVs to perform soft-clip analysis on", default=SVS_SOFTCLIP_SUPPORTED,
                           choices=SVS_SOFTCLIP_SUPPORTED)
    parser.add_argument("--overlap_ratio", help="Reciprocal overlap ratio", default=OVERLAP_RATIO, type=float,
                                required=False)
    parser.add_argument("--mean_read_length", type=float, default=MEAN_READ_LENGTH, help="Mean read length")
    parser.add_argument("--mean_read_coverage", type=float, default=MEAN_READ_COVERAGE, help="Mean read coverage")
    parser.add_argument("--min_ins_cov_frac", type=float, default=MIN_INS_COVERAGE_FRAC, help="Minimum read coverage around the insertion breakpoint.")
    parser.add_argument("--max_ins_cov_frac", type=float, default=MAX_INS_COVERAGE_FRAC, help="Maximum read coverage around the insertion breakpoint.")

    args = parser.parse_args()

    logger.info("Command-line: " + " ".join(sys.argv))

    parallel_generate_sc_intervals(args.bams, args.chromosomes, args.skip_bed, args.workdir,
                                   num_threads=args.num_threads, min_avg_base_qual=args.min_avg_base_qual,
                                   min_mapq=args.min_mapq, min_soft_clip=args.min_soft_clip,
                                   pad=args.pad, min_support=args.min_support,
                                   min_support_frac=args.min_support_frac, max_intervals=args.max_intervals,
                                   max_nm=args.max_nm, min_matches=args.min_matches, isize_mean=args.isize_mean, 
                                   isize_sd=args.isize_sd, svs_to_softclip=args.svs_to_softclip, 
                                   overlap_ratio=args.overlap_ratio, mean_read_length=args.mean_read_length,
                                   mean_read_coverage=args.mean_read_coverage, min_ins_cov_frac=args.min_ins_cov_frac,
                                   max_ins_cov_frac=args.max_ins_cov_frac)
