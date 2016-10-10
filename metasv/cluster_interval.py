import os
import logging
import re
from collections import defaultdict
import intervaltree
import pybedtools
import itertools
import argparse
import multiprocessing

def makedirs(dirlist):
    logger = logging.getLogger(makedirs.__name__)
    for dirname in dirlist:
        if not os.path.isdir(dirname):
            logger.info("Creating {}".format(dirname))
            os.makedirs(dirname)
        else:
            logger.warn("{} exists.".format(dirname))


def good_overlap(interval1, interval2, wiggle=100, overlap=0.5):
    best_shift = max(-wiggle, min(wiggle, interval2.begin - interval1.begin))
    overlap_length = min(interval1.end + best_shift, interval2.end) - max(interval1.begin + best_shift, interval2.begin)
    length1 = interval1.end - interval1.begin
    length2 = interval2.end - interval2.begin
    return overlap_length >= float(max(length1, length2)) * overlap


def check_range(s, value_type, value_range):
    value = value_type(s)
    if not (value_range[0] <= value <= value_range[1]):
        raise argparse.ArgumentTypeError("{} is outside the range [{}, {}].".format(s, value_range[0], value_range[1]))
    return value


def merge_clusters(interval1, interval2):
    return intervaltree.Interval(min(interval1.begin, interval2.begin), max(interval1.end, interval2.end), data=(min(interval1.data[0], interval2.data[0]), interval1.data[1] + interval2.data[1]))


def cluster_interval_list(intervals, clustered_intervals, wiggle=100, overlap=0.5, init_cluster_id=0, min_cluster_size=1, output_cluster_size=False):
   logger = logging.getLogger(cluster_interval_list.__name__ + str(multiprocessing.current_process()))

   if not intervals:
       return init_cluster_id

   if len(intervals) == 1:
       clustered_intervals.append(pybedtools.create_interval_from_list(intervals[0].fields + [str(init_cluster_id)] + (["1"] if output_cluster_size else [])))
       return init_cluster_id + 1

   intervals.sort(key=lambda x: (x.start, x.end))

   # Do a simple clustering to prune the interval set
   logger.info("Doing simple initial clustering")
   last_interval = intervals[0]
   current_cluster = []
   clusters = []
   for index, interval in enumerate(itertools.chain.from_iterable([intervals, [None]])):
       if last_interval and (not interval or last_interval.start < interval.start - 20 or abs(last_interval.end - interval.end) > 20):
           if current_cluster:
               start = min(map(lambda x: intervals[x].start, current_cluster))
               end = max(map(lambda x: intervals[x].end, current_cluster))
               clusters.append(intervaltree.Interval(start, end, data=(0, current_cluster)))
           current_cluster = []
           last_interval = interval
       if interval and interval.length:
           current_cluster.append(index)
   logger.info("Created {} initial clusters from {} intervals".format(len(clusters), len(intervals)))

   num_iterations = 0
   while True:
       clustered_indices = set()
       # Build interval tree using the current clustering
       interval_tree = intervaltree.IntervalTree()
       for index, cluster in enumerate(clusters):
           interval_tree.addi(cluster.begin, cluster.end, (index, cluster.data[1]))

       next_clusters = []
       for index, cluster in enumerate(clusters):
           if index in clustered_indices: continue

           clustered_indices.add(index)

           interval_tree.discard(cluster)

           overlapping_clusters = interval_tree[cluster.begin: cluster.end]
           merged_cluster = cluster
           for overlapping_cluster in overlapping_clusters:
               overlapping_index = overlapping_cluster.data[0]
               if overlapping_index in clustered_indices:
                   continue
               if good_overlap(merged_cluster, overlapping_cluster, wiggle=wiggle, overlap=overlap):
                   merged_cluster = merge_clusters(merged_cluster, overlapping_cluster)
                   clustered_indices.add(overlapping_index)
                   interval_tree.discard(overlapping_cluster)
           next_clusters.append(merged_cluster)
       num_iterations += 1
       if len(next_clusters) == len(clusters):
           break
       clusters = next_clusters

   clusters.sort(key=lambda x: (x.begin, x.end))

   logger.info("Clustering converged in {} iterations and created {} clusters".format(num_iterations, len(clusters)))
   for cluster in clusters:
       #logger.info("Cluster {} with interval {}:{}-{} and {} members".format(init_cluster_id, intervals[0].chrom, cluster.begin, cluster.end, len(cluster.data[1])))
       for interval_index in cluster.data[1]:
           clustered_intervals.append(pybedtools.create_interval_from_list(intervals[interval_index].fields + [str(init_cluster_id)] + ([str(len(cluster.data[1]))] if output_cluster_size else [])))
       init_cluster_id += 1

   return init_cluster_id


def cluster_intervals((bed, work_dir, init_cluster_id, wiggle, overlap, output_cluster_size)):
    logger = logging.getLogger(cluster_intervals.__name__ + str(multiprocessing.current_process()))

    bedtool = pybedtools.BedTool(bed)

    prev_interval = None
    clustered_intervals = []
    intervals = []
    current_cluster_id = init_cluster_id
    for interval in itertools.chain.from_iterable([bedtool, [None]]):
        if prev_interval and (not interval or prev_interval.fields[-1] != interval.fields[-1]):
            current_cluster_id = cluster_interval_list(intervals, clustered_intervals, wiggle=wiggle, overlap=overlap, init_cluster_id=current_cluster_id, output_cluster_size=output_cluster_size)
            intervals = []
        if interval:
            intervals.append(interval)
        prev_interval = interval

    logger.info("Clustering done. Generating the BED")
    clustered_bed = os.path.join(work_dir, "cluster.bed")
    pybedtools.BedTool(clustered_intervals).saveas(clustered_bed)

    return clustered_bed
    

def cluster_intervals_parallel(bedtool, work_dir, wiggle=100, overlap=0.5, output_cluster_size=False, nthreads=1):
    logger = logging.getLogger(cluster_intervals.__name__)

    makedirs([work_dir])

    logger.info("Simple merge based clustering using bedtools first")
    initial_cluster_bedtool = bedtool.cluster().saveas(os.path.join(work_dir, "bedtools_cluster.bed"))
    logger.info("Clustering intervals from {}".format(initial_cluster_bedtool.fn))

    per_thread = (len(initial_cluster_bedtool) + nthreads - 1) // nthreads

    # Next we split among multiple processes on a cluster level granularity
    prev_interval = None
    intervals = []
    partition_id = 0
    map_args = []
    init_cluster_id = 0
    for interval in itertools.chain.from_iterable([initial_cluster_bedtool, [None]]):
        if prev_interval and (not interval or prev_interval.fields[-1] != interval.fields[-1]) and (len(intervals) >= per_thread or not interval):
            # Emit the current list of intervals for processing in a single process
            process_workdir = os.path.join(work_dir, str(partition_id))
            process_bed = os.path.join(process_workdir, "intervals.bed")
            makedirs([process_workdir])
            pybedtools.BedTool(intervals).saveas(process_bed)
            map_args.append((process_bed, process_workdir, init_cluster_id, wiggle, overlap, output_cluster_size))
            init_cluster_id += len(intervals)
            intervals = []
            partition_id += 1
        intervals.append(interval)
        prev_interval = interval

    pool = multiprocessing.Pool(nthreads)
    cluster_output = pool.map_async(cluster_intervals, map_args).get()
    pool.close()

    # Now merge the per-process outputs
    logger.info("Merging clusters from {}".format(" ".join(cluster_output)))
    prev_interval = None
    cluster_id = 1
    final_intervals = []
    for interval in itertools.chain.from_iterable(map(pybedtools.BedTool, cluster_output)):
        if prev_interval and prev_interval.fields[-1] != interval.fields[-1]:
            cluster_id += 1
        final_intervals.append(pybedtools.create_interval_from_list(interval.fields[:-1] + [str(cluster_id)]))
        prev_interval = interval
    logger.info("Merged with {} final clusters".format(cluster_id))

    return pybedtools.BedTool(final_intervals)
            

def multisample_cluster_intervals(bedtools, samples, work_dir, reference, wiggle=100, overlap=0.9, nthreads=1):
    makedirs([work_dir])
    if not samples:
        samples = map(lambda x: x.fn, bedtools)
    reference_faidx = "{}.fai".format(reference)

    cat_intervals = []
    for sample, bedtool in zip(samples, bedtools):
        for interval in bedtool:
            cat_intervals.append(pybedtools.create_interval_from_list(interval.fields[:3] + ["{},{}".format(sample, interval.fields[3])] + interval.fields[4:]))

    cat_bedtool = pybedtools.BedTool(cat_intervals).saveas(os.path.join(work_dir, "cat.bed"))
    sorted_bedtool = cat_bedtool.sort(faidx=reference_faidx).saveas(os.path.join(work_dir, "sorted.bed"))
    return cluster_intervals_parallel(sorted_bedtool, work_dir, wiggle=wiggle, overlap=overlap, output_cluster_size=True, nthreads=nthreads)
