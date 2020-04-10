#!/usr/bin/env python
#-*- coding: utf-8 -*

from pyrosetta import *
from pyrosetta.rosetta.protocols.sewing.movers import AssemblyMover
from pyrosetta.rosetta.protocols.sewing.requirements import *
from pyrosetta.rosetta.protocols.sewing.scoring import *
from pyrosetta.rosetta.protocols.sewing.hashing import *
from pyrosetta.rosetta.protocols.sewing.data_storage import *

args = '''
-ignore_unrecognized_res
-detect_disulf false
-mh:match:ss1 true
-mh:match:ss2 true
-mh:match:aa1 false
-mh:match:aa2 false
-mh:score:use_ss1 true
-mh:score:use_ss2 true
-mh:score:use_aa1 false
-mh:score:use_aa2 false
-mh:path:motifs /usr/local/rosetta_src_2019.47.61047_bundle/main/database/additional_protocol_data/sewing/xsmax_bb_ss_AILV_resl0.8_msc0.3/xsmax_bb_ss_AILV_resl0.8_msc0.3.rpm.bin.gz
-mh:path:scores_BB_BB /usr/local/rosetta_src_2019.47.61047_bundle/main/database/additional_protocol_data/sewing/xsmax_bb_ss_AILV_resl0.8_msc0.3/xsmax_bb_ss_AILV_resl0.8_msc0.3
-mh:gen_reverse_motifs_on_load false
-mh::dump::max_rms 0.4
-mh::dump::max_per_res 20
'''

init(args)
pose = pose_from_sequence('A', "fa_standard")

ab = AssemblyMover()
ab.set_min_cycles(10000)
ab.set_max_cycles(100000)
ab.set_start_temperature(0.6)
ab.set_end_temperature(0.6)
ab.set_hashed(False)
ab.set_window_width(4)
ab.set_model_file_name('smotifs_H_5_40_L_2_6_H_5_40.segments')

# set scores;
ab.set_default_assembly_scorers()
ab.set_default_requirements()

#
ab.apply(pose)
