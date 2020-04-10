#!/usr/bin/env python
#-*- coding: utf-8 -*

from pyrosetta import *
from pyrosetta.rosetta.protocols.sewing.movers import AssemblyMover
from pyrosetta.rosetta.protocols.sewing.requirements import *
from pyrosetta.rosetta.protocols.sewing.scoring import *
from pyrosetta.rosetta.protocols import rosetta_scripts

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

xml = rosetta_scripts.XmlObjects.create_from_string('''
<MOVERS>
	<AssemblyMover
    	name="assemble"
    	minimum_cycles="10000"
    	maximum_cycles="100000"
    	start_temperature="0.6"
    	end_temperature="0.6"
    	hashed="false"
    	window_width="4"
    	model_file_name="smotifs_H_5_40_L_2_6_H_5_40.segments">
    <AssemblyScorers>
        <MotifScorer weight="1" />
        <InterModelMotifScorer weight="10" />
	</AssemblyScorers>
	<AssemblyRequirements>
        <ClashRequirement clash_radius="5" />
        <SizeInSegmentsRequirement minimum_size="4" maximum_size="5" />
        <DsspSpecificLengthRequirement dssp_code="H" minimum_length="12" maximum_length="1000" />
	</AssemblyRequirements>
</AssemblyMover>
</MOVERS>
''').get_mover('assemble')
#

xml.apply(pose)
