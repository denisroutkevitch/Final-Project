"""
This program folds protein 2reb using the
Monte Carlo with minimization method using
9-mer then 3-mer fragments for a low resolution
move, then small, shear, and minimization moves
for the high resolution move
"""

from pyrosetta import *
from pyrosetta import PyMOLMover
from pyrosetta.toolbox import mutate_residue
from pyrosetta.teaching import ScoreFunction
from rosetta.core.fragment import ConstantLengthFragSet
from rosetta.protocols.simple_moves import *
from rosetta.protocols.moves import *
import rosetta.protocols.rigid as rigid_moves
from rosetta.protocols.docking import *

init()

pmm = PyMOLMover()
pmm.keep_history(True)

pose = pose_from_pdb('3o2b.clean.pdb')
mutate_residue(pose, 85, 'B96')

poseinit = Pose()
poseinit.assign(pose)


# algorithm parameters
traj = 1
iterate = 25
jump_num = 1

"""
pose = Pose()
pose = pose_from_sequence(seq, 'centroid')
orig_pose = Pose()
lowest_pose = Pose()
low_pose = Pose()
orig_pose.assign(pose)
lowest_pose.assign(pose)
low_pose.assign(pose)
"""

# used mm_std score function
scorefxn = create_score_function('mm_std')


print 'Loading Movers'
movemap = MoveMap()
movemap.set_jump(jump_num, True)


pert_mover = rigid_moves.RigidBodyPerturbMover(jump_num, 8, 3)

# this was cuasing problems and i
# dont think its necessary for our small moves
# slide = DockingSlideIntoContact(jump_num) # for centroid mode
# slide = FaDockingSlideIntoContact(jump_num) # for fullatom mode


task_pack = standard_packer_task(pose)
task_pack.restrict_to_repacking()
# task_pack.temporarily_fix_everything()
# task_pack.temporarily_set_pack_residue(49, True)

pack_mover = PackRotamersMover(scorefxn, task_pack)


min_mover = MinMover()
min_mover.movemap(movemap)
min_mover.score_function(scorefxn) # use any scorefxn


mc = MonteCarlo(pose, scorefxn, 1.0)


dock_mover = SequenceMover()
dock_mover.add_mover(pert_mover)
# dock_mover.add_mover(slide)
dock_mover.add_mover(pack_mover)
dock_mover.add_mover(min_mover)


trial_dock = TrialMover(dock_mover, mc)
docker = RepeatMover(trial_dock, 10)


jd = PyJobDistributor("output", 1, scorefxn)
jd.native_pose = poseinit

print 'Starting Dock'

# low res folding algorithm with traj number of trajectories
while not jd.job_complete:
    print 'Starting Fold'
    try:
        pose.assign(poseinit)
        mc.reset(poseinit)
        docker.apply(pose)
        mc.recover_low(pose)

    except RuntimeError :
        print "ERROR:NAN occurred in H-bonding calculations!"
    
    # output
    # jd.additional_decoy_info = " LRMSD: " + str(scorefxn(pose))
    jd.output_decoy(pose)

print 'Finished'




