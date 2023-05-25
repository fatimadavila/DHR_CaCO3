#! /usr/bin/env python3

# Requires to run generate lattice constraints separatedly
############################ 190601 ####################################
# I think this is the third iteration of this protocol, it
# disallows design of surface residues on the tails of non-interface
# side.   -> op_resi_comb_sele_tot
# Also, it deletes the dumping of a commented pose.
########################################################################

# Python
import numpy as np
import argparse
import glob

# PyRosetta
from pyrosetta import *
from pyrosetta.rosetta import *
import pyrosetta.rosetta.core.pose as pose_obj
import pyrosetta.rosetta.core.select.residue_selector as residue_selector
import pyrosetta.toolbox.mutants as mutants

pyrosetta.init(
    extra_options='-symmetry_definition stoopid -old_sym_min true -beta -ex1 true -ex2aro true -lazy_ig False')

# Core Includes
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation

# Protocol Includes
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import task_operations as task_ops
from rosetta.protocols.rosetta_scripts import XmlObjects as xml_object
from rosetta.protocols import constraint_movers


def get_sec_struct_dict(some_pose):
    some_pose.display_secstruct()
    # Split pose in dictionary containing helix number and span of residues
    # corresponding to it
    key_number = 1
    d = {}
    for i in range(some_pose.total_residue()):
        j = i + 1
        if j == 1:
            if some_pose.secstruct(j) is 'H':
                new_key_name = key_number
                l = []
                l.append(j)
                d[new_key_name] = l
                key_number = key_number + 1
        else:
            if some_pose.secstruct(
                    j) is 'H' and some_pose.secstruct(i) is 'H':
                l.append(j)
                d[new_key_name] = l
            elif some_pose.secstruct(j) is 'H' and some_pose.secstruct(i) != 'H':
                new_key_name = key_number
                l = []
                l.append(j)
                d[new_key_name] = l
                key_number = key_number + 1
    return d


def select_correct_helix(ss_dict, side, rep_len):
    ss_dict_set = set(ss_dict[side])
    second_repeat_set = set(
        [i for i in range(rep_len + 1, 2 * rep_len + 1)])
    sorted_and_ssdict_secondrep = list(
        ss_dict_set.intersection(second_repeat_set))
    helix_selector = residue_selector.ResidueIndexSelector(
        ",".join(map(str, sorted_and_ssdict_secondrep)))
    return helix_selector


def select_surf_helix_nn(helix_selector):
    surface_selector_nn = residue_selector.LayerSelector()
    surface_selector_nn.set_use_sc_neighbors(1)
    surface_selector_nn.set_layers(0, 0, 1)
    surface_selector_nn.set_cutoffs(5.2, 2)  # Default 5.2 and 2.0
    # 0.3 selects boundary default 5, larger wider cone
    surface_selector_nn.set_angle_shift_factor(0.5)
    # 2 default smaller sharper and broader cone
    surface_selector_nn.set_angle_exponent(2)
    surface_selector_nn.set_sc_neighbor_dist_midpoint(9)  # Default 9
    surface_selector_nn.set_dist_exponent(1)  # Default 1
    surface_selector_nn.set_cache_selection(1)
    surface_helix_selector = residue_selector.AndResidueSelector(
        surface_selector_nn, helix_selector)
    return surface_helix_selector


def select_bound_helix_nn(helix_selector):
    bound_selector_nn = residue_selector.LayerSelector()
    bound_selector_nn.set_use_sc_neighbors(1)
    bound_selector_nn.set_layers(0, 1, 0)
    bound_selector_nn.set_cutoffs(5.2, 2)  # Default 5.2 and 2.0
    # 0.3 selects boundary default 5, larger wider cone
    bound_selector_nn.set_angle_shift_factor(0.5)
    # 2 default smaller sharper and broader cone
    bound_selector_nn.set_angle_exponent(2)
    bound_selector_nn.set_sc_neighbor_dist_midpoint(9)  # Default 9
    bound_selector_nn.set_dist_exponent(1)  # Default 1
    bound_selector_nn.set_cache_selection(0)
    bound_helix_selector = residue_selector.AndResidueSelector(
        bound_selector_nn, helix_selector)
    return bound_helix_selector


def select_surf_helix_sasa(helix_selector):
    surface_selector_sasa = residue_selector.LayerSelector()
    surface_selector_sasa.set_use_sc_neighbors(0)
    surface_selector_sasa.set_cache_selection(1)
    surface_selector_sasa.set_ball_radius(1.8)
    surface_selector_sasa.set_layers(0, 0, 1)
    surface_selector_sasa.set_cutoffs(20, 45)
    surface_helix_selector = residue_selector.AndResidueSelector(
        surface_selector_sasa, helix_selector)
    return surface_helix_selector


def get_rep_pos(rep_len, resi_selected):
    # Get repeat positions for a selection on the second inner repeat
    rep_pos_operations = [-rep_len, 0, rep_len, 2 * rep_len]
    resi_interface = []
    for i in resi_selected:
        for j in rep_pos_operations:
            resi_interface.append(int(i + j))
    return(resi_interface)


def default_tf():
    tf = TaskFactory()
    tf.push_back(operation.InitializeFromCommandline())  # IFC
    tf.push_back(operation.IncludeCurrent())  # IC
    aro_trp = task_ops.LimitAromaChi2Operation()
    aro_trp.include_trp(1)
    aro_trp.chi2min(45)
    aro_trp.chi2max(110)
    tf.push_back(aro_trp)  # AroTrp
    aro_chi = task_ops.LimitAromaChi2Operation()
    aro_chi.chi2min(70)
    aro_chi.chi2max(110)
    tf.push_back(aro_chi)  # AroChi
    return tf


def nn_surface_selector():
    surf_nn_sele = residue_selector.LayerSelector()
    surf_nn_sele.set_use_sc_neighbors(1)
    surf_nn_sele.set_layers(0, 0, 1)
    surf_nn_sele.set_cutoffs(5.2, 2)  # Default 5.2 and 2.0
    # 0.3 selects boundary default 5, larger wider cone
    surf_nn_sele.set_angle_shift_factor(0.8)
    # 2 default smaller sharper and broader cone
    surf_nn_sele.set_angle_exponent(2)
    surf_nn_sele.set_sc_neighbor_dist_midpoint(9)  # Default 9
    surf_nn_sele.set_dist_exponent(1)  # Default 1
    surf_nn_sele.set_cache_selection(0)
    return surf_nn_sele


def nn_boundary_selector():
    bound_nn_sele = residue_selector.LayerSelector()
    bound_nn_sele.set_use_sc_neighbors(1)
    bound_nn_sele.set_layers(0, 1, 0)
    bound_nn_sele.set_cutoffs(5.2, 2)  # Default 5.2 and 2.0
    # 0.3 selects boundary default 5, larger wider cone
    bound_nn_sele.set_angle_shift_factor(0.8)
    # 2 default smaller sharper and broader cone
    bound_nn_sele.set_angle_exponent(2)
    bound_nn_sele.set_sc_neighbor_dist_midpoint(9)  # Default 9
    bound_nn_sele.set_dist_exponent(1)  # Default 1
    bound_nn_sele.set_cache_selection(0)
    return bound_nn_sele


def main(pdb_name=None, input_side=None, allowed_aas_surface=None,
         allowed_aas_surf_tail=None, allowed_aas_bound_tail=None, nloop_sym_pack_mov=None, nloop_tail_pack_mov=None,
         nstruct_no=None):
    if input_side == 'A':
        cst_file_name = glob.glob('*_lattice_csts_A.cst')[0]
    elif input_side == 'B':
        cst_file_name = glob.glob('*_lattice_csts_B.cst')[0]
    print('Read in cst_file_name = ', cst_file_name)
    # Basic information 4-repeat-DHR
    output_str_temp = pdb_name[:-4] + '_side_' + input_side + \
        '_alloaas_' + allowed_aas_surface + \
        '_aasurftail_' + allowed_aas_surf_tail + \
        '_{:0>3d}.pdb'
    output_name = output_str_temp.format(nstruct_no)
    input_pose = pose_from_pdb(pdb_name)
    side_dict = {'A': 3, 'B': 4}
    opposite_side_dict = {'A': 4, 'B': 3}
    side = side_dict[input_side]  # 3rd Helix or side A
    opposite_side = opposite_side_dict[input_side]
    ss_dict = get_sec_struct_dict(input_pose)
    pos_len = input_pose.total_residue()
    tot_res_vector = np.arange(1, pos_len + 1)
    # Specific to DHR (2 H per repeat)
    rep_len = int(pos_len / (len(ss_dict) / 2))
    start_inner_rep = int(1 + rep_len)
    inner_repeat_sequence = str(
        input_pose.sequence())[
        start_inner_rep -
        1:int(
            start_inner_rep +
            rep_len -
            1)]
    # ss_dict Divided by two because each subunit has two elements
    tails = range(1, rep_len + 1), range(rep_len *
                                         int((len(ss_dict) / 2) - 1) + 1, pos_len + 1)
    ###

    print('Pose length : ', pos_len)
    print(
        'Dictionary of helical elements and their residue indices: ',
        ss_dict)
    print('Identified repeat length: ', rep_len)
    print('Inner repeat sequence :', inner_repeat_sequence)
    print('Tail ranges: ', tails)

    # Restore symmetry of the pose (thread inner repeat sequence in tails) <-
    # Make capless object
    for tail in tails:
        for counter, resi in enumerate(tail):
            mutants.mutate_residue(
                input_pose,
                resi,
                inner_repeat_sequence[counter],
                0.0,
                None)
            print('Tail residue index:', resi)
    print('Finished making capless pose.')
    ###

    # Mutate nn_surface residues to Ala
    helix_selector = select_correct_helix(ss_dict, side, rep_len)
    nn_selector_1hel = select_surf_helix_nn(helix_selector)
    resi_nn_selector_1hel = tot_res_vector[
        nn_selector_1hel.apply(input_pose)]
    resi_nn_selector_tot = get_rep_pos(rep_len, resi_nn_selector_1hel)
    print('Finished selecting residues.')
    for resi in resi_nn_selector_tot:
        mutants.mutate_residue(input_pose, resi, 'A', 0.0, None)
        print('Mutating residue {} to Alanine'.format(resi))
    ###

    # Create a combined selection with SASA and NN
    sasa_selector_1hel = select_surf_helix_sasa(helix_selector)
    comb_selector_1hel = residue_selector.OrResidueSelector(
        sasa_selector_1hel, nn_selector_1hel)
    resi_comb_selector_1hel = tot_res_vector[
        comb_selector_1hel.apply(input_pose)]
    resi_comb_selector_tot = get_rep_pos(
        rep_len, resi_comb_selector_1hel)
    print('Created combined SASA and NN selection for input side.')
    ###

    # Create a combined selection with SASA and NN for opposite side
    opposite_helix_selector = select_correct_helix(
        ss_dict, opposite_side, rep_len)
    op_nn_sele_1hel = select_surf_helix_nn(opposite_helix_selector)
    op_sasa_sele_1hel = select_surf_helix_sasa(opposite_helix_selector)
    op_comb_sele_1hel = residue_selector.OrResidueSelector(
        op_nn_sele_1hel, op_sasa_sele_1hel)
    op_resi_comb_sele_1hel = tot_res_vector[op_comb_sele_1hel.apply(
        input_pose)]
    op_resi_comb_sele_tot = get_rep_pos(rep_len, op_resi_comb_sele_1hel)
    print('Created combined SASA and NN selection for opposite side.')
    ###

    # Create XML object
    xml_string = '''
        <SCOREFXNS>
                <ScoreFunction name="beta_design_sym" weights="beta_nov16_cart.wts" symmetric="1">
                        <Reweight scoretype="hbond_sr_bb" weight="1.05"/>
                        <Reweight scoretype="hbond_lr_bb" weight="1.05"/>
                        <Reweight scoretype="hbond_bb_sc" weight="2.0" />
                        <Reweight scoretype="hbond_sc" weight="2"/>
                        <Reweight scoretype="fa_elec" weight="1.5"/>
                        <Set approximate_buried_unsat_penalty_burial_probe_radius="3.5" />
                        <Set approximate_buried_unsat_penalty_burial_atomic_depth="4.0" />
                        <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-0.2" />
                </ScoreFunction>
                <ScoreFunction name="beta_design_tails" weights="beta_nov16_cart.wts" symmetric="0">
                        <Reweight scoretype="hbond_sr_bb" weight="1.05"/>
                        <Reweight scoretype="hbond_lr_bb" weight="1.05"/>
                        <Reweight scoretype="hbond_bb_sc" weight="2.0" />
                        <Reweight scoretype="hbond_sc" weight="2"/>
                        <Reweight scoretype="fa_elec" weight="1.5"/>
                        <Set approximate_buried_unsat_penalty_burial_probe_radius="3.5" />
                        <Set approximate_buried_unsat_penalty_burial_atomic_depth="4.0" />
                        <Set approximate_buried_unsat_penalty_hbond_energy_threshold="-0.2" />
                </ScoreFunction>
                <ScoreFunction name="beta_relax" weights="beta_nov16_cart.wts" symmetric="0">
                        <Reweight scoretype="cart_bonded" weight="1.5"/>
                </ScoreFunction>
                <ScoreFunction name="beta_score" weights="beta_nov16_cart.wts" symmetric="0">
                        <Reweight scoretype="atom_pair_constraint" weight="1"/>
                </ScoreFunction>
        </SCOREFXNS>
        <RESIDUE_SELECTORS>
                <SecondaryStructure name="loop_select_no_termini"
                        overlap="1" minH="1" minE="1" include_terminal_loops="false"
                        use_dssp="true" ss="L" />
        </RESIDUE_SELECTORS>
        <TASKOPERATIONS>
                <InitializeFromCommandline name="IFC"/>
                <IncludeCurrent name="IC"/>
                <LimitAromaChi2 name="aroTrp" chi2min="45" chi2max="110" include_trp="true" />
                <LimitAromaChi2 name="aroChi" chi2min="70" chi2max="110" />
        </TASKOPERATIONS>
        <MOVERS>
                <RepeatProteinRelax name="setup_sym" numb_repeats="4" scorefxn="beta_design_sym"
                        loop_cutpoint_mode="false" minimize="false" relax_iterations="1" cartesian="true"
                        modify_symmetry_and_exit="true" remove_symmetry="false"/>
                <RepeatProteinRelax name="remove_sym" numb_repeats="4" scorefxn="beta_design_sym"
                        loop_cutpoint_mode="false" minimize="false" relax_iterations="1" cartesian="true"
                        modify_symmetry_and_exit="true" remove_symmetry="true"/>
                <FastRelax name="frelax"
                        scorefxn="beta_relax" min_type="lbfgs_armijo_nonmonotone" relaxscript="MonomerRelax2019"
                        repeats="1" ramp_down_constraints="0"
                        cartesian="1" bondangle="1" bondlength="1"
                        task_operations="IFC,IC,aroTrp,aroChi"/>
        </MOVERS>'''
    xmlobj = xml_object.create_from_string(xml_string)
    ###

    # Setup pose symmetry
    print('Setting up symmetry!')
    setup_sym_mover = xmlobj.get_mover('setup_sym')
    setup_sym_mover.apply(input_pose)
    print('DONE setting up symmetry.')
    ###

    # Design aas from one side of the protein
    print('Designing aas of side {}'.format(input_side))
    sfxn_design = xmlobj.get_score_function('beta_design_sym')

    not_surface_selector = residue_selector.NotResidueSelector(
        comb_selector_1hel)

    tf = default_tf()  # Create Task Factory

    repk_not_surface = operation.OperateOnResidueSubset()
    repk_not_surface.append_op(operation.RestrictToRepackingRLT())
    repk_not_surface.selector(not_surface_selector)
    tf.push_back(repk_not_surface)  # Restrict other to repacking
    des_aas_rlt = operation.RestrictAbsentCanonicalAASRLT()
    des_aas_rlt.aas_to_keep(allowed_aas_surface)
    des_aas_op = operation.OperateOnResidueSubset()
    des_aas_op.op(des_aas_rlt)
    des_aas_op.selector(comb_selector_1hel)
    tf.push_back(des_aas_op)  # Design constraints for surface

    sym_pack_mov = pack_min.symmetry.SymPackRotamersMover()
    sym_pack_mov.score_function(sfxn_design)
    sym_pack_mov.nloop(nloop_sym_pack_mov)
    sym_pack_mov.task_factory(tf)
    sym_pack_mov.apply(input_pose)
    print('DONE with interface design.')
    ###

    # Preparation removes symmetry
    print('Removing symmetry.')
    remove_sym_mover = xmlobj.get_mover('remove_sym')
    remove_sym_mover.apply(input_pose)
    no_tails_pose = input_pose.clone()

    # Fast relax
    fast_relax = xmlobj.get_mover('frelax')
    fast_relax.apply(no_tails_pose)

    # Dumping and scoring capless pose
    sfxn_score = xmlobj.get_score_function(
        'beta_score')  # with weighted atom_pair_csts
    apply_lattice_csts = constraint_movers.ConstraintSetMover()
    apply_lattice_csts.add_constraints(1)
    apply_lattice_csts.constraint_file(cst_file_name)
    apply_lattice_csts.apply(no_tails_pose)
    no_tails_pose.dump_scored_pdb(
        (output_name[:-4] + '_uncapped.pdb'), sfxn_score)

    # Re-design the tails
    sfxn_upehb_tails = xmlobj.get_score_function('beta_design_tails')
    tf_2 = default_tf()  # Create Task Factory

    # Restrict to repacking the inner loops
    loop_selector = xmlobj.get_residue_selector(
        'loop_select_no_termini')
    loop_resi_list = list(tot_res_vector[loop_selector.apply(
        input_pose)])
    keep_residues = list(
        set(loop_resi_list + resi_comb_selector_tot + op_resi_comb_sele_tot))
    print('Keep residues: ', keep_residues)
    tails_not_interface = [
        resi for tail in tails for resi in tail if resi not in keep_residues]
    str_tail_not_intf = ','.join(map(str, tails_not_interface))
    tail_not_interface_selector = residue_selector.ResidueIndexSelector(
        str_tail_not_intf)

    surf_nn_sele = nn_surface_selector()
    surf_nn_sele_tail_nitf = residue_selector.AndResidueSelector(
        surf_nn_sele, tail_not_interface_selector)
    des_aas_rlt_tail_surf = operation.RestrictAbsentCanonicalAASRLT()
    des_aas_rlt_tail_surf.aas_to_keep(allowed_aas_surf_tail)
    des_aas_op_tail_surf = operation.OperateOnResidueSubset()
    des_aas_op_tail_surf.op(des_aas_rlt_tail_surf)
    des_aas_op_tail_surf.selector(surf_nn_sele_tail_nitf)
    # Design constraints for surface tails
    tf_2.push_back(des_aas_op_tail_surf)

    bound_nn_sele = nn_boundary_selector()
    bound_nn_sele_tail_nitf = residue_selector.AndResidueSelector(
        bound_nn_sele, tail_not_interface_selector)
    des_aas_rlt_tail_bound = operation.RestrictAbsentCanonicalAASRLT()
    des_aas_rlt_tail_bound.aas_to_keep(allowed_aas_bound_tail)
    des_aas_op_tail_bound = operation.OperateOnResidueSubset()
    des_aas_op_tail_bound.op(des_aas_rlt_tail_bound)
    des_aas_op_tail_bound.selector(bound_nn_sele_tail_nitf)
    # Design constraints for boundary tails
    tf_2.push_back(des_aas_op_tail_bound)

    des_tails_sele = residue_selector.OrResidueSelector(
        surf_nn_sele_tail_nitf, bound_nn_sele_tail_nitf)
    ndes_tails_sele = residue_selector.NotResidueSelector(
        des_tails_sele)
    freeze_ndes_tails_op = operation.OperateOnResidueSubset()
    freeze_ndes_tails_op.append_op(operation.PreventRepackingRLT())
    # Select non-tails for freezing
    bound_h_nn_sele = select_bound_helix_nn(helix_selector)
    resi_bound_h_nn_sele_1hel = tot_res_vector[
        bound_h_nn_sele.apply(input_pose)]
    resi_bound_h_nn_sele_tot = get_rep_pos(
        rep_len, resi_bound_h_nn_sele_1hel)
    str_resi_bound_h_nn_sele_tot = ','.join(
        map(str, resi_bound_h_nn_sele_tot))
    bound_h_nn_sele_tot = residue_selector.ResidueIndexSelector(
        str_resi_bound_h_nn_sele_tot)
    freeze_doubleb_sele = residue_selector.AndResidueSelector(
        bound_nn_sele_tail_nitf, bound_h_nn_sele_tot)
    # Freeze residues at boundary between repeats
    tot_freeze_sele = residue_selector.OrResidueSelector(
        ndes_tails_sele, freeze_doubleb_sele)
    freeze_ndes_tails_op.selector(tot_freeze_sele)
    tf_2.push_back(freeze_ndes_tails_op)  # Freeze non-designable

    print("I'll start packing the tails!")

    tail_pack_mov = pack_min.PackRotamersMover()
    tail_pack_mov.nloop(nloop_tail_pack_mov)
    tail_pack_mov.score_function(sfxn_upehb_tails)
    tail_pack_mov.task_factory(tf_2)
    tail_pack_mov.apply(input_pose)
    ###

    fast_relax.apply(input_pose)
    apply_lattice_csts.apply(input_pose)
    input_pose.dump_scored_pdb(output_name, sfxn_score)

    print("I'm done with the {} iteration.".format(str(nstruct_no)))

if __name__ == "__main__":
    # Declare parser object for managing input options
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_pdb", type=str,
                        help="Name of the input pdb.")
    parser.add_argument("-s", "--side", type=str,
                        help="One letter side identifier for the protein.", choices=('A', 'B'), default='A')
    parser.add_argument('-aas', '--design_aas', type=str,
                        help='One letter aa string of allowed to design the interface.', default='EDA')
    parser.add_argument('-aast', '--design_aas_surf_tail', type=str,
                        help='One letter aa string of allowed identities on the surface of tails.',
                        default='DNSTQKREDH')
    parser.add_argument('-aabt', '--design_aas_bound_tail', type=str,
                        help='One letter aa string of allowed identities on the surface of tails.',
                        default='ALDNSTQKREDH')
    parser.add_argument('-nlpsym', '--nloop_sym', type=int,
                        help='Number of loops for the symmetric packing step.', default=5)  # default back to 5
    parser.add_argument('-nlptail', '--nloop_tail', type=int,
                        help='Number of loops for the tail packing step.', default=5)  # default back to 5
    parser.add_argument('-nstruct', '--number_of_runs', type=int,
                        help='Number of full runs.', default=1)
    args = parser.parse_args()

    nstruct = 0
    for n in range(args.number_of_runs):
        print('This is the {} iteration!'.format(n))
        # Run protocol
        main(pdb_name=args.input_pdb, input_side=args.side,
             allowed_aas_surface=args.design_aas,
             allowed_aas_surf_tail=args.design_aas_surf_tail,
             allowed_aas_bound_tail=args.design_aas_bound_tail,
             nloop_sym_pack_mov=args.nloop_sym,
             nloop_tail_pack_mov=args.nloop_tail, nstruct_no=n)
