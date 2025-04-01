import os
import logging
import re
from optimhc.feature_generator.ladder_peptide import (
    LadderPeptideFeatureGenerator,
    assign_brother_aggregated_feature
)
from optimhc.feature_generator.basic import BasicFeatureGenerator
from optimhc.feature_generator.PWM import PWMFeatureGenerator
from optimhc.feature_generator.mhcflurry import MHCflurryFeatureGenerator
from optimhc.feature_generator.netMHCpan import NetMHCpanFeatureGenerator
from optimhc.feature_generator.netMHCIIpan import NetMHCIIpanFeatureGenerator
from optimhc.feature_generator.DeepLC import DeepLCFeatureGenerator
from optimhc.feature_generator.spectra_similarity import SpectraSimilarityFeatureGenerator

logger = logging.getLogger(__name__)


# TODO: Add dynamic import of feature generators
def import_generators(config):
    """
    Dynamically import feature generator modules based on the provided configuration.
    This function is left as a placeholder for optional dynamic imports.

    Parameters:
        config (dict): Configuration dictionary loaded from YAML.
    """
    feature_generators = config['feature_generators']
    _ = [fg['type'] for fg in feature_generators if isinstance(fg, dict) and 'type' in fg]
    # Expand logic here if needed for dynamic imports in the future.


# TODO: Add a function to save original outputs of feature generators
def generate_features(psms, config):
    """
    Generate features from different generators according to the configuration.

    Parameters:
        psms (PsmContainer): A container object holding PSMs and relevant data.
        config (dict): Configuration dictionary loaded from YAML.
    """
    global_params = config['global_parameters']
    remove_modification = global_params['remove_modification']
    remove_pre_nxt_aa = global_params['remove_pre_nxt_aa']
    n_processes = global_params['n_processes']
    show_progress = global_params['show_progress']
    output_dir = config['output_dir']
    oxidation_tag = global_params['oxidation_tag']
    
    feature_generators = config['feature_generators']
    allele = config['allele']
    score = config['score']
    unique_peptides = list(set(psms.peptides))

    for generator_config in feature_generators:
        if not isinstance(generator_config, dict):
            logger.warning("Feature generator config is not a dictionary, skipping...")
            continue

        generator_type = generator_config.get('type')
        logger.info(f"Generating features with {generator_type}...")

        if generator_type == 'LadderPeptide':
            ladder_peptide = LadderPeptideFeatureGenerator(
                unique_peptides,
                min_overlap_length=generator_config.get('min_overlap_length', 8),
                min_entropy=generator_config.get('min_entropy', 0),
                min_length=generator_config.get('min_length', 8),
                max_length=generator_config.get('max_length', 25),
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                remove_modification=remove_modification
            )
            ladder_features = ladder_peptide.generate_features()
            full_data = ladder_peptide.get_full_data()

            psms.add_metadata(
                full_data[['Peptide', 'brother_count', 'ContigSequence']],
                psms_key=psms.peptide_column,
                metadata_key='Peptide',
                source='LadderPeptide'
            )
            psms.add_features(
                ladder_features,
                psms_key=psms.peptide_column,
                feature_key=ladder_peptide.id_column,
                source='LadderPeptide'
            )

            if score:
                assign_brother_aggregated_feature(psms, feature_columns=score, ladder_source='LadderPeptide')

        elif generator_type == 'Basic':
            basic_generator = BasicFeatureGenerator(
                psms.psms[psms.peptide_column].tolist(),
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                remove_modification=remove_modification
            )
            basic_features = basic_generator.generate_features()
            psms.add_features_by_index(basic_features[basic_generator.feature_columns], source='Basic')

        elif generator_type == 'PWM':
            pwm_generator = PWMFeatureGenerator(
                unique_peptides,
                alleles=allele,
                mhc_class=generator_config.get('mhc_class', 'I'),
                remove_modification=remove_modification,
                remove_pre_nxt_aa=remove_pre_nxt_aa
            )
            pwm_features = pwm_generator.generate_features()
            psms.add_features(
                pwm_features,
                psms_key=psms.peptide_column,
                feature_key=pwm_generator.id_column,
                source='PWM'
            )

        elif generator_type == 'MHCflurry':
            mhcflurry_generator = MHCflurryFeatureGenerator(
                unique_peptides,
                alleles=allele,
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                remove_modification=remove_modification
            )
            mhcflurry_features = mhcflurry_generator.generate_features()
            psms.add_features(
                mhcflurry_features,
                psms_key=psms.peptide_column,
                feature_key=mhcflurry_generator.id_column,
                source='MHCflurry'
            )

        elif generator_type == 'NetMHCpan':
            netmhcpan_generator = NetMHCpanFeatureGenerator(
                unique_peptides,
                alleles=allele,
                mode=generator_config.get('mode', 'best'),
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                remove_modification=remove_modification,
                n_processes=generator_config.get('n_processes', n_processes),
                show_progress=generator_config.get('show_progress', show_progress)
            )
            netmhcpan_features = netmhcpan_generator.generate_features()
            psms.add_features(
                netmhcpan_features,
                psms_key=psms.peptide_column,
                feature_key=netmhcpan_generator.id_column,
                source='NetMHCpan'
            )

        elif generator_type == 'NetMHCIIpan':
            netmhciipan_generator = NetMHCIIpanFeatureGenerator(
                unique_peptides,
                alleles=allele,
                mode=generator_config.get('mode', 'best'),
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                remove_modification=remove_modification,
                n_processes=generator_config.get('n_processes', n_processes),
                show_progress=generator_config.get('show_progress', show_progress)
            )
            netmhciipan_features = netmhciipan_generator.generate_features()
            psms.add_features(
                netmhciipan_features,
                psms_key=psms.peptide_column,
                feature_key=netmhciipan_generator.id_column,
                source='NetMHCIIpan'
            )

        elif generator_type == 'DeepLC':
            if oxidation_tag is not None:
                mod_dict = {
                    oxidation_tag: 'Oxidation'
                }
            deeplc_generator = DeepLCFeatureGenerator(
                psms,
                calibration_criteria_column=generator_config.get('calibration_criteria_column', 'spscore'),
                lower_score_is_better=generator_config.get('lower_score_is_better', False),
                calibration_set_size=generator_config.get('calibration_set_size', 0.1),
                processes=generator_config.get('processes', n_processes),
                model_path=generator_config.get('model_path', None),
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                mod_dict=mod_dict
            )
            deeplc_features = deeplc_generator.generate_features()
            psms.add_features_by_index(deeplc_features[deeplc_generator.feature_columns], source='DeepLC')

        elif generator_type == 'SpectraSimilarity':

            # Match PSMs with the spectra
            mzML_dir = generator_config.get('mzML_dir', None)
            if mzML_dir is None:
                logger.error("mzML_dir is not provided for SpectraSimilarity feature generator.")
                continue 

            pattern = generator_config.get('spectrum_id_pattern', None)
            mz_file_names = []
            spectrum_ids = psms.spectrum_ids

            if pattern:
                logger.info(f"Using pattern: {pattern} to extract mzML file names from spectrum IDs.")
                for spectrum_id in spectrum_ids:
                    mz_file_names.append(re.match(pattern, spectrum_id).group(1))
                logger.info(f"mzML file names: {list(set(mz_file_names))}")
            else:
                logger.info("Spectrum ID pattern is not provided.")
                if psms.ms_data_file_column is not None:
                    logger.info(f"Trying to extract mzML file names from {psms.ms_data_file_column}")
                    logger.info(f"MS data file format: {set(psms.psms[psms.ms_data_file_column])}")

                    for ms_data_file in psms.psms[psms.ms_data_file_column]:
                        mz_file_basename = os.path.basename(ms_data_file).split('.')[0]
                        if mz_file_basename.endswith('.mzML'):
                            mz_file_basename = mz_file_basename[:-5]
                        elif mz_file_basename.endswith('mzML'):
                            mz_file_basename = mz_file_basename[:-4]
                        mz_file_names.append(mz_file_basename)

                    logger.info(f"mzML file names: {list(set(mz_file_names))}")
                else:
                    logger.info("MS data file information is not provided.")
                    logger.info("Trying to use the default pattern: (.+?)\.\d+\.\d+\.\d+ to extract mzML file names from spectrum IDs.")
                    for spectrum_id in spectrum_ids:
                        mz_file_names.append(re.match('(.+?)\.\d+\.\d+\.\d+', spectrum_id).group(1))

            mz_file_paths = [os.path.join(mzML_dir, f'{mz_file}.mzML') for mz_file in mz_file_names]
            mz_file_paths_set = set(mz_file_paths)
            logger.info(f'mz_file_paths: {mz_file_paths_set}')
            
            for mz_file_path in mz_file_paths_set:
                if not os.path.exists(mz_file_path):
                    logger.error(f"mzML file not found: {mz_file_path}")
                    continue

            spectra_similarity_generator = SpectraSimilarityFeatureGenerator(
                spectrum_ids=psms.spectrum_ids,
                peptides=psms.peptides,
                charges=psms.charges,  
                scan_ids=psms.scan_ids,
                mz_file_paths=mz_file_paths,
                model_type=generator_config.get('model', 'CID'),
                collision_energies=[generator_config.get('collision_energy', None)] * len(psms.scan_ids),
                remove_modification=remove_modification,
                remove_pre_nxt_aa=remove_pre_nxt_aa,
                url=generator_config.get('url', 'koina.wilhelmlab.org:443'),
                top_n=generator_config.get('top_n', 36),
                tolerance_ppm=generator_config.get('tolerance_ppm', 20),
            )

            spectra_similarity_features = spectra_similarity_generator.generate_features()
            psms.add_features(
                spectra_similarity_features,
                psms_key=[psms.spectrum_column, psms.peptide_column, psms.charge_column],
                feature_key=spectra_similarity_generator.id_column,
                source='SpectraSimilarity'
            )

        else:
            logger.warning(f"Unknown feature generator type: {generator_type}")
