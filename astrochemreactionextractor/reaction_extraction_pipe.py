import os
import re
import sys
import time
import math
import json
from tqdm import tqdm
import pandas as pd
import matplotlib.pyplot as plt

from astrochemreactionextractor import logger
from astrochemreactionextractor.utils.utils_reaction import find_formulas_contexts, get_smiles_from_multipub, is_valid_formula, get_smiles_from_simple_species, extract_rate_from_page, completion_reaction

from astrochemreactionextractor.llm import deepseek_api_reactions_pdf, deepseek_api_formula, deepseek_api_species

sys.path.append('..')
sys.path.append('')
plt.style.use('default')

class ReactionPipeline:
    def __init__(self, input_path, output_path):
        self.input_path = input_path
        self.query_path = 'cache'
        self.output_path = output_path
        self.dataframe_simple_species = pd.DataFrame({
            'formula': [],
            'SMILES': []
        })
    def reaction_extract_pipe(self, reactions_paper, all_reactions_paper, file_content, reaction_content, shared_data, _timeout=None):
        if 'doi' in reactions_paper:
            return [], reactions_paper
        else:
            if os.path.exists(f'{self.query_path}/llm_df_reactions_{reactions_paper["title"]}.csv'):
                logger.info(f'{reactions_paper["title"]} has been processed, load the memory result')
                df_results1 = pd.read_csv(
                    f'{self.query_path}/llm_df_reactions_{reactions_paper["title"]}.csv')
                prompt_tokens_reactions = 0
                completion_tokens_reactions = 0
            else:
                logger.info(f'Start processing paper {reactions_paper["title"]}')
                df_results1, prompt_tokens_reactions, completion_tokens_reactions = deepseek_api_reactions_pdf(
                    reaction_content)
                df_results1.to_csv(
                    f'{self.query_path}/llm_df_reactions_{reactions_paper["title"]}.csv', index=False)
                logger.info(f'{reactions_paper["title"]} processing completed!!!')
            df_results1 = completion_reaction(df_results1)

        shared_data['prompt_tokens_reactions'] += prompt_tokens_reactions
        shared_data['completion_tokens_reactions'] += completion_tokens_reactions
        reaction_paper_list = []
        if df_results1.shape[0] > 0:
            logger.info(f'{reactions_paper["title"]} contain recations')
            shared_data['num_paper_contain_recation'] += 1
            num_reactions_temp = 0
            num_reactions_invalid_temp = 0
            num_retrieval_contexts = 0
            results_all = {}
            value_smiles_pubchem = {}
            num = 0

            for index, row in df_results1.iterrows():
                recation_normal = row['reactions_wx']
                formula_api, prompt_tokens_species, completion_tokens_species = deepseek_api_formula(
                    recation_normal)
                shared_data['prompt_tokens_species'] += prompt_tokens_species
                shared_data['completion_tokens_species'] += completion_tokens_species
                formula_list = [formula.strip().replace('\\text', '').replace(
                    '\\ce', '') for formula in formula_api.split('\n') if formula.strip() != '```']
                all_reaction_paper = {
                    "title": reactions_paper["title"],
                    "reaction": recation_normal,
                    "formulas": []
                }
                all_reaction_paper["formulas"].append(formula_list)
                all_reactions_paper["reactions"].append(all_reaction_paper)
                shared_data['num_reactions_all'] += 1
                if is_valid_formula(formula_list):
                    shared_data['num_reactions_valid'] += 1
                    num_reactions_temp += 1

                    reaction_paper = {
                        "title": reactions_paper["title"],
                        'coordinates': row['coordinates'],
                        'params': row['params'],
                        "reaction": '',
                        "formulas": []
                    }

                    cleaned_formula_list = [text.replace('·', '').replace('^', '').replace('_', '').replace('*', '').replace('(', '').replace(')', '') for text in formula_list]
                    zip_formula_pairs = [
                        (formula_ori, formula_cleaned)
                        for formula_ori, formula_cleaned in zip(formula_list, cleaned_formula_list)
                        if formula_cleaned not in self.dataframe_simple_species['formula'].values
                    ]
                    if zip_formula_pairs == []:
                        retrieval_contexts = ''
                        updated_formula_list = []
                        new_dict_exist = {}
                    else:
                        elements_formula_list, elements_formula_list_cleaned = zip(
                            *zip_formula_pairs)
                        new_dict_exist = {
                            key: results_all[key] for key in elements_formula_list if key in results_all}
                        updated_formula_list = [
                            key for key in elements_formula_list if key not in results_all]
                        retrieval_contexts = find_formulas_contexts(
                            file_content, updated_formula_list)

                    if retrieval_contexts == '':
                        num_retrieval_contexts += 1
                        results = {
                            molecule: 'None' for molecule in cleaned_formula_list}
                    else:
                        results, prompt_tokens_species1, completion_tokens_species1 = deepseek_api_species(
                            recation_normal, retrieval_contexts)

                    extracted_results = {
                        key: results[key] for key in updated_formula_list if key in results}
                    results_all.update(extracted_results)
                    num_value = 0
                    num_pub = 0
                    num_smiles_temp = 0
                    updated_results = {key: new_dict_exist.get(
                        key, results[key]) for key in results}
                    for key, value in updated_results.items():
                        smiles = 'None'
                        if not is_valid_formula([key]):
                            break
                        if get_smiles_from_simple_species(key) != 'None':
                            smiles = get_smiles_from_simple_species(key)
                            if smiles != 'None':
                                num_smiles_temp += 1
                        else:
                            if value != 'None':
                                num_value += 1
                                if value in value_smiles_pubchem:
                                    smiles = value_smiles_pubchem[value]
                                else:
                                    if '(' in value:
                                        value_list = re.split(
                                            r'\s*\(([^)]*)\)', value)
                                    else:
                                        value_list = value.split(' or ')
                                    if len(value_list) == 1:
                                        smiles = get_smiles_from_multipub(
                                            key, value)
                                    else:
                                        smiles_0 = get_smiles_from_multipub(key,
                                                                            value_list[0])
                                        smiles_1 = get_smiles_from_multipub(key,
                                                                            value_list[1].replace('or ', ''))
                                        if 'None' not in smiles_0:
                                            smiles = smiles_0
                                        else:
                                            smiles = smiles_1 if 'None' not in smiles_1 else 'None'
                                    if 'None' not in smiles:
                                        num_pub += 1
                                        num_smiles_temp += 1
                                    value_smiles_pubchem.update({value: smiles})

                        result_pair = f"{key}: {value}: {smiles}"

                        reaction_paper["formulas"].append(
                            result_pair)
                    if num_value != 0 and num_pub == 0:
                        num_reactions_invalid_temp += 1

                    reaction_paper["reaction"] = re.sub(
                        r' \(.*?\)', '', recation_normal)
                    reaction_paper_list.append(reaction_paper["reaction"])
                    reactions_paper["reactions"].append(reaction_paper)
                else:
                    num += 1

            if not os.path.exists(f'{self.query_path}/llm_result_species_{reactions_paper["title"]}.json'):
                with open(f'{self.query_path}/llm_result_species_{reactions_paper["title"]}.json', "w") as f:
                    json.dump(results_all, f, indent=4)

        return reaction_paper_list, reactions_paper


    def extract(self):
        shared_data = dict({'num_reactions_no_formulas_contexts': 0, 'num_papers_tar': 0, 'num_papers_tar_valid': 0,
                            'num_paper_contain_keywords': 0, 'num_paper_contain_recation': 0,
                            'num_paper_contain_recation_pass_validtest': 0, 'num_reactions_valid': 0,
                            'num_reactions_all': 0, 'num_reactions_no_chemical_name': 0,
                            'num_reactions_haschemical_name_nopub': 0, 'num_reactions_not_has_smiles': 0, 'num_formulas': 0,
                            'num_formulas_pubchem': 0, 'num_formulas_simple_species': 0, 'prompt_tokens_reactions': 0,
                            'completion_tokens_reactions': 0, 'prompt_tokens_species': 0, 'completion_tokens_species': 0,
                            'results_papers': [],
                            'results_papers_list': [],
                            'processed_papers_hasreaction': [],
                            'processed_papers_all': []})

        start_time = time.time()
        date_str = time.strftime("%Y-%m-%d")
        millis = int(start_time * 1000)
        input_file_list = []
        input_file_list = [os.path.join(
            self.input_path, filename) for filename in os.listdir(self.input_path)]

        count_pdf = 1
        for filepath in tqdm(input_file_list):
            if filepath.endswith('.pdf'):
                logger.info(f'Processing paper {count_pdf}: {os.path.basename(filepath)}!!!')
                count_pdf += 1
                shared_data['num_papers_tar'] += 1

                reactions_paper = {
                    "title": '',
                    "reactions": []
                }
                all_reactions_paper = {
                    "title": '',
                    "reactions": []
                }

                title = os.path.basename(filepath)
                if title != None:
                    reactions_paper['title'] = title[:200]
                    all_reactions_paper['title'] = title[:200]
                    if title in shared_data['processed_papers_all']:
                        print("the paper has been processed!")
                        continue

                try:
                    df_pdf, pdf_contents = extract_rate_from_page(filepath)
                except:
                    logger.error(f"PDF parsing failed!")
                    continue
                if pdf_contents == '':
                    continue 
                try:                
                    reactions_paper_list, reactions_paper = self.reaction_extract_pipe(
                        reactions_paper=reactions_paper, all_reactions_paper=all_reactions_paper, file_content=pdf_contents, reaction_content=df_pdf, shared_data=shared_data)
                except TimeoutError as e:
                    logger.error(f"Error: {e}，processing timeout")
                    continue
            else:
                logger.error("The file format is not supported, please upload a PDF.")
                continue

            if reactions_paper_list == []:
                continue
            if reactions_paper["reactions"] != []:
                shared_data['results_papers'].extend(reactions_paper['reactions'])
                shared_data['processed_papers_hasreaction'].append(title)
            shared_data['processed_papers_all'].append(title)
            shared_data['results_papers_list'].extend(reactions_paper_list)

            result_file = os.path.join(
                self.output_path, f'reactions_rate_{date_str}_{millis}.json')   #
            with open(result_file, "w", encoding="utf-8") as file:
                json.dump(list(shared_data['results_papers']),
                        file, indent=4, ensure_ascii=False)

            result_file2 = os.path.join(
                self.output_path, f'processed_papers_all_{date_str}_{millis}.json')
            with open(result_file2, "w") as file:
                json.dump(
                    list(shared_data['processed_papers_all']), file, indent=4)