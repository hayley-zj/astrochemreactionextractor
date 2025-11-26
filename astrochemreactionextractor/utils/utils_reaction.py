import re
import time
import math
import fitz  # PyMuPDF
import requests
import sympy as sp
import pandas as pd
from rdkit import Chem
from collections import defaultdict
from jellyfish import jaro_winkler_similarity

periodic_table_elements = {
    "J", "L", "c", "s", "E", "l", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs"
}

witelist_elements = ['e-', 'e−', 'e+', 'e', '𝑒−', '𝑒+', '𝑒', 'CR', 'CRP', 'Photon', 'photon', 'hν', 'hν(CR)', 'UV', 'CRUV', 'γ', 'CRPHOT', 'hu', 'hnu', 'h', 'gr', 'grain']

def has_chinese(text):
    return bool(re.search('[\u4e00-\u9fff]', text))

def contains_keywords(s):
    return re.search(r'reaction(s)?', s, re.IGNORECASE) is not None

def need_verify_validity_element(element):
    if len(element) >= 2:
        return False
    if any(char.isdigit() for char in element):
        return False
    return True

def standardize_smiles_list(smiles_list, addH=False):
    standardized_smiles_list = []
    for smiles in smiles_list:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                print(f"Invalid SMILES:{smiles}")
                continue
            if addH:
                mol_with_h = Chem.AddHs(mol)
                standard_smiles = Chem.MolToSmiles(mol_with_h, canonical=True)
            else:
                standard_smiles = Chem.MolToSmiles(mol, canonical=True)
            standardized_smiles_list.append(standard_smiles)
        except Exception as e:
            print(f"Error processing SMILES: {smiles}, error: {e}")
    return standardized_smiles_list


def is_valid_formula(formula_list):
    for formula in formula_list:
        if formula in witelist_elements:
            continue
        try:
            formula = formula.encode('latin1').decode('utf-8')
        except (UnicodeEncodeError, UnicodeDecodeError) as e:
            formula = formula.encode('utf-8').decode('utf-8')

        formula = re.sub(r'^[^-\s]+-', '', formula)   # 去掉如's-'开头的部分
        formula = re.sub(r'\(.*?\)', '', formula)
        formula = re.sub(r'[^a-zA-Z]', '', formula)
        formula = re.sub(r'\s*\([^)]*\)', '', formula)

        length = len(formula)
        if length <= 0:
            return False

        i = 0
        while i < length:
            if i == length - 1:
                if formula[i] not in periodic_table_elements:
                    return False
                i += 1
            else:
                if formula[i:i+2] in periodic_table_elements:
                    i += 2
                elif formula[i] in periodic_table_elements:
                    i += 1
                else:
                    return False
    return True


def check_elements_in_periodic_table(formula_list):
    for element in formula_list:
        if need_verify_validity_element(element) and element not in periodic_table_elements:
            return False
    return True


def get_smiles_from_simple_species(formula):
        return 'None'


def clean_text(text):
    text = re.sub(r'-\n', '', text)
    return text


def extract_doi_based_on_condition(input_string):
    if re.match(r'^\d', input_string):
        parts = input_string.split('.', 2)
        if len(parts) > 2:
            return f"{parts[0]}.{parts[1]}"
        else:
            return input_string
    else:
        parts = input_string.split('.', 1)
        return parts[0]


def get_smiles_from_pubchem(key, chemical_name, max_retries=2, retry_delay=1):
    molecule = re.sub(r'^\d+|\s|^[^-\s]+-', '', key).replace(" ", "")
    molecule_key = parse_molecule(molecule)
    retries = 0
    while retries < max_retries:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"

        try:
            response = requests.get(url)

            if response.status_code == 200:
                data = response.json()
                try:
                    smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                    standardize_non_empty_smiles = standardize_smiles_list(
                        [smiles], addH=True)

                    filtered_smiles_list = [
                        smiles for smiles in standardize_non_empty_smiles if molecule_key == parse_smiles(smiles)]
                    if len(filtered_smiles_list) != 0:
                        standardize_filtered_smiles = standardize_smiles_list(
                            filtered_smiles_list)
                        return standardize_filtered_smiles[0]
                    else:
                        return 'None'
                except (KeyError, IndexError):
                    print('Error: Unable to find SMILES for the given chemical name.')
        except (requests.exceptions.ConnectionError, requests.exceptions.ChunkedEncodingError) as e:
            print(f"ConnectionError occurred: {e}")
            print(f"Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)
        retries += 1
    return 'None'

def get_smiles_from_onepub(url, data_source, max_retries=2, retry_delay=1):
    retries = 0
    while retries < max_retries:
        try:
            response = requests.get(url)

            if response.status_code == 200:
                if data_source == 'CIR':
                    smiles = response.text
                    return smiles
                else:
                    data = response.json()
                    try:
                        if data_source == 'pubchem':
                            smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
                        else:
                            smiles = data['smiles']
                        return smiles
                    except (KeyError, IndexError):
                        print(
                            'Error: Unable to find SMILES for the given chemical name.')
        except Exception as e:
            print(f"Error occurred:{e}")
            print(f"{url} Retrying in {retry_delay} seconds...")
            time.sleep(retry_delay)
        retries += 1
    return 'None'


def extract_nonH(formula):
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    other_atoms = sorted({elem for elem, count in elements if elem != 'H'})
    return other_atoms


def extract_atom(molecule):
    reactant_elements = defaultdict(str)
    molecule = molecule.replace(" ", "")
    molecule = re.sub(r'^[^-\s]+-', '', molecule)
    molecule_elements = parse_molecule(molecule)
    for element, count in molecule_elements.items():
        if reactant_elements[element]:
            reactant_elements[element] += count
        else:
            reactant_elements[element] = count


def parse_smiles(formula):
    formula = formula.replace('c', 'C')
    atom_pattern = r'\[?([A-Z][a-z]*)\]?'
    atom_counts = defaultdict(int)

    stack = []
    i = 0
    while i < len(formula):
        if formula[i] == '(':
            stack.append(i)
        elif formula[i] == ')':
            start = stack.pop()
            i += 1
            continue
        elif formula[i] in '=#':
            i += 1
            continue

        match = re.match(atom_pattern, formula[i:])
        if match:
            atom = match.group(1)
            atom_counts[atom] += 1
            i += match.end()
        else:
            i += 1

    return dict(atom_counts)


def get_smiles_from_multipub(key, chemical_name, max_retries=4, retry_delay=2):
    molecule = re.sub(r'^\d+|\s|^[^-\s]+-', '', key).replace(" ", "")
    molecule_key = parse_molecule(molecule)
    retries = 0
    while retries < max_retries:
        url_1 = f"https://cactus.nci.nih.gov/chemical/structure/{chemical_name}/smiles"
        url_2 = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{chemical_name}/property/CanonicalSMILES/JSON"
        url_3 = f"https://opsin.ch.cam.ac.uk/opsin/{chemical_name}"

        smiles1 = get_smiles_from_onepub(url_1, 'CIR')
        smiles2 = get_smiles_from_onepub(url_2, 'pubchem')
        smiles3 = get_smiles_from_onepub(url_3, 'opsin')

        non_none_smiles = [s for s in [
            smiles1, smiles2, smiles3] if s != 'None']

        standardize_non_empty_smiles = standardize_smiles_list(
            non_none_smiles, addH=True)
        unique_non_empty_smiles = list(set(standardize_non_empty_smiles))

        filtered_smiles_list = [
            smiles for smiles in unique_non_empty_smiles if molecule_key == parse_smiles(smiles)]
        standardize_filtered_smiles = standardize_smiles_list(
            filtered_smiles_list)
        if standardize_filtered_smiles != []:
            output_smiles = '/'.join(standardize_filtered_smiles)
        else:
            output_smiles = 'None'
        return output_smiles

def read_and_split_txt(txt_content, max_tokens):
    # Split content by lines
    lines = txt_content.split('\n')

    chunks = []
    current_chunk = ""

    for line in lines:
        if len(current_chunk) + len(line) + 1 <= max_tokens:
            # Add the line to the current chunk
            if current_chunk:
                current_chunk += '\n'
            current_chunk += line
        else:
            # Add the current chunk to chunks and start a new chunk
            chunks.append(current_chunk)
            current_chunk = line

    # Add the last chunk if it's not empty
    if current_chunk:
        chunks.append(current_chunk)

    merged_chunks = []
    current_chunk = ""

    for chunk in chunks:
        if "→" not in chunk:
            current_chunk += chunk
        else:
            if current_chunk:
                merged_chunks.append(current_chunk)
            merged_chunks.append(chunk)
            current_chunk = ""

    if current_chunk:
        merged_chunks.append(current_chunk)

    return merged_chunks


def extract_formula_from_recation(reaction):
    pattern = re.compile(
        r'\\ce{([^}]+)}'
        r'|(\$?[a-z]?\$?-?[A-Z][a-z]?(?:\$_\d+\$)?(?:[A-Z][a-z]?(?:\$_\d+\$)?)*)'
        r'|([A-Za-z][A-Za-z0-9_{}^+-]*)'
    )

    formula_list = pattern.findall(reaction)
    formatted_formula_list = [
        next((part for part in formula if part), '') for formula in formula_list]

    return formatted_formula_list


def find_formulas_contexts(tex_content, chemical_list, context_size=8):
    contexts = []
    # Split the text into words
    words = tex_content.split()
    for chemical in chemical_list:
        matches = []
        try:
            pattern = re.compile(
                rf'\b(?:{re.escape(chemical)}\b|\({re.escape(chemical)}\b)\)')
        except re.error as e:
            print(f"Regular expression error: {e}")
            return ''
        matche_str = [match.group()
                      for match in pattern.finditer(tex_content)]
        if len(matche_str) != 0:
            matches = [m.start() for m in re.finditer(
                re.escape(matche_str[0]), tex_content)]
        else:
            matches = []

        for match in matches:
            # Find the word index where the match occurs
            word_index = len(re.findall(r'\S+', tex_content[:match]))
            # Extract context around the match
            start = max(word_index - context_size, 0)
            end = min(word_index + context_size, len(words))
            context = ' '.join(words[start:end])
            contexts.append(context)
    retrieval_contexts = '\n'.join([item for item in contexts])
    return retrieval_contexts


def extract_contents_reactions_from_pdf(pdf_document):
    document = fitz.open(pdf_document)

    pdf_reactions = ''
    pdf_contents = ''

    for page_num in range(document.page_count):
        page = document.load_page(page_num)
        text = page.get_text()
        pdf_contents += text
        print('page_num', page_num)

        pattern = re.compile(
            r'([\w\s\−\,\·\⋅\.\*\∗()\+\-]+\s*→\s*[\w\s\−\,\·\⋅\.\*()\+\-]+(?:\s*[+þ]\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+)*)'
            r'|(→\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+(?:\s*\+\s*[\w\s\−\,\·\⋅\.\*\∗()\+\-]+)*)'
        )

        matches = pattern.findall(text)
        if len(matches) != 0:
            formatted_equations = [
                next((part for part in reaction if part), '') for reaction in matches]

            pdf_reactions += '\n' + \
                '\n'.join([item for item in formatted_equations])
    return pdf_contents, pdf_reactions

def has_rate_coefficients(text):
    pattern = r'rate coef'
    match = re.search(pattern, text, re.IGNORECASE)
    if match:
        return True
    return False

def extract_rate_from_page(pdf_path):
    pattern = re.compile(
        r'([\w\s\−\,\·\⋅\.\*\∗()\+\þ\-]+\s*[→⇌])'
    )

    pdf_reactions = ''
    pdf_contents = ''
    try:
        doc = fitz.open(pdf_path)
    except Exception as e: 
        print(f"Exception: {e}")
        return pd.DataFrame(), pdf_contents
    data = []

    reactions_lastpage = ''
    lastpage_has_reaction = 0
    for page_num in range(len(doc)):
        reactions_lastpage_update = ''
        page = doc.load_page(page_num)
        width_page = page.rect.width
        height_page = page.rect.height
        page_contents = page.get_text()
        pdf_contents += page_contents

        if '→' not in page_contents and '⇌' not in page_contents:
            if has_rate_coefficients(page_contents) and lastpage_has_reaction:
                reactions_lastpage_update = reactions_lastpage + page_contents
                data.pop()
                data.append([reactions_lastpage_update, page_num - 1])
                lastpage_has_reaction = 0
            elif has_rate_coefficients(page_contents) and not lastpage_has_reaction:
                data.append([page_contents, page_num])
            continue
        matches = list(pattern.finditer(page_contents.replace(' (', '(')))
        if len(matches) != 0:
            text_before_reaction = page_contents[0: matches[0].start()]
            if has_rate_coefficients(text_before_reaction) and lastpage_has_reaction:
                reactions_lastpage_update = reactions_lastpage + text_before_reaction
                data.pop()
                data.append([reactions_lastpage_update, page_num - 1])
        if len(matches) == 1:
            text_after_reaction = page_contents[matches[0].start():]
            if has_rate_coefficients(page_contents):
                page_reactions =  page_contents
            else:
                page_reactions = page_contents[matches[0].start(
                ): matches[0].end() + 500]
            data.append([page_reactions, page_num])
            reactions_lastpage = page_reactions
            lastpage_has_reaction = 1
        elif len(matches) > 1:
            if 'Table' in page_contents:
                page_reactions = page_contents
            else:
                page_reactions = page_contents[matches[0].start():]
            data.append([page_reactions, page_num])
            reactions_lastpage = page_reactions
            lastpage_has_reaction = 1
        else:
            lastpage_has_reaction = 0
    df = pd.DataFrame(data, columns=['reactions_txt', 'coordinates'])
    return df, pdf_contents


def parse_molecule(molecule):
    pattern = r'(\([^)]*\))(\d*)'
    matches = re.findall(pattern, molecule)
    if matches:
        ddd = matches[0][0].replace('(', '').replace(')', '')
        numb = [i + matches[0][1] for i in ddd]
        bracket_str = ''.join(numb)

        remaining_str = re.sub(pattern, '', molecule)
        molecule = remaining_str + bracket_str

    molecule = molecule.replace('𝑛', 'n')
    elements = defaultdict(str)
    i = 0
    n = len(molecule)
    try:
        while i < n:
            if molecule[i] == '(' or molecule[i] == ')':
                i = i + 1
                continue
            else:
                element = molecule[i]
                i += 1
                while i < n and molecule[i].islower():
                    element += molecule[i]
                    i += 1
                count = ''
                while i < n and (molecule[i].isdigit() or molecule[i] in {'n', '+', '-'}):
                    count += molecule[i]
                    i += 1
                result_count = count if count and count != '+' and count != '-' else 1
                if len(element) == 1 and element.islower():
                    continue
                if re.search(r'n(\d)', count):
                    return elements
                if 'n' in count:
                    count_ex = count.replace('n', '*n')
                    try:
                        result_count = sp.sympify(
                            count_ex.replace('^', '**'), locals=dict(n=1))
                    except:
                        continue
                else:
                    if result_count != 1:
                        result_count = result_count.replace('+', '').replace('-', '')
                    if result_count == '':
                        continue
                try:
                    if isinstance(result_count, str):
                        result_count = ''.join(char for char in result_count if char not in '⁰¹²³⁴⁵⁶⁷⁸⁹')
                        result_count = result_count.translate(str.maketrans('₁₂₃₄₅₆₇₈₉', '123456789'))
                    if elements[element]:
                        elements[element] += int(result_count)
                    else:
                        elements[element] = int(result_count)
                except:
                    print("parse_molecule failure")
                    elements[element] = 1
    except:
        print("parse_molecule failure")
        elements[element] = 1
    return elements

def check_params_rate(params_str_list):
    keywords_to_check = ["alpha:", "α:", "rate coefficient", "EA:"]
    rate_str_list = []
    for param_str in params_str_list:
        result = {}
        if not isinstance(param_str['params'], str) and math.isnan(param_str['params']):
            continue
        for keyword in keywords_to_check:
            if re.search(keyword, param_str['params'], re.IGNORECASE):
                rate_str_list.append(param_str)
    return rate_str_list


def completion_reaction(df_results):
    last_reactant = None
    for i, row in df_results.iterrows():
        reaction = row['reactions_wx']
        if not isinstance(reaction, str) and math.isnan(reaction):
            continue
        reaction = reaction.replace('→', '->')
        reaction = re.sub(r'^\d+\.\s*', '', reaction)
        match = re.match(r'^->(.*->.*)$', reaction)
        if match:
            reaction = match.group(1).strip()
        if reaction.endswith('->') and reaction.count('->') == 1:
            reaction = '-> ' + reaction[:-2].strip()

        if reaction.startswith('->'):
            if last_reactant is not None:
                df_results.at[i,
                              'reactions_wx'] = f"{last_reactant} {reaction}"
        else:
            last_reactant = reaction.split('->')[0].strip()
            df_results.at[i, 'reactions_wx'] = reaction
    return df_results


def add_rate_coefficient(df_results, formula_dict):
    filtered_dict = {key.replace("Type ", "").replace("type ", ""): value for key, value in formula_dict.items() if '=' in value and len(value) >= 6}

    deduplicated_dict = {}
    for key, value in filtered_dict.items():
        if key.count(" and ") == 1 or key.count("/") == 1:
            if key.count(" and ") == 1:
                part1, part2 = key.split(" and ")
            if key.count("/") == 1:
                part1, part2 = key.split("/")
            deduplicated_dict[part1] = value
            deduplicated_dict[part2] = value
        else:
            deduplicated_dict[key] = value

    temp_dict = {}
    for key, value in deduplicated_dict.items():
        norm_key = re.sub(r'[\s-]reactions?', '', key, flags=re.IGNORECASE)
        if norm_key not in temp_dict or len(value) > len(temp_dict[norm_key][1]):
            temp_dict[norm_key] = (key, value)
    new_deduplicated_dict = {
        original_key: value for original_key, value in temp_dict.values()}
    if new_deduplicated_dict == {}:
        return df_results

    for i, row in df_results.iterrows():
        rate_coff = ''
        filtered_pairs = []
        data_params = {}
        for item in row['params'].split(';'):
            if ':' not in item:
                continue
            key, value = item.strip().split(':', 1)
            if value == '':
                continue
            if key.strip() not in ('reference', 'reaction ID', 'Width', 'Ref.'):
                data_params[key] = value
                filtered_pairs.append(f"{key.strip()}:{value.strip()}")
        filtered_string = '; '.join(filtered_pairs)
        if 'reaction type' in data_params or 'Reaction Type' in data_params:
            reaction_type = data_params.get(
                'reaction type') or data_params.get('Reaction Type')
            if reaction_type != None:
                reaction_type = re.sub(r' reactions?', '', reaction_type.strip(), flags=re.IGNORECASE)
                matches = [
                    (key, jaro_winkler_similarity(reaction_type.lower(), key.lower()))
                    for key in new_deduplicated_dict.keys()
                ]
                matches.sort(key=lambda x: -x[1])
                if matches[0][1] > 0.66:
                    rate_coff = new_deduplicated_dict[matches[0][0]]
        if rate_coff != '' and 'calculation formula' not in row['params']:
            df_results.at[i, 'params'] = filtered_string + \
                f'; calculation formula:{rate_coff}'
    return df_results

